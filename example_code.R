library(survival)
library(cmprsk)
library(caret)

#####################################################################################################################

find_id <- function( knots, t){
	knots = sort(knots, decreasing=FALSE)
	tid0 = which.min( abs( t - knots ) )
	if( t - knots[tid0] >= 0 ){
		return( tid0 )
	} else{
		return( tid0 - 1 )
	}
}

i_to_id <- function( i, ns ){
	i = i - 1
	m = length(ns)
	idxs = rep(1, m)
	for( j in 1:(m-1) ){
		idxs[j] = floor( i / prod( ns[ (j+1):m ] ) )
		i = i %% prod( ns[ (j+1):m ] )
	}
	idxs[m] = i
	idxs = idxs + 1
	return(idxs)
}

concord_ipw <- function( score_l, l, n_l, G_hat, types ){
	# G_hat is modified, so that for x > tau, G_hat is G(tau|Z)
	ns = rep(1, n_l)
	for( ll in 1:n_l ){
		tmp = score_l[ types[[ ll ]] ]
		types[[ ll ]] = types[[ ll ]][ order( tmp, decreasing=TRUE ) ]
		ns[ll] = length(types[[ ll ]])
	}
	inv_G_hat = 1 / G_hat
	ls_sum = rep(1, n_l)
	for( ll in 1:n_l ){
		ls_sum[ll] = sum( inv_G_hat[ types[[ ll ]] ] )
	}
	idxs = rep(1, n_l); cnt = 0
	eq_sum = rep(0, n_l); choice_seq = rep(1, n_l)
	while( idxs[l] <= ns[l] ){
		eq_sum[l] = inv_G_hat[ types[[ l ]][ idxs[l] ] ]
		ls_sum[l] = ls_sum[l] - inv_G_hat[ types[[ l ]][ idxs[l] ] ]
		while( idxs[l] + 1 <= ns[l] ){
			if( score_l[ types[[ l ]][ idxs[l] ] ] == score_l[ types[[ l ]][ idxs[l] + 1 ] ] ){
				idxs[l] = idxs[l] + 1
				eq_sum[l] = eq_sum[l] + inv_G_hat[ types[[ l ]][ idxs[l] ] ]
				ls_sum[l] = ls_sum[l] - inv_G_hat[ types[[ l ]][ idxs[l] ] ]
			}else{
				break
			}
		}
		for( ll in 1:n_l ){
			if( ll != l ){
				while( idxs[ll] <= ns[ll] & score_l[ types[[ ll ]][ idxs[ll] ] ] > score_l[ types[[ l ]][ idxs[l] ] ] ){
					ls_sum[ll] = ls_sum[ll] - inv_G_hat[ types[[ ll ]][ idxs[ll] ] ]
					idxs[ll] = idxs[ll] + 1
				}
				if( idxs[ll] > ns[ll] ){
					return(cnt)
				}
				eq_sum[ll] = 0
				while( idxs[ll] <= ns[ll] ){
					if( score_l[ types[[ l ]][ idxs[l] ] ] == score_l[ types[[ ll ]][ idxs[ll] ] ] ){
						eq_sum[ll] = eq_sum[ll] + inv_G_hat[ types[[ ll ]][ idxs[ll] ] ]; 
						ls_sum[ll] = ls_sum[ll] - inv_G_hat[ types[[ ll ]][ idxs[ll] ] ]; 
						choice_seq[ll] = 2
						idxs[ll] = idxs[ll] + 1
					}else{
						break
					}
				}
			}
		}
		tmp = 0;
		for( i in 1:prod( choice_seq ) ){
			neq = 1
			idl = i_to_id( i, choice_seq )
			tmp2 = eq_sum[l]
			for( ll in 1:n_l ){
				if( ll != l ){
					if( idl[ll] == 1 ){
						tmp2 = tmp2 * ls_sum[ll]
					}else{
						tmp2 = tmp2 * eq_sum[ll]
						neq = neq + 1
					}
				}
			}
			tmp = tmp + tmp2 / neq
		}
		cnt = cnt + tmp
		idxs[l] = idxs[l] + 1
	}
	return(cnt)
}

pdi <- function( tau, scores, x, delt_eps, G_hat, n_l ){
	# Compute PDI from given prognostic scores for n subjects, algorithm allows tie
	# scores: a n by n_l matrix of prognostic scores
	# x: a vector of observation times
	# delt_eps: a vector of event indicators
	# G_hat: a vector of estimated inverse weights
  
	res = rep(0, n_l);
	types = list();	ns = rep(1, n_l); denom = 1
	for( ll in 1:(n_l - 1) ){
		types[[ ll ]] = which( x <= tau & delt_eps == ll)
		ns[ll] = length( types[[ ll ]] )
		denom = denom * sum( 1 / G_hat[ types[[ ll ]] ] )
	}
	types[[ n_l ]] = which( x > tau )
	ns[n_l] = length( types[[ n_l ]] )
	denom = denom * sum( 1 / G_hat[ types[[ n_l ]] ] )
	if( any( ns == 0 ) ){
		return( rep(NA, n_l) )
	}else{
		for( l in 1:n_l ){
			numer = concord_ipw( scores[, l], l, n_l, G_hat, types )
			res[l] = numer / denom
		}
	}
	return( res )
}

comp_G_hat <- function(dat, tau){
  #function for calculating inverse weights using the KM method
	n = nrow(dat)
	km_fit <- survfit( Surv( x, delt_eps==0 ) ~ 1, type="kaplan-meier", data=dat )
	Gfunc <- stepfun(km_fit$time, c(1, km_fit$surv))
	G_hat <- Gfunc( dat$x )
	G_hat[ which( dat$x > tau ) ] = Gfunc( tau )
	return(G_hat)
}

comp_scores <- function( mdat, pdat, zcov, tau, n_l ){
	# fit Fine-Gray model on mdat and use the fitted model to generate prediction for pmat
	n = nrow(pdat)
	scores = matrix(1, nrow = n, ncol = n_l)
	for( l in 1:( n_l - 1 ) ){
		knots = unique( mdat$x[ mdat$delt_eps==l ] )
		tidx = find_id( knots, tau )
		if( tidx == 0 ){
			return( rep( 0, n ) )
		} else{
			fg <- crr(mdat$x, mdat$delt_eps, mdat[, zcov], failcode = l, cencode = 0)
			scores[, l] = sapply( 1:n, function(i) predict( fg, pdat[i, zcov] )[ tidx, 2 ] )
		}
		scores[, n_l] = scores[, n_l] - scores[, l]
	}	
	return(scores)
}

comp_pdi <- function(mdat, pdat, zcov, tau, n_l){
	# Compute inverse weights
	G_hat = comp_G_hat( dat=pdat, tau=tau )	

	# predict prognostic scores for each type out of outcome at time tau
	scores = comp_scores( mdat=mdat, pdat=pdat, zcov=zcov, tau=tau, n_l=n_l )

	# Compute PDI
	out = pdi( tau=tau, scores=scores, x=pdat$x, delt_eps=pdat$delt_eps, G_hat=G_hat, n_l=n_l )

	return(out)
}

fine_gray_pdi <- function( dat, tau, n_l, seed, cv_rep, zcov){
	# main function for computing PDI with competing risks data using Fine-Gray model
	# tau: the prespcified inspection time
	# n_l: number of outcome categories, including event-free
	# seed: random seed for cross validation
	# cv_rep: number of cross validation repetitions
        # zcov: covariates used in the FG model
  
	set.seed(seed)
	n = nrow(dat)
	pdis = rep(0, n_l)

	# Compute PDI
	for( i in 1:cv_rep ){
		tmp = createFolds( 1:n, k=2 )
		dat1 = dat[ tmp[[1]], ]; dat2 = dat[ tmp[[2]], ]	
		pdis = pdis + comp_pdi(mdat=dat1, pdat=dat2, zcov=zcov, tau=tau, n_l=n_l)
		pdis = pdis + comp_pdi(mdat=dat2, pdat=dat1, zcov=zcov, tau=tau, n_l=n_l)
	}
	pdis = pdis / ( cv_rep * 2 )
	
	set.seed(NULL)
	return( list(pdi_specific=pdis, pdi_overall=mean(pdis) ) )
}

##################################################################################################
################################      Example Codes    ###########################################
##################################################################################################

dat = read.csv( "...\\example_data.csv", header=TRUE )
n_l = 3  # number of outcome categories
tau = 1.4 # prespecified time for prediction

PDIs = fine_gray_pdi( dat=dat, tau=tau, n_l=n_l, seed=1, cv_rep=3, zcov=c("z1","z2") )
PDIs
