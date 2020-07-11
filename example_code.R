library(survival)

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
	# Compute PDI from given prognostic scores for n subjects
	# tau: the prespcified inspection time
	# n_l: number of outcome types, including event-free
	# scores: a n by n_l matrix of prognostic scores
	# x: a vector of observation times
	# delt_eps: a vector of failure types for each subject, with 0 for censored
	# G_hat: a vector of estimated survival time of censoring time for each subject, when observation time x > tau, replace G_hat by 1 ( KM esitmator ) or G_hat(tau | Z) ( covariate dependent estimator )
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

##################################################################################################
################################      Example Codes    ###########################################
##################################################################################################

dat = read.csv( "example_data.csv", header=TRUE )
n = nrow(dat); n_l = 3  # we have 3 types of outcome, with 1200 subjects
tau = 1.4

# Compute KM estimation of survival function of censoring time
km_fit <- survfit( Surv( x, delt==FALSE ) ~ 1, type="kaplan-meier", data=dat )
Gfunc <- stepfun(km_fit$time, c(1, km_fit$surv))
G_hat <- Gfunc( dat$x )
# when x > tau, G_hat(x) is replaced by G_hat(tau). They will appear both on the numerator and denominator, thus they will be canceling out. We can simply replace them by 1.
G_hat[ which( dat$x > tau ) ] = 1

# predict prognostic scores for each type out of outcome at time tau
zcov = c( "z1", "z2" ) # covariates used in the FG model
scores = matrix(1, nrow = n, ncol = n_l)
for( l in 1:( n_l - 1 ) ){
	knots = unique( dat$x[ dat$delt_eps==l ] )
	tidx = find_id( knots, tau )
	if( tidx == 0 ){
		return( rep( 0, n ) )
	} else{
		fg <- crr(dat$x, dat$delt_eps, dat[, zcov], failcode = l, cencode = 0)
		scores[, l] = sapply( 1:n, function(i) predict( fg, dat[i, zcov] )[ tidx, 2 ] )
	}
	scores[, n_l] = scores[, n_l] - scores[, l]
}

# Compute PDI
PDIs = pdi( tau=tau, scores=scores, x=dat$x, delt_eps=dat$delt_eps, G_hat=G_hat, n_l=n_l )

# > PDIs
# [1] 0.6303349 0.5842413 0.5630625

# incorporating incorrect covariate z3
zcov = c( "z1", "z3" ) # covariates used in the FG model
scores = matrix(1, nrow = n, ncol = n_l)
for( l in 1:( n_l - 1 ) ){
	knots = unique( dat$x[ dat$delt_eps==l ] )
	tidx = find_id( knots, tau )
	if( tidx == 0 ){
		return( rep( 0, n ) )
	} else{
		fg <- crr(dat$x, dat$delt_eps, dat[, zcov], failcode = l, cencode = 0)
		scores[, l] = sapply( 1:n, function(i) predict( fg, dat[i, zcov] )[ tidx, 2 ] )
	}
	scores[, n_l] = scores[, n_l] - scores[, l]
}

# Compute PDI
PDIs = pdi( tau=tau, scores=scores, x=dat$x, delt_eps=dat$delt_eps, G_hat=G_hat, n_l=n_l )

# > PDIs
# [1] 0.4966044 0.5914107 0.4068672