#check errors of input of tissue classi 
checkErrors <- function(mu, sigma, err=NULL){
    k <- length(mu)
    if (k != length(sigma))
        stop("The dimensions of 'mu' and 'sigma' do NOT match.")
    if(! is.null(err))
        if (any(err < 0)) stop("All 'err's have to be positive.")
}

#Compute the density of all vertices.
getDen <- function(yunique, n.yunique, ymatch, mu, sigma){
    k <- length(mu)
    dy <- rep(yunique, each=k)
    dmu <- rep(mu, n.yunique)
    dsigma <- rep(sigma, n.yunique)
    den <-  dnorm(dy, mean = dmu, sd = dsigma)
    den <- matrix(den, ncol=k, byrow=T)
    den <- den[ymatch,]
    den
}

#compute the relative errors
relerr <- function(x, y){
    max(abs(x - y)) / (1 + max(abs(x), abs(y)))
}

#Initialize indices
initialIndices <- function(y, nvert, mu, sigma, k){
    index <- max.col(matrix(dnorm(rep(y,k), mean=rep(mu, each=nvert),
                                      sd=rep(sigma, each=nvert)), ncol=k))
    indices <- do.call(cbind, lapply(1:k, function(x) indices==x))
    rbind(indices, rep(0L,k))
}
 
#check whether to stop the interation or not and whether ouput the curren state
checkStopVerbose <- function(muold, mu, sigmaold, sigma, err, it, maxit, verbose){
    flag <- 0
    remu <- relerr(muold, mu)
    resigma <- relerr(sigmaold, sigma)
    if (remu < err[1] && resigma < err[2] || it > maxit)
		flag <- 1
	if (verbose)
		cat(paste("Iteration ", it, ": relative error of mean = ", signif(remu, 1),
		"; sigma = ", signif(resigma, 1), "\n", sep=""))
    flag
}


#Initialize indices (multivariate), assume sub=FALSE, type="pure", subvox=NULL
initialIndicesMulti <- function(y, mu, sigma, k, dampFactor=NULL, forceDetectDamp=FALSE){
	nvert <- dim(y)[1]
	numdim <- dim(y)[2]

	eps <- eps()
	very_small <- mean(std(y)) * (1/k) * 0.0001
	
	den = matrix(0, nrow=nvert, ncol=k)

	for(j in 1:k){
		determinant <- det(sigma[,,j])
		invsigma <- array(0, dim(sigma[,,j]))
		#choice 1=============
		if(forceDetectDamp==FALSE & !is.null(dampFactor)){
			sigma[,,j]<-sigma[,,j] + diag(numdim) * dampFactor[j]
			invsigma <- solve(sigma[,,j])
			determinant <- det(sigma[,,j])
			
		}else if(forceDetectDamp==TRUE){
			damp<-findDampFactor(sigma[,,j], factor=1.05, d_cutoff=1e-60, startValue=0.0001)
			if(!is.null(damp)){
				sigma[,,j]<-sigma[,,j] + diag(numdim) * damp
				invsigma <- solve(sigma[,,j])
				determinant <- det(sigma[,,j])
			}
		}
		dist <- y - matrix(rep(t(mu[,j]),nvert),nrow=nvert,ncol=numdim,byrow=T)
		exponent <- -0.5 * (rowSums((dist %*% invsigma) * dist))
		const <- (numdim/2)*log(2*pi) + (1/2)*log(determinant)
		exp_diff <- exponent - const
		lim_exp_diff <- sapply(exp_diff, function(x) min(700, max(x, -700)))
		den[,j] <- exp(lim_exp_diff)
	}
	
	indices <- max.col(den)

    indices <- do.call(cbind, lapply(1:k, function(x) indices==x))
    rbind(indices, rep(0L,k))
}

relerr.multi <- function(x, y){
    max(abs(x - y)) / (1 + max(abs(x), abs(y)))
}

findDampFactor <- function(sigma, factor=1.05, d_cutoff=1e-60, startValue=0.0001){
	determinant <- det(sigma)
	#factor <- 1.05
	#d_cutoff<- 1e-60
	numdim<-dim(sigma)[1]
	damp<-NULL
	if(determinant<d_cutoff){
		print("Warning, need to add dampening factor to diagonal entries")
		damp <- startValue
		sigma_2 <- array(0, dim(sigma))
		sigma_2 <- sigma + diag(numdim) * damp
		determinant <- det(sigma_2)
		if(determinant<d_cutoff){
			repeat{
				damp <- damp * factor
				sigma_2 <- array(0, dim(sigma))
				sigma_2 <- sigma + diag(numdim) * damp
				determinant <- det(sigma_2)
				if(determinant>d_cutoff) break
			}
		}else{
			repeat{
				damp <- damp / factor
				sigma_2 <- array(0, dim(sigma))
				sigma_2 <- sigma + diag(numdim) * damp
				determinant <- det(sigma_2)
				if(determinant<d_cutoff){
					damp <- damp * factor
					break
				}
			}
		}
		print(damp)
	}
	damp
}

#check whether to stop the interation or not and whether ouput the curren state
checkStopVerboseMulti <- function(muold, mu, sigmaold, sigma, err, it, maxit, verbose){
    flag <- 0
    remu <- relerr.multi(muold, mu)
    resigma <- relerr.multi(sigmaold, sigma)
    if (remu < err[1] && resigma < err[2] || it > maxit)
        flag <- 1
    if (verbose)
        cat(paste("Iteration ", it, ": relative error of mean = ", signif(remu, 1),
        "; covariance = ", signif(resigma, 1), "\n", sep=""))
    flag
}
    
#check errors of input of tissue classi 
checkErrorsMulti <- function(mu, sigma, err=NULL){
	k <- dim(mu)[2]
	if (k != dim(sigma)[3])
		stop("The dimensions of 'mu' and 'sigma' do NOT match.")
	#if (!all(sigma >= 0))
	#	stop("All 'sigma's have to be positive")
	if(! is.null(err))
		if (any(err < 0)) stop("All 'err's have to be positive.")
}
