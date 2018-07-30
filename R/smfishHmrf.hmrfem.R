
smfishHmrf.hmrfem <- function(y, neighbors, numnei, blocks, beta=0.5, mu, sigma, err=1e-4, maxit=20, verbose){

    checkErrors(mu=mu, sigma=sigma, err=err)

    if (length(err) < 2) 
		err <- rep(err, length.out = 2)

    k <- length(mu)
    nvert <- length(y)

    maxnei <- ncol(neighbors)
    nblocks <- length(blocks)
    neighbors <- structure(as.integer(neighbors), dim = dim(neighbors))

    yunique <- sort(unique(y))
    n.yunique <- length(yunique)
    nvert <- length(y)
    ymatch <- match(y, yunique)

    muold <- rep(0,k)
    sigmaold <- rep(0,k)
    niter <- 10
    indices <- initialIndices(y, nvert, mu, sigma, k, sub=FALSE)

    it <- 0

    repeat{
        it <- it + 1
        muold <- mu
        sigmaold <- sigma
        for(j in 1:niter){
            den <- getDensity(yunique, n.yunique, ymatch, mu, sigma)
            #indices <- updateIndicesHMRFEM(neighbors, numnei, maxnei, blocks, nblocks, beta, k, indices, check, den)
            indices <- updateIndicesHMRFEM(neighbors, numnei, maxnei, blocks, nblocks, beta, k, indices, den)
        }
        prob  <- updateProbabilities(neighbors, numnei, indices, den, k, beta)
        mu <- updateMeans(prob, y)
        sigma <- updateStdevs(prob, y, mu, nvert, k)

        flag <- checkStopVerbose(muold, mu, sigmaold, sigma, err, it, maxit, verbose)
        if(flag==1) break
    }

    list(prob=prob, mu=mu, sigma=sigma)
}

smfishHmrf.hmrfem.multi <- function(y, neighbors, numnei, blocks, beta=0.5, mu, sigma, err=1e-4, maxit=20, verbose, dampFactor=NULL, forceDetectDamp=FALSE){
    checkErrorsMulti(mu=mu, sigma=sigma, err=err)
    if (length(err) < 2)
		err <- rep(err, length.out = 2)

	k <- dim(mu)[2] #mu is a m * k matrix
	nvert <- dim(y)[1] #number of data points, n
	numdim <- dim(y)[2] #number of dimensions, m
    maxnei <- ncol(neighbors)
    nblocks <- length(blocks)
    neighbors <- structure(as.integer(neighbors), dim = dim(neighbors))
    muold <- array(0, dim(mu)) # m * k matrix
    sigmaold <- array(0, dim(sigma)) # m * m * k matrix 
    niter <- 10

    indices <- initialIndicesMulti(y, mu, sigma, k, dampFactor=dampFactor, forceDetectDamp=forceDetectDamp)

    it <- 0
	eps <- eps()
	very_small <- mean(std(y)) * (1/k) * 0.0001

    repeat{
        it <- it + 1
        muold <- mu
        sigmaold <- sigma
		den = matrix(0, nrow=nvert, ncol=k)
		for(j in 1:k){
			determinant <- det(sigma[,,j])
			invsigma <- array(0, dim(sigma[,,j]))
			#choice 1======================
			if(forceDetectDamp==FALSE & !is.null(dampFactor)){
				sigma[,,j]<-sigma[,,j] + diag(numdim) * dampFactor[j]
				invsigma <- solve(sigma[,,j])
				determinant <- det(sigma[,,j])
			}
			else if(forceDetectDamp==TRUE){
				damp<-findDampFactor(sigma[,,j], factor=1.05, d_cutoff=1e-60, startValue=0.0001)
				if(!is.null(damp)){
					sigma[,,j]<-sigma[,,j] + diag(numdim) * damp
					invsigma <- solve(sigma[,,j])
					determinant <- det(sigma[,,j])
					#print(damp)
				}
			}
			dist <- y - matrix(rep(t(mu[,j]),nvert),nrow=nvert,ncol=numdim,byrow=T)
			exponent <- -0.5 * (rowSums((dist %*% invsigma) * dist))
			const <- (numdim/2)*log(2*pi) + (1/2)*log(determinant)
			exp_diff <- exponent - const
			lim_exp_diff <- sapply(exp_diff, function(x) min(700, max(x, -700)))
			den[,j] <- exp(lim_exp_diff)
		}
        for(j in 1:niter){
            indices <- updateIndicesHMRFEM(neighbors, numnei, maxnei, blocks, nblocks, beta, k, indices, den)
        }
		unnorm_prob <- updateUnnormProbabilitiesMulti(neighbors, numnei, indices, den, k, beta)
        prob  <- updateProbabilitiesMulti(neighbors, numnei, indices, den, k, beta)
        mu <- updateMeansMulti(prob, y)
        sigma <- updateCovariancesMulti(prob, y, mu, nvert, k)

		#need to check this!
        flag <- checkStopVerboseMulti(muold, mu, sigmaold, sigma, err, it, maxit, verbose)

        if(flag==1) break
    }
    list(prob=prob, mu=mu, sigma=sigma, unnormprob=unnorm_prob)
}

