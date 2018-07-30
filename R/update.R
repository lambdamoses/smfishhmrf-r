#Update probs
updateProbabilities <- function(neighbors, numnei, indices, den, k, beta){
    Nij <- matrix(0, nrow=nrow(neighbors), ncol=k)
	for(ii in 1:nrow(neighbors)){
		for(ij in 1:numnei[ii]){
			for(j in 1:k){
				Nij[ii,j] <- Nij[ii,j] + indices[neighbors[ii,ij],j]
			}
		}
	}

    #for(m in 1:nneigh){
    #    Nij <- Nij + indices[neighbors[,m],]
    #}
    #or <- Nij[,1]
    #for(i in 2:k)
    #    or <- or + Nij[,i]*(nneigh+1)^(i-1)
    #or <- or + 1

	col_sum <- colSums(indices) / nrow(neighbors)

    #prob <- den*check[or,]
	prob <- den * exp(beta * Nij)
	#prob <- den * exp(beta * Nij) * col_sum
	#sum_prob <- colSums(prob)

    prob <- prob/rowSums(prob)
	r_prob <- rowSums(prob)
	num_na <- length(r_prob[is.na(r_prob)])
	if(any(is.na(rowSums(prob)))){
		#print("is infinity!")
		print(paste(num_na, "is infinity", sep=" "))
	}
    prob[is.na(rowSums(prob)),] <- rep(1/k, k)
    prob[prob==Inf] <- 1/k
    prob
}

#Update probs
updateUnnormProbabilities <- function(neighbors, numnei, indices, den, k, beta){
    Nij <- matrix(0, nrow=nrow(neighbors), ncol=k)
	for(ii in 1:nrow(neighbors)){
		for(ij in 1:numnei[ii]){
			for(j in 1:k){
				Nij[ii,j] <- Nij[ii,j] + indices[neighbors[ii,ij],j]
			}
		}
	}

	prob <- den * exp(beta * Nij)
	#prob <- den * exp(beta * Nij) * col_sum

	#if(any(is.na(rowSums(prob)))){
	#	print("is infinity!")
	#}
    #prob[is.na(rowSums(prob)),] <- rep(1/k, k)
    #prob[prob==Inf] <- 1/k
    prob
}

#Update mus
updateMeans <- function(prob, y){
    inten <- y * prob
    mu <- colSums(inten)/colSums(prob)
    mu
}

#Update sds
updateStdevs <- function(prob, y, mu, nvert, k){
    mu <- matrix(mu, nrow=nvert, ncol=k, byrow=T)
    diff2 <- (y-mu)^2*prob
    sigma2 <- colSums(diff2)/colSums(prob)
    sqrt(sigma2)
}

updateIndicesHMRFEM <- function(neighbors, numnei, maxnei, blocks, nblocks, beta, k, indices, den){
    .Call("updateIndicesHMRFEM", blocks, neighbors, numnei, maxnei, beta, k, indices, den)
}


#Update probs
updateProbabilitiesMulti <- function(neighbors, numnei, indices, den, k, beta){
	updateProbabilities(neighbors, numnei, indices, den, k, beta)
}

updateUnnormProbabilitiesMulti <- function(neighbors, numnei, indices, den, k, beta){
	updateUnnormProbabilities(neighbors, numnei, indices, den, k, beta)
}

#Update mus
updateMeansMulti <- function(prob, y){
	#prob is nvert * k, y is nvert * numdim
	nvert = dim(y)[1]
	numdim = dim(y)[2]
	k = dim(prob)[2]
	mu = array(0, c(numdim, k))
	for(i in 1:k){
		mean_prob = mean(prob[,i])
		t_prob = matrix(rep(prob[,i],numdim), nrow=nvert, ncol=numdim, byrow=F)
		mu[,i] = colSums(t_prob * y) / (nvert * mean_prob)
	}
	mu
}
#Update sds
updateCovariancesMulti <- function(prob, y, mu, nvert, k){
	#prob is nvert * k, y is nvert * numdim
	nvert = dim(y)[1]
	numdim = dim(y)[2]
	k = dim(prob)[2]
	sigma2 = array(0, c(numdim, numdim, k))
	for(i in 1:k){
		dist <- y - matrix(rep(t(mu[,i]),nvert),nrow=nvert,ncol=numdim,byrow=T)
		mean_prob = mean(prob[,i])
		t_prob = matrix(rep(prob[,i],numdim), nrow=nvert, ncol=numdim, byrow=F)
		sigma2[,,i] <- t(t_prob * dist) %*% dist / (nvert * mean_prob)
	}
	sigma2
}

