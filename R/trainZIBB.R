trainZIBB <- function(y, X, ziMatrix, s, lambda, grpIdx, grpSize, K, rho, 
											para.init = NULL, tol.fista = 1e-4, imax.fista = 20, 
											tol.alm = 1e-4, imax.alm = 20, mode = "ZIBB", group) {
	n <- dim(X)[1]
	p <- dim(X)[2]
	q <- dim(ziMatrix)[2]
	m <- dim(y)[1]
	Y <- t(y)
	Y.c <- as.vector(apply(y, 2, sum))
	
	if (mode == "ZIBB") {
	  # Initial values: Beta.out, eta.out, and psi.out
	  if (is.null(para.init)) {
	    zibb.free <- fitZIBB(y, X, ziMatrix, mode = "free")
	    zibb.constrained <- fitZIBB(y, X, ziMatrix, mode = "constrained",
	                                betastart = zibb.free$betahat,
	                                psi.start = zibb.free$psi,
	                                eta.start = zibb.free$zeroCoef)
	    #*Re-parameterize over-dispersion
	    psi <- zibb.constrained$psi
	    # phi <- exp(psi) / (1 + exp(psi))
	    tmp.mat <- rbind(as.matrix(zibb.constrained$betahat),
	                     matrix(psi, 1, m),
	                     as.matrix(zibb.constrained$zeroCoef))
	    para.init <- c(as.vector(tmp.mat), as.vector(zibb.constrained$gamma))
	  }
	  
	  # print("Finish initialization!")
	  # print(paste("Number of NA's/NaN's is", sum(is.na(para.init))))
	  
	  # ALM algorithm
	  para.out <- alm.alg(para.init, s, lambda, grpIdx, grpSize, X, Y, Y.c, 
	                      ziMatrix, K, rho, tol.fista, imax.fista, 
	                      tol.alm, imax.alm)
	  
	  return(para.out)
	} else if (mode == "lasso") {
	  resp.train <- X[,p]
	  fit.lasso <- glmnet::glmnet(x=Y, y=resp.train, family="gaussian", alpha=1, 
	                              lambda=lambda, intercept = FALSE)
	  return(fit.lasso)
	} else if (mode == "glasso") {
	  resp.train <- X[,p]
	  fit.gl <- gglasso::gglasso(x=Y, y=resp.train, group=group, loss="ls", 
	                             lambda=lambda, intercept = FALSE)
	  return(fit.gl)
	} else {
	  print("Unsupported mode!")
	}
}