predictionZIBB <- function(y, X, ziMatrix, para, mode, case = 0, control = 1, 
                           mean = 0, sd = 1, lambda, fit.lasso, fit.glasso) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- dim(ziMatrix)[2]
  m <- dim(y)[1]
  Y <- t(y)
  Y.c <- as.vector(apply(y, 2, sum))
  
  if (mode == "discrete") {
    pred.out <- rep(0, n)
    X0 <- X
    X0[,p] <- case
    X1 <- X
    X1[,p] <- control
    tmp.vec0 <- loglikelihood(X0, Y, Y.c, ziMatrix, para)
    tmp.vec1 <- loglikelihood(X1, Y, Y.c, ziMatrix, para)
    tmp.val <- exp(tmp.vec1 - tmp.vec0)
    idx.large <- (tmp.vec1 - tmp.vec0) > 100
    idx.small <- (tmp.vec1 - tmp.vec0) < -100
    tmp.val[idx.large] <- exp(100)
    tmp.val[idx.small] <- exp(-100)
    pred.out <- (case + control * tmp.val) / (1 + tmp.val)
    return(pred.out)
  } else if (mode == "cont") {
    pred.out <- rep(0, n)
    pts <- seq(mean - 3 * sd, mean + 3 * sd, length.out = 100)
    num.pts <- length(pts)
    pts.aug <- matrix(rep(pts, n), n, num.pts, byrow = TRUE)
    probs <- matrix(0, n, num.pts)
    for (i in 1:num.pts) {
      X.tmp <- X
      X.tmp[,p] <- pts[i]
      probs[,i] <- loglikelihood(X.tmp, Y, Y.c, ziMatrix, para)
    }
    tmp.val <- exp(probs - probs[,1])
    idx.large <- (probs - probs[,1]) > 100
    idx.small <- (probs - probs[,1]) < -100
    tmp.val[idx.large] <- exp(100)
    tmp.val[idx.small] <- exp(-100)
    probs <- tmp.val * 
      matrix(rep(dnorm(pts, mean = mean, sd = sd), n), n, num.pts, byrow = TRUE)
    pred.out <- rowSums(pts.aug * probs) / rowSums(probs)
    return(pred.out)
  } else if (mode == "lasso") {
    pred.out <- predict(fit.lasso, newx=Y, s=lambda)
    return(as.vector(pred.out))
  } else if (mode == "glasso") {
    pred.out <- gglasso::predict.gglasso(fit.glasso, newx=Y, s=lambda, type = "link")
    return(as.vector(pred.out))
  } else {
    print("Unsupported mode!")
  }
}