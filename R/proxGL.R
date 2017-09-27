proxGL <- function(theta, s, lambda, grpSize) {
  output <- rep(0, length(theta))
  l2 <- sqrt(sum(theta^2))
  threshold <- s * lambda * sqrt(grpSize)
  if (l2 >= threshold) {
    output <- (1 - threshold / l2) * theta
  }
  return(output)
}