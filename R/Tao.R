Tao <-
function(times, delta, type.t, K, length) {
  t.unc <- sort(times[delta  ==  1])
  if (type.t == 1) {
    n <- length(t.unc)
    if (n > K) {
      quant <- c(t.unc, max(times))
      tao <- c(0, quantile(x = times, probs = (1:K) / K, names = FALSE))
      if (type.t == 1 && length(unique(tao)) != length(tao)) {
        warning("Too many repeated observations. Zero-length intervals may
                appear.")
      }
    } else {
      stop (paste("The partition length (", K,") must be smaller than the number 
                  of uncensored times (", n, ").", sep = ""))
    }
  }
  if (type.t == 2) {
    K.t2 <- ceiling(max(times))
      if(K.t2 %% length != 0){
        tao <- seq(0, K.t2 + length, by = length)
      } else{
        tao <- seq(0, K.t2 , by = length)
      }
  }
  if (type.t  ==  3) {
    K.t2 <- ceiling(max(times))
    tao <- seq(0, K.t2, K.t2 / K)
  }
  return(tao)
}
