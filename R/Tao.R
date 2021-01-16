Tao <-
function(times, delta, type.t, K, utao) {
  t.unc <- sort(times[delta  ==  1])
  if (type.t == 1) {
    n <- length(t.unc)
    if (n > K) {
      tao <- c(0, quantile(x = t.unc, probs = (1:(K-1)) / (K-1), names = FALSE), max(times))
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
    tao <- utao
  }
  if (type.t  ==  3) {
    K.t2 <- ceiling(max(times))
    tao <- seq(0, K.t2, K.t2 / K)
  }
  return(tao)
}
