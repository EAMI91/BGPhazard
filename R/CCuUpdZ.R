CCuUpdZ <- function(times, tao, lambda, W, delta.r, y, k_i){
  K <- length(lambda)
  z <- purrr::map_int(
    purrr::map2(
      purrr::map2(.x = purrr::map(W, ~-cumsum(.x*lambda)), 
                  .y = purrr::map(seq_len(length(times)), ~ -exp(sum(delta.r * y[.x,])) + (seq_len(K) - 1) * sum(delta.r * y[.x,]) - lgamma(seq_len(K))),
                  .f = ~exp(.x + .y)
      ),
      .y = k_i, ~(seq_len(K) >= .y) *.x
    ),
    ~sample(x = seq_len(K), size = 1, prob = .x))
  return(z)
}
