CCuUpdDelta <-
  function(delta, y , z, var.delta.str, var.delta.ini, acceptance.d) {
    delta.upd <- delta
    p <- length(delta)
    for(s in seq_len(p)){
      delta.str <- delta.upd
      delta.str[s] <- rnorm(1, mean = delta.upd[s], sd = sqrt(var.delta.str))
      pr <- dnorm(delta.str[s], mean = 0, sd = sqrt(var.delta.ini), log = T) -
        dnorm(delta.upd[s], mean = 0, sd = sqrt(var.delta.ini), log = T) +
        sum((delta.str[s] - delta.upd[s]) * (y[, s]) * (z - 1) - 
              (exp(unname(purrr::map_dbl(tibble::as_tibble(t(y),.name_repair = "minimal"),~sum(.x*delta.str))))) +
              (exp(unname(purrr::map_dbl(tibble::as_tibble(t(y),.name_repair = "minimal"),~sum(.x*delta.upd)))))
        )
      if(log(runif(1)) <= pr) {
        delta.upd <- delta.str
        acceptance.d[s] <- acceptance.d[s] + 1
      }
    }
    return(list(delta.upd,acceptance.d))
  }
