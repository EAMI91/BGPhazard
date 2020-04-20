CuLambdaSumm <-
  function(M, confidence = 0.95) {
    if (confidence <= 0 || confidence >= 1) {
      stop ("Invalid parameter: confidence must be between 0 and 1.")
    }
    v <- rlang::set_names(purrr::map(list("K",
                            "iterations",
                            "s",
                            "S",
                            c("simulations","Lambda"),
                            c("simulations","Pi"),
                            c("simulations","Z")
    ),
    ~extract(M,.x)),
    c("K","iterations","s","S","Lambda","Pi","Z"))
    K <- v$K
    iterations <- v$iterations
    pr <- (1 - confidence) / 2
    S <- v$S
    
    SUM.h <- rlang::set_names(tibble::tibble(a=seq_len(K),
                              b=purrr::map_dbl(v$Lambda, mean),
                              c=purrr::map_dbl(v$Lambda, quantile, probs = pr),
                              d=purrr::map_dbl(v$Lambda, quantile, probs = 0.5),
                              e=purrr::map_dbl(v$Lambda, quantile, probs = 1 - pr)
    ),
    c("k", "mean",  "lower", "median", "upper"))
    
    
    SUM.S <- rlang::set_names(tibble::tibble(a=v$s,
                       b=purrr::map_dbl(v$S, mean),
                       c=purrr::map_dbl(v$S, quantile, probs = pr), 
                       d=purrr::map_dbl(v$S, quantile, probs = 0.5), 
                       e=purrr::map_dbl(v$S, quantile, probs = 1-pr)),
                c("t", "S^(t)",  "lower", "median", "upper"))
    
    prop.pi <- v$Pi
    prop.pi <- dplyr::rename(tibble::as_tibble(t(c(mean(prop.pi), quantile(prop.pi, c(pr, 0.5, 1 - pr))))), "mean" = "V1")
    z <- v$Z
    z <- dplyr::rename(tibble::as_tibble(t(c(mean(z), quantile(z, c(pr, 0.5, 1 - pr))))), "mean" = "V1")
    out <- tibble::enframe(list(SUM.h = SUM.h, SUM.S = SUM.S, SUM.pi = prop.pi, SUM.z = z))
  }
