LambdaSumm <-
  function(M, confidence = 0.95) {
    if (confidence <= 0 || confidence >= 1) {
      stop ("Invalid parameter: confidence must be between 0 and 1.")
    }
    v <- list("K",
              "iterations",
              "s",
              "S",
              c("simulations","Lambda")
    ) %>% purrr::map(~extract(M,.x)) %>% rlang::set_names(c("K","iterations","s","S","Lambda"))
    K <- v$K
    iterations <- v$iterations
    pr <- (1 - confidence) / 2
    S <- v$S
    
    SUM.h <- tibble::tibble(a=seq_len(K),
                    b=v$Lambda %>% tibble::as_tibble() %>% purrr::map_dbl(mean),
                    c=v$Lambda %>% tibble::as_tibble() %>% purrr::map_dbl(quantile, probs = pr),
                    d=v$Lambda %>% tibble::as_tibble() %>% purrr::map_dbl(quantile, probs = 0.5),
                    e=v$Lambda %>% tibble::as_tibble() %>% purrr::map_dbl(quantile, probs = 1 - pr)
    ) %>% rlang::set_names(c("k", "mean",  "lower", "median", "upper"))
    
    
    SUM.S <- tibble::tibble(a=v$s,
                    b=v$S %>% purrr::map_dbl(mean),
                    c=v$S %>% purrr::map_dbl(quantile, probs = pr), 
                    d=v$S %>% purrr::map_dbl(quantile, probs = 0.5), 
                    e=v$S %>% purrr::map_dbl(quantile, probs = 1-pr)) %>%
      rlang::set_names(c("t", "S^(t)",  "lower", "median", "upper"))

    out <- list(SUM.h = SUM.h, SUM.S = SUM.S) %>% tibble::enframe()
  }
