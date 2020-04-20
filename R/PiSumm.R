PiSumm <-
  function(M, confidence = 0.95) {
    if (confidence <= 0 || confidence >= 1) {
      stop ("Invalid parameter: confidence must be between 0 and 1.")
    }
    v <- list("K",
              "tao",
              "iterations",
              "s",
              "S",
              c("simulations","PI")
    ) %>% purrr::map(~extract(M,.x)) %>% rlang::set_names(c("K","tao","iterations","s","S","PI"))
    K <- v$K
    iterations <- v$iterations
    pr <- (1 - confidence) / 2
    S <- v$S
    
    SUM.h <- tibble::tibble(a=seq_len(K),
                    b=v$PI %>% tibble::as_tibble() %>% purrr::map_dbl(mean),
                    c=v$PI %>% tibble::as_tibble() %>% purrr::map_dbl(quantile, probs = pr),
                    d=v$PI %>% tibble::as_tibble() %>% purrr::map_dbl(quantile, probs = 0.5),
                    e=v$PI %>% tibble::as_tibble() %>% purrr::map_dbl(quantile, probs = 1 - pr)
    ) %>% rlang::set_names(c("k", "mean",  "lower", "median", "upper"))
    
    
    SUM.S <- tibble::tibble(a=c(0,v$s),
                    b=c(1,v$S %>% purrr::map_dbl(mean)),
                    c=c(1,v$S %>% purrr::map_dbl(quantile, probs = pr)), 
                    d=c(1,v$S %>% purrr::map_dbl(quantile, probs = 0.5)), 
                    e=c(1,v$S %>% purrr::map_dbl(quantile, probs = 1-pr))) %>%
      rlang::set_names(c("t", "S^(t)",  "lower", "median", "upper"))
    
    out <- list(SUM.h = SUM.h, SUM.S = SUM.S) %>% tibble::enframe()
  }


