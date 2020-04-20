CGaLambdaSumm<-
  function(M, new=NULL, confidence = 0.95) {
    if (confidence <= 0 || confidence >= 1) {
      stop ("Invalid parameter: confidence must be between 0 and 1.")
    }
    if(!is.null(new)){
      if(ncol(new) != ncol(M %>% extract("data"))){
        stop("Covariables doesn't match.")
      }
      if(!(match(names(new),M %>% extract("data") %>% names) %>% is.na %>% sum() == 0)){
        stop(paste("Invalid colnames, should be", paste(M %>% extract("data") %>% names, collapse = ", ")))
      } 
    }
    v <- list("K",
              "iterations",
              "tao",
              "data",
              "s",
              "S",
              "S.m",
              c("simulations","Lambda"),
              c("simulations","Lambda.m"),
              c("simulations","Theta")
    ) %>% purrr::map(~extract(M,.x)) %>% rlang::set_names(c("K","iterations","tao","data",
                                                                                "s","S","S.m","Lambda","Lambda.m",
                                                                                "Theta"))
    tao <- v$tao
    Lambda.b <- v$Lambda
    Lambda.m <- v$Lambda.m
    Theta <- v$Theta
    Z <- v$Z.m %>% tibble::enframe() %>% dplyr::select(mean.obs = value) 
    Pi <- v$Pi.m %>% tibble::enframe() %>% dplyr::select(mean.obs  = value) 
    S.b <- v$S
    S.m <- v$S.m
    Lambda.obs <- NULL
    S.obs <- NULL
    if(!is.null(new)) {
      names <- names(v$data)
      new <- tibble::as_tibble(new)
      new <- dplyr::select(new, !!! colnames(data[,c(-1,-2)]))
      Lambda.obs <- v$Lambda
      
      Lambda.obs <- purrr::map2(.x = Lambda.obs, .y =seq_len(nrow(new)),function(a,b){
        eff <- as.numeric(exp(Theta%*%as.numeric(new[b,])))
        a <- dplyr::mutate_all(a,.f = ~.x*eff)
        return(a)
      })
      
      writeLines("Generating survival function estimates for new observations.")
      pb <- dplyr::progress_estimated(length(v$s))
      
      S.obs <-  do.call(dplyr::bind_cols, purrr::map(v$s, function(s = .x){
        pb$tick()$print()
        tibble::as_tibble(matrix(data = purrr::map_dbl(purrr::map(purrr::cross2(seq_len(v$iterations),purrr::map2(Lambda.obs, tibble::as_tibble(t(new)),~list(.x,.y))), 
                                                                  .f= ~(s > tao[-1]) * tao[-1] * as.numeric(.x[[2]][[1]][.x[[1]],]) + 
                                                                    (s > tao[-length(tao)] & s <= tao[-1]) * s * as.numeric(.x[[2]][[1]][.x[[1]],])  -
                                                                    (s > tao[-length(tao)]) * tao[-(length(tao))] * as.numeric(.x[[2]][[1]][.x[[1]],]) 
        ), ~exp(-sum(.x))),
        ncol = nrow(new), byrow = F))
      }))
      
      
      S.obs <- purrr::map(seq_len(nrow(new)),
                          ~S.obs[,seq(.x, (nrow(new))*length(v$s), nrow(new))])
      cat("\n Done.")
      
    } 
    
    pr <- (1 - confidence) / 2
    Lambda <- c(Lambda.b,Lambda.m,Lambda.obs)
    SUM.h <- purrr::map(Lambda, ~rlang::set_names(tibble::tibble(a=seq_len(v$K),
                                                         b=purrr::map_dbl(.x, mean,na.rm = T),
                                                         c=purrr::map_dbl(.x, quantile, probs = pr, na.rm = T),
                                                         d=purrr::map_dbl(.x, quantile, probs = 0.5, na.rm = T),
                                                         e=purrr::map_dbl(.x, quantile, probs = 1 - pr, na.rm = T)
    ), 
    c("k", "mean",  "lower", "median", "upper")))
    
    S <- c(S.b,S.m,S.obs)
    SUM.S <- purrr::map(S, ~ rlang::set_names(tibble::tibble(a=v$s,
                                                             b=purrr::map_dbl(.x, mean, na.rm = T),
                                                             c=purrr::map_dbl(.x, quantile, probs = pr, na.rm = T), 
                                                             d=purrr::map_dbl(.x, quantile, probs = 0.5, na.rm = T), 
                                                             e=purrr::map_dbl(.x, quantile, probs = 1-pr, na.rm = T)),
                                              c("t", "S^(t)",  "lower", "median", "upper")))
    
    out <- tibble::enframe(list(SUM.h = SUM.h, SUM.S = SUM.S, K = nrow(new)))
  }
