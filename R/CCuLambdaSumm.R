CCuLambdaSumm<-
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
              "covs.x",
              "covs.y",
              "s",
              "S",
              c("simulations","Lambda"),
              c("simulations","Lambda.m"),
              c("simulations","Delta"),
              c("simulations","Theta"),
              c("simulations","Pi.m"),
              c("simulations","Z.m")
    ) %>% purrr::map(~extract(M,.x)) %>% rlang::set_names(c("K","iterations","tao","data","covs.x","covs.y",
                                                                              "s","S","Lambda","Lambda.m","Delta",
                                                                              "Theta","Pi.m","Z.m"))
    tao <- v$tao
    Lambda <- v$Lambda.m
    Theta <- v$Theta
    Z <- v$Z.m %>% tibble::enframe() %>% dplyr::select(mean.obs = value) 
    Pi <- v$Pi.m %>% tibble::enframe() %>% dplyr::select(mean.obs  = value) 
    S <- v$S
    if(!is.null(new)) {
      names <- names(v$data)
      new <- tibble::as_tibble(new)
      Lambda <- v$Lambda
      new.x <- dplyr::select(new, !!! colnames(v$covs.x))
      new.y <- dplyr::select(new, !!! colnames(v$covs.y))
      
      Z <- matrix(data = (purrr::map_int(purrr::map(purrr::map(purrr::cross2(tibble::as_tibble(t(new.y)),
                                                                      tibble::as_tibble(t(v$Delta))),
                                                               purrr::lift(`*`)),
                                                    ~exp(sum(.x))),
                                         ~rpois(n = 1, lambda = .x)) + 1 ),
                  nrow = v$iterations, ncol = nrow(new.y), byrow = T)
      Z[Z>v$K] <- v$K
      Z <- rlang::set_names(tibble::as_tibble(Z),
                            paste0("new_obs_", seq_len(nrow(new.y)))) 
      
      
      Lambda <-  purrr::map_dfc(purrr::cross2(seq_len(v$K), Z),~ dplyr::pull((.x[[1]] <= .x[[2]])*Lambda[,.x[[1]]],1))
      Lambda <-  purrr::map(seq_len(ncol(Lambda)/v$K), ~tibble::as_tibble(Lambda[,seq_len(v$K) + (.x-1)*v$K]))
      
      
      Pi <- rlang::set_names(tibble::as_tibble(matrix(purrr::map2_dbl(purrr::map_dfc(purrr::cross2(seq_len(v$iterations),tibble::as_tibble(t(new.x))),
                                                                                     .f = ~ exp(sum(Theta[.x[[1]],] * .x[[2]])) * (tao[-1] - tao[-length(tao)])),
                                                                      .y =  tibble::as_tibble(t(purrr::reduce(purrr::map(Lambda, ~rlang::set_names(.x,paste0("V",seq_len(v$K)))), dplyr::bind_rows))), 
                                                                      .f = ~exp(-sum(.x*.y))),
                                                      nrow = v$iterations, ncol = nrow(new.x))),
                             paste0("new_obs_", seq_len(nrow(new.x))))
      writeLines("Generating survival function estimates for new observations.")
      pb <- dplyr::progress_estimated(length(v$s))
      
      Lambda <- purrr::map2(.x = Lambda, .y =seq_len(nrow(new)),function(a,b){
        eff <- as.numeric(exp(Theta%*%as.numeric(new[b,])))
        a <- dplyr::mutate_all(a,.f = ~.x*eff)
        return(a)
      })
      
      S <-  do.call(dplyr::bind_cols, purrr::map(v$s, function(s = .x){
        pb$tick()$print()
        tibble::as_tibble(matrix(data = purrr::map_dbl(purrr::map(purrr::cross2(seq_len(v$iterations),purrr::map2(Lambda, tibble::as_tibble(t(new.x)),~list(.x,.y))), 
                                                                  .f= ~(s > tao[-1]) * tao[-1] * as.numeric(.x[[2]][[1]][.x[[1]],])+ 
                                                                    (s > tao[-length(tao)] & s <= tao[-1]) * s * as.numeric(.x[[2]][[1]][.x[[1]],]) -
                                                                    (s > tao[-length(tao)]) * tao[-(length(tao))] * as.numeric(.x[[2]][[1]][.x[[1]],]) 
        ), ~exp(-sum(.x))),
        ncol = nrow(new.x), byrow = F))
      }))
      
      
      S <- purrr::map(seq_len(nrow(new.x)),
                      ~S[,seq(.x, (nrow(new.x))*length(v$s), nrow(new.x))])
      cat("\n Done.")
      
    } 
    
    pr <- (1 - confidence) / 2
    
    SUM.h <- purrr::map(Lambda, ~rlang::set_names(tibble::tibble(seq_len(v$K),
                                                                 purrr::map_dbl(.x, mean,na.rm = T),
                                                                 purrr::map_dbl(.x, quantile, probs = pr, na.rm = T),
                                                                 purrr::map_dbl(.x, quantile, probs = 0.5, na.rm = T),
                                                                 purrr::map_dbl(.x, quantile, probs = 1 - pr, na.rm = T)
    ), 
    c("k", "mean",  "lower", "median", "upper")))
    
    
    SUM.S <- purrr::map(S, ~ rlang::set_names(tibble::tibble(v$s,
                                                      purrr::map_dbl(.x, mean, na.rm = T),
                                                      purrr::map_dbl(.x, quantile, probs = pr, na.rm = T), 
                                                      purrr::map_dbl(.x, quantile, probs = 0.5, na.rm = T), 
                                                      purrr::map_dbl(.x, quantile, probs = 1-pr, na.rm = T)),
                                       c("t", "S^(t)",  "lower", "median", "upper")))
    
    SUM.Pi <-  purrr::map(Pi, ~rlang::set_names(tibble::tibble(mean(.x, na.rm=T), 
                                                               quantile(.x, probs = pr,na.rm=T),
                                                               quantile(.x, probs = 0.5,na.rm=T),
                                                               quantile(.x, probs = 1 - pr,na.rm=T)),
                                                c("mean",  "lower", "median", "upper")))
    
    z <- purrr::map(Z, ~rlang::set_names(tibble::tibble(mean(.x,na.rm=T), 
                                                        quantile(.x, probs = pr,na.rm=T),
                                                        quantile(.x, probs = 0.5,na.rm=T),
                                                        quantile(.x, probs = 1 - pr,na.rm=T)),
                                         c("mean",  "lower", "median", "upper")))
    
    out <- tibble::enframe(list(SUM.h = SUM.h, SUM.S = SUM.S, SUM.z = z, SUM.pi = SUM.Pi, K = nrow(new),
                                simulations = tibble::enframe(list(Z = Z, Pi = Pi))))
  }
