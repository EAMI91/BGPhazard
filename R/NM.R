NM <-
  function(times, delta, type.t, K, utao) {
    tao <- Tao(times, delta, type.t, K, utao)
    t.unc <- sort(times[delta == 1])
    n <- readr::parse_integer(as.character(table(cut(t.unc,breaks = tao))))
    w <- purrr::map(times, .f= ~ (.x > tao[-1]) * tao[-1] + 
               (.x > tao[-length(tao)] & .x <= tao[-1]) * .x -
               (.x > tao[-length(tao)]) * tao[-(length(tao))])
    
    m <- purrr::reduce(w,`+`)
    out <- list(n = n, m = m, tao = tao, t.unc = t.unc)
    return(out)
  }
