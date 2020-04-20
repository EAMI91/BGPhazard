BeNM <-
  function(times, delta) {
    K <- max(times)
    tao <- seq(0, K)
    t.unc <- sort(times[delta == 1])
    n <- cut(t.unc,tao) %>% table %>% as.character %>% readr::parse_integer()
    m <- tao[-1] %>% purrr::map_int(~length(times[times > .x]))
    out <- list(n = n, m = m, tao = tao, K = K, t.unc = t.unc)
    out
  }

