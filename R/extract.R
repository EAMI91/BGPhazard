extract <- function(tb, dir){
  for(i in seq_along(dir)){
    tb %<>% dplyr::filter(name == !!rlang::parse_expr("dir[i]")) %>% dplyr::pull(value) %>% .[[1]]  
  }
  return(tb)
}
