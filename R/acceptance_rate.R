# Not exported
acceptance_rate <- function(.x) {
  len_x <- length(.x)
  
  s <- sum(.x[1:(len_x-1)] == .x[2:len_x])
  
  return(1 - s/(len_x - 1))
}
