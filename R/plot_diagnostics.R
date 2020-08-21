plot_traceplots <- function(matrix) {

  suppressMessages(
  p <- matrix %>%
    tibble::as_tibble(.name_repair = "unique") %>%
    dplyr::mutate(individual = 1:dplyr::n()) %>%
    tidyr::pivot_longer(cols = -.data$individual) %>%
    dplyr::mutate(name = as.double(stringr::str_remove_all(.data$name, "\\."))) %>%
    ggplot2::ggplot() +
    ggplot2::aes(x = .data$name, y = .data$value) +
    ggplot2::geom_line(color = "steelblue1") +
    ggplot2::facet_wrap(~.data$individual, scales = "free_y") +
    ggplot2::labs(x = "Iteration", y = "") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  )

  return(p)

}

plot_ergodic_means <- function(matrix) {

  suppressMessages(
  p <- matrix %>%
    tibble::as_tibble(.name_repair = "unique") %>%
    dplyr::mutate(individual = 1:dplyr::n()) %>%
    tidyr::pivot_longer(cols = -.data$individual) %>%
    dplyr::mutate(
      name = as.double(stringr::str_remove_all(.data$name, "\\."))
      ) %>%
    dplyr::group_by(.data$individual) %>%
    dplyr::mutate(value = cumsum(.data$value)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(value = .data$value / .data$name) %>%
    ggplot2::ggplot() +
    ggplot2::aes(x = .data$name, y = .data$value) +
    ggplot2::geom_line(color = "steelblue1") +
    ggplot2::facet_wrap(~.data$individual, scales = "free_y") +
    ggplot2::labs(x = "Iteration", y = "") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  )

  return(p)

}
