plot_hazards <- function(df, int_len, y_lab) {

  suppressMessages(
    p <- df %>%
      dplyr::select(
        int = .data$Interval,
        mean = .data$Mean,
        low = .data$`Prob. Low 95%`,
        high = .data$`Prob. High 95%`
      ) %>%
      dplyr::mutate(
        int_start = seq(from = 0,  by = int_len, length.out = NROW(df)),
        int_end = seq(from = 0 + int_len, by = int_len, length.out = NROW(df))
        ) %>%
      ggplot2::ggplot(ggplot2::aes(x = .data$int_start)) +
      ggplot2::geom_segment(
        ggplot2::aes(y = .data$mean, xend = .data$int_end, yend = .data$mean),
        color = "steelblue3",
        size = 0.7
        ) +
      ggplot2::geom_segment(
        ggplot2::aes(y = .data$low, xend = .data$int_end, yend = .data$low),
        lty = 2,
        color = "lightcoral",
        size = 0.5
        ) +
      ggplot2::geom_segment(
        ggplot2::aes(y = .data$high, xend = .data$int_end, yend = .data$high),
        lty = 2,
        color = "lightcoral",
        size = 0.5
        ) +
      ggplot2::scale_x_continuous(
        breaks = seq(from = 0,  by = 4*int_len, length.out = NROW(df)/4)
      ) +
      ggplot2::labs(x = "t", y = y_lab) +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid = ggplot2::element_blank())

  )

  return(p)

}

plot_survival <- function(df, int_len, y_lab) {

  suppressMessages(
    p <- df %>%
      dplyr::select(
        .data$t,
        mean = .data$`S(t)`,
        low = .data$`Prob. Low 95%`,
        high = .data$`Prob. High 95%`
      ) %>%
      dplyr::mutate(
        t = as.double(.data$t),
        t_end = c(.data$t[2:length(.data$t)], .data$t[[length(.data$t)]]),
        y_start = .data$mean,
        y_end = c(.data$mean[1:(length(.data$mean) - 1)], 0)
        ) %>%
      dplyr::rename(t_start = .data$t) %>%
      ggplot2::ggplot(ggplot2::aes(x = .data$t_start)) +
      ggplot2::geom_segment(
        ggplot2::aes(y = .data$y_start, xend = .data$t_end, yend = .data$y_end),
        color = "steelblue1",
        size = 0.7
        ) +
      ggplot2::geom_segment(
        ggplot2::aes(y = .data$low, xend = .data$t_end, yend = .data$low),
        lty = 2,
        color = "lightcoral",
        size = 0.5
        ) +
      ggplot2::geom_segment(
        ggplot2::aes(y = .data$high, xend = .data$t_end, yend = .data$high),
        lty = 2,
        color = "lightcoral",
        size = 0.5
        ) +
      ggplot2::scale_x_continuous(
        breaks = seq(from = 0,  by = 4 * int_len, length.out = NROW(df)/4 + 1)
      ) +
      ggplot2::labs(x = "t", y = y_lab) +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid = ggplot2::element_blank())
    )

  return(p)

}
