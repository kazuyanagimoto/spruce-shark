tar_figure <- tar_plan(
  fn_figure = list(
    theme_proj = theme_proj,
    color_accent = color_accent,
    color_accent2 = color_accent2,
    color_accent3 = color_accent3
  )
)

theme_proj <- function(
  font_title = "Roboto Condensed",
  font_text = "Roboto Condensed Light",
  size_base = 11
) {
  ggplot2::theme_classic(base_family = font_text, base_size = size_base) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = size_base,
        face = "bold",
        family = font_title
      ),
      plot.subtitle = ggplot2::element_text(
        size = size_base,
        face = "plain",
        family = font_text
      ),
      plot.caption = ggplot2::element_text(
        size = size_base * 0.6,
        color = "grey50",
        face = "plain",
        family = font_text,
        margin = ggplot2::margin(t = 10)
      ),
      legend.position = "none",
      legend.title = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(
        size = size_base * 0.75,
        family = font_text,
        face = "plain"
      ),
      axis.line = ggplot2::element_line(linewidth = 0.3),
      axis.ticks = ggplot2::element_line(linewidth = 0.3),
      axis.title = ggplot2::element_text(
        family = font_text,
        face = "plain",
        size = size_base * 0.8
      ),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 5)),
      axis.text = ggplot2::element_text(family = font_text, face = "plain"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        size = size_base,
        hjust = 0,
        family = font_title,
        face = "bold"
      )
    )
}

color_accent <- "#107895"
color_accent2 <- "#9a2515"
color_accent3 <- "#e64173"
