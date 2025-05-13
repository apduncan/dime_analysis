#' Plot styling constants or functions
library(ggplot2)

BIOACTIVE_COLORS <- c(
  "Low Bioactive" = "#619CFF", 
  "High Bioactive" = "#00BA38"
)
BIOACTIVE_COLORS3 <- setNames(
  c("#00BA38", "#28a9e0", "grey"),
  c("High Bioactive", "Low Bioactive", "Baseline")
)
THEME_DIME <- theme(
  panel.background = element_rect(fill = "white"),
  panel.grid = element_line(colour = "lightgrey")
)
BIO_TEXT_COLOR <- setNames(
  c("black", "#007824", "#20789e"),
  names(BIOACTIVE_COLORS3)[order(BIOACTIVE_COLORS3)]
)

FIGURE_DIR <- c("output", "figures")

write_figure <- function(
  plt,
  pth,
  ...
) {
  #' Write a figure in multiple formats
  #'
  #' Writes in png, svg and pdf. Path can includes subdirectories.
  #' Return final path. Returns the path just to the PDF, for ease of using
  #' in targets.
  #'
  #' @param plt ggplot2 figure
  #' @param pth Path to save figure, without extention
  #' @param ... Passed to ggsave
  #' @returns Path to PDF version of figure

  for (xtn in c(".png", ".svg", ".pdf")) {
    full_pth <- do.call(file.path, as.list(c(FIGURE_DIR, glue("{pth}{xtn}"))))
    print(full_pth)
    ggsave(full_pth, plt, create.dir = TRUE, dpi = 300, ...)
  }
  return(full_pth)
}