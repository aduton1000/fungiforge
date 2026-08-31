# Shared theme + palette for the fungiforge comparative layer (mirrors the bacterial
# r-analysis conventions so the two studies look identical side by side).
suppressMessages({ library(ggplot2) })

# fixed compartment palette (same as the bacterial ASSARM-PHI figures)
COMPARTMENT_COLS <- c(AIR = "#2C7FB8", HUMAN = "#D95F0E", SURFACE = "#31A354")
COMP_LEVELS      <- c("AIR", "HUMAN", "SURFACE")

theme_pub <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title  = element_text(face = "bold", size = 12),
        legend.position = "bottom",
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold"))
theme_set(theme_pub)

set.seed(20260729)

# self-resolving paths: every objective script derives its own dir, then ../output + ../figures
ff_dirs <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  sd <- dirname(sub("^--file=", "", args[grep("^--file=", args)]))
  if (length(sd) == 0 || sd == "") sd <- "."
  root <- normalizePath(file.path(sd, ".."))
  list(root = root,
       output = file.path(root, "output"),
       figures = file.path(root, "figures"))
}

ff_save <- function(plot, name, w = 8, h = 6, dir = NULL) {
  if (is.null(dir)) dir <- ff_dirs()$figures
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(dir, paste0(name, ".png")), plot, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(dir, paste0(name, ".pdf")), plot, width = w, height = h, bg = "white")
}
