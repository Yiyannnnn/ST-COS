#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(lubridate)
  library(scales)
  library(ggpubr)
  library(cowplot)
})

out_dir <- "./Reference_based"

df <- data.frame(
  Methods = c(rep("ST-COS", 8), rep("SRTsim", 8), rep("scDesign3", 8)),
  N_genes = c(50, 100, 500, 1000, 3000, 5000, 8000, 10000,
              50, 100, 500, 1000, 3000, 5000, 8000, 10000,
              50, 100, 500, 1000, 3000, 5000, 8000, 10000),
  Time_str = c(
    "0:00:04", "0:01:25", "0:08:47", "0:24:58", "0:38:45", "0:19:02", "0:09:44", "0:05:59",
    "0:00:02", "0:00:34", "0:01:21", "0:04:46", "0:18:06", "0:11:43", "0:07:46", "0:04:26",
    "0:02:54", "0:05:31", "1:00:17", "4:36:34", "17:44:07", "17:59:14", "18:13:43", "19:46:35"
  )
)

df$Time_sec <- period_to_seconds(hms(df$Time_str))
df$Methods  <- factor(df$Methods, levels = c("ST-COS", "SRTsim", "scDesign3"))

method_colors <- c("ST-COS" = "#1F77B4", "SRTsim" = "#D62728", "scDesign3" = "#2CA02C")
method_shapes <- c("ST-COS" = 16, "SRTsim" = 17, "scDesign3" = 15)

base_aes <- list(
  geom_line(linewidth = 0.9),
  geom_point(size = 3),
  scale_color_manual(values = method_colors),
  scale_shape_manual(values = method_shapes),
  scale_x_continuous(
    breaks = c(0, 2000, 4000, 6000, 8000, 10000),
    labels = label_comma()
  ),
  labs(
    x = "Number of Genes",
    color = NULL, shape = NULL
  ),
  theme_pubr(base_size = 12),
  theme(
    legend.position  = "top",
    legend.text      = element_text(size = 11),
    axis.title       = element_text(face = "bold"),
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 9, color = "gray40")
  )
)

## ---- Panel A: regular scale ------------------------------------------
p_linear <- ggplot(df, aes(x = N_genes, y = Time_sec / 3600,
                           color = Methods, shape = Methods)) +
  base_aes +
  labs(y = "Runtime (hours)",
       title = "Linear Scale") +
  scale_y_continuous(labels = function(x) sprintf("%.1fh", x))

## ---- Panel B: log-log scale ------------------------------------------
time_fmt <- function(x) {
  h <- x %/% 3600; m <- (x %% 3600) %/% 60
  ifelse(h > 0, paste0(h, "h", m, "m"), paste0(m, "m"))
}

p_loglog <- ggplot(df, aes(x = N_genes, y = Time_sec,
                           color = Methods, shape = Methods)) +
  base_aes +
  scale_y_log10(
    breaks = c(5, 30, 120, 600, 1800, 7200, 36000, 72000),
    labels = time_fmt
  ) +
  scale_x_log10(
    breaks = c(50, 100, 500, 1000, 3000, 5000, 10000),
    labels = label_comma()
  ) +
  annotation_logticks(sides = "bl", size = 0.3) +
  labs(y = "Runtime (log scale)",
       title = "Log-Log Scale")

## ---- Combine ----------------------------------------------------------
combined <- plot_grid(
  p_linear, p_loglog,
  ncol   = 2,
  labels = c("A", "B"),
  label_size = 13,
  align  = "h", axis = "tb"
)

ggsave(file.path(out_dir, "benchmark_runtime_combined.pdf"),
       combined, width = 12, height = 5, bg = "white")
message("Saved: benchmark_runtime_combined.pdf")

## individual files
ggsave(file.path(out_dir, "benchmark_runtime_linear.pdf"),
       p_linear, width = 6, height = 5, bg = "white")
ggsave(file.path(out_dir, "benchmark_runtime_loglog.pdf"),
       p_loglog, width = 6, height = 5, bg = "white")
message("Saved individual PDFs.")
