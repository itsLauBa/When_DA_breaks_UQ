library(ggplot2)

output_dir <- file.path("results")
image_dir <- "plots"
dir.create(image_dir, recursive = TRUE, showWarnings = FALSE)


# Coverage plots (Figure 2 and Figure 3 in Appendix D) #########################
classical <- read.csv(file.path(output_dir, "classical_results.csv"), stringsAsFactors = FALSE)
high_dim <- read.csv(file.path(output_dir, "high_dim_results.csv"), stringsAsFactors = FALSE)

method_colors <- c(
  Credible = "#D55E00",
  Oracle = "#0072B2",
  Split = "#009E73",
  "Split-Lasso-cert" = "#009E73"
)
method_shapes <- c(Credible = 16, Oracle = 17, Split = 15, "Split-Lasso-cert" = 15)
method_labels <- c(
  Credible = "Credible",
  Oracle = "Oracle",
  Split = "Split",
  "Split-Lasso-cert" = "Certified split Lasso"
)

coverage_plot <- function(data, methods, title_text) {
  data <- data[data$method %in% methods, ]
  data$method <- factor(data$method, levels = methods)
  data$design <- factor(data$design, levels = unique(data$design))
  data$augmentation <- factor(data$augmentation, levels = c("Dropout", "Gaussian"))
  data$coverage_lower <- pmax(0, data$coverage - 1.96 * data$coverage_se)
  data$coverage_upper <- pmin(1, data$coverage + 1.96 * data$coverage_se)
  strength_levels <- sprintf("%.2f", sort(unique(data$strength)))
  data$strength_label <- factor(
    sprintf("%.2f", data$strength),
    levels = strength_levels
  )

  ggplot(
    data,
    aes(
      x = strength_label,
      y = coverage,
      color = method,
      shape = method,
      group = method
    )
  ) +
    geom_hline(
      yintercept = 0.975,
      linewidth = 0.6,
      linetype = "dashed",
      color = "grey35"
    ) +
    geom_line(linewidth = 0.6, alpha = 0.8) +
    geom_linerange(
      aes(ymin = coverage_lower, ymax = coverage_upper),
      linewidth = 0.65,
      show.legend = FALSE
    ) +
    geom_point(size = 2.4) +
    facet_grid(
      rows = vars(design),
      cols = vars(augmentation),
      scales = "free_x",
      space = "free_x",
      switch = "y"
    ) +
    scale_y_continuous(
      limits = c(0, 1.02),
      breaks = seq(0, 1, by = 0.25),
      expand = expansion(mult = c(0, 0.01))
    ) +
    scale_color_manual(
      values = method_colors[methods],
      labels = method_labels[methods],
      drop = FALSE
    ) +
    scale_shape_manual(
      values = method_shapes[methods],
      labels = method_labels[methods],
      drop = FALSE
    ) +
    labs(
      # title = title_text,
      # subtitle = "Point estimates with 95% Monte Carlo intervals; dashed horizontal line: 97.5% target",
      x = "Augmentation strength",
      y = "Empirical coverage",
      color = NULL,
      shape = NULL
    ) +
    theme_minimal(base_size = 10.5) +
    theme(
      plot.title.position = "plot",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.35),
      panel.spacing = grid::unit(0.7, "lines"),
      strip.background = element_rect(fill = "grey94", color = NA),
      strip.background.y = element_blank(),
      strip.placement = "outside",
      strip.text = element_text(face = "bold", color = "grey20", size = 15),
      strip.text.y.left = element_text(angle = 0, hjust = 1),
      axis.text = element_text(size = 15),
      legend.position = "bottom",
      axis.title.x = element_text(margin = margin(t = 8), size = 15),
      axis.title.y = element_text(margin = margin(r = 8), size = 15),
      legend.text = element_text(size = 12)
    )
}

classical_plot <- coverage_plot(
  classical,
  c("Credible", "Oracle", "Split"),
  "Classical regime"
)

ggsave(
  file.path(image_dir, "coverage_classical_regime.png"),
  classical_plot,
  width = 10.5,
  height = 8.5,
  units = "in",
  dpi = 220,
  bg = "white"
)


high_dim_plot <- coverage_plot(
  high_dim,
  c("Credible", "Oracle", "Split-Lasso-cert"),
  "High-dimensional regime"
)

ggsave(
  file.path(image_dir, "coverage_high_dimensional_regime.png"),
  high_dim_plot,
  width = 10.5,
  height = 9.5,
  units = "in",
  dpi = 220,
  bg = "white"
)

# Width comparison plot (Figure 4 in Appendix D) ###############################

zero <- high_dim[high_dim$method == "Split-Lasso-0", ]
cert <- high_dim[high_dim$method == "Split-Lasso-cert", ]
key_columns <- c("design", "augmentation", "strength")
paired <- merge(zero, cert, by = key_columns, suffixes = c("_zero", "_cert"), sort = FALSE)
paired$cell_id <- seq_len(nrow(paired))

width_data <- rbind(
  data.frame(
    cell_id = paired$cell_id,
    certificate = "Zero certificate",
    rel_width = paired$rel_width_zero
  ),
  data.frame(
    cell_id = paired$cell_id,
    certificate = "Ridge-grid certificate",
    rel_width = paired$rel_width_cert
  )
)
width_data$certificate <- factor(
  width_data$certificate,
  levels = c("Zero certificate", "Ridge-grid certificate")
)

certificate_colors <- c(
  "Zero certificate" = "#CC79A7",
  "Ridge-grid certificate" = "#009E73"
)

certificate_plot <- ggplot(
  width_data,
  aes(x = certificate, y = rel_width, group = cell_id)
) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  geom_line(color = "grey45", alpha = 0.45, linewidth = 0.5) +
  geom_point(aes(color = certificate), size = 2) +
  scale_y_log10() +
  scale_color_manual(values = certificate_colors, guide = "none") +
  labs(
    y = "Relative interval width (log scale)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title.position = "plot",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text = element_text(size = 15),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 15),
    legend.text = element_text(size = 12)
  )

ggsave(
  file.path(image_dir, "certificate_width_diagnostic.png"),
  certificate_plot,
  width = 9.5,
  height = 6,
  units = "in",
  dpi = 220,
  bg = "white"
)
