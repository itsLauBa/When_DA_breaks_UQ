library(tidyverse)
library(latex2exp)

# initializtions as in Section 4.3
n <- 10
theta <- 2
lambda <- seq(0, 12, length.out = 50)
lambda_highlight <- 5
alpha <- 0.05
sigma <- 1

b_v <- lambda / (lambda + n) * theta
s2_v <- n / (n + lambda)^2
sigma_post <- 1 / (n + lambda)
z <- qnorm(1 - alpha / 2)


# calculate upper and lower interval bounds as in equation (3)
upper <- (z * sqrt(sigma_post) - b_v) / (sqrt(sigma * s2_v))
lower <- (-z * sqrt(sigma_post) - b_v) / (sqrt(sigma * s2_v))

coverage <- pnorm(upper) - pnorm(lower)
nominal_cov <- 1 - alpha


# calculate coverage for lambda = 5
b_v_highlight <- lambda_highlight / (lambda_highlight + n) * theta
s2_v_highlight <- n / (n + lambda_highlight)^2
sigma_post_highlight <- 1 / (n + lambda_highlight)
coverage_highlight <- pnorm(
  (z * sqrt(sigma_post_highlight) - b_v_highlight) /
    (sqrt(sigma * s2_v_highlight))
) -
  pnorm(
    (-z * sqrt(sigma_post_highlight) - b_v_highlight) /
      (sqrt(sigma * s2_v_highlight))
  )


# create plot
plot_data <- tibble(lambda, coverage)

curve_color <- "#0072B2"
highlight_color <- "#D55E00"
reference_color <- "#4D4D4D"

p <- ggplot(plot_data, aes(lambda, coverage)) +
  geom_ribbon(
    aes(ymin = pmin(coverage, nominal_cov), ymax = nominal_cov),
    fill = highlight_color,
    alpha = 0.08
  ) +
  geom_hline(
    yintercept = nominal_cov,
    color = reference_color,
    linewidth = 0.65,
    linetype = "22"
  ) +
  geom_line(color = curve_color, linewidth = 1.25, lineend = "round") +
  geom_segment(
    x = lambda_highlight,
    xend = lambda_highlight,
    y = 0,
    yend = coverage_highlight,
    color = highlight_color,
    linewidth = 0.7,
    linetype = "22"
  ) +
  geom_point(
    x = lambda_highlight,
    y = coverage_highlight,
    shape = 21,
    size = 3.5,
    stroke = 1,
    color = highlight_color,
    fill = "white"
  ) +
  annotate(
    "label",
    x = lambda_highlight + 0.55,
    y = coverage_highlight + 0.08,
    hjust = 0,
    label = sprintf(
      "atop(lambda == %g, Coverage == \"%.1f%%\")",
      lambda_highlight,
      100 * coverage_highlight
    ),
    parse = TRUE,
    color = highlight_color,
    fill = "white",
    label.size = 0,
    label.padding = unit(0.18, "lines"),
    size = 5,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 11.8,
    y = nominal_cov + 0.025,
    hjust = 1,
    label = "Nominal 95% coverage",
    color = reference_color,
    size = 5
  ) +
  scale_x_continuous(
    name = expression(lambda),
    breaks = seq(0, 12, 2),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    name = "Frequentist coverage",
    labels = scales::label_percent(accuracy = 1),
    breaks = seq(0, 1, 0.25),
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.01))
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 12, base_family = "sans") +
  theme(
    plot.background = element_rect(colour = "#FFFFFF"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#E6E6E6", linewidth = 0.45),
    axis.line.x = element_line(color = "#333333", linewidth = 0.55),
    axis.ticks.x = element_line(color = "#333333", linewidth = 0.55),
    axis.ticks.length = unit(3, "pt"),
    axis.title.x = element_text(margin = margin(t = 8), size = 20),
    axis.title.y = element_text(margin = margin(r = 8), size = 20),
    axis.text = element_text(color = "#333333", size = 20),
    plot.margin = margin(8, 14, 8, 8)
  )

dir.create("plots", recursive = TRUE, showWarnings = FALSE)

ggsave("plots/coverage_example.png", p)
