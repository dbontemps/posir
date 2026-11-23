# Table 1
library(posir)

alphas <- c(0.5, 0.2, 0.1, 0.05, 0.01, 0.005, 0.001)
deltas <- 10:1/10
dimension <- 1

quantiles <- sapply(deltas, FUN = function(delta) qposir(1 - alphas, delta,
                                                         d = dimension))
quantiles <- round(quantiles, 3)
# dimnames(quantiles) <- list(alpha = alphas, delta = deltas)
colnames(quantiles) <- deltas
rownames(quantiles) <- alphas

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

q <- as.matrix(quantiles)
cols <- "Reds"
# Exact Gaussian two-sided quantiles, used as a reference
alpha_num <- as.numeric(rownames(q))
q_gauss   <- qnorm(1 - alpha_num / 2)

# Each row divided by its exact Gaussian quantile
ratio <- sweep(q, 1, q_gauss, "/")

dat_r <- as.data.frame(ratio) %>%
  rownames_to_column("alpha") %>%
  pivot_longer(!alpha, names_to = "delta", values_to = "ratio") %>%
  mutate(alpha = factor(alpha, levels = rownames(q)),
         delta = factor(delta, levels = colnames(q)))

# -----------------------------------------------
# heatmap of ratios to standard Gaussian quantile
# -----------------------------------------------
figname <- sprintf("heatmap_quantile-ratios_%sd", dimension)
filename <- sprintf("figures/%s.png", figname)
p_heatmap <- ggplot(dat_r, aes(delta, alpha, fill = ratio)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", ratio),
                colour = ratio > 0), size = 3.5) +
  scale_fill_distiller(palette = cols, direction = 1, breaks = 1:4) +
  scale_colour_manual(values = c("black", "white"), guide = "none") +
  scale_y_discrete(limits = rev) +
  labs(x = expression(delta), y = expression(alpha),
       fill = expression(q[delta] / q[Gauss])) +
  theme_minimal() +
  guides(fill = guide_colourbar(title.position = "left", title.vjust = 0.8)) +
  # theme(legend.position = "bottom",
  #       legend.box.spacing = unit(0, "pt"),      # b/w axis label and legend
  #       legend.margin = margin(0, 0, 0, 0),      # internal legend margin
  #       axis.title.x = element_text(margin = margin(t = 10, b = 0)),
  #       axis.text.x  = element_text(margin = margin(t = 4)))
  theme(legend.position = "none",
        legend.box.spacing = unit(0, "pt"),      # b/w axis label and legend
        legend.margin = margin(0, 0, 0, 0),      # internal legend margin
        axis.title.x = element_text(margin = margin(t = 10, b = 0)),
        axis.text.x  = element_text(margin = margin(t = 4)))
ggsave(filename, width = 5, height = 4)

## raw quantiles (no ratio)
dat_q <- as.data.frame(q) %>%
  rownames_to_column("alpha") %>%
  pivot_longer(!alpha, names_to = "delta", values_to = "quantile") %>%
  mutate(alpha     = factor(alpha, levels = rownames(q)),
         delta_num = as.numeric(delta))

# Exact Gaussian quantiles at delta = 1 (reference points)
dat_g <- tibble(alpha = factor(rownames(q), levels = rownames(q)),
                delta_num = 1,
                quantile  = q_gauss)

# -------------------
# curves of quantiles
# -------------------
figname <- sprintf("curves_quantiles_%sd", dimension)
filename <- sprintf("figures/%s.png", figname)
ylab <- sprintf("Empirical quantiles")
ggplot(dat_q, aes(delta_num, quantile, colour = alpha, group = alpha)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  geom_point(data = dat_g, shape = 4, size = 2.5, stroke = 1,
             show.legend = FALSE) +
  scale_colour_manual(values = RColorBrewer::brewer.pal(9, "Greens")[3:9]) +
  guides(colour = guide_legend(ncol = 1, reverse = TRUE)) +
  labs(x = expression(delta),
       y = ylab, colour = expression(alpha)) +
  theme_minimal() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.07, 0.35),   # (x, y) in [0, 1], left side
        legend.justification = c(0, 1),          # anchor the legend's top-left corner
        legend.margin = margin(0, 0, 0, 0),
        legend.key.size = unit(0.4, "cm"),
        axis.title.x = element_text(margin = margin(t = 10, b = 0)),
        axis.text.x  = element_text(margin = margin(t = 4)))
ggsave(filename, width = 5, height = 4)
