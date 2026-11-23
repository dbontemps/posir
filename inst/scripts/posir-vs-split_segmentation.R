# Segmentation: Comparision between POSIR and basic data splitting data
# splitting:
# - observations with even index are used to obtain segmentation
# - observations with odd index are use to build confidence intervals

#  If $k$ is the (potentially random) number of regions in the final
#  segmentation, the confidence interval for the average value of $\mu$ in the
#  region $I$ will have as radius $q_\hat{\sigma}/\sqrt{\ell(I)/2}$, where
#  $\ell(I)$ is the length of $I$, the "/2" is due to the fact that only half of
#  the observations are used, and $q$ is the Gaussian quantile of level
#  $1-\alpha/(2k)$ (Bonferroni’s rule). We compare with the radius of POSIR
#  confidence intervals $K_{1-\alpha,\delta} \hat{\sigma}/\sqrt{\ell(I)}$:

library("posir")
CI_width_ratio <- function(c, L, nb_seg = 1, alpha = 0.05) {
  # \frac{\sqrt{2} q_{1-\frac{\alpha}{2 L}}}{K_{1-\alpha,\delta}}.
  delta = c / L
  ratio <- sqrt(2) * qnorm(1 - alpha/(2 * L * nb_seg))/posir:::posir_quantiles(alpha, delta)
  dim(ratio) <- NULL
  ratio
}
seq_L <- c(5, 10, 20, 50)

config_list <- lapply(seq_L, function(L) {
  cbind(L, c = seq(from = L*.005, to = 1, length.out = 101))
})
configs <- as.data.frame(Reduce(rbind, config_list))
configs <- rbind(cbind(configs, nb_seg = 1),
                 cbind(configs, nb_seg = 3),
                 cbind(configs, nb_seg = 5),
                 cbind(configs, nb_seg = 10))
alpha <- 0.05

res <- dplyr::mutate(configs, ratio =
                        CI_width_ratio(c, L, nb_seg, alpha = 0.05))
res$L <- as.factor(res$L)
seg_labels <- c("1" = "1 segmentation",
                "3" = "3 segmentations",
                "5" = "5 segmentations",
                "10" = "10 segmentations")

library("ggplot2")
p <- ggplot(res, aes(x = c, y = ratio, colour = L, group = L)) +
  geom_line() +
  facet_wrap(vars(nb_seg), labeller = labeller(nb_seg = seg_labels)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkgray") +
  theme(legend.position = c(0.9, 0.3)) +
  labs(title = expression(paste("Quantile ratios for ", delta, " = c/L"))) +
  guides(color = guide_legend(reverse = TRUE)) +
  theme_minimal()
p

ggsave(p, file = "posir-vs-split_segmentation.pdf",
       width = 6, height = 4)
