# Segmentation: Comparision between POSIR and basic data splitting data
# splitting:
# - observations with even index are used to obtain segmentation
# - observations with odd index are use to build confidence intervals
# see POSIR paper for additional details.


seq_k <- c(5, 10, 20, 50)

config_list <- lapply(seq_k, function(k) {
  cbind(k, c = seq(from = k*.005, to = 1, length.out = 101))
})
configs <- as.data.frame(Reduce(rbind, config_list))
res <- dplyr::mutate(configs, ratio = CI_width_ratio(c,k))
res$L <- as.factor(res$k)

library("ggplot2")
p <- ggplot(res, aes(x = c, y = ratio, colour = L, group = L)) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkgray") +
  theme(legend.position = c(0.9, 0.3)) +
  labs(title = expression(paste("Quantile ratios for ", delta, " = c/L"))) +
  guides(color = guide_legend(reverse = TRUE))
p

ggsave(p, file = "posir-vs-split_segmentation.pdf",
       width = 6, height = 4)
