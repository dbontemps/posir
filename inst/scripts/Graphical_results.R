
library(rlist)
library("ggplot2")
ggplot2::theme_set(theme_minimal())
library(xtable)

print_quantiles = function(filename) {
  alphagrid <- c(0.5, 0.2, 0.1, 0.05, 0.01, 0.005, 0.001)
  deltagrid <- seq(10, 1, -1)/10
  myrownames = sapply(alphagrid, toString)
  mycolnames = sapply(deltagrid, function(x){paste("X",toString(x),sep="")})
  X = read.table(filename)
  X = X[myrownames, mycolnames]
  colnames(X) = deltagrid
  print(xtable(X,digits=3))
}

print_quantiles("Table_quantiles_1D.txt")
print_quantiles("Table_quantiles_2D.txt")

read_levels = function(prename, sufname, criterion, alphagrid, deltagrid) {
  myY = list()
  for(k in criterion) {
    X = read.table(paste(prename, toString(k), sufname, sep=""))
    myrownames = sapply(alphagrid, toString)
    mycolnames = sapply(deltagrid, function(x){paste("X",toString(x),sep="")})
    X = X[myrownames, mycolnames]
    colnames(X) = deltagrid
    X$alpha <- as.factor(alphagrid)
    df <- tidyr::pivot_longer(X, cols = !contains("alpha"),
                              names_to = "delta", values_to = "level")
    df$delta <- as.numeric(df$delta)
    df$crit <- k
    myY=list.append(myY, df)
  }
  return(myY)
}

plot_levels = function(Y, alphagrid, critname, title, filename) {
  level <- Reduce(rbind, Y)
  level$crit <- as.factor(level$crit)
  targets <- data.frame(y = alphagrid, alpha = factor(alphagrid))

  p <- ggplot(level, aes(x = delta, y = level, group = crit, color = crit)) +
    geom_line() +
    geom_hline(aes(yintercept = y), targets, linetype = "dotted") +
    facet_wrap(vars(alpha), labeller = label_both, scales = "free") +
    theme(legend.position = c(0.8, 0.2)) +
    scale_colour_brewer(palette = "Purples") +
    scale_y_log10() +
    labs(title = title,
         x = "delta",
         y = "Effective error level",
         color = critname)
  p
  ggsave(p, filename = filename, height = 4, width = 6) #, height = 4, width = 6
}

alphagrid <- c(0.5, 0.1, 0.05, 0.01, 0.001)
deltagrid <- (1:10)/10
prename = "1D_Discretization_"
sufname = "_Laplace_grid_1to.1by.1/Table_niveaux_effectifs.txt"
preIname = "effective-confidence-level_"
ExtI = ".png"

filename = paste(preIname, "Laplace", ExtI, sep="")
Y = read_levels(prename, sufname, c(30, 50, 100, 400), alphagrid, deltagrid)
plot_levels(Y, alphagrid, "n",
            "Laplace for various n",
            filename)

sufname = "_Pareto2.1_grid_1to.1by.1/Table_niveaux_effectifs.txt"
filename = paste(preIname, "centered-Pareto_2.1", ExtI, sep="")
Y = read_levels(prename, sufname,
                c(30, 50, 100, 400, 1000, 5000, 10000, 20000, 50000),
                alphagrid, deltagrid)
plot_levels(Y, alphagrid, "n",
            "Centered Pareto of shape 2.1 for various n",
            filename)

sufname = "_Pareto3_grid_1to.1by.1/Table_niveaux_effectifs.txt"
filename = paste(preIname, "centered-Pareto_3", ExtI, sep="")
Y = read_levels(prename, sufname, c(30, 50, 100, 400, 1000, 5000),
                alphagrid, deltagrid)
plot_levels(Y, alphagrid, "n",
            "Centered Pareto of shape 3 for various n",
            filename)

prename = "1D_Discretization_100_Pareto"
sufname = "_grid_1to.1by.1/Table_niveaux_effectifs.txt"
filename = paste(preIname, "n100_centered-Paretos", ExtI, sep="")
Y = read_levels(prename, sufname, c(2.1, 2.4, 3, 4, 10),
                alphagrid, deltagrid)
plot_levels(Y, alphagrid, "Pareto shape",
            "Centered Pareto with n = 100 for various shapes",
            filename)

prename = "2D_Discretization_"
sufname = "_Pareto3_grid_1to.1by.1/Table_niveaux_effectifs.txt"
filename = "2D_effective-confidence-level_centered-Pareto_3.png"
Y = read_levels(prename, sufname, c(50, 100, 200, 400),
                alphagrid, deltagrid)
plot_levels(Y, alphagrid, "n",
            "2D centered Pareto of shape 3 for various n",
            filename)

