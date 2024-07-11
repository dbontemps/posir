
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/v1/examples>. -->

# posir

<!-- badges: start -->
<!-- badges: end -->

## A Research project

posir is a help package to compute simulations in the framework defined
in the preprint “Inference post region selection” by Bontemps, Bachoc,
and Neuvial, 2024.

## What is it exactly for?

The goal of posir is to simulate trajectories of 1D or 2D POSIR
processes, as defined in the paper. This allows to estimate quantiles as
well as effective error levels of simutaneous confidence intervals.

## Installation

You can install the development version of posir like so:

``` r
# posir can be compilated from Rstudio (or even simple R):
# 1) go to the posir directory
# 2) open the Rstudio project "posir.Rproj"
# 3) run the following:
#    > library(devtools)
#    > install()
```

## Available results

The results of the simulations presented in the paper, as well as the R
scripts it used to exploit the package, are also included in the
installation.
