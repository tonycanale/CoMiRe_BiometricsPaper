# CoMiRe: Convex Mixture Regression for Quantitative Risk Assessment 

This repository contains the R functions and scripts to implement the **Convex Mixture Regression model** of

- Canale, A., Durante, D., and Dunson, D. B., (2018), [Convex Mixture Regression for Quantitative Risk Assessment](http://arxiv.org/abs/1701.02950), __Biometrics__, to appear

Install the R package `CoMiRe_VERSION.tar.gz` with `R CMD INSTALL CoMiRe_VERSION.tar.gz` or

```R
install.packages("devtools")
devtools::install_github("tonycanale/CoMiRe/CoMiRe")
```

and then follow the tutorial [Tutorial.md](Tutorial.md) to reproduce the analysis of the US Collaborative Perinatal Project (CPP) data on the effect of Dichlorodiphenyldichloroethylene (DDE) on premature delivery.

The folder `./CoMiRe` contains the package source files.

The output of [Tutorial.md](Tutorial.md) are obtained using R version 3.4.0
