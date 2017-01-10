# ---
# Analysis of the US Collaborative Perinatal Project (CPP) data on
# the effect of Dichlorodiphenyldichloroethylene (DDE) on premature delivery 
# [Longnecker et al., 2001, "Association between maternal serum concentration of the DDT metabolite 
#  DDE and preterm and small-for- gestational-age babies at birth.", The Lancet, 358]
# as discussed in 
# [Canale, Durante, Dunson, 2017, "Convex Mixture Regression for Quantitative Risk Assessment"]
# ---

# load the CoMiRe set of functions of the 'CoMiRe' package
# install the package with R CMD INSTALL CoMiRe_*.tar.gz before running this script
require(CoMiRe)

# load the CPP data
load("cpp.rda")
x = cpp$dde    # DDE 
y = cpp$gest   # gestational week at delivery
n <- NROW(cpp)

# Bayesian preprocess
H <-10
prior <- list(mu0=mean(y), dirpar=rep(1, 6), kappa=1, a=2, b=2, H=10, J=6, alpha=1)
mcmc <- list(nrep=5000, nb=4000, thin=5, ndisplay=8)
fit.bayes.start <- mix.gibbs(y, mcmc=mcmc, prior=prior, seed=1)
  
# start state + prior centering in the point estimates of the plain mixture
start.state = list(mu0=sort(fit.bayes.start$post.means$mu), tau0=fit.bayes.start$post.means$tau[order(fit.bayes.start$post.means$mu)], 
                   mu1=mean(fit.bayes.start$post.means$mu), tau1=mean(fit.bayes.start$post.means$tau), 
                   pi0=fit.bayes.start$post.means$pi[order(fit.bayes.start$post.means$mu)], pi1=1, 
                   w=prior$dirpar/sum(prior$dirpar))

# compute the basis 
source("cpp_basis.R")

# Bayesian estimation 
mcmc <- list(nrep=50000, nb=4000, thin=5, ndisplay=4)
prior <- list(mu0=mean(y), dirpar=rep(1, 6), kappa=1, a=2, b=2, H=10, J=6, alpha=1)
fit.bayes <- comire.gibbs(y, x, basistype="empirical", basisX=basisX, 
                          mcmc=mcmc, prior=prior, state=start.state, seed=1)

# Posterior inference on relevant functionals and figures in cpp_figures.R
