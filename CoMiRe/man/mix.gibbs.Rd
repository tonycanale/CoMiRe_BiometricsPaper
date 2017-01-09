\name{mix.gibbs}
\alias{mix.gibbs}
\title{
Gibbs sampler for standard mixture of Gaussians}
\description{
Posterior inference via Gibbs sampler for standard mixture of Gaussians}
\usage{
mix.gibbs(y, grid = NULL, mcmc, prior, state = NULL, seed)
}
\arguments{
  \item{y}{
Numeric vector for the response}
  \item{grid}{
a list giving the parameters for plotting the posterior mean density over finite grid. It must include the following values:\code{grids}, logical value (if \code{TRUE} the provided grid is used, otherwise standard grid is automatically used); \code{xgrid} a numerical vector with the actual values of the grid for y.}
  \item{mcmc}{
a list giving the MCMC parameters. It must include the
          following integers: \code{nb} giving the number of burn-in iterations, \code{nrep} giving
          the total number of iterations, \code{thin} giving the thinning interval, \code{ndisplay} giving the multiple of iterations to be displayed on screen while the algorithm is running (a message will be printed every \code{ndisplay} iterations)}
  \item{prior}{
a list containing the values of the hyperparameters. It must include the following values: \code{mu0}, the prior mean for each location paramter, \code{kappa} the prior precision for each location paramter, \code{a} and \code{b} prior scale and shape parameter for the gamma distribution of each precision parameter, \code{H} total number of components, and \code{alpha} prior for the mixture weights.}
  \item{state}{
a list giving the current value of the parameters. This list is used if the current analysis is the continuation of a
          previous analysis or if we want to start the MCMC algorithm from some particular value of the parameters.}
  \item{seed}{
Seed for random initialization}
}
\value{
A list containing 
\item{post.means}{A list containing the posterior mean of all the mixture parameters (it does not check or account for label switching!)}
\item{ci}{A list containing the 95 credible intervals for all the quantities stored in \code{post.means}}
\item{mcmc}{A list containing all the MCMC chains}
}
\references{
Canale, A., Durante, D., and Dunson, D. (2017), Convex Mixture Regression for Quantitative Risk Assessment }
\author{
Antonio Canale
}
\keyword{mixture}
