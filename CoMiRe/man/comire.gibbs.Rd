\name{comire.gibbs}
\alias{comire.gibbs}
\title{
Gibbs sampler for CoMiRe model}
\description{
Posterior inference via Gibbs sampler for CoMiRe model}
\usage{
comire.gibbs(y, x, grid = NULL, mcmc, prior, state = NULL, seed, max.x = max(x))
}
\arguments{
  \item{y}{
Numeric vector for the response}
  \item{x}{
Numeric vector for the covariate
}
  \item{grid}{
a list giving the parameters for plotting the posterior mean density and the posterior mean beta(x) over finite grids. It must include the following values:\code{grids}, logical value (if \code{TRUE} the provided grids are used, otherwise standard grids are automatically used); \code{xgrid} and \code{ygrid}, numerical vectors with the actual values of the grid for y and x.}
  \item{mcmc}{a list giving the MCMC parameters. It must include the
          following integers: \code{nb} giving the number of burn-in iterations, \code{nrep} giving
          the total number of iterations, \code{thin} giving the thinning interval, \code{ndisplay} giving the multiple of iterations to be displayed on screen while the algorithm is running (a message will be printed every \code{ndisplay} iterations)}
  \item{prior}{
a list containing the values of the hyperparameters. It must include the following values: \code{mu0}, the prior mean for each location paramter, \code{kappa} the prior variance for each location paramter, \code{a} and \code{b} prior scale and shape parameter for the gamma distribution of each precision parameter, \code{dirpar}, numeric vector of size J for the Dirichlet prior on the beta basis weights, \code{J}, parameter controlling the number of Beta basis (used only if \code{basistype="BetaCDF"}),\code{H} total number of components in the mixture at x_0, and \code{alpha} prior for the mixture weights.
}
  \item{state}{
a list giving the current value of the parameters. This list is used if the current analysis is the continuation of a
          previous analysis or if we want to start the MCMC algorithm from some particular value of the parameters.}
  \item{seed}{
Seed for random initialization}
  \item{max.x}{
Maximum value allowed for \code{x}. Used only when \code{basistype="BetaCDF"} to transform the covariate in [0,1]}
}
\value{
A list containing 
\item{post.means}{A list containing the posterior mean density beta over the grid, of all the mixture parameters and of f_{x_0} and f_{x_1} over the \code{y.grid}}
\item{ci}{A list containing the 95 credible intervals for all the quantities stored in \code{post.means}}
\item{mcmc}{A list containing all the MCMC chains}
}
\references{
Canale, A., Durante, D., and Dunson, D. (2017), Convex Mixture Regression for Quantitative Risk Assessment 
}
\author{
Antonio Canale}
\keyword{CoMiRe}% use one of  RShowDoc("KEYWORDS")
