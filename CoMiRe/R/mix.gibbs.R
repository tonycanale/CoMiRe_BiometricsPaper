mix.gibbs <-
function(y, grid=NULL, mcmc, prior, state=NULL, seed)
{
  #internal working variables
  n <- length(y)
  H <- prior$H
  print_now <- c(mcmc$nb + 1:mcmc$ndisplay*(mcmc$nrep)/mcmc$ndisplay)
  
  #create the objects to store the MCMC samples
  pi <- matrix(NA, mcmc$nrep+mcmc$nb, H)
  mu <- tau <- matrix(NA, mcmc$nrep+mcmc$nb, H)
  
  #initialize each quantity:
  #paramters of the model
  if(is.null(state))
  {
    set.seed(seed)
    pi[1,] = rep(1/H, H)
    tau[1,] = rgamma(H, prior$a, prior$b)
    mu[1,] = sort(rnorm(H, prior$mu0, prior$kappa/tau[1,]))
  }
  else
  {
    pi[1,] = state$pi
    mu[1,] =state$mu
    tau[1,] = state$tau
  }
 occ <- rep(1, (mcmc$nrep+mcmc$nb))
  #start the MCMC simulation 
  set.seed(seed)
  for(ite in 2:(mcmc$nrep+mcmc$nb))
  {
    #0. Print the iteration
    if(ite==mcmc$nb) cat("Burn in done\n")
    if(ite %in% print_now) cat(ite, "iterations over", mcmc$nrep+mcmc$nb, "\n")

    # Update the cluster labels 
    c = sapply(1:n, labelling_mix, y=y, pi=pi[ite-1,], mu=mu[ite-1,], tau=tau[ite-1,])

    # Update the mixture weights sampling from dirichlet
    n_h = table(factor(c, levels=1:H))
    occ[ite] <- length(table(c))
    pi[ite,] = as.double(rdirichlet2(as.double(prior$alpha + n_h)))
        
    #6. Updated mu and tau from the usual normal inverse gamma
    hat_a = prior$a + n_h/2
    mean_h = tapply(y, factor(c, levels=1:H), mean)
    deviance_h = (n_h-1)*tapply(y, factor(c, levels=1:H), var)
    mean_h[n_h==0]  = 0
    deviance_h[n_h<=1] = 0
    hat_b = prior$b + 0.5*(deviance_h + n_h/(1+prior$kappa*n_h)*(mean_h-prior$mu0)^2)
    hat_kappa = 1/(1/prior$kappa + n_h)
    hat_mu = hat_kappa*(1/prior$kappa*prior$mu0 + n_h*mean_h)
    tau[ite, ] = rgamma(H, hat_a, hat_b)
    mu[ite, ] = rnorm(H, hat_mu, sqrt(hat_kappa/tau[ite,]))
  }
  
  post.mean.mu <- colMeans(mu[-c(1:mcmc$nb),][1:((mcmc$nrep)/mcmc$thin)*mcmc$thin,])
  post.mean.tau <- colMeans(tau[-c(1:mcmc$nb),][1:((mcmc$nrep)/mcmc$thin)*mcmc$thin,])
  post.mean.pi <- colMeans(pi[-c(1:mcmc$nb),][1:((mcmc$nrep)/mcmc$thin)*mcmc$thin,])
  
  ci.mu <- apply(mu[-c(1:mcmc$nb),][1:((mcmc$nrep)/mcmc$thin)*mcmc$thin,], 2, quantile, probs=c(0.025, 0.975))
  ci.tau <- apply(tau[-c(1:mcmc$nb),][1:((mcmc$nrep)/mcmc$thin)*mcmc$thin,], 2, quantile, probs=c(0.025, 0.975))
  ci.pi <- apply(pi[-c(1:mcmc$nb),][1:((mcmc$nrep)/mcmc$thin)*mcmc$thin,], 2, quantile, probs=c(0.025, 0.975))
  list(
    post.means=list(mu=post.mean.mu, tau=post.mean.tau, pi=post.mean.pi),
    ci = list(mu=ci.mu, tau=ci.tau, pi=ci.pi),
    mcmc = list(mu=mu, tau=tau, pi=pi), occ=occ)
}
