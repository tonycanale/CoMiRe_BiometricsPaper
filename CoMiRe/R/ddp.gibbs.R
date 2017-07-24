ddp.gibbs <-
function(y, x, grid=NULL, mcmc, prior, state=NULL, seed)
{
  #internal working variables
  n <- length(y)
  H <- prior$H
  print_now <- c(mcmc$nb + 1:mcmc$ndisplay*(mcmc$nrep)/mcmc$ndisplay)
  
  #create the objects to store the MCMC samples
  pi <- matrix(NA, mcmc$nrep+mcmc$nb, H)
  beta1 <- beta2 <- tau <- matrix(NA, mcmc$nrep+mcmc$nb, H)
  S0inv <- solve(prior$S0)
  
  # solo una cazzo di covariata
  X <- cbind(1, x)
  
  #initialize each quantity:
  #paramters of the model
  if(is.null(state))
  {
    set.seed(seed)
    pi[1,] = rep(1/H, H)
    tau[1,] = rgamma(H, prior$a, prior$b)
    beta1[1,] = sort(rnorm(H, mean(y), 1/tau[1,]))
    beta2[1,] = 0
  }
  else
  {
    pi[1,] = state$pi
    beta1[1,] =state$mu
    beta2[1,] =state$mu
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
    c = sapply(1:n, labelling_ddp, y=y, Xreg=X, pi=pi[ite-1,], beta1=beta1[ite-1,],
               beta2=beta2[ite-1,], tau=tau[ite-1,])
    
    # Update the mixture weights sampling from dirichlet
    n_h = table(factor(c, levels=1:H))
    occ[ite] <- length(table(c))
    pi[ite,] = as.double(rdirichlet2(as.double(prior$alpha + n_h)))
        
    #6. Updated mu and tau from the usual normal inverse gamma
    hat_a = prior$a + n_h/2
    pars = sapply(X=1:H, FUN=sample_beta, c=c, Xreg=X, y=y, 
                     tau=tau[ite-1, ], S0inv=S0inv, beta0=prior$beta0, prior$a, prior$b)
    beta1[ite, ] <- pars[1,]
    beta2[ite, ] <- pars[2,]
    tau[ite, ] <- pars[3,]
  }
  
  post.mean.beta1 <- colMeans(beta1[-c(1:mcmc$nb),][1:((mcmc$nrep)/mcmc$thin)*mcmc$thin,])
  post.mean.beta2 <- colMeans(beta2[-c(1:mcmc$nb),][1:((mcmc$nrep)/mcmc$thin)*mcmc$thin,])
  post.mean.tau <- colMeans(tau[-c(1:mcmc$nb),][1:((mcmc$nrep)/mcmc$thin)*mcmc$thin,])
  post.mean.pi <- colMeans(pi[-c(1:mcmc$nb),][1:((mcmc$nrep)/mcmc$thin)*mcmc$thin,])
  
  ci.beta1 <- apply(beta1[-c(1:mcmc$nb),][1:((mcmc$nrep)/mcmc$thin)*mcmc$thin,], 2, quantile, probs=c(0.025, 0.975))
  ci.beta2 <- apply(beta2[-c(1:mcmc$nb),][1:((mcmc$nrep)/mcmc$thin)*mcmc$thin,], 2, quantile, probs=c(0.025, 0.975))
  ci.tau <- apply(tau[-c(1:mcmc$nb),][1:((mcmc$nrep)/mcmc$thin)*mcmc$thin,], 2, quantile, probs=c(0.025, 0.975))
  ci.pi <- apply(pi[-c(1:mcmc$nb),][1:((mcmc$nrep)/mcmc$thin)*mcmc$thin,], 2, quantile, probs=c(0.025, 0.975))
  list(
    post.means=list(beta1=post.mean.beta1, beta2=post.mean.beta2, 
                    tau=post.mean.tau, pi=post.mean.pi),
    ci = list(beta1=ci.beta1, beta2=ci.beta2, tau=ci.tau, pi=ci.pi),
    mcmc = list(beta1=beta1, beta2= beta2,  tau=tau, pi=pi), occ=occ)
}
