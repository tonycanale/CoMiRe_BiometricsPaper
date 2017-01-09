comire.gibbs <-
function(y, x, basistype="BetaCDF", basisX, grid=NULL, mcmc, prior, state=NULL, seed, max.x=1)
{
  #internal working variables
  n <- length(y)
  H <- prior$H
  print_now <- c(mcmc$nb + 1:mcmc$ndisplay*(mcmc$nrep)/mcmc$ndisplay)
  
  #create the objects to store the MCMC samples
  pi0 <- matrix(NA, mcmc$nrep+mcmc$nb, H)
  mu0 <- tau0 <- matrix(NA, mcmc$nrep+mcmc$nb, H)
  mu1 <- tau1 <- pi1 <- rep(NA, mcmc$nrep+mcmc$nb)
  w <- matrix(NA, mcmc$nrep+mcmc$nb, length(prior$dirpar))
  
  #initialize each quantity:
  #paramters of the model
  if(is.null(state))
  {
    set.seed(seed)
    pi0[1,] = rep(1/H, H)
    pi1[1] = 1
    tau0[1,] = rgamma(H, prior$a, prior$b)
    mu0[1,] = rnorm(H, prior$mu0, prior$kappa/tau0[1,])
    tau1[1] = rgamma(1, prior$a, prior$b)
    mu1[1] = rnorm(1, prior$mu0, prior$kappa/tau1[1])
    w[1,] = prior$dirpar/sum(prior$dirpar)
  }
  else
  {
    pi0[1,] = state$pi0
    pi1[1] = 1
    mu0[1,] =state$mu0
    tau0[1,] = state$tau0
    mu1[1] =state$mu1
    tau1[1] = state$tau1
    w[1,] = state$w
  }
  if(is.null(grid$grids))
  {
    x.grid <- seq(0,250, length=250)
    y.grid <- seq(min(y)-sqrt(var(y)), max(y) + sqrt(var(y)), length = 100)	
    beta_x = matrix(NA, mcmc$nrep+mcmc$nb, length(x.grid))	
    f0 = matrix(NA, mcmc$nrep+mcmc$nb, length(y.grid))
    f1 = matrix(NA, mcmc$nrep+mcmc$nb, length(y.grid))
  }
  else
  {
    x.grid = grid$xgrid
    y.grid = grid$ygrid
    beta_x = matrix(NA, mcmc$nrep+mcmc$nb, length(x.grid))	
    f0 = matrix(NA, mcmc$nrep+mcmc$nb, length(y.grid))
    f1 = matrix(NA, mcmc$nrep+mcmc$nb, length(y.grid))
  }
  
  # basis of X
  if(basistype=="BetaCDF")
  {
    J <- prior$J
    x.std <- x/max.x
    Fl <- matrix(NA,length(x), 2^J)
    for(l in 1:(2^J))
    {
      Fl[,l] <- pbeta(x.std, l, 2^J-l+1)
    }
    Fl.grid <- matrix(NA,length(x.grid), 2^J)
    for(l in 1:(2^J))
    {
      Fl.grid[,l] <- pbeta(x.grid, l, 2^J-l+1)
    }
  }
  if(basistype=="empirical")
  {
    Fl <- basisX(x)
    J <- NCOL(Fl)
    Fl.grid <- basisX(x.grid)
  }
  if((basistype!="BetaCDF") & (basistype!="empirical")) 
    stop("basistype need to be 'empirical' or 'BetaCDF'\n")
    
  #beta_i is the interpolating function evaluated at x_i
  beta_i =  as.double(Fl %*% w[1, ])
  
  # f0i and f1i
  f0i = sapply(1:n, mixdensity, y=y, pi=pi0[1,], mu=mu0[1,], tau=tau0[1,])
  f1i = dnorm(y, mu1[1], sqrt(1/tau1[1]))
  
  #start the MCMC simulation 
  set.seed(seed)
  for(ite in 2:(mcmc$nrep+mcmc$nb))
  {
    #0. Print the iteration
    if(ite==mcmc$nb) cat("Burn in done\n")
    if(ite %in% print_now) cat(ite, "iterations over", mcmc$nrep+mcmc$nb, "\n")
    
    #1. Update d_i marginalising out b_i from
    d = rbinom(n, 1, prob=(beta_i*f1i)/((1-beta_i)*f0i + beta_i*f1i))
    
    #2. Update b_i from the multinomial 
    b = sapply(1:n, labelling_b, w[ite-1,], Fl=Fl, f0i=f0i, f1i=f1i)
    
    #3. Update c_i, marginalizing over b_i and d_i, from the multinomial 
    ind0 <- c(1:n)[d==0]
    ind1 <- c(1:n)[d==1]
    c0 = sapply(ind0, labelling_c, y=y, pi=pi0[ite-1,], mu=mu0[ite-1,], tau=tau0[ite-1,])

    #4. Update the mixture weights sampling from dirichlet
    n_0h = table(factor(c0, levels=1:H))
    pi0[ite,] = as.double(rdirichlet(1, as.double(prior$alpha + n_0h)))
    pi1[ite] = 1
    
    #5. Update w from the Dirichlet and obtain an updated function beta_i
    dirpar.post = as.double(prior$dirpar + table(factor(b, levels=1:length(w[ite-1,]))))
    w[ite, ] = as.double(rdirichlet(1, dirpar.post))
    beta_i = as.numeric(Fl %*% w[ite, ])
    beta_i[beta_i>1] <- 1
    beta_i[beta_i<0] <- 0
    
    #6. Updated mu and tau from the usual normal inverse gamma
    # in cluster 0
    n_h = table(factor(c0, levels=1:H))
    hat_a = prior$a + n_h/2
    mean_h = tapply(y[ind0], factor(c0, levels=1:H), mean)
    deviance_h = (n_h-1)*tapply(y[ind0], factor(c0, levels=1:H), var)
    mean_h[n_h==0]  = 0
    deviance_h[n_h<=1] = 0
    hat_b = prior$b + 0.5*(deviance_h + n_h/(1+prior$kappa*n_h)*(mean_h-prior$mu0)^2)
    hat_kappa = 1/(1/prior$kappa + n_h)
    hat_mu = hat_kappa*(1/prior$kappa*prior$mu0 + n_h*mean_h)
    tau0[ite, ] = rgamma(H, hat_a, hat_b)
    mu0[ite, ] = rtruncnorm(H, a=mu1[ite-1], b=Inf, hat_mu, sqrt(hat_kappa/tau0[ite,]))
    
    # in cluster 1
    n_h = sum(d==1)
    hat_a = prior$a + n_h/2
    mean_h = mean(y[ind1])
    deviance_h = (n_h-1)*var(y[ind1])
    mean_h[n_h==0]  = 0
    deviance_h[n_h<=1] = 0
    hat_b = prior$b + 0.5*(deviance_h + n_h/(1+prior$kappa*n_h)*(mean_h-prior$mu0)^2)
    hat_kappa = 1/(1/prior$kappa + n_h)
    hat_mu = hat_kappa*(1/prior$kappa*prior$mu0 + n_h*mean_h)
    tau1[ite] = rgamma(1, hat_a, hat_b)
    mu1[ite] = rtruncnorm(1, a=-Inf, b=min(mu0[ite, ]), hat_mu, sqrt(hat_kappa/tau1[ite]))# c(30, 35, 40)#
    
    # update the values of the densities in the observed points
    f0i = sapply(1:n, mixdensity, y=y, pi=pi0[ite,], mu=mu0[ite,], tau=tau0[ite,])
    f1i = dnorm(y, mu1[ite], sqrt(1/tau1[ite]))
    
    #7. compute some posterior quanities of interest
    beta_x[ite, ] = Fl.grid %*% w[ite, ]
    f0[ite, ] = sapply(1:length(y.grid), mixdensity, y=y.grid, pi=pi0[ite,], mu=mu0[ite,], tau=tau0[ite,])
    f1[ite, ] = dnorm(y.grid, mu1[ite], sqrt(1/tau1[ite]))
    
  }
  post.mean.beta <- colMeans(beta_x[-c(1:mcmc$nb),])
  post.mean.mu0 <- colMeans(mu0[-c(1:mcmc$nb),])
  post.mean.tau0 <- colMeans(tau0[-c(1:mcmc$nb),])
  post.mean.mu1 <- mean(mu1[-c(1:mcmc$nb)])
  post.mean.tau1 <- mean(tau1[-c(1:mcmc$nb)])
  post.mean.f0 <- colMeans(f0[-c(1:mcmc$nb),])
  post.mean.f1 <- colMeans(f1[-c(1:mcmc$nb),])
  post.mean.pi0 <- colMeans(pi0[-c(1:mcmc$nb),])
  post.mean.pi1 <- mean(pi1[-c(1:mcmc$nb)])
  post.mean.theta <- colMeans(w[-c(1:mcmc$nb),])
  
  
  ci.beta <- apply(beta_x[-c(1:mcmc$nb),], 2, quantile, probs=c(0.025, 0.975))
  ci.mu0 <- apply(mu0[-c(1:mcmc$nb),], 2, quantile, probs=c(0.025, 0.975))
  ci.tau0 <- apply(tau0[-c(1:mcmc$nb),], 2, quantile, probs=c(0.025, 0.975))
  ci.mu1 <- quantile(mu1[-c(1:mcmc$nb)], probs=c(0.025, 0.975))
  ci.tau1 <- quantile(tau1[-c(1:mcmc$nb)], probs=c(0.025, 0.975))
  ci.f0 <- apply(f0[-c(1:mcmc$nb),], 2, quantile, probs=c(0.025, 0.975))
  ci.f1 <- apply(f1[-c(1:mcmc$nb),], 2, quantile, probs=c(0.025, 0.975))
  ci.pi0 <- apply(pi0[-c(1:mcmc$nb),], 2, quantile, probs=c(0.025, 0.975))
  ci.pi1 <- quantile(pi1[-c(1:mcmc$nb)], probs=c(0.025, 0.975))
  list(
    post.means=list( beta=post.mean.beta, theta=post.mean.theta, 
                     mu0=post.mean.mu0, tau0=post.mean.tau0, 
                     mu1=post.mean.mu1, tau1=post.mean.tau1, 
                     f0=post.mean.f0, f1=post.mean.f1, pi0=post.mean.pi0, pi1=post.mean.pi1),
    ci = list( beta=ci.beta, mu0=ci.mu0, tau0=ci.tau0, 
               mu1=ci.mu1, tau1=ci.tau1, 
               f0=ci.f0,f1=ci.f1, pi0=ci.pi0, pi1=ci.pi1),
    mcmc = list(beta=beta_x, theta=w, mu0=mu0, tau0=tau0, mu1=mu1, tau1=tau1, f0=f0, f1=f1, pi0=pi0, pi1=pi1))
}
