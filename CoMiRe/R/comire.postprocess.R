# postetrior predictive check
post.pred.check <- function(x.cpoints, fit, mcmc, max.x=max(x), J)
{
  index <- c((mcmc$nb+1):(mcmc$nrep+mcmc$nb))[1:((mcmc$nrep)/mcmc$thin)*mcmc$thin]
  x.grid <- seq(0, max.x, length=100)
  y.grid <- seq(min(y)-sqrt(var(y)), max(y) + sqrt(var(y)), length = 100)	
  y.pred <- matrix(NA, length(index), length(x.cpoints))
  knots <- seq(0, max.x, length=J-3)
  phi <- function(x) iSpline(x, df=3, knots = knots, Boundary.knots=c(0,max.x))
  ite <- 1
  for(i in index)
  {
    for(j in 1:length(x.cpoints))
    {
      beta <-  as.double(phi(x.cpoints[j]) %*% fit$mcmc$theta[i,])
      pii <- round(c(beta * fit$mcmc$pi1[i], (1-beta)*fit$mcmc$pi0[i,]),8)
      if(pii[1]==1) ind <- 1
      else {
        ind <- sample(1:length(pii), 1, prob=pii)
      }
      y.pred[ite,j] <- rnorm(1, c(fit$mcmc$mu1[i], fit$mcmc$mu0[i,])[ind], 
                             1/sqrt(c(fit$mcmc$tau1[i], fit$mcmc$tau0[i,])[ind]))
    }
    ite <- ite + 1
  }
  return(y.pred)
}
# posterior mean density
fit.cdf.mcmc <- function(x, y.grid, fit, mcmc, H=10, max.x)
{
  index <- c((mcmc$nb+1):(mcmc$nrep+mcmc$nb))[1:((mcmc$nrep)/mcmc$thin)*mcmc$thin]
  knots <- seq(0, max.x, length=J-3)
  phi <- function(x) iSpline(x, df=3, knots = knots, Boundary.knots=c(0,max.x))
  res <- rep(NA, length(index))
  dens.mcmc <- function(i, pi0, mu0, tau0, mu1, tau1, theta, phi_x, y.grid)
  {
    beta_x <-  as.double(phi_x %*% theta[i,])
    res <-  beta_x * dnorm(y.grid, mu1[i], 1/sqrt(tau1[i]))
    for(h in 1:ncol(pi0))
    {
      res <- res + (1-beta_x)*pi0[i,h] * dnorm(y.grid, mu0[i,h], 1/sqrt(tau0[i,h]))     
    }
    res
  }
  res <- sapply(index, dens.mcmc, fit$mcmc$pi0, fit$mcmc$mu0, fit$mcmc$tau0,
                fit$mcmc$mu1, fit$mcmc$tau1, fit$mcmc$theta, phi(x), y.grid)
  cbind( rowMeans(res), t(apply(res,1,quantile, probs=c(0.025, 0.975))))
}
# additional risk plot
add.risk <- function(a=37, fit, mcmc, xgrid, y, alpha=0.05)
{
  index <- c((mcmc$nb+1):(mcmc$nrep+mcmc$nb))[1:((mcmc$nrep)/mcmc$thin)*mcmc$thin]
  r_a <- matrix(NA, length(index), length(xgrid))
  y.grid <- seq(min(y)-sqrt(var(y)), max(y) + sqrt(var(y)), length = 100)
  low.int <- floor(min(y))
  ite <- 1
  for(i in index)
  {
    f0 <- approxfun(x=y.grid, y=fit$mcmc$f0[i,])
    f1 <- approxfun(x=y.grid, y=fit$mcmc$f1[i,])
    r_a[ite,] <- fit$mcmc$beta[i,]*(integrate(function(x)f1(x), low.int, a)$value-
                                            integrate(function(x)f0(x), low.int, a)$value )
    
    ite <- ite + 1
  }
  risk.data <- data.frame(apply(r_a,2,mean), 
                        apply(r_a,2,quantile, probs=alpha/2), 
                        apply(r_a,2,quantile, probs=1-alpha/2), 
                        xgrid)
colnames(risk.data) <- c("risk","low","upp", "x")
list(mcmc.risk=r_a, summary.risk= risk.data)
}
# benchmark dose
BMD <- function(level, risk, x, alpha=0.05)
{
  bmd <- function(r,q,range=c(0,max(x))) 
  {
    ris <- splinefun(y=r,x=x)
    if(ris(range[2])-q<0) return(NA)
    else
      uniroot(function(x)ris(x)-q,range)$root
  }
  bmd.apply <- function(q) apply(risk, 1, bmd, q=q, range=c(0,180))
  if(length(level)==1){
    bmd.data <- bmd.apply(level)
    return(c(mean(bmd.data), quantile(bmd.data, probs=c(alpha/2,1-alpha/2,alpha))))
  }
  else
  {
  bmd.mcmc.comire <-  sapply(level, bmd.apply)
  bmd.data <- data.frame(level,
                         colMeans(bmd.mcmc.comire, na.rm=TRUE), 
                         t(apply(bmd.mcmc.comire,2,quantile, prob=c(alpha/2,1-alpha/2), na.rm=TRUE)),
                         apply(bmd.mcmc.comire,2,quantile, prob=alpha, na.rm=TRUE))
  colnames(bmd.data) <- c("q", "BMD","low","upp", "BMDL")
  return(bmd.data)
  }
}
# risk plot
riskplot <- function(risk.data, xlabel="x", x=NULL)
  {
  out <- ggplot(risk.data, aes(x,risk)) + geom_line(lty=1, col=4) + 
  geom_ribbon(aes(ymax=upp, ymin=low), fill=4,alpha=.1) +
  labs(y=expression(R[A](x, a)), x=xlabel)+ theme_bw() + 
  theme(plot.margin=unit(c(1,0,0,0),"lines")) + 
  coord_cartesian(ylim=c(0,.8), xlim=c(0,180)) 
  if(is.null(x)) return(out)
  else{
    onlyx <- data.frame(x, zero=rep(0,length(x)))
    out <- out + geom_point(data=onlyx, aes(x, zero), alpha=1, cex=.5, pch="|")
  }
  return(out)
}
# bmd plot
bmd.plot <- function(bmd.data)
{
  ggplot(bmd.data, aes(q,BMD)) + 
  geom_line(aes(q, BMD), lty=1, col=4) + geom_ribbon(aes(ymax=upp, ymin=low), fill=4,alpha=.1) +
  labs(y=expression(BMD[q]), x="q")+ theme_bw() + 
  theme(plot.margin=unit(c(1,0,0,0),"lines"))  
}