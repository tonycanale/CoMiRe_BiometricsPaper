labelling_b <- function(i, w, phi, f0i, f1i)
{
  probs <- w*((1-phi[i,])*f0i[i] + phi[i,]*f1i[i])
  if(any(probs<0)) probs[probs<0]=0
  sample(1:ncol(phi), 1, prob=probs)
}
#
labelling_c <- function(i, y, pi, mu, tau)
{
  probs <- pi*dnorm(y[i], mu, sqrt(1/tau))
  if(any(probs<0)) probs[probs<0]=0
  sample(1:length(mu), 1, prob=probs)
}
#
labelling_mix <- function(i, y, pi, mu, tau)
{
  probs <- pi*dnorm(y[i], mu, sqrt(1/tau))
  if(any(probs<0)) probs[probs<0]=0
  sample(1:length(mu), 1, prob=probs)
}
#
labelling_ddp <- function(i, y, Xreg, pi, beta1, beta2, tau)
{
  mu <- as.double(Xreg[i,] %*% rbind(beta1,beta2))
  probs <- pi*dnorm(y[i], mu, sqrt(1/tau))
  sample(1:length(mu), 1, prob=probs)
}
#
sample_beta <- function(ind, c, Xreg, y, tau, S0inv, beta0, atau, btau)
{
  n_h <- sum(c==ind)
  if(n_h==0)
  {
    var <- ginv(S0inv)
    coeffs <- rmnorm(1, beta0, var)
    precision <- rgamma(1, atau, btau)    
  }
  if(n_h==1)
  {
    var <- ginv(tau[ind]*matrix(t(Xreg[c==ind,])) %*% Xreg[c==ind,] + S0inv)
    mean <- as.double(var %*% (tau[ind]*matrix(t(Xreg[c==ind,])) %*% y[c==ind] + S0inv %*% beta0))
    coeffs <- rmnorm(1, mean, var)
    precision <- rgamma(1, atau + n_h/2, btau + 0.5*sum((y[c==ind] - Xreg[c==ind,] %*% coeffs)^2))  
  }
  if(n_h>1)
  {
  var <- ginv(tau[ind]*t(Xreg[c==ind,]) %*% Xreg[c==ind,] + S0inv)
  mean <- as.double(var %*% (tau[ind]*t(Xreg[c==ind,]) %*% y[c==ind] + S0inv %*% beta0))
  coeffs <- rmnorm(1, mean, var)
  precision <- rgamma(1, atau + n_h/2, btau + 0.5*sum((y[c==ind] - Xreg[c==ind,] %*% coeffs)^2))
  }
  c(coeffs, precision)
  }
#
mixdensity <- function(i, y, pi, mu, tau)
{
	kernels <- dnorm(y[i], mu, sqrt(1/tau))
	fji <- sum(pi * kernels)
	fji
}
#
rdirichlet2 <- function(alpha) 
{
    x <- rgamma(length(alpha), alpha, 1)
    x <- x/sum(x)
	if(any(x<=0 | x>=1))
		rdirichlet2(alpha)
	else return(x)
}