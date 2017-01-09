labelling_b <- function(i, w, Fl, f0i, f1i)
{
  probs <- w*((1-Fl[i,])*f0i[i] + Fl[i,]*f1i[i])
  sample(1:ncol(Fl), 1, prob=probs)
}
#
labelling_c <- function(i, y, pi, mu, tau)
{
  probs <- pi*dnorm(y[i], mu, sqrt(1/tau))
  sample(1:length(mu), 1, prob=probs)
}
#
labelling_mix <- function(i, y, pi, mu, tau)
{
  probs <- pi*dnorm(y[i], mu, sqrt(1/tau))
  sample(1:length(mu), 1, prob=probs)
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