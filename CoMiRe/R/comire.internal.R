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
mixdensity <- function(i, y, pi, mu, tau)
{
	kernels <- dnorm(y[i], mu, sqrt(1/tau))
	fji <- sum(pi * kernels)
	fji
}
#
pssq <- function(index, data, cluster, locations) sum((data[cluster==index] - locations[index])^2)