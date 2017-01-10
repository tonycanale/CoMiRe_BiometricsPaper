#-----
# Code to produce the Figures of the paper
#-----

require(ggplot2)
require(gridExtra)

#-----
# Figure 1
#-----
ggplot(data=cpp) + geom_point(aes(x=dde, y=gest), alpha=.5, cex=.5) + 
  geom_smooth(aes(x=dde, y=gest), method="loess", span = 1, col=1) + 
  xlab("Dichlorodiphenyldichloroethylene (DDE)") + ylab("Gestational age at delivery") + theme_bw()


#-----
# Figure 2
#-----

# binning of x into groups
x.cpoints <- c(7.5, 22.5,37.5,52.5,67.5, 125)
groups <- as.numeric(cut(x, breaks = c(0, 15,30,45,60, 75, 180)))

# binning of y
y1 <- cpp[cpp$dde<15,"gest"]
y2 <- cpp[cpp$dde>=15 & cpp$dde<30,"gest"]
y3 <- cpp[cpp$dde>=30 & cpp$dde<45,"gest"]
y4 <- cpp[cpp$dde>=45 & cpp$dde<60,"gest"]
y5 <- cpp[cpp$dde>=60 & cpp$dde<75,"gest"]
y6 <- cpp[cpp$dde>=75,"gest"]

# fit the density for given valus of beta
fit.dens.mcmc <- function(y, fit, ind, H=10)
{
  index <- c((mcmc$nb+1):(mcmc$nrep+mcmc$nb))[1:((mcmc$nrep)/mcmc$thin)*mcmc$thin]
  ite <- 1
  res <- rep(NA, 10000)
  for(i in index)
  {
    res[ite] <-  (fit$mcmc$beta[i,ind]) * dnorm(y, fit$mcmc$mu1[i], 1/sqrt(fit$mcmc$tau1[i]))
    for(h in 1:H)
    {
      res[ite] <- res[ite] + (1-fit$mcmc$beta[i,ind])*fit$mcmc$pi0[i,h] * dnorm(y, fit$mcmc$mu0[i,h], 1/sqrt(fit$mcmc$tau0[i,h]))     }
    ite <- ite + 1
  }
  mean(res)
}

xnow <- seq(25, 48, length=50)

curvefit1 <- rep(NA,50); for(i in 1:50) curvefit1[i] <- fit.dens.mcmc(xnow[i], fit.bayes, 8)
f1 <- function(x)approxfun(x=xnow, y=curvefit1)(x)
p1 <-ggplot(data.frame(y1), aes(y1)) + 
  geom_histogram(aes(y =..density..), fill="grey", alpha=.3, binwidth=0.5) + 
  geom_density(lty=2) + labs(x="Gestational age at delivery (DDE<15)", y="") + 
  theme_bw() + stat_function(fun=f1) + xlim(c(25,48)) + coord_cartesian(ylim=c(0, 0.25)) +
  theme(plot.margin=unit(c(1,0,0,0),"lines"), axis.title=element_text(size=10))

curvefit2 <- rep(NA,50); for(i in 1:50) curvefit2[i] <- fit.dens.mcmc(xnow[i], fit.bayes, 23)
f2 <- function(x)approxfun(x=xnow, y=curvefit2)(x)
p2 <-ggplot(data.frame(y2), aes(y2)) + geom_histogram(aes(y =..density..), fill="grey", alpha=.3, binwidth=0.5) + 
  geom_density(lty=2)  + labs(x="Gestational age at delivery (15<DDE<30)", y="") + theme_bw() + xlim(c(25,48)) + ylim(c(0, 0.25)) +
  stat_function(fun=f2) +theme(plot.margin=unit(c(1,0,0,0),"lines"), axis.title=element_text(size=10))

curvefit3 <- rep(NA,50); for(i in 1:50) curvefit3[i] <- fit.dens.mcmc(xnow[i], fit.bayes, 38)
f3 <- function(x)approxfun(x=xnow, y=curvefit3)(x)
p3 <- ggplot(data.frame(y3), aes(y3)) + geom_histogram(aes(y =..density..), fill="grey", alpha=.3, binwidth=0.5) + 
  geom_density(lty=2) + labs(x=expression("Gestational age at delivery (30<DDE<45)"), y="") + theme_bw() + xlim(c(25,48))+ ylim(c(0, 0.25)) +
  stat_function(fun=f3) +theme(plot.margin=unit(c(1,0,0,0),"lines"), axis.title=element_text(size=10))

curvefit4 <- rep(NA,50); for(i in 1:50) curvefit4[i] <- fit.dens.mcmc(xnow[i], fit.bayes, 53)
f4 <- function(x)approxfun(x=xnow, y=curvefit4)(x)
p4 <- ggplot(data.frame(y4), aes(y4)) + geom_histogram(aes(y =..density..), fill="grey", alpha=.3, binwidth=0.5) + 
  geom_density(lty=2) + labs(x="Gestational age at delivery (45<DDE<60)", y="") + theme_bw() + xlim(c(25,48)) +ylim(c(0, 0.25)) +
  stat_function(fun=f4) +theme(plot.margin=unit(c(1,0,0,0),"lines"), axis.title=element_text(size=10))

curvefit5 <- rep(NA,50); for(i in 1:50) curvefit5[i] <- fit.dens.mcmc(xnow[i], fit.bayes, 68)
f5 <- function(x)approxfun(x=xnow, y=curvefit5)(x)
p5 <- ggplot(data.frame(y5), aes(y5)) + geom_histogram(aes(y =..density..), fill="grey", alpha=.3, binwidth=0.5) + 
  geom_density(lty=2) + labs(x="Gestational age at delivery  (60<DDE<75)", y="") + theme_bw()+ xlim(c(25,48)) +coord_cartesian(ylim=c(0, 0.25))+
  stat_function(fun=f5) +theme(plot.margin=unit(c(1,0,0,0),"lines"), axis.title=element_text(size=10))

curvefit6 <- rep(NA,50); for(i in 1:50) curvefit6[i] <- fit.dens.mcmc(xnow[i], fit.bayes, 126)
f6 <- function(x)approxfun(x=xnow, y=curvefit6)(x)
p6 <- ggplot(data.frame(y6), aes(y6)) + geom_histogram(aes(y =..density..), fill="grey", alpha=.3, binwidth=0.5) + 
  geom_density(lty=2) + labs(x="Gestational age at delivery  (DDE>75)", y="") + theme_bw()  + xlim(c(25,48)) + ylim(c(0, 0.25)) +
  stat_function(fun=f6) +theme(plot.margin=unit(c(1,0,0,0),"lines"), axis.title=element_text(size=10))

grid.arrange(p1, p2, p3, p4, p5, p6, ncol=3, nrow=2)


#---------
# Figure 5
#----------

# ecdf of the binned data
y.grid <- seq(min(y)-sqrt(var(y)), max(y) + sqrt(var(y)), length = 1000)	
ecdf1 <- ecdf(y1)(y.grid)
ecdf2 <- ecdf(y2)(y.grid)
ecdf3 <- ecdf(y3)(y.grid)
ecdf4 <- ecdf(y4)(y.grid)
ecdf5 <- ecdf(y5)(y.grid)
ecdf6 <- ecdf(y6)(y.grid)

# fit the density for given valus of beta
fit.cdf.mcmc <- function(y, fit, ind, H=10)
{
  index <- c((mcmc$nb+1):(mcmc$nrep+mcmc$nb))[1:((mcmc$nrep)/mcmc$thin)*mcmc$thin]
  ite <- 1
  res <- rep(NA, length(index))
  for(i in index)
  {
    res[ite] <-  (fit$mcmc$beta[i,ind]) * pnorm(y, fit$mcmc$mu1[i], 1/sqrt(fit$mcmc$tau1[i]))
    for(h in 1:H)
    {
      res[ite] <- res[ite] + (1-fit$mcmc$beta[i,ind])*fit$mcmc$pi0[i,h] * pnorm(y, fit$mcmc$mu0[i,h], 1/sqrt(fit$mcmc$tau0[i,h]))     }
    ite <- ite + 1
  }
  c( mean(res), as.double(quantile(res, probs = c(0.025, 0.975))))
}

cdf_fit_1 <- matrix(NA,50,3); for(i in 1:50) cdf_fit_1[i,] <- fit.cdf.mcmc(xnow[i], fit.bayes, 8)
data1 <- data.frame(cdf_fit_1, xnow, ecdf1,y.grid)
names(data1)[1:3] <- c("cdf_fit_1","low","upp")

cdf1 <- ggplot(data1) +  
  geom_step(aes(x=y.grid, y=ecdf1)) + 
  geom_line(aes(x=xnow, y=cdf_fit_1), col="blue") +
  geom_ribbon(aes(ymax=upp, ymin=low, x=xnow), fill=4,alpha=.1) + 
  coord_cartesian(xlim=c(25,48)) +  labs(x="Gestational age at delivery (DDE<15)", y="") + 
  theme_bw() + coord_cartesian(ylim=c(0, 1)) +
  theme(plot.margin=unit(c(1,0,0,0),"lines"), axis.title=element_text(size=10))

cdf_fit_2 <- matrix(NA,50,3); for(i in 1:50) cdf_fit_2[i,] <- fit.cdf.mcmc(xnow[i], fit.bayes, 23)
data2 <- data.frame(cdf_fit_2, xnow, ecdf2,y.grid)
names(data2)[1:3] <- c("cdf_fit","low","upp")

cdf2 <- ggplot(data2) +  
  geom_step(aes(x=y.grid, y=ecdf2)) + 
  geom_line(aes(x=xnow, y=cdf_fit), col="blue") +
  geom_ribbon(aes(ymax=upp, ymin=low, x=xnow), fill=4,alpha=.1) + 
  coord_cartesian(xlim=c(25,48)) +  labs(x="Gestational age at delivery (15<DDE<30)", y="") + 
  theme_bw() + coord_cartesian(ylim=c(0, 1)) +
  theme(plot.margin=unit(c(1,0,0,0),"lines"), axis.title=element_text(size=10))

cdf_fit_3 <- matrix(NA,50,3); for(i in 1:50) cdf_fit_3[i,] <- fit.cdf.mcmc(xnow[i], fit.bayes, 38)
data3 <- data.frame(cdf_fit_3, xnow, ecdf3, y.grid)
names(data3)[1:3] <- c("cdf_fit","low","upp")

cdf3 <- ggplot(data3) +  
  geom_step(aes(x=y.grid, y=ecdf3)) + 
  geom_line(aes(x=xnow, y=cdf_fit), col="blue") +
  geom_ribbon(aes(ymax=upp, ymin=low, x=xnow), fill=4,alpha=.1) + 
  coord_cartesian(xlim=c(25,48)) +  labs(x="Gestational age at delivery (30<DDE<45)", y="") + 
  theme_bw() + coord_cartesian(ylim=c(0, 1)) +
  theme(plot.margin=unit(c(1,0,0,0),"lines"), axis.title=element_text(size=10))

cdf_fit_4 <- matrix(NA,50,3); for(i in 1:50) cdf_fit_4[i,] <- fit.cdf.mcmc(xnow[i], fit.bayes, 53)
data4 <- data.frame(cdf_fit_4, xnow, ecdf4, y.grid)
names(data4)[1:3] <- c("cdf_fit","low","upp")

cdf4 <- ggplot(data4) +  
  geom_step(aes(x=y.grid, y=ecdf4)) + 
  geom_line(aes(x=xnow, y=cdf_fit), col="blue") +
  geom_ribbon(aes(ymax=upp, ymin=low, x=xnow), fill=4,alpha=.1) + 
  coord_cartesian(xlim=c(25,48)) +  labs(x="Gestational age at delivery (45<DDE<60)", y="") + 
  theme_bw() + coord_cartesian(ylim=c(0, 1)) +
  theme(plot.margin=unit(c(1,0,0,0),"lines"), axis.title=element_text(size=10))

cdf_fit_5 <- matrix(NA,50,3); for(i in 1:50) cdf_fit_5[i,] <- fit.cdf.mcmc(xnow[i], fit.bayes, 68)
data5 <- data.frame(cdf_fit_5, xnow, ecdf5, y.grid)
names(data5)[1:3] <- c("cdf_fit","low","upp")

cdf5 <- ggplot(data5) +  
  geom_step(aes(x=y.grid, y=ecdf5)) + 
  geom_line(aes(x=xnow, y=cdf_fit), col="blue") +
  geom_ribbon(aes(ymax=upp, ymin=low, x=xnow), fill=4,alpha=.1) + 
  coord_cartesian(xlim=c(25,48)) +  labs(x="Gestational age at delivery (60<DDE<75)", y="") + 
  theme_bw() + coord_cartesian(ylim=c(0, 1)) +
  theme(plot.margin=unit(c(1,0,0,0),"lines"), axis.title=element_text(size=10))

cdf_fit_6 <- matrix(NA,50,3); for(i in 1:50) cdf_fit_6[i,] <- fit.cdf.mcmc(xnow[i], fit.bayes, 126)
data6 <- data.frame(cdf_fit_6, xnow, ecdf6, y.grid)
names(data6)[1:3] <- c("cdf_fit","low","upp")

cdf6 <- ggplot(data6) +  
  geom_step(aes(x=y.grid, y=ecdf6)) + 
  geom_line(aes(x=xnow, y=cdf_fit), col="blue") +
  geom_ribbon(aes(ymax=upp, ymin=low, x=xnow), fill=4,alpha=.1) + 
  coord_cartesian(xlim=c(25,48)) +  labs(x="Gestational age at delivery (DDE>75)", y="") + 
  theme_bw() + coord_cartesian(ylim=c(0, 1)) +
  theme(plot.margin=unit(c(1,0,0,0),"lines"), axis.title=element_text(size=10))

grid.arrange(cdf1, cdf2, cdf3, cdf4, cdf5, cdf6, ncol=3, nrow=2)


#-----
# Figure 6
#-----

#-----
# Histogram of  x
#-----
xdomain <- ggplot(cpp, aes(dde)) + geom_histogram(aes(y = ..density..), fill="grey", alpha=.7) + theme_bw() +
  labs(x="Dichlorodiphenyldichloroethylene (DDE)")
xdomain
#-----
# probability y < threshold (e.g. 37)
#-----
a <- 37
y.grid <- seq(min(y)-sqrt(var(y)), max(y) + sqrt(var(y)), length = 100)	
index <- c((mcmc$nb+1):(mcmc$nrep+mcmc$nb))[1:((mcmc$nrep)/mcmc$thin)*mcmc$thin]
pi_mean <- matrix(NA, length(index), 250)
ite <- 1
for(i in index)
{
  f0 <- approxfun(x=y.grid, y=fit.bayes$mcmc$f0[i,])
  f1 <- approxfun(x=y.grid, y=fit.bayes$mcmc$f1[i,])
  pi_mean[ite,] <- integrate(function(x)f0(x), 25.07, a)$value*(1-fit.bayes$mcmc$beta[i,]) + 
    fit.bayes$mcmc$beta[i,]*integrate(function(x)f1(x), 25.07, a)$value
  ite <- ite + 1
}

# comparison with GAM
library(gam)
gam1 <- gam(y<=a ~ s(x), family="binomial")
p37.gam <- predict(gam1, type="response", newdata = data.frame(x=seq(0,250, length=250)))

comire.data <- data.frame(apply(pi_mean,2,mean), apply(pi_mean,2,quantile, probs=0.025),apply(pi_mean,2,quantile, probs=0.975), 
                          x= c(seq(0,1, length=250)*250), p37.gam)
names(comire.data)[1:3] <- c("mean","low","upp")
comire.plot <- ggplot(comire.data, aes(x)) + 
  geom_line(aes(y=mean), lty=1) + geom_ribbon(aes(ymax=upp, ymin=low), fill=1,alpha=.1) +
  geom_line(aes(y=p37.gam), lty=2) + 
  labs(y=expression("pr( y "<=" 37 )"), x="Dichlorodiphenyldichloroethylene (DDE)")+ 
  theme_bw() + theme(plot.margin=unit(c(1,0,0,0),"lines")) + 
  coord_cartesian(ylim=c(0,.8), xlim=c(0,180)) 
comire.plot

#-----
# beta(x) function
#-----
beta.data <- data.frame(fit.bayes$post.means$beta, t(fit.bayes$ci$beta), x=c(seq(0,250,length=250)))
colnames(beta.data) <- c("mean","low","upp","x")
betaplot <- ggplot(beta.data, aes(x,mean)) + geom_line(lty=1, col=1) + geom_ribbon(aes(ymax=upp, ymin=low), fill=1,alpha=.1) +
  labs(y=expression(beta(x)), x="Dichlorodiphenyldichloroethylene (DDE)")+ theme_bw() + theme(plot.margin=unit(c(1,0,0,0),"lines")) + 
  coord_cartesian(ylim=c(0,.8), xlim=c(0,200)) 
betaplot

#-----
# extra risk function
#-----
a <- 37
index <- c((mcmc$nb+1):(mcmc$nrep+mcmc$nb))[1:((mcmc$nrep)/mcmc$thin)*mcmc$thin]
r_a <- matrix(NA, length(index), 250)
y.grid <- seq(min(y)-sqrt(var(y)), max(y) + sqrt(var(y)), length = 100)
ite <- 1
for(i in index)
{
  f0 <- approxfun(x=y.grid, y=fit.bayes$mcmc$f0[i,])
  f1 <- approxfun(x=y.grid, y=fit.bayes$mcmc$f1[i,])
  r_a[ite,] <- fit.bayes$mcmc$beta[i,]*(integrate(function(x)f1(x), 25.2, a)$value-integrate(function(x)f0(x), 25.1, a)$value )/(1-integrate(function(x)f0(x), 25.1, a)$value)
  ite <- ite + 1
}
risk.data <- data.frame(apply(r_a,2,mean), apply(r_a,2,quantile, probs=0.025), apply(r_a,2,quantile, probs=0.975), seq(0,250, length=250))
colnames(risk.data) <- c("risk","low","upp", "x")

riskplot <- ggplot(risk.data, aes(x,risk)) + geom_line(lty=1) + geom_ribbon(aes(ymax=upp, ymin=low), fill=1,alpha=.1) +
  labs(y=expression(r[37](x)), x="Dichlorodiphenyldichloroethylene (DDE)")+ theme_bw() + 
  theme(plot.margin=unit(c(1,0,0,0),"lines")) + 
  coord_cartesian(ylim=c(0,.8), xlim=c(0,200)) 
riskplot

# benchmark dose
bmd <- function(r, q, x) x[which.min(abs(r-q))]
bmd.mcmc <- apply(r_a, 1, bmd, q=0.1, x=risk.data$x)
quantile(bmd.mcmc, probs = c(0.025, 0.975))
mean(bmd.mcmc)