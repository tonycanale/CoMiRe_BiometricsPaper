#--
# Compute the values for the basis
#--

#----
library("drc")
#----

#----
# Weibull
fit.drm <- drm(gest~dde, data=cpp, fct = weibull2(),type="continuous")
b1<-coefficients(fit.drm)[1]
c1<-coefficients(fit.drm)[2]
d1<-coefficients(fit.drm)[3]
e1<-coefficients(fit.drm)[4]
weibull<-function(x){
	y<-(exp(-exp(b1*(log(x)-log(e1)))))
	return(y)}

#----
#Log-logistic
fit.drm <- drm(gest~dde, data=cpp, fct = llogistic(),type="continuous")
b2<-coefficients(fit.drm)[1]
c2<-coefficients(fit.drm)[2]
d2<-coefficients(fit.drm)[3]
e2<-coefficients(fit.drm)[4]
f2<-coefficients(fit.drm)[5]

log_logi<-function(x){
	y<-1-(1+exp(b2*(log(x)-log(e2))))^(-f2)
	return(y)}

#----
#Quantal linear
fit.drm <- drm(gest~dde, data=cpp, fct = weibull2(fixed=c(1,NA,NA,NA)),type="continuous")
c3<-coefficients(fit.drm)[1]
d3<-coefficients(fit.drm)[2]
e3<-coefficients(fit.drm)[3]
linear<-function(x){
	y<-1-(exp(-exp((log(x)-log(e3)))))
	return(y)}

#----
#Quantal quadratic
fit.drm <- drm(gest~dde, data=cpp, fct = weibull2(fixed=c(2,NA,NA,NA)),type="continuous")
c4<-coefficients(fit.drm)[1]
d4<-coefficients(fit.drm)[2]
e4<-coefficients(fit.drm)[3]
quadratic<-function(x){
	y<-1-(exp(-exp(2*(log(x)-log(e4)))))
	return(y)}

#----
# Quantal cubic
fit.drm <- drm(gest~dde, data=cpp, fct = weibull2(fixed=c(3,NA,NA,NA)),type="continuous")
c5<-coefficients(fit.drm)[1]
d5<-coefficients(fit.drm)[2]
e5<-coefficients(fit.drm)[3]
cubic<-function(x){
	y<-1-(exp(-exp(3*(log(x)-log(e5)))))
	return(y)}

#--
#Log-normal
fit.drm <- drm(gest~dde, data=cpp, fct = lnormal(),type="continuous")
b6<-coefficients(fit.drm)[1]
c6<-coefficients(fit.drm)[2]
d6<-coefficients(fit.drm)[3]
e6<-coefficients(fit.drm)[4]
log_probit<-function(x){
	y<-1-pnorm(b6*(log(x)-log(e6)))
	return(y)}
	
# --
# Create final function merging the J* = 6 basis
# --

basisX <- function(x)
{
  cbind(
    weibull(x),
    log_logi(x),
    linear(x),
    quadratic(x),
    cubic(x),
    log_probit(x)
  )
}