#pscl package scripts for hurdle model and zinb
##################

library(pscl)
# count data conditional on the zero hurdle
#count model is usually poisson or neg binomial with log link

#zero hurdle model is binomial 
#get an output for 2 different models: (1) positive counts, (2) zero counts
#outcome of hurdle component is the occurrence of a non-zero (positive) count
#"positive coefficients in the hurdle component indicate that an increase in the regressor increases the probability of a non-zero count" - package & function description
#can include all covariates in each equation or make a simpler model


#left-truncated count component + right-censored hurdle component
#try both a poisson fit and a neg-binom fit
hurdle1 <- hurdle(CFUs ~ Temperature*Age | Temperature*Age, data=II_Data,dist="poisson")
summary(hurdle1)
hurdle2 <- hurdle(CFUs ~ Temperature+Age+Temperature*Age | Temperature+Age+Temperature*Age, data=II_Data,dist="negbin") #logit-negbin
summary(hurdle2)

#compare the likelihood ratios for the hurdle models based on poisson and neg-binom fits
lrtest(hurdle1,hurdle2) #hurdle 2 better

#use negative binomial fit because log likelihood decreases compared to poisson

#does the interaction term matter?
hurdle4 <- hurdle(CFUs ~ Temperature*Age | Temperature+Age, data=II_Data,dist="negbin")
summary(hurdle4)
lrtest(hurdle2,hurdle4) #no significant diff - remove the interaction term

hurdle3 <- hurdle(CFUs ~ Temperature+Age | Temperature+Age, data=II_Data,dist="negbin")
summary(hurdle3)
lrtest(hurdle3,hurdle4) #there is a significant difference - need to include the interaction term in the count model
#so far, hurdle 4 is best

#for the hurdle component, check the importance of the main effects:
hurdle6 <- hurdle(CFUs ~ Temperature*Age | Temperature, data=II_Data,dist="negbin")
summary(hurdle6)
lrtest(hurdle4,hurdle6) #hurdle 4 is better

hurdle7 <- hurdle(CFUs ~ Temperature*Age | Age, data=II_Data,dist="negbin")
summary(hurdle7)
lrtest(hurdle4,hurdle7) #ns - but LogLik does slightly increase with removal of Temp

#hurdle5 <- hurdle(CFUs ~ Temperature*Age | 1, data=II_Data,dist="negbin")
#summary(hurdle5)

lrtest(hurdle7,hurdle5) #including age is better than the intercept alone

#going with hurdle4 because of increase in logLik but lack of significance -- rename to hurdle_model:
hurdle_model <- hurdle4
summary(hurdle_model)

###
#try zero-inflated model instead: 
#can include all covariates in each equation or make a simpler model

z1 <- zeroinfl(CFUs ~ Temperature*Age | ## Predictor for the Poisson process
                 Temperature*Age, ## Predictor for the Bernoulli process;
               dist = 'poisson',
               data = II_Data)
#z1 leads to overfitting

z2 <- zeroinfl(CFUs ~ Temperature*Age | ## Predictor for the Poisson process
                 1, ## Predictor for the Bernoulli process;
               dist = 'poisson',
               data = II_Data) #simple model
#z2 is the null poisson zero-inflated model
summary(z2)

#z3 is the null neg bin. zero-inflated model
z3 <- zeroinfl(CFUs ~ Temperature*Age | ## Predictor for the Poisson process
                 1, ## Predictor for the Bernoulli process;
               dist = 'negbin',
               data = II_Data)
summary(z3)

z4 <- zeroinfl(CFUs ~ Temperature*Age | Temperature*Age ## Predictor for the Poisson process
               , ## Predictor for the Bernoulli process;
               dist = 'negbin',
               data = II_Data)
summary(z4)
#z4 leads to overfitting; Std errors of ZI model are too large

z5 <- zeroinfl(CFUs ~ Temperature*Age | Temperature+Age ## Predictor for the Poisson process
               , ## Predictor for the Bernoulli process;
               dist = 'negbin',
               data = II_Data)
summary(z5)
#z5 leads to overfitting; Std errors of ZI model  coeff for Intercept and Age are too large

z6 <- zeroinfl(CFUs ~ Temperature*Age | Temperature ## Predictor for the Poisson process
               , ## Predictor for the Bernoulli process;
               dist = 'negbin',
               data = II_Data)
summary(z6)

z7 <- zeroinfl(CFUs ~ Temperature*Age | Age ## Predictor for the Poisson process
               , ## Predictor for the Bernoulli process;
               dist = 'negbin',
               data = II_Data)
summary(z7)

#z7 leads to overfitting; Std errors of ZI model Int and Age are too large

#are predictors important to the model?
lrtest(z3,z4) #yes, the complex model is better (z4)

#do we need to include the interaction term for the second part of the model?
lrtest(z4,z5) #NS - go with z5 - do not include interaction in 2nd half
lrtest(z5,z6) #include age in 2nd half
lrtest(z5,z7) #NS - can exclude temperature

#null model only for both parts:
z0 <- update(z5, . ~ 1)
summary(z0)

lrtest(z0,z5)
lrtest(z0,z7) #complex model is better than null model


z8 <- zeroinfl(CFUs ~ Temperature+Age | Age ## Predictor for the Poisson process
               , ## Predictor for the Bernoulli process;
               dist = 'negbin',
               data = II_Data)
summary(z8) #overfitting

z9 <- zeroinfl(CFUs ~ Temperature+Age | Temperature+Age ## Predictor for the Poisson process
               , ## Predictor for the Bernoulli process;
               dist = 'negbin',
               data = II_Data) #overfitting
summary(z9)

lrtest(z5,z8) 
lrtest(z7,z8) 
lrtest(z9,z5) #better to include interaction in count model

############
#try alternate methods for model selection:
require(pscl)
require(mpath)

#chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://cran.r-project.org/web/packages/mpath/vignettes/static_german.pdf

#stepwise backward selection:
#be.zeroinfl()
fitbe_z4 <- be.zeroinfl(z4, data=II_Data, dist="negbin", alpha=0.01, trace=FALSE) 
summary(fitbe_z4)

cat("loglikofzero-inflatedmodelwithbackwardselection",logLik(fitbe_z4))
cat("BICofzero-inflatedmodelwithbackwardselection",AIC(fitbe_z4,k=log(dim(II_Data)[1])))
fitbe_z4$formula


################

#decide between the hurdle model and zero-inflated negative binomial model:
#Hurdle: the zero hurdle component = the probability of observing a positive count
#hurdle model = two sources of zeros
#ZINB model: the zero-inflation component = probability of observing a zero count from the point mass component
#one source of zeros
AIC(z5, hurdle_model)
BIC(z5, hurdle_model)

summary(z5)
summary(hurdle_model)


fm <- list("ZINB" = z7, "Hurdle-NB" = hurdle_model)
t(sapply(fm[1:2], function(x) round(x$coefficients$zero, digits = 3)))
#the coefficients for the zero component are opposite signs (makes sense because models predict probability of positive/negative)
t(sapply(fm[1:2], function(x) round(exp(x$coefficients$zero), digits = 3)))

#install.packages("countreg", repos="http://R-Forge.R-project.org")
library(countreg)
graphics.off()
#rootogram(z3,max=15000)
#rootogram(glmNegBinomial,max=10000)

rootogram(z7, main = "Zero-Inflated Negative Binomial", ylim = c(-5, 5), max = 10000)
rootogram(hurdle_model, main = "Negative Binomial Hurdle", ylim = c(-10, 5), max = 10000)
qqrplot(z7, main = "Zero-Inflated Negative Binomial") #very intensive computationally
#qqrplot(hurdle_model, main = "Negative Binomial Hurdle") #intensive computationally

gc() #clear memory in R first
#hurdle_predictions <- predict(hurdle_model,type="prob") #this is computationally intensive and may take awhile

NB_Hurdle_zeros <- round(sum(hurdle_predictions[,1]))
NB_Hurdle_ones <- round(sum(hurdle_predictions[,2]))

#expected mean counts:
mu_hurdle <- predict(hurdle_model,type="response")

plot(hurdle_model$fitted.values~hurdle_model$residuals)

hurdle_model$coefficients


predhudlenb <-predict(hurdle_model,newdata=II_Data[II_Data$CFUs!=0,],type="response")
plot(II_Data$CFUs[II_Data$CFUs!=0],type="b",col="red") #actual
lines(round(predhudlenb),col="blue") #predicted values


#https://www.r-bloggers.com/2020/01/count-data-models/
require("ModelMetrics")

predhnb<- predict(hurdle_model,newdata=II_Data,type = "response")
plot((II_Data$CFUs),type="p")
lines(round(predhnb),col="blue")
str(predhnb)

rmsemodelhnb<-ModelMetrics::rmse(II_Data$CFUs,round(predhnb)) #root mean square error
maemodelhnb<-mae(II_Data$CFUs,round(predhnb)) #mean absolute error
rmsemodelhnb
maemodelhnb

############################

####
#https://stats.oarc.ucla.edu/r/dae/zinb/
library(boot)

#get start values for each part of the model:
dput(round(coef(z7, "count"), 4))
dput(round(coef(z7, "zero"), 4))

f <- function(data, i) {
  require(pscl)
  z7 <- zeroinfl(CFUs ~ Temperature*Age | Age,
                 dist = 'negbin',
                 data = data[i,],
                start = list(count = c(-4.3642, 0.5336, -0.5563,0.024), 
                             zero = c(2.453, -4.5153)))
  as.vector(t(do.call(rbind, coef(summary(z7)))[, 1:2]))
}

set.seed(10)
(res <- boot(II_Data, f, R = 1200, parallel = "snow", ncpus = 4))

## basic parameter estimates with percentile and bias adjusted CIs
parms <- t(sapply(c(1, 3, 5, 9, 11,13), function(i) {
  out <- boot.ci(res, index = c(i, i + 1), type = c("perc", "bca"))
  with(out, c(Est = t0, pLL = percent[4], pUL = percent[5],
              bcaLL = bca[4], bcaUL = bca[5]))
}))

## add row names
as.data.frame(parms)
row.names(parms) <- names(coef(z7))

## print results
parms

## compare with normal based approximation
confint(z7)
# bootstrapped errors should be larger

###
#get incident risk ratio (IRR) for negative binomial model or odds ratio (OR) for logit model (zero inflation model)

## exponentiated parameter estimates with percentile and bias adjusted CIs
expparms <- t(sapply(c(1, 3, 5, 7, 9,11), function(i) {
  out <- boot.ci(res, index = c(i, i + 1), type = c("perc", "bca"), h = exp)
  with(out, c(Est = t0, pLL = percent[4], pUL = percent[5],
              bcaLL = bca[4], bcaUL = bca[5]))
}))


## add row names
as.data.frame(expparms)
row.names(expparms) <- names(coef(z7))
## print results
expparms


newdata1 <- expand.grid(c(27,30,32), rep(c(1,5,10,15)))
#newdata1 <- expand.grid(factor(27:32), rep(c(1,5,10,15)))

colnames(newdata1) <- c("Temperature", "Age")
newdata1$phat <- predict(z7, newdata1)
newdata1$phatlog <- log(newdata1$phat)

ggplot(newdata1, aes(x = Temperature, y = phat, colour = factor(Age))) +
  geom_point() +
  geom_line() +
  facet_grid(~Age,labeller = labeller(Age=agelabs, Temperature=templabs)) +
  labs(x = "Temperature", y = "Predicted CFUs vs. Observed CFUs")+
  geom_point(aes(x=Temperature,y=CFUs),size=(0.5),color="dark gray",data=II_Data)+
  geom_point(aes(x=Temperature,y=mean_CFUs),color="red",data=II_Summary)+theme_bw()


########
#check the dispersion statistic:
E2 <- resid(z7, type = "pearson")
N  <- nrow(II_Data)
p  <- length(coef(z5)) + 1  # '+1' is for variance in NB model
sum(E2^2) / (N - p) 
# = 1.005472, close to 1 = good
#closer to 1 than the negative binomial model


#############


Zero_inflated_neg_binom_model <- z7 #rename for ease
summary(Zero_inflated_neg_binom_model)

plot(density(II_Data$CFUs))

plot(Zero_inflated_neg_binom_model$fitted.values~Zero_inflated_neg_binom_model$residuals)
Zero_inflated_neg_binom_model$coefficients

predzinb <- predict(Zero_inflated_neg_binom_model,
                    newdata=II_Data[II_Data$CFUs!=0,],
                    type="response")
plot(II_Data$CFUs[II_Data$CFUs!=0],type="b",col="red")
lines(round(predzinb),col="blue") #predicted values
