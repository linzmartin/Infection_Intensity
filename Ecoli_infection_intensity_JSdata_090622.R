#analysis of Jordyn's infection intensity data
#last updated: 08/24/22
#available on github
###########################################################

####################################
# clear existing workspace
rm(list = ls(all = TRUE))
graphics.off()
shell("cls")

#set wd to your project folder
#setwd("C:/Users/linzm/'OneDrive - Vanderbilt'/Hillyer_Lab/Infection_Intensity") #modfiy this as need be
getwd() #check working directory
###########################

#####################
#load libraries needed:
library(readxl)
library(writexl)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(car) 
library(ggpubr)

library(lattice)
library(MASS)

require(pscl) 
library(lmtest)
library(sandwich)
require(foreign)
require(fitdistrplus)
##############################

#import the data:
II_Data <- read_xlsx("Infection_Intensity_Data_JS_May2022.xlsx",
                     sheet = "R",
                     col_types = c("numeric","numeric","numeric",
                      "numeric","numeric","numeric","numeric"))

str(II_Data)
II_Data <- as.data.frame(II_Data)

#format factors:
#II_Data$Temperature <- as.factor(II_Data$Temperature)
#II_Data$Age <- as.factor(II_Data$Age)
#Temperature and Age are continuous variables
II_Data$Block <- as.factor(II_Data$Block) #categorical

#inspect data:
str(II_Data)
#levels(II_Data$Age)
#levels(II_Data$Temperature)


##########################################################
########### Graphing the Raw Data ########################
##########################################################
#create labels for axes
templabs <- c("27˚C","30˚C","32˚C")
names(templabs)<- c("27","30","32")

agelabs <- c("1 day","5 days", "10 days", "15 days")
names(agelabs)<-c("1","5","10","15")

#Graph 1a: raw data, temp vs. CFUs
library(dplyr)
library(magrittr)
library(plyr)
?count()
II_Data %>%
  count(Temperature,CFUs,Age)%>%
  ggplot(aes(x=Temperature,y=CFUs))+ 
  scale_shape_identity(guide="legend")+
  geom_point()+
  facet_grid(Temperature, labeller = labeller(Age=agelabs, Temperature=templabs))+
  labs(title=~italic("E. coli")~" Infection Intensity",
       x="Temperature (˚C)",parse=TRUE)+
  ylab(expression(~italic("E. coli")~" CFUs"))+theme_bw()+
  theme(legend.position = "none")

II_Data %>%
  count(Temperature,CFUs,Age)%>%
  ggplot(aes(x=Temperature,y=CFUs,color=Age))+ 
  scale_shape_identity(guide="legend")+
  geom_point()+
  facet_grid(~Age, labeller = labeller(Age=agelabs, Temperature=templabs))+
  labs(title=~italic("E. coli")~" Infection Intensity",
       x="Temperature (˚C)",parse=TRUE)+
  ylab(expression(~italic("E. coli")~" CFUs"))+theme_bw()+
  theme(legend.position = "none")

II_Data %>%
  count(Temperature,CFUs,Age)%>%
  ggplot(aes(x=Age,y=CFUs,color=Temperature))+ 
  scale_shape_identity(guide="legend")+
  geom_point()+
  facet_grid(~Temperature, labeller = labeller(Age=agelabs, Temperature=templabs))+
  labs(title=~italic("E. coli")~" Infection Intensity",
       x="Age (days post eclosion)",parse=TRUE)+
  ylab(expression(~italic("E. coli")~" CFUs"))+theme_bw()+
  theme(legend.position = "none")

#Graph 1b: raw data, age vs. CFUs
II_Data %>%
  count(Temperature,CFUs,Age)%>%
  ggplot(aes(x=Age,y=CFUs, color=Temperature))+ 
  scale_shape_identity(guide="legend")+
  geom_point()+
  facet_grid(Temperature~Age, labeller = labeller(Age=agelabs, Temperature=templabs))+
  labs(title=~italic("E. coli")~" Infection Intensity",
       x="Adult Age (days post eclosion)",parse=TRUE)+
  ylab(expression(~italic("E. coli")~" CFUs"))+theme_bw()+
  theme(legend.position = "none")

#Graph 2a: box plot of raw data, temp vs. CFUs
II_Data %>%
  ggplot(aes(x=Temperature,y=CFUs,color=Age))+
  facet_grid(Temperature~Age, labeller = labeller(Age=agelabs, Temperature=templabs))+
  labs(title=~italic("E. coli")~" Infection Intensity",
       x="Temperature (˚C)")+
  ylab(expression(~italic("E. coli")~" CFUs"))+
  #stat_summary(fun=mean,geom="point",shape=23,size=4)+
  geom_boxplot(notch = TRUE)+theme_bw()+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1.3)+
  theme(legend.position = "none")

II_Data %>%
  ggplot(aes(x=Temperature,y=CFUs,color=Age,group=Temperature))+
  facet_grid(~Age, labeller = labeller(Age=agelabs, Temperature=templabs))+
  labs(title=~italic("E. coli")~" Infection Intensity",
       x="Temperature (˚C)")+
  ylab(expression(~italic("E. coli")~" CFUs"))+
  geom_boxplot(notch = FALSE)+theme_bw()+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  theme(legend.position = "none")

#Graph 2b: box plot of raw data, temp vs. CFUs
II_Data %>%
  ggplot(aes(x=Age,y=CFUs,color=Temperature))+
  facet_grid(Temperature~Age, labeller = labeller(Age=agelabs, Temperature=templabs))+
  labs(title=~italic("E. coli")~" Infection Intensity",
       x="Adult Age (days post eclosion)")+
  ylab(expression(~italic("E. coli")~" CFUs"))+
  #stat_summary(fun=mean,geom="point",shape=23,size=4)+
  geom_boxplot(notch = TRUE)+theme_bw()+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1.3)+
  theme(legend.position = "none")

II_Data %>%
  ggplot(aes(x=Age,y=CFUs,color=Temperature))+
  facet_grid(~Temperature, labeller = labeller(Age=agelabs, Temperature=templabs))+
  labs(title=~italic("E. coli")~" Infection Intensity",
       x="Adult Age (days post eclosion)")+
  ylab(expression(~italic("E. coli")~" CFUs"))+
  geom_boxplot(notch = TRUE)+theme_bw()+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  theme(legend.position = "none")

#calculate the mean, median, SE, sd, etc.:
II_Summary <- II_Data %>%
  group_by(Temperature,Age) %>%
  summarise(mean_CFUs = mean(CFUs),
            median_CFUs = median(CFUs),
            sd_CFUs = sd(CFUs),
            n_CFUs = n(),
            SE_CFUs = sd(CFUs)/sqrt(n()))

head(II_Summary)

#arithmetic means graphed:
II_Summary %>%
  count(Temperature,mean_CFUs,Age,SE_CFUs)%>%
  ggplot(aes(x=Temperature,y=mean_CFUs,color=Age,
             ymin=mean_CFUs-SE_CFUs, ymax=mean_CFUs+SE_CFUs))+
  geom_point()+
  facet_grid(Temperature~Age, labeller = labeller(Age=agelabs, Temperature=templabs))+
  labs(title="Mean "~italic("E. coli")~" Infection Intensity",
       x="Temperature (˚C)")+
  ylab(expression(~italic("E. coli")~" CFUs"))+
  geom_errorbar()+ theme_bw()+
  scale_shape(guide="none") + theme(legend.position = "none")


II_Summary %>%
  count(Temperature,mean_CFUs,Age,SE_CFUs)%>%
  ggplot(aes(x=Temperature,y=mean_CFUs,color=Age,
             ymin=mean_CFUs-SE_CFUs, ymax=mean_CFUs+SE_CFUs))+
  geom_point()+
  facet_grid(~Age, labeller = labeller(Age=agelabs, Temperature=templabs))+
  labs(title="Mean "~italic("E. coli")~" Infection Intensity",
       x="Temperature (˚C)")+
  ylab(expression(~italic("E. coli")~" CFUs"))+
  geom_errorbar()+ theme_bw()+
  scale_shape(guide="none") + theme(legend.position = "none")

II_Summary %>%
  count(Temperature,mean_CFUs,Age,SE_CFUs)%>%
  ggplot(aes(x=Age,y=mean_CFUs,color=Temperature,
             ymin=mean_CFUs-SE_CFUs, ymax=mean_CFUs+SE_CFUs))+
  geom_point()+
  facet_grid(~Temperature, labeller = labeller(Age=agelabs, Temperature=templabs))+
  labs(title="Mean "~italic("E. coli")~" Infection Intensity",
       x="Adult Age (days post eclosion)")+
  ylab(expression(~italic("E. coli")~" CFUs"))+
  geom_errorbar()+ theme_bw()+
  scale_shape(guide="none") + theme(legend.position = "none")

#################################

#check for correlations
require(GGally)
ggpairs(II_Data[,c('Temperature','Age','CFUs')])
ggpairs(II_Data[,c('Temperature','Age')])

##################################
#data analysis info:
#Alexander, N., 2012. Review: analysis of parasite and other skewed counts. Tropical Medicine & International Health 17, 684–693.. doi:10.1111/j.1365-3156.2012.02987.x

#we are not log-transforming the data because of the zero values at 27C

#possible models: negative binomial regression, zero-inflated regression, Poisson


####################################################
#Distribution of data and model fitting:
#what distribution does the data follow?
#are the count data overdispersed?

#citation:
#Delignette-Muller, M. L., & Dutang, C. (2015). fitdistrplus: An R Package for Fitting Distributions. Journal of Statistical Software, 64(4), 1–34. https://doi.org/10.18637/jss.v064.i04
#https://data.princeton.edu/wws509/r/overdispersion
#https://fukamilab.github.io/BIO202/04-C-zero-data.html 
#https://stats.oarc.ucla.edu/r/dae/negative-binomial-regression/
#https://www.jstatsoft.org/article/view/v027i08

#plot a histogram of the data:
plot(table(II_Data$CFUs),xlab = "CFUs",ylab="count",data=II_Data)
plot(CFUs ~ Temperature, data = II_Data)
plot(CFUs ~ Age, data = II_Data)
plot(CFUs ~ Infectious_Dose, data = II_Data)

#corrected log equation:
c_log <- function(x) log(x + 0.5) #add 0.5 so values can be log transformed

#take log to visualize:
plot(c_log(CFUs) ~ Temperature, data = II_Data)
plot(c_log(CFUs) ~ Age, data = II_Data)


#plot the empirical distribution and density (CDF)
plotdist(II_Data$CFUs, histo=TRUE,demp=TRUE)

plotdist(c_log(II_Data$CFUs),histo=TRUE,demp=TRUE)

#descriptive statistics of the data and graph
#estimating skewness and kurtosis, using bootstrap samples 
  #(constructed by random sampling with replacement from original data set) 
  # have high variance regardless of bootstrapping

descdist(II_Data$CFUs, boot = 1000)
#skewness: 1.431 left skewed
#kurtosis: 4.26

#what percent of the data is made up of zeros?
100*sum(II_Data$CFUs ==0)/nrow(II_Data) #3.16%


#check mean and variance of CFU data
ggplot(II_Data, aes(x=CFUs)) + 
  geom_histogram(bins = 100)

r <- c(mean(II_Data$CFUs), var(II_Data$CFUs))
c(mean=r[1],var=r[2],ratio=r[2]/r[1]) #ratio is 1.2e6 (greater than 1.1 and considered large / overdispersed)
#Poisson assumption of equal mean and variance does not hold true
#use Neg Binom (mean does not equal variance, which is assumed by Poisson dist.


#######################################################
#assess fits:
#check to see if the data follow a poisson distribution
fit_poisson_CFUs <- fitdist((II_Data$CFUs),"pois")
summary(fit_poisson_CFUs)
plot(fit_poisson_CFUs)

#check to see if the data follow a negative binomial distribution
fit_negbinom_CFUs <- fitdist((II_Data$CFUs), "nbinom")
summary(fit_negbinom_CFUs)
plot(fit_negbinom_CFUs)

#compare the two fits:

par(mfrow = c(2, 2))
#graphics.off()
plot.legend <- c("pois", "nbinom")
denscomp(list(fit_poisson_CFUs,fit_negbinom_CFUs),legendtext = plot.legend)
qqcomp(list(fit_poisson_CFUs,fit_negbinom_CFUs),legendtext = plot.legend)
cdfcomp(list(fit_poisson_CFUs,fit_negbinom_CFUs),legendtext = plot.legend)
ppcomp(list(fit_poisson_CFUs,fit_negbinom_CFUs),legendtext = plot.legend)

#check goodness of fit statistics - chi square for discrete data
gofstat(list(fit_negbinom_CFUs, fit_poisson_CFUs),
        fitnames = c("nbinom", "pois"))


#negative binomial model is a much better fit than poisson
##############

#model fitting:
#fit Poisson, quasiPoisson, and negative binomial distributions with log link (modeling the mean not the actual data using log link)
glmPoisson <- glm(CFUs ~ Temperature*Age, data=II_Data, family = poisson(link = "log"))
glmQuasiPoisson <- glm(CFUs ~ Temperature*Age, data=II_Data, family = quasipoisson(link = "log"))
glmNegBinomial <- glm.nb(CFUs ~ Temperature*Age, data=II_Data,link="log")
glmNegBinomial_additive <- glm.nb(CFUs ~ Temperature+Age, data=II_Data,link="log") #leaves out interaction term

#check models:
summary(glmPoisson) 
#Residual deviance: 342540627  on 692  degrees of freedom
#AIC: 342549916

summary(glmQuasiPoisson)
#Residual deviance: 342540627  on 692  degrees of freedom
#AIC: NA

summary(glmNegBinomial)
#Residual deviance:  867.08  on 692  degrees of freedom
#AIC: 18787
# 2 x log-likelihood:  -18776.9520  
#Theta:  0.5549 
#Std. Err.:  0.0255


summary(glmNegBinomial_additive)
#Residual deviance:  868.84  on 693  degrees of freedom
#AIC: 18801
#Theta:  0.5455 
#Std. Err.:  0.0250
# 2 x log-likelihood:  -18793.1080 


###
#checking for overdispersion in Poisson model:
E2 <- resid(glmPoisson, type = "pearson")
N  <- nrow(II_Data)
p  <- length(coef(glmPoisson))   
sum(E2^2) / (N - p) # =  550206.2 #overdispersion present!
#show it another way:
#how much larger/smaller is the variance relative to the mean for the poisson?
pr <- residuals(glmPoisson,"pearson")
phi <- sum(pr^2)/df.residual(glmPoisson)
round(c(phi,sqrt(phi)),4) #variance is 550206x the mean

qchisq(0.95, df.residual(glmNegBinomial)) #this gives the threshold value
#threshold value for this model is 754.3079
deviance(glmNegBinomial) #deviance = 867
pr <- residuals(glmNegBinomial,"pearson")
sum(pr^2) #pearson chi square value - compare to the threshold, if greater, not signif
#chi square value is less than threshold at 527
#suggests model is a good fit

#for negative binomial model, need to use log likelihood to compare to poisson model
-2*(logLik(glmPoisson)-logLik(glmNegBinomial))
#'log Lik.' 342531131 (df=4) when temp and age are continuous vars.
#'log Lik.' 193430370 (df=12) for categorical predictor vars.
lrtest(glmPoisson,glmNegBinomial)
#chisq=342531131;  pr(>chisq) = < 2.2e-16 ***
#Neg binom is best fit compared to Poisson model

#also should compare diff neg binomial model options (should the interaction term be included in the model?)
glmNegBinomial <- glm.nb(CFUs ~ Temperature*Age, data=II_Data,link="log")
glmNegBinomial_additive <- glm.nb(CFUs ~ Temperature+Age, data=II_Data,link="log") #leaves out interaction term
anova(glmNegBinomial,glmNegBinomial_additive)
#this indicates that the interaction term is a significant predictor of CFUs response 
#(LR stat = 16.1559, pr(chi)=5.833665e-05)

#theta is the precision of the multiplicative random effect = 1/σ2 for the neg binom model
glmNegBinomial$theta #theta = 0.5548731
1/glmNegBinomial$theta #estimated variance = 1.802214

###
# Dispersion statistic for neg binomial model:
#how much larger/smaller is the variance relative to the mean for the Neg Binomial model?
phi <- sum(pr^2)/df.residual(glmNegBinomial)
round(c(phi,sqrt(phi)),4) #variance is 0.7623 times smaller than the mean
#show another way too:
E2 <- resid(glmNegBinomial, type = "pearson")
N  <- nrow(II_Data)
p  <- length(coef(glmNegBinomial)) + 1  # '+1' is for variance in NB model
sum(E2^2) / (N - p) #= 0.763 ==== under-dispersion due to neg binomial model
#may have over-corrected for the zero inflation by using a negative binomial fit --
#check to see if other fits may be better


###
#can use a zero-inflated model or a hurdle model instead to account for excess zeros

#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5628625/#:~:text=Structural%20zeros%20refer%20to%20zero,zero%20due%20to%20sampling%20variability.
#https://lanzaroark.org/docs/zero.pdf
#https://www.jstatsoft.org/article/view/v027i08 
#https://m-clark.github.io/models-by-example/hurdle.html
#https://search.r-project.org/CRAN/refmans/pscl/html/hurdle.html
#https://data.library.virginia.edu/getting-started-with-hurdle-models/
#https://www.r-bloggers.com/2013/05/veterinary-epidemiologic-research-count-and-rate-data-zero-counts/
#https://stats.oarc.ucla.edu/r/dae/negative-binomial-regression/
#https://stats.oarc.ucla.edu/r/dae/zinb/


require(MuMIn)
require(glmmTMB)

#build full zero-inflated neg binom model with glmmTMB
#nbinom2 = negative binomial where the variance increases quadratically with the mean as σ2 = µ(1 + µ/θ), with θ > 0
ZINB_full <- glmmTMB(CFUs ~ Temperature*Age, 
                ziformula = ~ Temperature*Age,
                data=II_Data, family=nbinom2(link='log'), 
                na.action = 'na.fail')

ms <- dredge(ZINB_full, rank='AIC')
#convergence problems with full model

head(ms) #this agrees with my stepwise selection: 
#model should include only age in the ZI component

ZINB_complex <- glmmTMB(CFUs ~ Temperature*Age, 
                        ziformula = ~ Temperature+Age, 
                        #REML=TRUE,
                        data=II_Data, 
                        family=nbinom2(link='log'), 
                        na.action = 'na.fail')

ZINB_nullzi <- glmmTMB(CFUs ~ Temperature*Age, 
                       ziformula = ~ 1, 
                       #REML=TRUE,
                       data=II_Data, 
                       family=nbinom2(link='log'), 
                       na.action = 'na.fail')

ZINB_simple <- glmmTMB(CFUs ~ Temperature*Age, 
                       ziformula = ~ Age, 
                       #REML=TRUE,
                       data=II_Data, 
                       family=nbinom2(link='log'), 
                       na.action = 'na.fail')

lrtest(ZINB_complex,ZINB_simple) #ZINB_simple is better
lrtest(ZINB_simple,ZINB_nullzi) #ZINB_simple is better

summary(ZINB_simple)
AIC(ZINB_simple)
logLik(ZINB_simple)

confint(ZINB_simple) ## Wald/delta-method CIs 


#check model specs:
library(DHARMa)

diagnose(ZINB_simple)
#calculate residuals for model:
#calculates randomized quantile residuals according to the algorithm
simulationOutput <- simulateResiduals(fittedModel = ZINB_simple, plot = F)

# the plot function runs 4 tests
# i) KS test i) Dispersion test iii) Outlier test iv) quantile test
plot(simulationOutput, quantreg = TRUE)

plotResiduals(simulationOutput, form = II_Data$Temperature)
plotResiduals(simulationOutput, form = II_Data$Age)

#simulationOutput = recalculateResiduals(simulationOutput, group = II_Data$Age)


output <- data.frame(resid = resid(ZINB_simple), fitted = fitted(ZINB_simple))
ggplot(output, aes(fitted, resid)) + geom_jitter(position = position_jitter(width = 0.25), 
                                                 alpha = 0.5) + stat_smooth(method = "loess")


ggplot(output, aes(fitted, resid)) +
  geom_jitter(position=position_jitter(width=.25), alpha=.5) +
  stat_quantile(method="rq")

output <- within(output, {
    broken <- cut(fitted, hist(fitted, plot=FALSE)$breaks)
  })

ggplot(output, aes(broken, resid)) +
  geom_boxplot() +
  geom_jitter(alpha=.25)




#################
######
#https://www.jstatsoft.org/article/view/v027i08 source for model estimations
#compare ML estimators for different models: 
model_list <- list("ML-Pois" = glmPoisson, "Quasi-Pois" = glmQuasiPoisson, "NB" = glmNegBinomial,
                   #"Hurdle-NB" = hurdle_model, 
                   "ZINB" = ZINB_simple)

#inspect estimated regression coefficients for each model:
sapply(model_list, function(x) coef(x)[1:8])

#library(sandwich)
#cbind("ML-Pois" = sqrt(diag(vcov(glmPoisson))),
 #     "Adj-Pois" = sqrt(diag(sandwich(glmPoisson))),
  #    sapply(model_list[1], function(x) sqrt(diag(vcov(x)))[1:8]))

rbind(logLik = sapply(model_list, function(x) round(logLik(x), digits = 0)),
      Df = sapply(model_list, function(x) attr(logLik(x), "df")))
#ML-Pois and Quasi-Pois can be eliminated
#better models are NB, Hurdle-NB, and ZINB
#NB is improved by addition of hurdle or zinb (nearly identical)

######

#how are the zeros captured by the models? 
gc() #clear memory to make space for computation
zinb_predictions <- predict(ZINB_simple, type = "zprob") #even more computationally intensive - too large to process in R alone
#=22
str(zinb_predictions)
summary(zinb_predictions)

as.data.frame(zinb_predictions)

observed_zeros <- round(sum(II_Data$CFUs < 1)) 
Pois_zeros <- round(sum(dpois(0, fitted(glmPoisson))))
NegBin_zeros <- round(sum(dnbinom(0, mu = fitted(glmNegBinomial), size = glmNegBinomial$theta)))
zinb_zeros <- round(sum(zinb_predictions[,1]))


comparing_zeros <- c("Obs"=observed_zeros,"ML-Pois"=Pois_zeros,"NB"=NegBin_zeros,"ZINB"=zinb_zeros)
comparing_zeros

#expected mean counts:
mu_zinb <- predict(ZINB_simple, type="response")

#compare AIC values now:
Model_AICs <- c(PoisAIC = AIC(glmPoisson),NegBinAIC=AIC(glmNegBinomial),ZeroInfl_AIC = AIC(ZINB_simple))
Model_AICs



#########################
#plot model predictions:
newdata <- unique(II_Data[,c("Temperature","Age")]) 

predszi <- predict(ZINB_simple, newdata, se.fit=TRUE, type="response") 
predszi <- predict(ZINB_nullzi, newdata, se.fit=TRUE, type="response")
predsnb <- predict(glmNegBinomial,newdata,se.fit=TRUE,type="response")

newdata$predFE = predszi$fit 
newdata$predFE = predsnb$fit

newdata$predFE.min = predszi$fit-1.98*predszi$se.fit 
newdata$predFE.max = predszi$fit+1.98*predszi$se.fit 

newdata$predFE.min = predsnb$fit-1.98*predsnb$se.fit 
newdata$predFE.max = predsnb$fit+1.98*predsnb$se.fit 
#negbinpreds <- predict(glmNegBinomial, newdata, se.fit=TRUE, type="response") 
#newdata$predFE <- negbinpreds$fit


require(plyr)
real <- ddply(II_Data, ~Temperature+Age, summarize, m=mean(CFUs)) 

ggplot(newdata, aes(x=Temperature, y=predFE, colour=factor(Age)))+
  geom_point() + facet_grid(~Age,labeller = labeller(Age=agelabs, Temperature=templabs))+
  geom_errorbar(aes(ymin=predFE.min, ymax=predFE.max)) + 
  geom_point(data=real, aes(x=Temperature, y=m,color="red"))+ 
  geom_point(aes(x=Temperature,y=CFUs),size=(0.5),color="dark gray",data=II_Data)+
  theme_bw()+
  ylab("Average CFUs") + xlab("Temp")



ggplot(newdata1, aes(x = Temperature, y = phat, colour = factor(Age))) +
  geom_point() +
  geom_line() +
  facet_grid(~Age,labeller = labeller(Age=agelabs, Temperature=templabs)) +
  labs(x = "Temperature", y = "Predicted CFUs vs. Observed CFUs")+
  geom_point(aes(x=Temperature,y=CFUs),size=(0.5),color="dark gray",data=II_Data)+
  geom_point(aes(x=Temperature,y=mean_CFUs),color="red",data=II_Summary)+
  theme_bw()






##############



#Neg binomial with log link regression coefficients for each of the variables along 
#with standard errors, z-scores, and p-values 
#A second block follows that corresponds to the inflation model. 
#This includes binomial logit coefficients for predicting excess zeros along 
#with their standard errors, z-scores, and p-values.

#is this model better than the negative binomial model alone (non-hurdle)?
#?vuong() Vuong non-nested test compares models that don't nest with likelihood ratios
#vuong(hurdle_model, glmNegBinomial) #too computationally intensive for R


