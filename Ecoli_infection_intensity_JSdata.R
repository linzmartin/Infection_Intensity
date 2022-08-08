#analysis of Jordyn's infection intensity data
#last updated: 08/08/22
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
                     sheet = "R")
                     #col_types = c()) #import all PO Data

#format factors:
II_Data$Temperature <- as.factor(II_Data$Temperature)
II_Data$Age <- as.factor(II_Data$Age)
II_Data$Block <- as.factor(II_Data$Block)

str(II_Data)

##########################################################
########### Graphing the Raw Data ########################
##########################################################
#create labels for axes
templabs <- c("27˚C","30˚C","32˚C")
names(templabs)<- c("27","30","32")

agelabs <- c("1 day","5 days", "10 days", "15 days")
names(agelabs)<-c("1","5","10","15")

#Graph 1a: raw data, temp vs. CFUs
II_Data %>%
  count(Temperature,CFUs,Age)%>%
  ggplot(aes(x=Temperature,y=CFUs,color=Age))+ 
  scale_shape_identity(guide="legend")+
  geom_point()+
  facet_grid(Temperature~Age, labeller = labeller(Age=agelabs, Temperature=templabs))+
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
  ggplot(aes(x=Temperature,y=CFUs,color=Age))+
  facet_grid(~Age, labeller = labeller(Age=agelabs, Temperature=templabs))+
  labs(title=~italic("E. coli")~" Infection Intensity",
       x="Temperature (˚C)")+
  ylab(expression(~italic("E. coli")~" CFUs"))+
  geom_boxplot(notch = TRUE)+theme_bw()+
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

head(II_Summary)

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
plot(table(II_Data$CFUs))
plot(CFUs ~ Temperature, data = II_Data)
plot(CFUs ~ Age, data = II_Data)
plot(CFUs ~ Infectious_Dose, data = II_Data)

#corrected log equation:
c_log <- function(x) log(x + 0.5) #add 0.5 so values can be log transformed

#take log to visualize:
plot(c_log(CFUs) ~ Temperature, data = II_Data)
plot(c_log(CFUs) ~ Age, data = II_Data)

###########
#plot the empirical distribution and density (CDF)
plotdist(II_Data$CFUs, histo=TRUE,demp=TRUE)

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
#Residual deviance: 193439731  on 684  degrees of freedom
#AIC: 193449036

summary(glmQuasiPoisson)
#Residual deviance: 193439731  on 684  degrees of freedom
#AIC: NA

summary(glmNegBinomial)
#Residual deviance:  853.02  on 684  degrees of freedom
#AIC: 18668
# 2 x log-likelihood:  -18641.9380 

summary(glmNegBinomial_additive)
#Residual deviance:  864.47  on 690  degrees of freedom
#AIC: 18767
# 2 x log-likelihood:  -18752.8340

###
#checking for overdispersion in Poisson model:
E2 <- resid(glmPoisson, type = "pearson")
N  <- nrow(II_Data)
p  <- length(coef(glmPoisson))   
sum(E2^2) / (N - p) # =  291737.3 #overdispersion present!
#show it another way:
#how much larger/smaller is the variance relative to the mean for the poisson?
pr <- residuals(glmPoisson,"pearson")
phi <- sum(pr^2)/df.residual(glmPoisson)
round(c(phi,sqrt(phi)),4) #variance is 291737x the mean

qchisq(0.95, df.residual(glmNegBinomial)) #this gives the threshold value
#threshold value for this model is 745.953
deviance(glmNegBinomial) #deviance = 853
pr <- residuals(glmNegBinomial,"pearson")
sum(pr^2) #pearson chi square value - compare to the threshold, if greater, not signif
#chi square value is less than threshold at 482
#suggests model is a good fit

#for negative binomial model, need to use log likelihood to compare to poisson model
-2*(logLik(glmPoisson)-logLik(glmNegBinomial))
#'log Lik.' 193430370 (df=12)
lrtest(glmPoisson,glmNegBinomial)
#chisq=193430370;  pr(>chisq) = < 2.2e-16 ***
#Neg binom is best fit

#also should compare diff neg binomial model options (should the interaction term be included in the model?)
glmNegBinomial <- glm.nb(CFUs ~ Temperature*Age, data=II_Data,link="log")
glmNegBinomial_additive <- glm.nb(CFUs ~ Temperature+Age, data=II_Data,link="log") #leaves out interaction term
anova(glmNegBinomial,glmNegBinomial_additive)
#this indicates that the interaction term is a significant predictor of CFUs response 
#(LR stat = 110.8968, pr(chi)=0)

#theta is the precision of the multiplicative random effect = 1/σ2 for the neg binom model
glmNegBinomial$theta #theta = 0.6415172
1/glmNegBinomial$theta #estimated variance = 1.558805

###
# Dispersion statistic for neg binomial model:
#how much larger/smaller is the variance relative to the mean for the Neg Binomial model?
phi <- sum(pr^2)/df.residual(glmNegBinomial)
round(c(phi,sqrt(phi)),4) #variance is 0.7048 times smaller than the mean

E2 <- resid(glmNegBinomial, type = "pearson")
N  <- nrow(II_Data)
p  <- length(coef(glmNegBinomial)) + 1  # '+1' is for variance in NB model
sum(E2^2) / (N - p) #= 0.705 ==== under-dispersion due to neg binomial model

###
#can use a zero-inflated model or a hurdle model instead to account for excess zeros
#zero-inflated model has two sources of zeros: sampling zeros and structural zeros
#hurdle model only has 1 source: sampling zeros (there are no structural zeros)
#we will try a hurdle model - the zero values that occur are due to sampling (and could be non-zero, so they are not structural zeros)

#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5628625/#:~:text=Structural%20zeros%20refer%20to%20zero,zero%20due%20to%20sampling%20variability.
#https://lanzaroark.org/docs/zero.pdf
#https://www.jstatsoft.org/article/view/v027i08 
#https://m-clark.github.io/models-by-example/hurdle.html
#https://search.r-project.org/CRAN/refmans/pscl/html/hurdle.html
#https://data.library.virginia.edu/getting-started-with-hurdle-models/

library(pscl)
# count data conditional on the zero hurdle
#count model is usually poisson or neg binomial with log link
#zero hurdle model is binomial 
#get an output for 2 different models: (1) positive counts, (2) zero counts
#outcome of hurdle component is the occurrence of a non-zero (positive) count
#"positive coefficients in the hurdle component indicate that an increase in the regressor increases the probability of a non-zero count" - package & function description
#can include all covariates in each equation or make a simpler model

hurdle1 <- hurdle(CFUs ~ Temperature*Age | Temperature*Age, data=II_Data,dist="poisson")
summary(hurdle1)

hurdle2 <- hurdle(CFUs ~ Temperature+Age+Temperature*Age | Temperature+Age+Temperature*Age, data=II_Data,dist="negbin") #logit-negbin
summary(hurdle2)

#use negative binomial fit because log likelihood decreases compared to poisson

#does the interaction term matter?
hurdle3 <- hurdle(CFUs ~ Temperature+Age | Temperature+Age, data=II_Data,dist="negbin")
summary(hurdle3)

lrtest(hurdle3,hurdle2) #significantly different - need to include interaction term

#going with hurdle2 -- rename to hurdle_model:
hurdle_model <- hurdle2

hurdle4 <- hurdle(CFUs ~ Temperature*Age | Temperature+Age, data=II_Data,dist="negbin")
summary(hurdle4)
lrtest(hurdle2,hurdle4) #no significant diff - remove the interaction term
#go with hurdle4


plot(factor(CFUs == 0) ~ Temperature*Age, data = II_Data, main = "Zero")
plot(log(sightings) ~ depth, data = hdata1, subset = sightings > 0,
     main = "Count")

###
#try zero-inflated model instead: 
#can include all covariates in each equation or make a simpler model

z1 <- zeroinfl(CFUs ~ Temperature*Age | ## Predictor for the Poisson process
                Temperature*Age, ## Predictor for the Bernoulli process;
               dist = 'poisson',
               data = II_Data)
z2 <- zeroinfl(CFUs ~ Temperature*Age | ## Predictor for the Poisson process
                 1, ## Predictor for the Bernoulli process;
               dist = 'poisson',
               data = II_Data) #simple model

z3 <- zeroinfl(CFUs ~ Temperature*Age | ## Predictor for the Poisson process
                       1, ## Predictor for the Bernoulli process;
                     dist = 'negbin',
                     data = II_Data)
z4 <- zeroinfl(CFUs ~ Temperature*Age | Temperature*Age ## Predictor for the Poisson process
                , ## Predictor for the Bernoulli process;
               dist = 'negbin',
               data = II_Data)

summary(z1) #not applicable
summary(z2) #poisson with log link
summary(z3)
summary(z4) #go with z4

pchisq(2 * (logLik(z4) - logLik(z3)), df = 3, lower.tail = FALSE)
lrtest(z4,z3)


######
#https://www.jstatsoft.org/article/view/v027i08 source for model estimations
#compare ML estimators for different models: 
model_list <- list("ML-Pois" = glmPoisson, "Quasi-Pois" = glmQuasiPoisson, "NB" = glmNegBinomial,
                   "Hurdle-NB" = hurdle_model, "ZINB" = z4)

#inspect estimated regression coefficients for each model:
sapply(model_list, function(x) coef(x)[1:8])

library(sandwich)
cbind("ML-Pois" = sqrt(diag(vcov(glmPoisson))),
      "Adj-Pois" = sqrt(diag(sandwich(glmPoisson))),
      sapply(model_list[1], function(x) sqrt(diag(vcov(x)))[1:8]))

rbind(logLik = sapply(model_list, function(x) round(logLik(x), digits = 0)),
      Df = sapply(model_list, function(x) attr(logLik(x), "df")))
#ML-Pois and Quasi-Pois can be eliminated
#better models are NB, Hurdle-NB, and ZINB
#NB is improved by addition of hurdle or zinb (nearly identical)
#go with hurdle model because zeros are due to sampling, not structural zeros


#how are the zeros captured by the models? 

hurdle_predictions <- predict(hurdle_model,type="prob") #this is computationally intensive and may take awhile
#zinb_predictions <- predict(z4, type = "prob") #even more computationally intensive - too large to process in R alone

observed_zeros <- round(sum(II_Data$CFUs < 1)) 
Pois_zeros <- round(sum(dpois(0, fitted(glmPoisson))))
NegBin_zeros <- round(sum(dnbinom(0, mu = fitted(glmNegBinomial), size = glmNegBinomial$theta)))
NB_Hurdle_zeros <- round(sum(hurdle_predictions[,1]))

comparing_zeros <- c("Obs"=observed_zeros,"ML-Pois"=Pois_zeros,"NB"=NegBin_zeros,"NB-Hurdle"=NB_Hurdle_zeros)
comparing_zeros

#we can see that the hurdle model prediction matches the observed # of zeros

#compare AIC values now:
Model_AICs <- c(PoisAIC = AIC(glmPoisson),NegBinAIC=AIC(glmNegBinomial), NB_HurdleAIC = AIC(hurdle_model),ZeroInfl_AIC = AIC(z4))
Model_AICs

plot(hurdle_model)
################################
##################
#don't run this#
"""
#can also use family = neg binomial in glm, but standard errors differ due to theta estimation vs. need to provide theta
mnbg <- glm(CFUs ~ Temperature*Age,
            family=negative.binomial(glmNegBinomial$theta), data=II_Data)
all(abs(coef(glmNegBinomial)==coef(mnbg)) < 1e-12)

deviance(mnbg)
"""
###################

#####
#data exploration:

head(II_Data)

#convert CFU data into z-scores - standardization to mean zero, StDev = 1
#https://www.r-bloggers.com/2021/12/how-to-use-the-scale-function-in-r/#:~:text=Scale()%20is%20a%20built,center%E2%80%9D%20value%20from%20the%20argument.&text=center%3A%20When%20scaling%2C%20whether%20the%20mean%20should%20be%20subtracted.
II_Data$scaled_CFUs <- scale(II_Data$CFUs)

#plot z-scores
II_Data %>%
  ggplot(aes(x=Age,y=scaled_CFUs,color=Temperature))+
  facet_grid(~Temperature, labeller = labeller(Age=agelabs, Temperature=templabs))+
  labs(title=~italic("E. coli")~" Infection Intensity, normalized",
       x="Adult Age (days post eclosion)")+
  ylab(expression(~italic("E. coli")~" CFUs normalized"))+
  geom_boxplot(notch = TRUE)+theme_bw()+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  theme(legend.position = "none")


#marginal means: 
#https://easystats.github.io/modelbased/articles/estimate_means.html

#raw means may be biased due to different number of observations in each group
require(modelbased)

model_mm_1 <- lm(CFUs ~ Temperature*Age, data = II_Data)
Anova(model_mm_1)
means1 <- estimate_means(model_mm_1)
means1

II_Data %>%
  ggplot(aes(x=Age,y=CFUs,color=Temperature))+
  facet_grid(~Temperature, labeller = labeller(Age=agelabs, Temperature=templabs))+
  geom_violin()+
  geom_jitter(width=0.2,alpha=0.2)+
  geom_pointrange(
    data=means1, 
    aes(y=Mean, ymin=CI_low,ymax=CI_high), 
    size=0.5,color="black")+
  labs(title="Marginal Means of"~italic("E. coli")~" Infection Intensity",
       x="Adult Age (days post eclosion)")+
  ylab(expression("Marginal Means of"~italic("E. coli")~" CFUs"))+
  theme_bw()+ theme(legend.position = "none")

#find marginal means on normalized z-scores of CFUs and plot:
model_mm_2 <- lm(scaled_CFUs ~ Temperature*Age, data = II_Data)
Anova(model_mm_2)
means2 <- estimate_means(model_mm_2)
means2

II_Data %>%
  ggplot(aes(x=Age,y=scaled_CFUs,color=Temperature))+
  facet_grid(~Temperature, labeller = labeller(Age=agelabs, Temperature=templabs))+
  geom_violin()+
  geom_jitter(width=0.2,alpha=0.2)+
  geom_pointrange(
    data=means2, 
    aes(y=Mean, ymin=CI_low,ymax=CI_high), 
    size=0.5,color="black")+
  labs(title="Marginal Means of Normalized"~italic(" E. coli")~" Infection Intensity",
       x="Adult Age (days post eclosion)")+
  ylab(expression("Marginal Means of"~italic("E. coli")~" CFUs"))+
  theme_bw()+ theme(legend.position = "none")

#add infectious dose data as a model factor:
require(lme4)
model_mm_3 <- lmer(CFUs ~ Temperature*Age + (1|Infectious_Dose), REML = FALSE, data = II_Data)
Anova(model_mm_3)
model_mm_3
means3 <- estimate_means(model_mm_3)
means3

?anova()
anova(model_mm_1,model_mm_3)
AIC(model_mm_1,model_mm_3)

II_Data %>%
  ggplot(aes(x=Age,y=CFUs,color=Temperature))+
  facet_grid(~Temperature, labeller = labeller(Age=agelabs, Temperature=templabs))+
  geom_violin()+
  geom_jitter(width=0.2,alpha=0.2)+
  geom_pointrange(
    data=means3, 
    aes(y=Mean, ymin=CI_low,ymax=CI_high), 
    size=0.5,color="black")+
  labs(title="Marginal Means of"~italic("E. coli")~" Infection Intensity",
       x="Adult Age (days post eclosion)")+
  ylab(expression("Marginal Means of"~italic("E. coli")~" CFUs"))+
  theme_bw()+ theme(legend.position = "none")
############

#plot both means1 and means3 to see difference of including infectious dose as random effect:
II_Data %>%
  ggplot(aes(x=Age,y=CFUs,color=Temperature))+
  facet_grid(~Temperature, labeller = labeller(Age=agelabs, Temperature=templabs))+
  geom_violin()+
  geom_jitter(width=0.2,alpha=0.2)+
  geom_pointrange(
    data=means1, 
    aes(y=Mean, ymin=CI_low,ymax=CI_high), 
    size=0.5,color="black")+
  geom_pointrange(
    data = means3,
    aes(y = Mean, ymin = CI_low, ymax = CI_high),
    size = 0.5,
    color = "orange"
  )+
  labs(title="Marginal Means of"~italic("E. coli")~" Infection Intensity",
       x="Adult Age (days post eclosion)",
       subtitle="Black dots: lm(CFUs ~ Temperature*Age)\nOrange dots: lmer(CFUs ~ Temperature*Age + (1|Infectious_Dose))")+
  ylab(expression("Marginal Means of"~italic("E. coli")~" CFUs"))+
  theme_bw() +theme(legend.position = "none")
