#zinb with glmmTMB package



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

ZINBcoef <- as.data.frame(confint(ZINB_simple))
head(ZINBcoef)

###
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


###################################################

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


