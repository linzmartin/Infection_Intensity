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
