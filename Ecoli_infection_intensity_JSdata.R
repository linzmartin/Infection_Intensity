#analysis of Jordyn's infection intensity data

#05/20/22

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
  ylab(expression(~italic("E. coli")~" CFUs"))+
  theme(legend.position = "none") +theme_bw()

#Graph 1b: raw data, age vs. CFUs
II_Data %>%
  count(Temperature,CFUs,Age)%>%
  ggplot(aes(x=Age,y=CFUs,color=Temperature))+ 
  scale_shape_identity(guide="legend")+
  geom_point()+
  facet_grid(Temperature~Age, labeller = labeller(Age=agelabs, Temperature=templabs))+
  labs(title=~italic("E. coli")~" Infection Intensity",
       x="Adult Age (days post eclosion)",parse=TRUE)+
  ylab(expression(~italic("E. coli")~" CFUs"))+
  theme(legend.position = "none") +theme_bw()

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


##################################
#data analysis info:
#Alexander, N., 2012. Review: analysis of parasite and other skewed counts. Tropical Medicine & International Health 17, 684–693.. doi:10.1111/j.1365-3156.2012.02987.x

#we are not log-transforming the data because of the zero values

#possible models: negative binomial regression, zero-inflated regression, Poisson


###############

#negative binomial distribution for count data:
require(foreign)
require(ggplot2)
require(MASS)


require(fitdistrplus)
?fitdist()
# log normal distribution - can only be used on positive values (any zeros and it will fail)

#fit_lognormal <- fitdist((II_Data$CFUs), "lnorm")
#fit_gamma <- fitdist((II_Data$CFUs), "gamma")

fit_poisson_CFUs <- fitdist((II_Data$CFUs),"pois")
summary(fit_poisson_CFUs)
plot(fit_poisson_CFUs)

fit_negbinom_CFUs <- fitdist((II_Data$CFUs), "nbinom")
summary(fit_negbinom_CFUs)
plot(fit_negbinom_CFUs)

#compare the two fits:
cdfcomp(list(fit_poisson_CFUs,fit_negbinom_CFUs)) 
#negative binomial model is a much better fit than poisson
