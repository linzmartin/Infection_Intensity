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
