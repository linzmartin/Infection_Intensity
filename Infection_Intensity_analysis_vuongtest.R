#can try subsetting the dataset to make the vuong test faster
set.seed(123)
index<-sample(2,nrow(II_Data),replace = TRUE,p=c(.6,.4))
train<-II_Data[index==1,]
test<-II_Data[index==2,]
summary(test)
summary(train)

z5_train <- zeroinfl(CFUs ~ Temperature*Age | Temperature+Age ## Predictor for the Poisson process
                     , ## Predictor for the Bernoulli process;
                     dist = 'negbin',
                     data = train)

negbinom_train <- glm.nb(CFUs ~ Temperature*Age, data=train,link="log")

summary(z5_train)
summary(negbinom_train)

#compare zero-inflated model with regular negative binomial model using vuong test
gc() #clear memory to make room for this process
#vuong(z5, glmNegBinomial)
vuong(z5_train, negbinom_train)

