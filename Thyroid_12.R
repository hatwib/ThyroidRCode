library(Boruta)
library(glmnet)
library(stringr)
library(randomForest)
library(MASS)
library(geoR)
library(car)
library(caret)
library(DMwR)
library(ROSE)
library(pROC)
library(ROCR)
library(plyr)
library(mlbench)
library(smotefamily)
library(ROCR)
rm(list=ls())

setwd("C:/Class/Dissertation/Data/Meka/Thyroid/CSV")

fileName = matrix("",nrow = 56, ncol = 1)
filterName = matrix("",nrow = 56, ncol = 1)
k <- 1
for(i in seq(5, 17, by = 2)){
  for(j in seq(5, 12, by = 1)){
    fileName[k] = str_c("data_BSIF",i,i,j,"_thyroid.csv",collapse ="")
    filterName[k] = paste0(paste0(paste0("Filter ",i),", Bit Size "),j)
    k = k + 1
  }
}

mydata <- read.csv(file = fileName[1], header=TRUE)

# ROC.results <- matrix(0, nrow = (56), ncol = 3)
# for(f in 1:56){
#   mydata <- read.csv(file = fileName[f], header=TRUE)
# }

mycolors = c('red','blue')



boxplot(split(mydata$BSIF1,mydata$Class_Label))

with(mydata,plot(mydata$BSIF2,mydata$BSIF1,col=mycolors[mydata$Class_Label+1]))

new_mydata <- Boruta(Class_Label~., data = mydata, doTrace = 0)
mydata <- subset(mydata, select = c(c(getSelectedAttributes(new_mydata, withTentative = F)), 'Class_Label'))

updata <- smotefamily::ADAS(mydata[,-ncol(mydata)],mydata$Class_Label)$data
names(updata)[ncol(updata)] <- 'Class_Label'
updata$Class_Label <- as.factor(updata$Class_Label)
mydata <- updata


table(mydata$Class_Label)

with(mydata,plot(mydata$BSIF1,mydata$BSIF2,col=mycolors[mydata$Class_Label]))

min(mydata$BSIF1)
max(mydata$BSIF1)
mean(mydata$BSIF1)
median(mydata$BSIF1)
mode(mydata$BSIF1)

showMethods(min)
neuralnet::neuralnet()
require(neuralnet)

edit(getAnywhere('neuralnet'), file='source_rfcv.r')

getAnywhere('ADAS')
trace("neuralnet",edit=FALSE)
install.packages("installr"); library(installr)
updateR()
install.packages("devtools")
library("devtools")
