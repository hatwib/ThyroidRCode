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
install.packages("dbscan")
#library(formattable)
rm(list=ls())

setwd("C:/Class/Dissertation/Data/Meka/Thyroid/CSV")
mydata <- read.csv(file = "allData_BorutaSelection.csv", header=TRUE)

tmpdata <- read.csv(file = "data_BSIF555_thyroid.csv", header=TRUE)

mydata <- mydata[,-1] 



mydata <- data.frame(mydata,Class_Label= tmpdata$Class_Label)

mydata$Class_Label <- ifelse(mydata$Class_Label == 1, 'Malignant', 'Benign')
mydata$Class_Label <- as.factor(mydata$Class_Label)
mydata <- mydata[sample(nrow(mydata)),]

###Option for Resolve Class imbalance in Training Set Only
# samples <- sample(NROW(mydata), NROW(mydata) * .632)
# data.train <- mydata[samples, ]
# updata <- smotefamily::ADAS(data.train[,-ncol(data.train)],data.train$Class_Label)$data
# names(updata)[ncol(updata)] <- 'Class_Label'
# updata$Class_Label <- as.factor(updata$Class_Label)
# data.train <- updata
# data.test <- mydata[-samples, ]
###---


###Option for Resolve Class imbalance in both Training and Testing Set
updata <- smotefamily::ADAS(mydata[,-ncol(mydata)],mydata$Class_Label)$data
names(updata)[ncol(updata)] <- 'Class_Label'
updata$Class_Label <- as.factor(updata$Class_Label)
samples <- sample(NROW(updata), NROW(updata) * .632)
data.train <- updata[samples, ]
data.test <- updata[-samples, ]
###---
class(updata$Class_Label)
table(updata$Class_Label)
table(data.train$Class_Label)
table(data.test$Class_Label)

myTrainingControl <- trainControl(method = "cv", 
                                  number = 10, 
                                  savePredictions = TRUE, 
                                  classProbs = TRUE, 
                                  verboseIter = FALSE)


data.model <- train(Class_Label ~., data.train, method = "rf", trControl = myTrainingControl , ntree = 500, importance=TRUE, nodesize = 1)
#data.model$finalModel$err.rate
data.pred <- predict(data.model, data.test, type="prob") # Prediction
data.result <- roc(data.test$Class_Label, data.pred$Malignant)
data.result$auc

coords(data.result, "best", best.method="closest.topleft", ret=c("threshold"))
coords(data.result, "best", best.method="closest.topleft", ret=c("accuracy"))

mymodel.rf <- train(Class_Label ~., data.test, method = "rf", #preProcess = c("scale","center"), 
                                                        trControl = myTrainingControl, 
                                                        ntree = 500, importance=TRUE, nodesize = 1)
pred.rf <- predict(mymodel.rf, data.test, type="prob")
result.rf <- roc(data.test$Class_Label, pred.rf$Malignant)
result.rf$auc


plot(data.result, col = gray(0.4), print.thres="best", print.thres.best.method="closest.topleft", main =  "ROC of Random Forest with Adaptive-SMOTE Sampling Technique")
lines(result.rf, col = "steelblue2")

plot(result.rf, col = gray(0.4), print.thres="best", print.thres.best.method="closest.topleft", main =  "ROC of Random Forest with Adaptive-SMOTE Sampling Technique")
lines(data.result, col = "steelblue2")
legend(x=0.19,y=0.5, legend = c("Testing", "Training"), 
       lty = c(1, 1), lwd = rep(2, 3), col = c(gray(0.4), "steelblue2"), 
       text.width = 0.4, cex = 0.85)


library(ROCR)
data.model <- randomForest(x = data.train[,-ncol(data.train)], y = data.train$Class_Label, ntree = 1500, importance=TRUE, proximity=TRUE)

# pred <- predict(data.model, data.test)
# predictions=as.vector(data.model$votes[,2])
# pred = prediction(predictions,data.test$Class_Label)
pred.test <- predict(data.model, newdata=data.test, type='prob')
data.result.test <- roc(data.test$Class_Label, pred.test[,2])



pred.train <- predict(data.model, newdata=data.train, type='prob')
data.result.train <- roc(data.train$Class_Label, pred.train[,2])

#predictions=as.vector(data.model$votes[,2])
#pred=prediction(predictions,target)


data.result$auc

plot(data.result.test, col = "steelblue2", print.thres="best", print.thres.best.method="closest.topleft", 
     main =  "ROC of Random Forest with Adaptive-SMOTE Sampling Technique")
lines(data.result.train, col = gray(0.4))
legend(x=0.19,y=0.5, legend = c("Testing", "Training"), 
       lty = c(1, 1), lwd = rep(2, 3), col = c(gray(0.4), "steelblue2"), 
       text.width = 0.4, cex = 0.85)

coords(data.result.test, "best", best.method="closest.topleft", ret=c("threshold"))
coords(data.result.test, "best", best.method="closest.topleft", ret=c("accuracy"))

coords(data.result.train, "best", best.method="closest.topleft", ret=c("threshold"))
coords(data.result.train, "best", best.method="closest.topleft", ret=c("accuracy"))


perf_AUC=performance(pred,"auc") #Calculate the AUC value
AUC=perf_AUC@y.values[[1]]

perf_ROC=performance(pred.train,"tpr","fpr") #plot the actual ROC curve
plot(perf_ROC, main="ROC plot")
text(0.5,0.5,paste("AUC = ",format(AUC, digits=5, scientific=FALSE)))