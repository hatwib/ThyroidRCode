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
# mydata <- read.csv(file = "allData_BorutaSelection.csv", header=TRUE)
#myTrainingControl <- trainControl(method = "cv", number = 10, savePredictions = TRUE, classProbs = TRUE, verboseIter = TRUE)
#myTrainingControl <- trainControl(method="repeatedcv", number=10, repeats=3, savePredictions=TRUE, classProbs=TRUE, verbose=TRUE)
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


progress.bar <- create_progress_bar("text")
progress.bar$init(56)
ROC.results <- matrix(0, nrow = (56), ncol = 3)
f <- 1
for(f in 1:56){
  mydata <- read.csv(file = fileName[f], header=TRUE)


  mydata$Class_Label <- ifelse(mydata$Class_Label == 1, 'Malignant', 'Benign')
  mydata$Class_Label <- as.factor(mydata$Class_Label)
  mydata <- mydata[sample(nrow(mydata)),]

  n_cols = ncol(mydata)
  for (h in 1:(n_cols-1)){
    mydata[,h] = (mydata[,h]-min(mydata[,h]))/(max(mydata[,h])-min(mydata[,h]))
  }

  rm_cols <- nearZeroVar(mydata)
  rm_idx <- length(rm_cols)


  if(rm_idx > 0)
    mydata <- mydata[,-c(rm_cols)]

  n_cols = ncol(mydata)
  if(n_cols <= 3)
    next

  new_mydata <- Boruta(Class_Label~., data = mydata, doTrace = 0)
  mydata <- subset(mydata, select = c(c(getSelectedAttributes(new_mydata, withTentative = F)), 'Class_Label'))
  n_cols = ncol(mydata)
  if(n_cols <= 3)
    next
# ###names(mydata)[ncol(mydata)] <- 'Class_Label'

###Option for Resolve Class imbalance in Training Set Only
  samples <- sample(NROW(mydata), NROW(mydata) * .632)

  data.train <- mydata[samples, ]
  updata <- smotefamily::ADAS(data.train[,-ncol(data.train)],data.train$Class_Label)$data
  names(updata)[ncol(updata)] <- 'Class_Label'
  updata$Class_Label <- as.factor(updata$Class_Label)
  data.train <- updata

  data.test <- mydata[-samples, ]
  updata <- smotefamily::ADAS(data.test[,-ncol(data.test)],data.test$Class_Label)$data
  names(updata)[ncol(updata)] <- 'Class_Label'
  updata$Class_Label <- as.factor(updata$Class_Label)
  data.test <- updata
###---


###Option for Resolve Class imbalance in both Training and Testing Set
# updata <- smotefamily::ADAS(mydata[,-ncol(mydata)],mydata$Class_Label)$data
# names(updata)[ncol(updata)] <- 'Class_Label'
# updata$Class_Label <- as.factor(updata$Class_Label)
# samples <- sample(NROW(updata), NROW(updata) * .632)
# data.train <- updata[samples, ]
# data.test <- updata[-samples, ]
###---

  class(updata$Class_Label)
  table(updata$Class_Label)
  table(data.train$Class_Label)
  table(data.test$Class_Label)
  
  data.model <- randomForest(x = data.train[,-ncol(data.train)], y = data.train$Class_Label, ntree = 1500, importance=TRUE, proximity=TRUE, nodesize = 1)
  #data.model2 <- train(Class_Label ~., data.train, method = "rf", ntree = 1500, importance=TRUE, proximity=TRUE, nodesize = 1)
  
  # pred <- predict(data.model, data.test)
  # predictions=as.vector(data.model$votes[,2])
  # pred = prediction(predictions,data.test$Class_Label)
  pred.test <- predict(data.model, newdata=data.test, type='prob')
  data.result.test <- roc(data.test$Class_Label, pred.test[,2])
  
  
  pred.train <- predict(data.model, newdata=data.train, type='prob')
  data.result.train <- roc(data.train$Class_Label, pred.train[,2])
  
  #predictions=as.vector(data.model$votes[,2])
  #pred=prediction(predictions,target)
  
  auc <- data.result.test$auc[1]
  thr <- coords(data.result.test, "best", best.method="closest.topleft", ret=c("threshold"))
  acc <- coords(data.result.test, "best", best.method="closest.topleft", ret=c("accuracy"))
  ROC.results[f,] <- c(auc,thr,acc)
  
  png(paste0(paste0("ROC_",f),".png"),width = 512, height = 384, units = "px",
      pointsize = 12, bg = "white", res = NA, family = "", restoreConsole = TRUE,
      type = c("windows", "cairo", "cairo-png"))
  
  plot(data.result.test, col = "steelblue2", print.thres="best", print.thres.best.method="closest.topleft", 
       main = filterName[f],
       xaxt="n",yaxt="n",xlab=paste0(paste0("Specificity (1 - False Positive Rates)\nAUC =",(acc*100)),"%"), ylab = "Sensitivity (True Positive Rate)")
  
  axis(2, at = seq(0,1,by=.2), labels = seq(0,1, by=.2), tick = TRUE)
  axis(1, at = seq(1,0,by=-.2), labels = seq(1,0, by=-.2), tick = TRUE)
  axis(1, at=0:3, las = 1,xaxp = c(0, 0.5, 1), cex.axis = 1)
  dev.off()
  
  Sys.sleep(10)
  progress.bar$step()

}
write.csv(ROC.results, file = "ROC_results.csv")