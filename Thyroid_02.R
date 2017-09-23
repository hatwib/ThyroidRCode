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

data.auc <- matrix(0, nrow = 56, ncol = 1)
data.threshold <- matrix(0, nrow = 56, ncol = 1)
data.accuracy <- matrix(0, nrow = 56, ncol = 1)

results.auc <- matrix(0, nrow = 56, ncol = 4)
results.threshold <- matrix(0, nrow = 56, ncol = 4)
results.accuracy <- matrix(0, nrow = 56, ncol = 4)

progress.bar <- create_progress_bar("text")
progress.bar$init(56)

for(f in 1:56){
  mydata <- read.csv(file = fileName[f], header=TRUE)
  
  mydata <- mydata[sample(nrow(mydata)),]
  n_mydata <- mydata

  n_cols = ncol(n_mydata)
  
  rm_idx <- 0
  rm_cols <- c(0)
  for (h in 1:(n_cols-1))
  {
    if(sd(mydata[,h]) == 0){
      rm_idx <- rm_idx + 1
      rm_cols[rm_idx] = h
    }else{
      mn_value = min(n_mydata[,h])
      mx_value = max(n_mydata[,h])
      n_mydata[,h] = (n_mydata[,h]-mn_value)/(mx_value-mn_value)
      coefY <- powerTransform(n_mydata$Class_Label ~ n_mydata[,h], family="yjPower")
      transformX <- yjPower(n_mydata[,h], lambda=coef(coefY), jacobian.adjusted = FALSE)
      n_mydata[,h] =  transformX
    }
  }

  if(rm_idx > 0)
    n_mydata <- n_mydata[,-c(rm_cols)]
  
  prin_comp <- prcomp(n_mydata, scale. = T)
  
  std_dev <- prin_comp$sdev
  pr_var <- std_dev^2
  prop_varex <- pr_var/sum(pr_var)
  
  PCA_row <- 1
  PCA_limit <- 0.95
  
  while(sum(prop_varex[1:PCA_row]) < PCA_limit)
    PCA_row  <- PCA_row + 1
  
  mydata <- data.frame(prin_comp$x[,1:PCA_row], Class_Label = n_mydata$Class_Label)
  mydata$Class_Label <- ifelse(mydata$Class_Label == 1, 'Malignant', 'Benign')
  mydata$Class_Label <- as.factor(mydata$Class_Label)

  samples <- sample(NROW(mydata), NROW(mydata) * .60)
  data.train <- mydata[samples, ]
  data.test <- mydata[-samples, ]
  sampling.list <- list("SMOTE","ROSE","Up Sampling","Down Sampling")
  data.train.list <- list(SMOTE(Class_Label ~ .,data  = data.train),
                          ROSE(Class_Label ~ .,data  = data.train)$data,
                          upSample(x = data.train[, -ncol(data.train)],y = data.train$Class_Label, yname = 'Class_Label'),
                          downSample(x = data.train[, -ncol(data.train)],y = data.train$Class_Label, yname = 'Class_Label'))

  myTrainingControl <- trainControl(method = "cv", 
                                    number = 10, 
                                    savePredictions = TRUE, 
                                    classProbs = TRUE, 
                                    verboseIter = FALSE)
  
  data.model <- train(Class_Label ~., data.train, method = "rf", trControl = myTrainingControl , ntree = 500, importance=TRUE, nodesize = 1)
  data.pred <- predict(data.model, data.test, type="prob") # Prediction
  data.result <- roc(data.test$Class_Label, data.pred$Malignant)
  data.auc[f,] <- data.result$auc
  data.threshold[f,] <- coords(data.result, "best", best.method="closest.topleft", ret=c("threshold"))
  data.accuracy[f,] <- coords(data.result, "best", best.method="closest.topleft", ret=c("accuracy"))
  
  mymodel.rf <- lapply(data.train.list, function(x) train(Class_Label ~., x, method = "rf", 
                                                          #preProcess = c("scale","center"), 
                                                          trControl = myTrainingControl, 
                                                          ntree = 500, importance=TRUE, nodesize = 1))
  
  pred.rf <- lapply(mymodel.rf, function(x) predict(x, data.test, type="prob"))
  result.rf <- lapply(pred.rf, function(x) roc(data.test$Class_Label, x$Malignant))
  result.rf
      png(paste0(paste0("ROC_",f),".png"),width = 1024, height = 699, units = "px", 
      pointsize = 12, bg = "white", res = NA, family = "", restoreConsole = TRUE,
      type = c("windows", "cairo", "cairo-png"))
      par(mfrow=c(2,2))
      mapply(function(x,y) {plot(x, print.thres="best", print.thres.best.method="closest.topleft", main = paste0(paste0(paste0("ROC for ",filterName[f])," with "),y))
      lines(data.result, col = "red")
  
    }, result.rf, sampling.list)
  dev.off()
  
  results.threshold[f,] <- mapply(function(x) coords(x, "best", best.method="closest.topleft", ret=c("threshold")), result.rf)
  results.accuracy[f,] <- mapply(function(x) coords(x, "best", best.method="closest.topleft", ret=c("accuracy")), result.rf)
  results.auc[f,] <- mapply(function(x) x$auc, result.rf)
  progress.bar$step()
}
