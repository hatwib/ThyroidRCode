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
rm(list=ls())

setwd("C:/Class/Dissertation/Data/Meka/Thyroid/CSV")

myTrainingControl <- trainControl(method = "cv", 
                                  number = 10, 
                                  savePredictions = TRUE, 
                                  classProbs = TRUE, 
                                  verboseIter = FALSE)
#sampling.list <- list("SMOTE","ROSE","Up Sampling","Down Sampling")
algorithms.list <- list("lvq","gbm","svmRadial","rf")
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


Ensemble.results.auc <- matrix(0, nrow = 56, ncol = 4)
Ensemble.results.threshold <- matrix(0, nrow = 56, ncol = 4)
Ensemble.results.accuracy <- matrix(0, nrow = 56, ncol = 4)

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
    if(wilcox.test(mydata[,h],mydata$Class_Label, data=mydata)$p.value > 0.05){
      rm_idx <- rm_idx + 1
      rm_cols[rm_idx] = h
    }
  }
  
  if(rm_idx > 0)
    n_mydata <- n_mydata[,-c(rm_cols)]
  
  n_cols = ncol(n_mydata)
  rm_idx <- 0
  rm_cols <- c(0)
  for (h in 1:(n_cols-1))
  {
    if(sd(n_mydata[,h]) == 0){
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
  n_cols = ncol(n_mydata)
  rm_idx <- 0
  rm_cols <- c(0)
  for (h in 1:(n_cols-1))
  {
    if(sd(n_mydata[,h]) == 0){
      rm_idx <- rm_idx + 1
      rm_cols[rm_idx] = h
    }
  }
  
  if(rm_idx > 0)
    n_mydata <- n_mydata[,-c(rm_cols)]
  n_cols = ncol(n_mydata)
  
  model <- randomForest(Class_Label ~., n_mydata, ntree = 50, importance=TRUE, nodesize = 1)
  imp_attributes <- sort(data.frame(t(varImp(model, scale=TRUE))),decreasing = TRUE)

  imp_row <- 1
  imp_limit <-  imp_attributes[1,1]
  while(imp_attributes[imp_row] >= imp_limit*0.5)
    imp_row  <- imp_row + 1


  mydata <- data.frame(subset(n_mydata, select = c(names(imp_attributes[1:imp_row]))), Class_Label = n_mydata$Class_Label)
  
  mydata$Class_Label <- ifelse(mydata$Class_Label == 1, 'Malignant', 'Benign')
  mydata$Class_Label <- as.factor(mydata$Class_Label)
  
  model.list <- lapply(algorithms.list, function(fx) train(Class_Label~., data=mydata, method=fx, trControl=myTrainingControl))
  
  results <- resamples(list(LVQ=model.list[[1]], GBM=model.list[[2]], SVM=model.list[[3]], RF=model.list[[4]]))
  
  
  Ensemble.results.accuracy[f,] <- summary(results)$statistics$Accuracy[,4]
  
      png(paste0(paste0("BoxPlot",f),".png"),width = 512, height = 350, units = "px", 
      pointsize = 12, bg = "white", res = NA, family = "", restoreConsole = TRUE,
      type = c("windows", "cairo", "cairo-png"))
      bwplot(results, main = "Box Plot")
      dev.off()
      
      # png(paste0(paste0("Dotplot_",f),".png"),width = 512, height = 350, units = "px", 
      #     pointsize = 12, bg = "white", res = NA, family = "", restoreConsole = TRUE,
      #     type = c("windows", "cairo", "cairo-png"))
      # dotplot(results)
      # dev.off()
  
  
  progress.bar$step()
  

}
