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


getFeatures <- function(dt, withT) {
  new_dt <- Boruta(Class_Label~., data = dt, doTrace = 0)
  if(withT)
    new_dt <- subset(dt, select = c(c(getSelectedAttributes(new_dt, withTentative = F), 'Class_Label')))
  else{
    
    new_dt <- subset(dt, select = c(c(getSelectedAttributes(new_dt, withTentative = F), 'Class_Label')))
  }
  return(new_dt)
}

myTrainingControl <- trainControl(method = "cv", number = 10, savePredictions = TRUE, classProbs = TRUE, verboseIter = FALSE)
myTrainingControl <- trainControl(method="repeatedcv", number=10, repeats=3, savePredictions=TRUE, classProbs=TRUE, verbose=FALSE)

sampling.list <- list("SMOTE","Borderline-SMOTE","Density based-SMOTE","Adaptive-SMOTE","ROSE","Up","Down")
algorithms.list <- list("C5.0","rf","gbm") #"treebag", "lda","lvq","gbm","svmRadial", "glm"

fileName = matrix("",nrow = 56, ncol = 1)

filterName = matrix("",nrow = 56, ncol = 1)
k <- 1
for(i in seq(5, 17, by = 2)){
  for(j in seq(5, 12, by = 1)){
    #fileName[k] = str_c("data_BSIF",i,i,j,"_thyroid.csv",collapse ="")
    filterName[k] = paste0(paste0(paste0("_F",i),"_B"),j)
    k = k + 1
  }
}


Ensemble.results.auc <- matrix(0, nrow = 56, ncol = 3)
Ensemble.results.threshold <- matrix(0, nrow = 56, ncol = 3)
Ensemble.results.accuracy <- matrix(0, nrow = (56*4), ncol = 3)

f_idx <- 1

rm(allData)
rm(ldata)
rm(mydata)
f <- 1
options(digits=22)
allData <- read.csv(file = fileName[1], header=TRUE)

# Histogram of example data and normalized data
n_cols = ncol(allData)
for (h in 2:(n_cols-1)){
  allData[,h] = (allData[,h]-min(allData[,h]))/(max(allData[,h])-min(allData[,h]))
}

rm_cols <- nearZeroVar(allData)
rm_idx <- length(rm_cols)

new_mydata <- Boruta(Class_Label~., data = allData, doTrace = 0)
allData <- subset(allData, select = c(getSelectedAttributes(new_mydata, withTentative = F)))



for(i in 1:ncol(allData)){
  names(allData)[i] <- paste0(paste0(names(allData)[i],"_"),1)
}

# allData <- data.frame(allData, Class_Label = mydata$Class_Label)
# 
# allData$Class_Label <- ifelse(allData$Class_Label == 1, 'Malignant', 'Benign')
# allData$Class_Label <- as.factor(allData$Class_Label)
# 
# allData$Class_Label <- as.factor(allData$Class_Label)
# table(allData$Class_Label)
# bm_mydata2 <- BLSMOTE(allData[,-16],allData$Class_Label, K=2)
# class(bm_mydata2$data$class)
#allData <- cbind(allData, mydata[-1,-c(ncol(mydata))])

f_idx <- 1
progress.bar <- create_progress_bar("text")
progress.bar$init(55*3)
for(f in 1:56){
  mydata <- read.csv(file = fileName[f], header=TRUE)
  
  n_cols = ncol(mydata)
  for (h in 1:(n_cols-1))
  {
    mydata[,h] = (mydata[,h]-min(mydata[,h]))/(max(mydata[,h])-min(mydata[,h]))
  }
  rm_cols <- nearZeroVar(mydata)
  rm_idx <- length(rm_cols)
  
  if(rm_idx > 0)
    mydata <- mydata[,-c(rm_cols)]
  
  n_cols = ncol(mydata)
  
  mydata$Class_Label <- as.factor(mydata$Class_Label)
  
  new_mydata <- Boruta(Class_Label~., data = mydata, doTrace = 0)
  
  if(length(getSelectedAttributes(new_mydata, withTentative = F)) < 5)
    next()
  mydata <- subset(mydata, select = c(c(getSelectedAttributes(new_mydata, withTentative = F)), 'Class_Label'))
  
  # n_cols = ncol(mydata)
  # if(n_cols > 0){
  #   for(i in 1:n_cols){
  #     names(mydata)[i] <- paste0(paste0(names(mydata)[i],"_"),f)
  #   }
  #   allData <- cbind(allData, mydata)
  # }
  progress.bar$step()

  mydata$Class_Label <- ifelse(mydata$Class_Label == 1, 'Malignant', 'Benign')
  mydata$Class_Label <- as.factor(mydata$Class_Label)

  mydata <- mydata[sample(nrow(mydata)),]



  data.train.list <- list(smotefamily::SMOTE(mydata[,-ncol(mydata)],mydata$Class_Label)$data,
                        smotefamily::BLSMOTE(mydata[,-ncol(mydata)],mydata$Class_Label)$data,
                        smotefamily::DBSMOTE(mydata[,-ncol(mydata)],mydata$Class_Label)$data,
                        smotefamily::ADAS(mydata[,-ncol(mydata)],mydata$Class_Label)$data,
                        ROSE(Class_Label ~ .,data  = mydata)$data,
                        upSample(x = mydata[, -ncol(mydata)],y = mydata$Class_Label, yname = 'Class_Label'),
                        downSample(x = mydata[, -ncol(mydata)],y = mydata$Class_Label, yname = 'Class_Label'))


# data.train.list <- list(downSample(x = mydata[, -ncol(mydata)],y = mydata$Class_Label, yname = 'Class_Label'))
ADAS

  png(paste0(paste0("allData",filterName[f]),".png"),width = 2048, height = 768, units = "px", 
      pointsize = 12, bg = "white", res = NA, family = "", restoreConsole = TRUE,
      type = c("windows", "cairo", "cairo-png"))
  par(mfrow=c(1,7))
  
  for(d in 1 : 7){
    tryCatch({
      data.trained <- data.frame(data.train.list[[d]])
      # if(d  >= 6){
      #   
      #   data.trained <- cbind(data.trained$x,data.trained$y)
      # }
      names(data.trained)[ncol(data.trained)] <- 'Class_Label'
      #data.trained$Class_Label <- ifelse(data.trained$Class_Label == 1, 'Malignant', 'Benign')
      data.trained$Class_Label <- as.factor(data.trained$Class_Label)
      
      model.list <- lapply(algorithms.list, function(fx) train(Class_Label~., data=data.trained, method=fx, trControl=myTrainingControl)) # , verbose=FALSE
      progress.bar$step()
      results <- resamples(list('C5.0 Decision Tree'= model.list[[1]], 'Random Forest'=model.list[[2]], 'Stocastic Gradient Boosting'=model.list[[3]]))
      
      Ensemble.results.accuracy[f_idx,] <- summary(results)$statistics$Accuracy[,4] #f_idx <- 1
      f_idx <- f_idx + 1
      
      main_data <- data.frame(results$value[,-c(1,3,5,7)], check.names = FALSE)
      Names <- c(rep("C5.0",30),  rep( "RM",30), rep("SGB",30))
      
      Values <- matrix(c(0), ncol = 1, nrow = (30*3))
      j <- 1
      for( i in 1:3){
        Values[j:(j+29),1] <-  main_data[,i]
        j <- j + 30
      }
      
      data= data.frame(Names,Values)
      boxplot(c((data$Values*100),(data$Values*100)) ~ c(data$Names,data$Names), col=c('blue','grey','orange'), ylab = "Model Accuracy (%)" ,xaxt="n" , xlab = "", 
              main = paste0(paste0(paste0("Boosting with 10 Fold CV and ",sampling.list[[d]])," Sampling "),filterName[f]))
      axis(1, at=1:3, las = 1, labels = c("C5.0 Decision Tree","Random Forest","Stochastic Gradient Boosting"), cex.axis = 1)
      
      mylevels <-levels(data$Names)
      levelProportions <- summary(data$Names)/nrow(data)
      
      for(i in 1:length(mylevels)){
        
        thislevel <- mylevels[i]
        thisvalues <- data[data$Names==thislevel, "Values"]
        
        # take the x-axis indices and add a jitter, proportional to the N in each level
        myjitter<-jitter(rep(i, length(thisvalues)), amount=levelProportions[i]/2)
        points(myjitter, thisvalues, pch=20, col=rgb(0,0,0,.2)) 
        
      }
      Sys.sleep(20)
      progress.bar$step()
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  Sys.sleep(20)
  dev.off()

}
Ensemble.results.accuracy <- data.frame(Ensemble.results.accuracy)
names(Ensemble.results.accuracy) <- c("C5.0","RM","GSB")
write.csv(Ensemble.results.accuracy, file = "Ensemble_Filter_Bit.csv")


#====================================================================


