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
library(readxl)
#install.packages("dbscan")
#library(formattable)
rm(list=ls())

setwd("C:/Class/Dissertation/Data/Meka/Thyroid/CSV/Public/")

fileName = matrix("",nrow = 56, ncol = 1)

filterName = matrix("",nrow = 56, ncol = 1)
k <- 1
for(i in seq(5, 17, by = 2)){
  for(j in seq(5, 12, by = 1)){
    fileName[k] = str_c("data_BSIF",i,i,j,"_thyroid.xlsx",collapse ="")
    filterName[k] = paste0(paste0(paste0("F",i),"-B"),j)
    k = k + 1
  }
}

f<- 44

png(paste0(paste0("Pub_PCA_Bo_Re_",f),".png"),width = 1024, height = 699, units = "px", 
    pointsize = 12, bg = "white", res = NA, family = "", restoreConsole = TRUE,
    type = c("windows", "cairo", "cairo-png"))

par(mfrow=c(2,2))

progress.bar <- create_progress_bar("text")
progress.bar$init(6)

for(f in c(44,45,48,52,53,56)){
  mydata <- read_xlsx(fileName[f], sheet = "Data")
  names(mydata)[ncol(mydata)] <- 'Class_Label'
  mydata$Class_Label <- as.factor(mydata$Class_Label)
  
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
  if(n_cols < 2){
    progress.bar$step()
    next
  }
  prin_comp <- prcomp(mydata[,-c(n_cols)], scale. = T)
  mydata <- data.frame(prin_comp$x, Class_Label = mydata$Class_Label)
  
  
  new_mydata <- Boruta(Class_Label~., data = mydata, doTrace = 0)
  
  #plotImpHistory(new_mydata)
  
  l <- length(new_mydata$finalDecision)
  m <- as.numeric(max(as.matrix(attStats(new_mydata))[,4]))
  
  plot(new_mydata,xlim=c(l-30,(l+3)), las="2")
  
  text((l-30),(m),"Importance Levels",adj=0, cex=.9)
  legend((l-30),m, fill=c("green", "yellow","red", "blue"),
         legend = c("Important", "Tentative", "Unimportant","Shadow Attributes"),
         bty = "n", text.width = .3, cex = .85)
  

  
  #Boruta_data <- cbind(Boruta_data,new_mydata)
  Boruta_data <- subset(mydata, select = c(getSelectedAttributes(new_mydata, withTentative = T)))
  Finaldata <- data.frame(Boruta_data, Class_Label = mydata$Class_Label)
  
  #mydata$Class_Label <- ifelse(mydata$Class_Label == 1, 'Malignant', 'Benign')
  #mydata$Class_Label <- ifelse(mydata$Class_Label == 'Malignant', 1 , 0)
  Finaldata$Class_Label <- as.factor(Finaldata$Class_Label)
  n_cols = ncol(Finaldata)
  
  if(n_cols < 2){
    progress.bar$step()
    next
  }
  
  # for(i in 1:(n_cols)){
  #   names(mydata)[i] <- paste0(paste0(names(mydata)[i],"_"),filterName[f])
  # }
  SMOTE_data <- smotefamily::ADAS(Finaldata[,-ncol(Finaldata)],Finaldata$Class_Label)$data
  names(SMOTE_data)[ncol(SMOTE_data)] <- 'Class_Label'
  SMOTE_data$Class_Label <- as.factor(SMOTE_data$Class_Label)
  
  SMOTE_data <- SMOTE_data[sample(nrow(SMOTE_data)),]
  samples <- sample(NROW(SMOTE_data), NROW(SMOTE_data) * .632)
  data.train <- SMOTE_data[samples, ]
  data.test <- SMOTE_data[-samples, ]
  
  data.model <- randomForest(x = data.train[,-ncol(data.train)], y = data.train$Class_Label, importance=TRUE, proximity=TRUE)
  pred.test <- predict(data.model, newdata=data.test, type='prob')
  data.result.test <- roc(data.test$Class_Label, pred.test[,1])
  
  a <- data.result.test$auc
  
  plot(data.result.test, col = "steelblue2", print.thres="best", print.thres.best.method="closest.topleft", 
       main =  paste0("Public Thyroid Cancer MRI Data",filterName[f]),
       xaxt="n",yaxt="n",xlab=paste0(paste0("Specificity (1 - False Positive Rates)\nAUC = ",(a*100)),"%"), ylab = "Sensitivity (True Positive Rate)")
  
  axis(2, at = seq(0,1,by=.2), labels = seq(0,1, by=.2), tick = TRUE)
  axis(1, at = seq(1,0,by=-.2), labels = seq(1,0, by=-.2), tick = TRUE)
  
  if(f %% 2 == 0){
    dev.off()
    png(paste0(paste0("Pub_PCA_Bo_Re_",f),".png"),width = 1024, height = 699, units = "px", 
        pointsize = 12, bg = "white", res = NA, family = "", restoreConsole = TRUE,
        type = c("windows", "cairo", "cairo-png"))
    
    par(mfrow=c(2,2))
  }
  
  progress.bar$step()
}
dev.off()