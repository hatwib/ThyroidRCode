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
install.packages("dbscan")
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

f <- 1

options(digits=22)
allData <- read_xlsx(fileName[1], sheet = "Data", col_names = TRUE)

names(allData)[ncol(allData)] <- 'Class_Label'
allData$Class_Label <- as.factor(allData$Class_Label)

n_cols = ncol(allData)
ldata <- allData[,c(n_cols)]
for (h in 1:(n_cols-1))
{
  allData[,h] = (allData[,h]-min(allData[,h]))/(max(allData[,h])-min(allData[,h]))
}
rm_cols <- nearZeroVar(allData)
rm_idx <- length(rm_cols)

if(rm_idx > 0)
  allData <- allData[,-c(rm_cols)]
n_cols = ncol(allData)
prin_comp <- prcomp(allData[,-c(n_cols)], scale. = T)

allData <- data.frame(prin_comp$x, Class_Label = allData$Class_Label)

Boruta_data <- Boruta(Class_Label~., data = allData, doTrace = 0)
# png(paste0(paste0("PublicImportance_",f),".png"),width = 1024, height = 699, units = "px", 
#     pointsize = 12, bg = "white", res = NA, family = "", restoreConsole = TRUE,
#     type = c("windows", "cairo", "cairo-png"))
# 
# par(mfrow=c(2,2))
# 
# plotImpHistory(Boruta_data)
# plot(Boruta_data, las="2")
# text(0,11.1,"Importance Levels",adj=0, cex=.9)
# legend(0,11, fill=c("green", "yellow","red", "blue"),
#        legend = c("Important", "Tentative", "Unimportant","Shadow Attributes"),
#        bty = "n", text.width = .3, cex = .85)


allData <- subset(allData, select = c(getSelectedAttributes(Boruta_data, withTentative = F)))
n_cols = ncol(allData)
for(i in 1:ncol(allData)){
  names(allData)[i] <- paste0(paste0(names(allData)[i],"_"),filterName[1])
}


f <- 2
progress.bar <- create_progress_bar("text")
progress.bar$init(55)
for(f in 2:56){
  #mydata <- read.csv(file = fileName[f], header=TRUE)
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
  
  # plotImpHistory(new_mydata)
  # plot(new_mydata, las="2")
  # text(0,11.1,"Importance Levels",adj=0, cex=.9)
  # legend(0,11, fill=c("green", "yellow","red", "blue"),
  #        legend = c("Important", "Tentative", "Unimportant","Shadow Attributes"),
  #        bty = "n", text.width = .3, cex = .85)
  # 
  # if(f %% 2 == 0){
  #   dev.off()
  #   png(paste0(paste0("PublicImportance_",f),".png"),width = 1024, height = 699, units = "px", 
  #       pointsize = 12, bg = "white", res = NA, family = "", restoreConsole = TRUE,
  #       type = c("windows", "cairo", "cairo-png"))
  #   
  #   par(mfrow=c(2,2))
  # }
  # dev.off()
  
  # Boruta_data <- cbind(Boruta_data,new_mydata)
  mydata <- subset(mydata, select = c(getSelectedAttributes(new_mydata, withTentative = F)))
  
  n_cols = ncol(mydata)
  
  if(n_cols < 1){
    progress.bar$step()
    next
  }
  
  for(i in 1:(n_cols)){
    names(mydata)[i] <- paste0(paste0(names(mydata)[i],"_"),filterName[f])
  }
  allData <- cbind(allData, mydata)
  
  progress.bar$step()
}
dev.off()
mydata <- data.frame(allData, Class_Label = ldata)


SMOTE_data <- smotefamily::ADAS(mydata[,-ncol(mydata)],mydata$Class_Label)$data
names(SMOTE_data)[ncol(SMOTE_data)] <- 'Class_Label'
SMOTE_data$Class_Label <- as.factor(SMOTE_data$Class_Label)

table(SMOTE_data$Class_Label)

#Boruta_SMOTE_data <- Boruta(Class_Label~., data = SMOTE_data, doTrace = 0)
#,srt=90
#attStats(Boruta_SMOTE_data)
# Final_data <- subset(SMOTE_data, select = c(getSelectedAttributes(Boruta_SMOTE_data, withTentative = F)))
# Final_data <- data.frame(Final_data, Class_Label = SMOTE_data$Class_Label)

Final_data <- SMOTE_data

Final_data <- Final_data[sample(nrow(Final_data)),]
samples <- sample(NROW(Final_data), NROW(Final_data) * .632)
data.train <- Final_data[samples, ]
data.test <- Final_data[-samples, ]

data.model <- randomForest(x = data.train[,-ncol(data.train)], y = data.train$Class_Label, importance=TRUE, proximity=TRUE, trControl=myTrainingControl)
# data.model <- train(Class_Label ~., data.train, method = "rf", trControl = myTrainingControl , ntree = 1500, importance=TRUE, proximity=TRUE, nodesize = 1)
varImpPlot(data.model, type =1)
# pred <- predict(data.model, data.test)
# predictions=as.vector(data.model$votes[,2])
# pred = prediction(predictions,data.test$Class_Label)
pred.test <- predict(data.model, newdata=data.test, type='prob')
data.result.test <- roc(data.test$Class_Label, pred.test[,1])



plot(data.result.test, col = "steelblue2", print.thres="best", print.thres.best.method="closest.topleft", 
     xaxt="n",yaxt="n",xlab="Specificity (1 - False Positive Rates)\nAUC = 99.86%", ylab = "Sensitivity (True Positive Rate)")

axis(2, at = seq(0,1,by=.2), labels = seq(0,1, by=.2), tick = TRUE)
axis(1, at = seq(1,0,by=-.2), labels = seq(1,0, by=-.2), tick = TRUE)



for (h in 1 : 56){
  Boruta_SMOTE_data <- Boruta_data[1]
  class(Boruta_SMOTE_data)
  attr(Boruta_SMOTE_data, "class") <- "Boruta"
  plot(Boruta_SMOTE_data, xlim=range(136:146))
  attStats(Boruta_SMOTE_data)
  
  
}

#write.csv(allData, file = "allData_BorutaSelection.csv")
plot(new_mydata)
#for(f in 1:1){
ldata <- read.csv(file = fileName[1], header=TRUE)

# ldata.train.list <- list(smotefamily::SMOTE(ldata[,-ncol(ldata)],ldata$Class_Label, K=5)$data,
#                          ds <- smotefamily::BLSMOTE(ldata[,-ncol(ldata)],ldata$Class_Label)$data,
#                          smotefamily::DBSMOTE(ldata[,-ncol(ldata)],ldata$Class_Label)$data,
#                          smotefamily::ADAS(ldata[,-ncol(ldata)],ldata$Class_Label)$data,
#                         ROSE(Class_Label ~ .,data  = ldata)$data,
#                         upSample(x = ldata[, -ncol(ldata)],y = ldata$Class_Label, yname = 'Class_Label'),
#                         ds <- downSample(x = ldata[, -ncol(ldata)],y = ldata$Class_Label, yname = 'Class_Label'))


# ds <- cbind(ds$x,ds$y)
# names(ds)[ncol(ds)] <- 'Class_Label'
# 
# ds$Class_Label <- ifelse(ds$Class_Label == 1, 'Malignant', 'Benign')
# ds[,ncol(ds)] <- as.factor(ds[,ncol(ds)])
# levels(ds[,ncol(ds)])
# make.names(levels(ds[,ncol(ds)]))
# levels(ds[,ncol(ds)])
# str(ds)
# md <- train(Class_Label~., data=ds, method='rf', trControl=myTrainingControl)


mydata <- data.frame(allData, Class_Label = ldata[,33])


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

progress.bar <- create_progress_bar("text")
progress.bar$init(2*7)
# png("allData2_1.png",width = 1024, height = 768, units = "px", 
#     pointsize = 12, bg = "white", res = NA, family = "", restoreConsole = TRUE,
#     type = c("windows", "cairo", "cairo-png"))

par(mfrow=c(2,2))
for(d in 5 : 7){
  tryCatch({
    if(d == 5){
      # dev.off()
      png("allData2_2.png",width = 1024, height = 768, units = "px", 
          pointsize = 12, bg = "white", res = NA, family = "", restoreConsole = TRUE,
          type = c("windows", "cairo", "cairo-png"))
      par(mfrow=c(2,2))
    }
    data.trained <- data.frame(data.train.list[[d]])
    # if(d >= 6){
    #   data.trained <- cbind(data.trained$x,data.trained$y)
    # }
    names(data.trained)[ncol(data.trained)] <- 'Class_Label'
    #data.trained$Class_Label <- ifelse(data.trained$Class_Label == 1, 'Malignant', 'Benign')
    data.trained$Class_Label <- as.factor(data.trained$Class_Label) # class(data.trained$Class_Label)
    
    model.list <- lapply(algorithms.list, function(fx) train(Class_Label~., data=data.trained, method=fx, trControl=myTrainingControl)) # , verbose=FALSE
    progress.bar$step()
    results <- resamples(list('C5.0 Decision Tree'= model.list[[1]], 'Random Forest'=model.list[[2]], 'Stocastic Gradient Boosting'=model.list[[3]]))
    
    Ensemble.results.accuracy[f_idx,] <- summary(results)$statistics$Accuracy[,4]
    f_idx <- f_idx + 1
    
    main_data <- data.frame(results$value[,-c(1,3,5,7)], check.names = FALSE)
    Names <- c(rep("C5.0",30),  rep( "RM",30), rep("SGB",30))
    
    Values <- matrix(c(0), ncol = 1, nrow = (30*3))
    j <- 1
    for( i in 1:3){
      Values[j:(j+29),1] <-  main_data[,i]
      j <- j + 30
    }
    
    data=data.frame(Names,Values)
    boxplot((data$Values*100) ~ data$Names, col=c('blue','grey','orange'), ylab = "Model Accuracy (%)" ,ylim = c(50, 100),
            xaxt="n" , xlab = "", 
            main = paste0(paste0("Boosting with 10 Fold CV and ",sampling.list[[d]])," Sampling"))
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


Ensemble.results.accuracy <- data.frame(Ensemble.results.accuracy)
names(Ensemble.results.accuracy) <- c("C5.0","RM","GSB")
write.csv(Ensemble.results.accuracy, file = "Ensemble2.csv")


#====================================================================


