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
options(digits=22)
setwd("C:/Class/Dissertation/Data/Meka/Thyroid/CSV")

myTrainingControl <- trainControl(method = "cv", 
                                  number = 10, 
                                  savePredictions = TRUE, 
                                  classProbs = TRUE, 
                                  verboseIter = TRUE)
myTrainingControl <- trainControl(method="repeatedcv", number=10, repeats=3, savePredictions=TRUE, classProbs=TRUE, verbose=TRUE)
sampling.list <- list("SMOTE","ROSE","Up Sampling","Down Sampling")
algorithms.list <- list("lvq","rpart","gbm","svmRadial","rf","C5.0", "glm") #"treebag", "lda",
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


Ensemble.results.auc <- matrix(0, nrow = 56, ncol = 7)
Ensemble.results.threshold <- matrix(0, nrow = 56, ncol = 7)
Ensemble.results.accuracy <- matrix(0, nrow = (56*4), ncol = 7)

f_idx <- 1
progress.bar <- create_progress_bar("text")
progress.bar$init(224)
for(f in 1:56){
  mydata <- read.csv(file = fileName[f], header=TRUE)
  n_cols = ncol(mydata)
  
  # n_mydata <- preProcess(mydata[, -n_cols], method = c("center", "scale", "YeoJohnson", "nzv"))  
  # n_mydata$ica
  # comboInfo <- nearZeroVar(mydata, saveMetrics= TRUE)
  # dim(filteredDescr)
  # comboInfo$linearCombos
  # mydata <- mydata[sample(nrow(mydata)),]
  n_mydata <- mydata

  n_cols = ncol(n_mydata)
  rm_idx <- 0
  rm_cols <- c(0)

  for (h in 1:(n_cols-1))
  {
    if(wilcox.test(mydata[,h],mydata$Class_Label, data=mydata)$p.value >= 0.05){
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
  n_cols = ncol(mydata)
  rm_idx <- 0
  rm_cols <- c(0)
  for (h in 1:(n_cols-1))
  {
    if(sd(mydata[,h]) == 0){
      rm_idx <- rm_idx + 1
      rm_cols[rm_idx] = h
    }
  }
  
  if(rm_idx > 0)
    mydata <- mydata[,-c(rm_cols)]
  n_cols = ncol(mydata)
  
  # model <- randomForest(Class_Label ~., n_mydata, ntree = 50, importance=TRUE, nodesize = 1)
  # imp_attributes <- sort(data.frame(t(varImp(model, scale=TRUE))),decreasing = TRUE)
  # 
  # imp_row <- 1
  # imp_limit <-  imp_attributes[1,1]
  # while(imp_attributes[imp_row] >= imp_limit*0.5)
  #   imp_row  <- imp_row + 1
  # 
  # 
  # mydata <- data.frame(subset(n_mydata, select = c(names(imp_attributes[1:imp_row]))), Class_Label = n_mydata$Class_Label)
  # 
  mydata$Class_Label <- ifelse(mydata$Class_Label == 1, 'Malignant', 'Benign')
  mydata$Class_Label <- as.factor(mydata$Class_Label)
  
  samples <- sample(NROW(mydata), NROW(mydata) * .60)
  data.train <- mydata[samples, ]
  data.test <- mydata[-samples, ]
  
  data.train.list <- list(SMOTE(Class_Label ~ .,data  = data.train),
                          ROSE(Class_Label ~ .,data  = data.train)$data,
                          upSample(x = data.train[, -ncol(data.train)],y = data.train$Class_Label, yname = 'Class_Label'),
                          downSample(x = data.train[, -ncol(data.train)],y = data.train$Class_Label, yname = 'Class_Label'))
  #class(data.train.list[[4]]$Class_Label)
  #data.train.list <- lapply(data.train.list, function(x) {x$Class_Label <- (as.numeric(x$Class_Label)-1); return(x)})
  #data.test$Class_Label <- as.numeric(data.test$Class_Label)-1
  #model.list <- lapply(algorithms.list, function(fx) train(Class_Label~., data=mydata, method=fx, trControl=myTrainingControl))
  png(paste0(paste0("BoxPlotGraph_",f),".png"),width = 1024, height = 768, units = "px", 
       pointsize = 12, bg = "white", res = NA, family = "", restoreConsole = TRUE,
       type = c("windows", "cairo", "cairo-png"))
  par(mfrow=c(2,2))
  for(d in 1 : 4){
    tryCatch({
        data.trained <- data.frame(data.train.list[[d]])
        model.list <- lapply(algorithms.list, function(fx) train(Class_Label~., data=data.trained, method=fx, trControl=myTrainingControl)) # , verbose=FALSE
        results <- resamples(list('Learning Vector Quantization'=model.list[[1]], 'SimpleCART'=model.list[[2]], 'Stochastic Gradient Boosting'=model.list[[3]], 'SVM'=model.list[[4]], 'Random Forest'=model.list[[5]], C5.0=model.list[[6]], 'Logistic Regression(via GLM)'=model.list[[7]]))
    
        Ensemble.results.accuracy[f_idx,] <- summary(results)$statistics$Accuracy[,4]
        f_idx <- f_idx + 1
     
        main_data <- data.frame(results$value[,-c(1,3,5,7,9,11,13,15)], check.names = FALSE)
        Names <- c(rep("LVQ",30), rep( "CART",30), rep( "SGD",30), rep( "SVM",30), rep( "RM",30), rep("C5.0",30), rep( "LR(via GLM)",30))
    
        Values <- matrix(c(0), ncol = 1, nrow = (30*7))
        j <- 1
        for( i in 1:7){
          Values[j:(j+29),1] <-  main_data[,i]
            j <- j + 30
        }
        
        data=data.frame(Names,Values)
        boxplot(data$Values ~ data$Names, col=topo.colors(7) ,xaxt="n" , xlab = "", main = paste0(paste0(filterName[f]," with "),sampling.list[[d]]))
        axis(1, at=1:7, las = 1, labels = c("LVQ","CART","SGD","SVM","RF","C5.0","GLM"), cex.axis = 1)
    
        mylevels <-levels(data$Names)
        levelProportions <- summary(data$Names)/nrow(data)
        
        for(i in 1:length(mylevels)){
          
          thislevel <- mylevels[i]
          thisvalues <- data[data$Names==thislevel, "Values"]
          
          # take the x-axis indices and add a jitter, proportional to the N in each level
          myjitter<-jitter(rep(i, length(thisvalues)), amount=levelProportions[i]/2)
          points(myjitter, thisvalues, pch=20, col=rgb(0,0,0,.2)) 
          
        }
        Sys.sleep(10)
        progress.bar$step()
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  Sys.sleep(10)
  dev.off()
  
  
  
  
}
Ensemble.results.accuracy <- data.frame(Ensemble.results.accuracy)
names(Ensemble.results.accuracy) <- c("LVQ","CART","SGD","SVM","RM","C5.0","LR(via GLM)")
write.csv(Ensemble.results.accuracy, file = "EnsembleAccuracies.csv")