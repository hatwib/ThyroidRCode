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
rm(list=ls())

setwd("C:/Class/Dissertation/Data/Meka/Thyroid/CSV")

fileName = matrix("",nrow = 56, ncol = 1)
k = 1
filterName = matrix("",nrow = 56, ncol = 1)

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
  
  h <- c(1)
  n_cols = ncol(n_mydata)
  
  for (h in 1:(n_cols-1))
  {
    mn_value = min(n_mydata[,h])
    mx_value = max(n_mydata[,h])
    n_mydata[,h] = (n_mydata[,h]-mn_value)/(mx_value-mn_value)
    if(sd(n_mydata[,h]) > 0){
      coefY <- powerTransform(n_mydata$Class_Label ~ n_mydata[,h], family="yjPower")
      transformX <- yjPower(n_mydata[,h], lambda=coef(coefY), jacobian.adjusted = FALSE)
      n_mydata[,h] =  transformX
    }else
      print(h)
  }
  
  
  prin_comp <- prcomp(n_mydata, scale. = T)
  
  std_dev <- prin_comp$sdev
  pr_var <- std_dev^2
  prop_varex <- pr_var/sum(pr_var)
  
  PCA_row <- 1
  PCA_limit <- 0.85
  
  while(sum(prop_varex[1:PCA_row]) < PCA_limit)
    PCA_row  <- PCA_row + 1
  
  
  mydata <- data.frame(prin_comp$x[,1:PCA_row], Class_Label = n_mydata$Class_Label)
  #mydata <- n_mydata
  
  #mydata <- read.csv(file = "data_BSIF5510_ovarianMRI.csv", header=TRUE)
  mydata$Class_Label <- ifelse(mydata$Class_Label == 1, 'Malignant', 'Benign')
  mydata$Class_Label <- as.factor(mydata$Class_Label)
  #table(mydata$Class_Label)
  
  #testidx <- sample(nrow(mydata), nrow(mydata)*0.75)
  #testidx <- which(1:nrow(mydata)%%4==0)
  
  #mydata_train1 <- mydata[-testidx,]
  #mydata_test <- mydata[testidx,]
  
  #mydata_train <- ROSE(Class_Label ~ ., data  = mydata_train)
  #mydata_train <- SMOTE(Class_Label ~ ., data  = mydata_train)
  #mydata_train <- upSample(x = mydata_train[, -ncol(mydata_train)],y = mydata_train$Class_Label, yname = 'Class_Label')
  #mydata_test <- downSample(x = mydata_test[, -ncol(mydata_test)],y = mydata_test$Class_Label, yname = 'Class_Label')
  #table(mydata_train$Class_Label)
  #table(mydata_train1$Class_Label)
  #table(mydata_test$Class_Label)
  
  #table(data.train$Class_Label)
  #mydata$Class_Label <- factor(mydata$Class_Label)  # setosa should be removed from factor
  samples <- sample(NROW(mydata), NROW(mydata) * .60)
  data.train <- mydata[samples, ]
  data.test <- mydata[-samples, ]
  sampling.list <- list("SMOTE","ROSE","Up Sampling","Down Sampling")
  data.train.list <- list(SMOTE(Class_Label ~ .,data  = data.train),
                          ROSE(Class_Label ~ .,data  = data.train)$data,
                          upSample(x = data.train[, -ncol(data.train)],y = data.train$Class_Label, yname = 'Class_Label'),
                          downSample(x = data.train[, -ncol(data.train)],y = data.train$Class_Label, yname = 'Class_Label'))
  
  #data.train.list <- list(data.train,data.train.SMOTE,data.train.ROSE,data.train.UP,data.train.DOWN)
  
  myTrainingControl <- trainControl(method = "cv", 
                                    number = 10, 
                                    savePredictions = TRUE, 
                                    classProbs = TRUE, 
                                    verboseIter = FALSE)
  
  data.model <- train(Class_Label ~., data.train, method = "rf", trControl = myTrainingControl , ntree = 500, importance=TRUE, nodesize = 1)
  data.pred <- predict(data.model, data.test, type="prob") # Prediction
  data.result <- roc(data.test$Class_Label, data.pred$Benign)
  data.auc[f,] <- data.result$auc
  data.threshold[f,] <- coords(data.result, "best", best.method="closest.topleft", ret=c("threshold"))
  data.accuracy[f,] <- coords(data.result, "best", best.method="closest.topleft", ret=c("accuracy"))
  
  mymodel.rf <- lapply(data.train.list, function(x) train(Class_Label ~., x, method = "rf", 
                                                          #preProcess = c("scale","center"), 
                                                          trControl = myTrainingControl, 
                                                          ntree = 500, importance=TRUE, nodesize = 1))
  
  # forest.model <- train(Class_Label ~., data.train)
  # forest.model.SMOTE <- train(Class_Label ~., data.train.SMOTE)
  # forest.model.ROSE <- train(Class_Label ~., data.train.ROSE)
  # forest.model.UP <- train(Class_Label ~., data.train.UP)
  # forest.model.DOWN <- train(Class_Label ~., data.train.DOWN)
  
  
  pred.rf <- lapply(mymodel.rf, function(x) predict(x, data.test, type="prob"))
  result.rf <- lapply(pred.rf, function(x) roc(data.test$Class_Label, x$Benign))
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
#result.predicted.prob <- predict(forest.model.DOWN, data.test, type="prob") # Prediction
#result.predicted.class <- predict(forest.model.DOWN, data.test, type="raw") # Class
#pre_results <- table(result.predicted.class, data.test$Class_Label)
result.roc <- roc(data.test$Class_Label, result.predicted.prob$Malignant) # Draw ROC curve.
plot(result.roc, print.thres="best", print.thres.best.method="closest.topleft")
result.coords <- coords(result.roc, "best", best.method="closest.topleft", ret=c("threshold", "accuracy"))
print(result.coords)


model <- randomForest(Class_Label~., data=mydata_train, importance=TRUE, ntree=500)

pred <- predict(model, newdata=mydata_test, type='prob')
performance(pred,"tpr","fpr")
result.roc <- roc(mydata_test$Class_Label, pred[,1])
plot(result.roc, print.thres="best", print.thres.best.method="closest.topleft")
result.coords <- coords(result.roc, "best", best.method="closest.topleft", ret=c("threshold", "accuracy"))
print(result.coords)

pred <- predict(model, newdata=mydata_test, type='class')
pre_results <- table(pred, mydata_test$Class_Label)
my_results <- (pre_results[1]+pre_results[4])/ nrow(mydata_test)


confusionMatrix(pred, mydata_test$Class_Label)

plot(perf,col='red',lwd=3)
abline(a=0,b=1,lwd=2,lty=2,col="gray")

plot_pred_type_distribution(prediction)
myRoc <- roc(Class_Label., auc.polygon=TRUE, grid=TRUE, plot=FALSE)


pca_mydata <- mydata[,-(n_cols)]
prin_comp.bc_norm <- prcomp(pca_mydata, scale. = T)

std_dev.bc_norm <- prin_comp.bc_norm$sdev
pr_var.bc_norm <- std_dev.bc_norm^2
prop_varex.bc_norm <- pr_var.bc_norm/sum(pr_var.bc_norm)

#dev.off() 
plot(prop_varex.bc_norm[1:100], xlab = "Principal Component", ylab = "Variance",type = "b")
PCA_row = 10
lines(c(PCA_row,PCA_row),c(-1,prop_varex.bc_norm[PCA_row]), col='blue')


#dev.off() 
plot(cumsum(prop_varex.bc_norm[1:20]), xlab = "Principal Component", ylab = "Cumulative Variance", type = "b")
#lines(c(19,19),c(-1,sum(prop_varex.z_norm[1:19])), col='blue')
lines(c(-10,PCA_row),c(sum(prop_varex.bc_norm[1:PCA_row]),sum(prop_varex.bc_norm[1:PCA_row])), col='blue')

pca_mydata <- data.frame(Class_Label = mydata$Class_Label, dim(prin_comp.bc_norm$x)


pca_mydata <- pca_mydata[,-c(22:129)]


pca_mydata <- mydata[sample(nrow(pca_mydata)),]

testidx <- sample(nrow(pca_mydata), nrow(pca_mydata)*0.5)

mydata_train <- pca_mydata[testidx,]
mydata_test <- pca_mydata[-testidx,]

mydata_train <- mydata[sample(which(mydata_train$Class_Label == 'Malignant'),  replace = TRUE),]
#table(mydata_train$Class_Label)

model <- randomForest(Class_Label~., data=mydata_train, importance=TRUE, ntree=2000)

prediction <- predict(model, newdata=mydata_test, type='class')
pre_results <- table(prediction, mydata_test$Class_Label)
my_results <- (pre_results[1]+pre_results[4])/ nrow(mydata_test)