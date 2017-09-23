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
library(devtools)
library(ggbiplot)
rm(list=ls())
options(digits=22)
setwd("C:/Class/Dissertation/Data/Meka/Thyroid/CSV")
mydata <- read.csv(file = "allData_BorutaSelection.csv", header=TRUE)

#myTrainingControl <- trainControl(method = "cv", number = 10, savePredictions = TRUE, classProbs = TRUE, verboseIter = TRUE)
myTrainingControl <- trainControl(method="repeatedcv", number=10, repeats=3, savePredictions=TRUE, classProbs=TRUE, verbose=TRUE)

tmpdata <- read.csv(file = "data_BSIF555_thyroid.csv", header=TRUE)
mydata <- mydata[,-1] 
mydata <- data.frame(mydata,Class_Label= tmpdata$Class_Label)
mydata$Class_Label <- as.factor(mydata$Class_Label)
mydata <- mydata[sample(nrow(mydata)),]

#mydata$Class_Label <- ifelse(mydata$Class_Label == 1, 'Malignant', 'Benign')


######### Start PCA

mydata <- read.csv(file = "data_BSIF171712_thyroid.csv", header=TRUE)
pca_mydata <- mydata[,1:(ncol(mydata)-1)]
prin_comp <- prcomp(pca_mydata, center = TRUE, scale. = TRUE)
std_dev <- prin_comp$sdev
pr_var <- std_dev^2
prop_varex <- pr_var/sum(pr_var)

tmp_var <- cumsum(prop_varex[1:2])[2]
if(tmp_var > max_var)
  max_var <- tmp_var

g <- ggbiplot(prin_comp, obs.scale = 1, var.scale = 1, groups = mydata$Class_Label, ellipse = TRUE, circle = TRUE, var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', legend.position = 'top')
print(g)


pca_log_mydata <- log(mydata[,1:(ncol(mydata)-1)]+1)
prin_comp_log <- prcomp(pca_log_mydata, center = TRUE, scale. = TRUE)
plot(prin_comp_log, type = "l")

g <- ggbiplot(prin_comp_log,pc.biplot = FALSE, choices = c(1,2), obs.scale = 1, var.scale = 1, groups = mydata$Class_Label,labels =NULL, ellipse = TRUE, circle = TRUE
              , var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', legend.position = 'top')
print(g)
summary(prin_comp_log)
summary(prin_comp_log)


std_dev <- prin_comp_log$sdev
pr_var <- std_dev^2
prop_varex <- pr_var/sum(pr_var)
PCA_row <- 1
PCA_limit <- 0.95

while(sum(prop_varex[1:PCA_row]) < PCA_limit)
  PCA_row  <- PCA_row + 1

plot(prop_varex[1:10], xlab = "Principal Component", ylab = "Variance",type = "b")
lines(c(PCA_row,PCA_row),c(-1,prop_varex.bc_norm[PCA_row]), col='blue')
plot(cumsum(prop_varex[1:10]), xlab = "Principal Component", ylab = "Cumulative Variance", type = "b")
lines(c(-10,PCA_row),c(sum(prop_varex[1:PCA_row]),sum(prop_varex[1:PCA_row])), col='blue')



mydata <- data.frame(prin_comp$x[,1:5], Class_Label = mydata$Class_Label)


######## End PCA

# mydata$Class_Label
# 

###Option for Resolve Class imbalance in Training Set Only
mydata$Class_Label <- ifelse(mydata$Class_Label == 1, 'Malignant', 'Benign')
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
# updata <- updata[sample(nrow(updata)),]
# 
# samples <- sample(NROW(updata), NROW(updata) * .632)
# data.train <- updata[samples, ]
# data.test <- updata[-samples, ]
###---

class(updata$Class_Label)
table(updata$Class_Label)
table(data.train$Class_Label)
table(data.test$Class_Label)

class(data.train$Class_Label)
class(data.test$Class_Label)

data.test$Class_Label <- ifelse(data.test$Class_Label == 1, 'Malignant', 'Benign')
data.train$Class_Label <- ifelse(data.train$Class_Label == 1, 'Malignant', 'Benign')

which(data.test$Class_Label == 'Benign')  
which(pred.test[,2] < 0.553)
View(data.test)
data.result.test$cases



data.model <- randomForest(x = data.train[,-ncol(data.train)], y = data.train$Class_Label, ntree = 1500, importance=TRUE, proximity=TRUE
                           , trControl=myTrainingControl)
# data.model <- train(Class_Label ~., data.train, method = "rf", trControl = myTrainingControl , ntree = 1500, importance=TRUE, proximity=TRUE, nodesize = 1)

# pred <- predict(data.model, data.test)
# predictions=as.vector(data.model$votes[,2])
# pred = prediction(predictions,data.test$Class_Label)
pred.test <- predict(data.model, newdata=data.test, type='prob')
data.result.test <- roc(data.test$Class_Label, pred.test[,1])



plot(data.result.test, col = "steelblue2", print.thres="best", print.thres.best.method="closest.topleft", 
     main =  "ROC of Random Forest with Adaptive-SMOTE Sampling Technique",
     xaxt="n",yaxt="n",xlab="Specificity (1 - False Positive Rates)\nAUC = 99.14%", ylab = "Sensitivity (True Positive Rate)")

axis(2, at = seq(0,1,by=.2), labels = seq(0,1, by=.2), tick = TRUE)
axis(1, at = seq(1,0,by=-.2), labels = seq(1,0, by=-.2), tick = TRUE)


 
#  
# # pred.train <- predict(data.model, newdata=data.train, type='prob')
# # data.result.train <- roc(data.train$Class_Label, pred.train[,2])
# # 
# # #predictions=as.vector(data.model$votes[,2])
# # #pred=prediction(predictions,target)
# # 
# # 
# # data.result.test$auc
# # 
# # plot(data.result.test, col = "steelblue2", print.thres="best", print.thres.best.method="closest.topleft", 
# #      main =  "ROC of Random Forest with Adaptive-SMOTE Sampling Technique")
# # lines(data.result.train, col = gray(0.4))
# # legend(x=0.19,y=0.5, legend = c("Testing", "Training"), 
# #        lty = c(1, 1), lwd = rep(2, 3), col = c(gray(0.4), "steelblue2"), 
# #        text.width = 0.4, cex = 0.85)
# # 
# # 
# # data.result.test
# # 
# # coords(data.result.test, "best", best.method="closest.topleft", ret=c("threshold"))
# # coords(data.result.test, "best", best.method="closest.topleft", ret=c("accuracy"))
# # 
# # coords(data.result.train, "best", best.method="closest.topleft", ret=c("threshold"))
# # coords(data.result.train, "best", best.method="closest.topleft", ret=c("accuracy"))
# # 
# # 
# # perf_AUC=performance(pred,"auc") #Calculate the AUC value
# # AUC=perf_AUC@y.values[[1]]
# # 
# # perf_ROC=performance(pred.train,"tpr","fpr") #plot the actual ROC curve
# # plot(perf_ROC, main="ROC plot")
# # text(0.5,0.5,paste("AUC = ",format(AUC, digits=5, scientific=FALSE)))
# 



###### With PCA
prin_comp <- prcomp(mydata, scale. = T)

std_dev <- prin_comp$sdev
pr_var <- std_dev^2
prop_varex <- pr_var/sum(pr_var)
PCA_row <- 1
PCA_limit <- 0.85

while(sum(prop_varex[1:PCA_row]) < PCA_limit)
  PCA_row  <- PCA_row + 1

plot(prop_varex[1:100], xlab = "Principal Component", ylab = "Variance",type = "b")
lines(c(PCA_row,PCA_row),c(-1,prop_varex.bc_norm[PCA_row]), col='blue')
plot(cumsum(prop_varex[1:100]), xlab = "Principal Component", ylab = "Cumulative Variance", type = "b")
lines(c(-10,PCA_row),c(sum(prop_varex[1:PCA_row]),sum(prop_varex[1:PCA_row])), col='blue')



mydata <- data.frame(prin_comp$x[,1:PCA_row], Class_Label = mydata$Class_Label)

mydata <- sample(mydata[,1:3])

###Option for Resolve Class imbalance in Training Set Only
mydata <- mydata[sample(nrow(mydata)),]
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
# updata <- updata[sample(nrow(updata)),]
# 
# samples <- sample(NROW(updata), NROW(updata) * .632)
# data.train <- updata[samples, ]
# data.test <- updata[-samples, ]
###---

data.model <- randomForest(x = data.train[,-ncol(data.train)], y = data.train$Class_Label, ntree = 1500, importance=TRUE, proximity=TRUE
                           , trControl=myTrainingControl)
# data.model <- train(Class_Label ~., data.train, method = "rf", trControl = myTrainingControl , ntree = 1500, importance=TRUE, proximity=TRUE, nodesize = 1)

# pred <- predict(data.model, data.test)
# predictions=as.vector(data.model$votes[,2])
# pred = prediction(predictions,data.test$Class_Label)
pred.test <- predict(data.model, newdata=data.test, type='prob')
data.result.test <- roc(data.test$Class_Label, pred.test[,1])



plot(data.result.test, col = "steelblue2", print.thres="best", print.thres.best.method="closest.topleft", 
     main =  "ROC of Random Forest with Adaptive-SMOTE Sampling Technique",
     xaxt="n",yaxt="n",xlab="Specificity (1 - False Positive Rates)\nAUC = 99.14%", ylab = "Sensitivity (True Positive Rate)")

axis(2, at = seq(0,1,by=.2), labels = seq(0,1, by=.2), tick = TRUE)
axis(1, at = seq(1,0,by=-.2), labels = seq(1,0, by=-.2), tick = TRUE)


