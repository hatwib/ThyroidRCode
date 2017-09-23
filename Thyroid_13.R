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

######### Start PCA

fileName = matrix("",nrow = 56, ncol = 1)
filterName = matrix(0,nrow = 56, ncol = 2)
k <- 1
for(i in seq(5, 17, by = 2)){
  for(j in seq(5, 12, by = 1)){
    fileName[k] = str_c("data_BSIF",i,i,j,"_thyroid.csv",collapse ="")
    filterName[k,] <- c(i,j)
    k = k + 1
  }
}


progress.bar <- create_progress_bar("text")
progress.bar$init(56)
ROC.results <- matrix(0, nrow = (56), ncol = 13)
max_var <- 0
for(f in 1:56){
  mydata <- read.csv(file = fileName[f], header=TRUE)
  
  rm_cols <- nearZeroVar(mydata)
  rm_idx <- length(rm_cols)
  
  
  if(rm_idx > 0)
    mydata <- mydata[,-c(rm_cols)]
  
  n_cols = ncol(mydata)
  if(n_cols <= 3)
    next
  
  pca_mydata <- mydata[,1:(ncol(mydata)-1)]
  prin_comp <- prcomp(pca_mydata, center = TRUE, scale. = TRUE)
  std_dev <- prin_comp$sdev
  pr_var <- std_dev^2
  prop_varex <- pr_var/sum(pr_var)
  tmp_var <- cumsum(prop_varex[1:2])[2]
  if(tmp_var > max_var){
    max_var <- tmp_var
    print(max_var)
  }
    

  ROC.results[f,] <- c(filterName[f,],tmp_var,cumsum(prop_varex[1:10]))
  
  progress.bar$step()
}

axis(2, at = seq(0,1,by=.2), labels = seq(0,1, by=.2), tick = TRUE)
axis(1, at = seq(1,0,by=-.2), labels = seq(1,0, by=-.2), tick = TRUE)
mycol <- c("red","green","blue")
plot(ROC.results[,13], type = 'l', xaxt="n", ylim = c(0.1,1) )
axis(1, at=1:56,labels = paste0(paste0(ROC.results[,1],"-"),ROC.results[,2]), tick = TRUE, cex.axis = 0.5)

for(t in 4:12)
  lines(ROC.results[,t], col=t)

paste0(paste0(ROC.results[,1],"-"),ROC.results[,2])

7-5 -> 9
7-8 -> 12
mydata <- read.csv(file = fileName[12], header=TRUE)
pca_mydata <- mydata[,1:(ncol(mydata)-1)]
mydata$Class_Label <- as.factor(mydata$Class_Label)

prin_comp <- prcomp(pca_mydata, center = TRUE, scale. = TRUE)
std_dev <- prin_comp$sdev
pr_var <- std_dev^2
prop_varex <- pr_var/sum(pr_var)


g <- ggbiplot(prin_comp, choices = c(2,1) , obs.scale = 1, var.scale = 1, groups = mydata$Class_Label, ellipse = TRUE, circle = TRUE, var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', legend.position = 'top')
print(g)


