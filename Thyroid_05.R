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
  
  mydata <- n_mydata
  
  mydata$Class_Label <- ifelse(mydata$Class_Label == 1, 'Malignant', 'Benign')
  mydata$Class_Label <- as.factor(mydata$Class_Label)
  
  algorithms.list <- list("lvq","gbm","svmRadial")
  control <- trainControl(method="repeatedcv", number=10, repeats=3)
  ls()
  model.list <- mapply(function(fx) train(Class_Label~., data=mydata, method=fx, trControl=control), algorithms.list)
  
  results <- resamples(list(LVQ=model.list[[1]], GBM=model.list[[2]], SVM=model.list[[3]], RF=model.list[[4]]))
  summary(results)
  
}