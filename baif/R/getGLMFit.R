getGLMFit <- function(x,y,Cv=TRUE,Ld,Alp,dist){
  rt_fit <- if(Cv){
    glmnet::cv.glmnet(x = x,y = y,family=dist, lambda=Ld, alpha=Alp);
  }else{
    glmnet::glmnet(x = x,y = y,family=dist, lambda=Ld, alpha=Alp);
  }
  return(rt_fit)
}
