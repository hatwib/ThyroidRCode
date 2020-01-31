znorm <- function(x){
  if(sd(x)==0){
    return(c(rep(0,length(x))))
  }else{
      return((x- mean(x)) /sd(x))
    }
  }
