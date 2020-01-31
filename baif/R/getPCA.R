getPCA <- function(ft){
  tmp_ang <- NULL;
  ft <- as.data.frame(lapply(ft, znorm));
  for(level in 0:5){
    if(level==0){
      tmp <- ft %>%
        base::subset(select = c(stringr::str_subset(colnames(ft),".s[1-5]", negate = T) )) %>%
        base::t();
    }else{
      tmp <- ft %>%
        base::subset(select = c(stringr::str_subset(colnames(ft),paste(".s",level,sep="")) )) %>%
        base::t();
    }
    #rownames(tmp) <- c("h.asm","h.con","h.cor","h.var","h.idm","h.sav","h.sva","h.sen","h.ent","h.dva","h.den","h.f12","h.f13")
    rm_cols <- caret::nearZeroVar(tmp)
    if(length(rm_cols) > 0){
      tmp <- tmp[,-c(rm_cols)]
    }
    prin_comp <- stats::prcomp(tmp, scale. = T);
    tmp <- data.frame(t(prin_comp$x[,1]))
    if(is.null(tmp_ang)){
      tmp_ang <- tmp;
    }else{
      tmp_ang <- cbind(tmp_ang,tmp);
    }
  }
  return(tmp_ang);
}
