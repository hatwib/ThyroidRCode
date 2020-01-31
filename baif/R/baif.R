estimateKernel <-  function(filename, n_lasso = 10, n_alpha = 10,
                            filter = seq(5, 17, by = 2)){ #bitSize = matrix(seq(5, 12, by = 1), nrow = 1)
  all_docs <- c();
  objects = list(seq.int(9.1e-2,1, by = 1e-2),seq.int(9e-2,0, by = -1e-3))
  #registerDoParallel();
  tryCatch({
    if(file.exists(filename)){
      img = EBImage::readImage(filename);
      img_gray <- EBImage::Image(img,dim(img)[-3],'Grayscale');
      img_nmask = EBImage::bwlabel(img_gray);
      img_nmask = EBImage::rmObjects(img_nmask, objects, reenumerate = F)
      img_close <- EBImage::opening(img_nmask, EBImage::makeBrush(3, "box"));
      img_nmask = EBImage::bwlabel(img_close);
      sts <- base::order(EBImage::computeFeatures.shape(img_nmask, img,
                                         methods.ref=c("computeFeatures.haralick"))[,'s.area'], decreasing = T)[1];
      n_cols <- dim(img)[1]; n_rows <- dim(img)[2];
      t_b <- sort(which(img_nmask[as.integer(n_cols/2),]==sts));
      l_r <- sort(which(img_nmask[,as.integer(n_rows/2)]==sts));
      img_trim <- img_gray[c(min(l_r):max(l_r)),c(min(t_b):max(t_b))];
      img_target <- EBImage::clahe(img_gray, nx = 2, ny = 2, limit=16, bins = 128, keep.range = T);
      img_target <- img_target[c(min(l_r):max(l_r)),c(min(t_b):max(t_b))];
      #display(img_target, method = "raster");
      img_target <- max(img_target) - img_target;
      img_znorm <- znorm(img_target);

      H = dim(img_trim)[1]; W = dim(img_trim)[2];
      all_docs <- foreach(i_f = c(1:length(filter)), .combine = cbind) %do% {
        tryCatch({
          f <- filter[i_f]; n <- ((f-1)/2);
          n <- ((f-1)/2);
          c_sample <- sample(c(n:(H-n)));
          r_sample <- sample(c(n:(W-n)));
          m_sample <- foreach(i =1:200, .combine = rbind) %do% {
            cbind(sample(r_sample,1),sample(c_sample,1))
          }
          d <-c();
          d <- foreach(i_w = m_sample[,1], i_h = m_sample[,2], .combine = rbind) %do% {
            c(t(EBImage::imageData(img_trim)[c((i_h-n):(i_h+n)),c((i_w-n):(i_w+n))]),
              EBImage::imageData(img_znorm)[i_h,i_w]);
          }
          d <- data.frame(d);
          colnames(d) <- c(paste("Coef",c(1:(ncol(d)-1)), sep=""), "Class_Label"); rownames(d) <- NULL;
          tr_x <- model.matrix(Class_Label ~ ., d)[, -1];
          tr_y <- d$Class_Label;

          lasso = seq(0, 1, length.out = n_lasso);
          alpha = seq(1e-2, 1e2, length.out = n_alpha);
          # docs <- foreach(i_n = c(1:2), .combine = c) %dopar% {
          docs <- foreach(i_lasso = 1:n_lasso, .combine = c) %do% {
          #doc <- foreach(i_lasso = 1:n_lasso, .combine = c) %do% {
              #if(i_n == 1){
              fit_x <- getGLMFit(tr_x,tr_y,T,NULL,lasso[i_lasso], "gaussian");

              cf <- matrix(coef(fit_x,s=seq(min(fit_x$lambda),max(fit_x$lambda), length.out = n_lasso)), ncol = n_lasso)[-1,];
              return(list(cf));
              #}else{
              #  fit_x <- getGLMFit(tr_x,as.factor(ifelse(tr_y > 0,1,0)),F,alpha[i_alpha],lasso[i_lasso], "binomial");
              #}
              #cf <- coef(fit_x);
              #cf <- scale(cf,center = 0, scale=F)[-1];
              #cf <- coef(fit_x,s=seq(min(fit_x$lambda),max(fit_x$lambda), length.out = n_lasso))[-1];
            #}
            #return(list(doc));
          }
        return(list(docs));
      },error = function(err) {
        print(err);
      })
      }
    }

    #######
    # tot_docs <- foreach(i_1 = c(1:2),.combine = cbind) %do% {
    #   ctbl <- foreach(i_2 = c(1:7),  .combine = cbind) %do% {
    #     ctb <- foreach(i_3 = c(1:nl),.combine = cbind) %do% {
    #       ct <- foreach(n_n = c(1:15) , .combine = cbind) %do% {
    #         return(as.data.frame(acc_tot[[n_n]][[i_1]][[i_2]][[i_3]]));
    #       }
    #       colnames(ct) <- NULL;
    #       bp <- graphics::boxplot(t(ct), plot = F);
    #       for(l in bp$group)
    #         ct[l,which(ct[l,] %in% bp$out[bp$group==l])] = mean(t(ct[l,]));
    #       return(as.data.frame(rowMeans(ct)));
    #     }
    #     colnames(ctb) <- NULL;
    #     return(list(ctb));
    #   }
    #   return(list(ctbl));
    # }
    #######

    },error = function(err){
    print(err);
  })
  return(all_docs);
}


