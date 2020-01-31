fileDIR ="D:/Dissertation/Data/NewData/Cimalab Data/thyroid/";
NL <- 11;
i_i <- 1; j_j <- 1;
fname <- paste(fileDIR,stringr::str_c(i_i,"_",j_j,".jpg",collapse =""),sep="");

baif::estimateKernel(fname, 5,5, c(5));
