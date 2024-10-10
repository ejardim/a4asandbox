rm(list = ls())
library(FLa4a)


stk <- readRDS("github/danai/STECF 23-09 - Annex_I/STECF 23-09 - Annex I/Stock_Objs_2309/HKE_1_5_6_7_stk.rds")

idx <- readRDS("github/danai/STECF 23-09 - Annex_I/STECF 23-09 - Annex I/Stock_Objs_2309/HKE_1_5_6_7_idx.rds")
idx <- FLIndices(idx)


fmodsk <- list()
for(i in 3:12) {
  fmodsk[[paste0(i)]] <- as.formula(paste0("~s(age, k =3)+s(year, k=",i,")"))
}

fitsk <- .sca(stk, idx, fmodel = fmodsk)

plotFitStats <- function(fits){
  gcv = lapply(fits,function(x) fitSumm(x)['gcv',])
  bic = lapply(fits, function(x) BIC(x))
  
  df <- data.frame(unlist(gcv), unlist(bic))
  
  df$fit <- as.numeric(gsub("fit", "",names(gcv)))
  names(df) <- c("GCV","BIC","fit")
  df <- df[complete.cases(df),]
  
  plot(df$fit, df$GCV, type = "b", col = "blue", 
       ylim = c(0.75*min(df$GCV), 1.25*max(df$GCV)), ylab = "", xlab = "fit")
  par(new = TRUE)
  plot(df$fit, df$BIC, type = "b", col = "red", 
       ylim = c(0.75*min(df$BIC), 1.25*max(df$BIC)), axes = FALSE, xlab = "", ylab = "")
  axis(4)                
  # mtext("y2 values", side = 4, line = 3)
  abline(v=df[min(df$GCV)==df$GCV,]$fit, col = "blue",lty = 2)
  abline(v=df[min(df$BIC)==df$BIC,]$fit, col = "red",lty = 2)
  legend("topleft", legend = c("GCV", "BIC"), col = c("blue", "red"), lty = 1)
  
  
}

plotFitStats(fitsk)
