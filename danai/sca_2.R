rm(list = ls())
library(FLCore)
library(FLa4a)
library(ggplotFL)

stk <- readRDS("github/a4asandbox/danai/STECF 23-09 - Annex_I/STECF 23-09 - Annex I/Stock_Objs_2309/HKE_1_5_6_7_stk.rds")

idx <- readRDS("github/a4asandbox/danai/STECF 23-09 - Annex_I/STECF 23-09 - Annex I/Stock_Objs_2309/HKE_1_5_6_7_idx.rds")
idx <- FLIndices(idx)

fit <- readRDS("github/a4asandbox/danai/STECF 23-09 - Annex_I/STECF 23-09 - Annex I/Stock_Objs_2309/HKE_1_5_6_7_fit.rds")

fmods = list(~s(year, k =10) + s(age, k= 3),
             ~s(year, k =9) + s(age, k= 3)
)

qmods = list(list(~factor(age)),
             list(~s(age, k =3)))



.sca <- function(stock, indices, fmodel = missing, qmodel = missing, 
                  srmodel = missing, n1model = missing, vmodel = missing)
{
  #-----------------------------------------------------------------
  # set models if missing
  if(missing(fmodel)) fmodel <- list(defaultFmod(stock))
  if(missing(qmodel)) qmodel <- list(defaultQmod(indices))
  if(missing(n1model)) n1model <- list(defaultN1mod(stock))
  if(missing(vmodel)) vmodel <- list(defaultVmod(stock, indices))
  if(missing(srmodel)) srmodel <- list(defaultSRmod(stock))
  
  dm <- expand.grid(1:length(fmodel), 
                    1:length(qmodel),
                    1:length(n1model),
                    1:length(vmodel),
                    1:length(srmodel))
  
  dm <- as.data.frame(t(as.matrix(dm)))
  
 fits <- lapply(dm, function(x){
    sca(stock, indices, 
        fmodel = fmodel[[x[1]]], 
        qmodel = qmodel[[x[2]]],
        n1model = n1model[[x[3]]],
        vmodel = vmodel[[x[4]]], 
        srmodel = srmodel[[x[5]]])
  })
 
 names(fits) <- paste0("fit", c(1:length(fits)))
 return(a4aFitSAs(fits))
}

myFits <- .sca(stk, idx, fmodel = fmods, qmodel = qmods)


