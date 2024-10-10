library(R2jags);

library(ggplot2);library(reshape2)

library(FLCore);library(FLa4a)

library(mgcv)

library(foreach)

library(doParallel)

library(MCMCpack)

library(moments)
library(ggplotFL)

library(iterators) 
library(doParallel) 
library(tcltk)
#aux function to create subsample for bootstraping
#if there is an error in a4afit a new sample is created to reach sample size of 1000
boot_samp<-function(a4a.normal,BBs.sim,BBi.sim,dir){
  error<-0
  boot_stk<-BBs.sim
  #varC<-apply((log(catch.n(a4a.normal))-log(catch.n(BBs.sim)))^2,2,sum)/7
  varC<-predict(a4a.normal)$vmodel$catch
  #varC<-expand(varC,age=0:6)
  #boot_stk@catch.n<-catch.n(a4a.normal)*exp(rnorm(1,0,(varC)))
  boot_stk@catch.n<-rlnorm(1,log(catch.n(a4a.normal)),varC)
  
  
  boot_ind<-BBi.sim
  #varI1<-apply((log(index(a4a.normal)[[1]])-log(index(BBi.sim[[1]])))^2,2,sum)/6
  varI1<-predict(a4a.normal)$vmodel$pelgas
  #varI1<-expand(varI1,age=1:6)
  #boot_ind[[1]]@index<-index(a4a.normal)[[1]]*exp(rnorm(1,0,(varI1)))
  boot_ind[[1]]@index<-rlnorm(1,log(index(a4a.normal)[[1]]),varI1)
  
  #varI2<-apply((log(index(a4a.normal)[[2]])-log(index(BBi.sim[[2]])))^2,2,sum)
  varI2<-predict(a4a.normal)$vmodel$dpm
  #boot_ind[[2]]@index<-index(a4a.normal)[[2]]*exp(rnorm(1,0,sqrt(varI2)))
  boot_ind[[2]]@index<-rlnorm(1,log(index(a4a.normal)[[2]]),varI2)
  
  boot_fit <- a4aSCA(boot_stk,boot_ind,fmodel=fmodel,qmodel=qmodel,vmodel=vmodel,wkdir=dir)
  
  if(is.na(boot_fit@catch.n[1,1])){
    error<-1
    rm(boot_fit)
    boot_fit<-boot_samp(a4a.normal,BBs.sim,BBi.sim)[[1]]
  }
  
  return(list(boot_fit,error))
  
}


# simulate mvnorm with empirical T (fixes bug in mvrnorm)

mvrEmpT <- function(n, mu, Sigma, tol = 1e-6, empirical=TRUE){
  if(empirical){
    if(n>length(mu)){
      mm <- mvrnorm(n, mu, Sigma, tol=tol, empirical=T)
    } else {
      mm <- mvrnorm(length(mu)+1, mu, Sigma, tol=tol, empirical=T)
      mm <- mm[1:n,]	
    }
  } else {
    mm <- mvrnorm(n, mu, Sigma, tol=tol, empirical=FALSE)
  }
  
  # output with right dims for FLPar
  if(is(mm, "matrix")) t(mm) else (t(t(mm)))
  
}

#==================================================================== 
#    simulate  methods
#==================================================================== 



simulate_L<-  function(object, nsim = 1, seed = NULL, empirical=TRUE) {
            out <- object
            out @ pars <- simulate_L0(pars(object), nsim = nsim, seed = seed, empirical=empirical)
            
            # now get catch.n, stock.n, harvest and index
            preds <- predict(out)
            out @ harvest <- preds $ stkmodel $  harvest
            out @ stock.n <- out @ catch.n <- out @ harvest
            out @ stock.n[1,] <- preds $ stkmodel $  rec
            out @ stock.n[-1,1] <- preds $ stkmodel $ ny1[-1,]
            
            # plusgroup?
            dms <- dims(object)
            plusgrp <- !is.na(dms $ plusgroup) && dms $ plusgroup >= dms $ max
            
            # fill stock.n (waste space save time)
            stkn <- stock.n(out)
            Zs <- harvest(out) + m(out)
            for (a in 2:dms $ age) {
              stkn[a,-1] <- stkn[a-1, 1:(dms $ year-1)] * exp( - Zs[a-1, 1:(dms $ year-1)] )
            }
            # if plus group
            if (plusgrp) {
              for (y in 1:(dms $ year-1)) 
                stkn[a,y+1,] <- stkn[a,y+1,] + stkn[a, y,] * exp( - Zs[a, y,] )
            } 
            
            out@stock.n <- stkn
            
            # calculate catch
            zfrac <- harvest(out) / Zs * (1 - exp(-Zs))
            out @ catch.n <- zfrac * stkn
            
            # work out indices
            out @ index <- idxs <- preds $ qmodel
            
            for (i in seq(idxs)) {
              idx <- idxs[[i]]
              dnms <- dimnames(idx)	
              iages <- dnms$age
              iyears <- dnms$year
              when <- mean(range(qmodel(pars(out))[[i]])[c("startf", "endf")])
              # is it a biomass idx ?
              bioidx <- FALSE
              if(attr(index(object)[[i]], "FLIndexBiomass")) bioidx <- TRUE 
              # if biomass index use ages in range or all ages have to be accounted
              # WARNING: spagheti code
              if(bioidx){
                #		if(missing(stock)){
                #		    warning("Can't simulate the biomass index. Please provide FLStock to get stock weights.")
                #		    out @ index[[i]][] <- object@index[[i]][]
                #	    } else {
                rng <- attr(index(object)[[i]], "range")
                if(is.na(rng["min"])) iages <- dimnames(stkn)[[1]] else iages <- ac(rng["min"]:rng["max"])
                stk <- stkn[iages]*exp(-Zs[iages] * when)
                stk <- mcf(list(e1=stk, e2=wt(out)[iages]))
                stk$e2[] <- stk$e2[,,,,,1]
                stk <- quantSums(do.call("*", stk))[,iyears]
                out @ index[[i]] <- stk * out @ index[[i]]
                attr(out@index[[i]], "FLIndexBiomass") <- attr(object@index[[i]], "FLIndexBiomass")
                attr(out@index[[i]], "range") <- attr(object@index[[i]], "range")
                
                #	    }
                # or else it's a age based index
              } else {
                stk <- (stkn * exp(-Zs * when))[iages, iyears]
                out @ index[[i]] <- stk * out @ index[[i]]
                attr(out@index[[i]], "FLIndexBiomass") <- attr(object@index[[i]], "FLIndexBiomass")
                attr(out@index[[i]], "range") <- attr(object@index[[i]], "range")
              }
            }
            out
          }

simulate_L0<- function(object, nsim = 1, seed=NULL, empirical=TRUE) {    
            out <- object
            
            # get parameter estimates
            b <- coef(object)
            
            # get parameter variance matrices
            #V <- vcov(object)
            Vtemp<-getADMBCovariance('C:\\probak\\sim_sto')$cov
            dv<-dim(Vtemp)[1]
            l<-length(object@stkmodel@params)
            #19=Y+A=13+6
            reorder<-c(1:(l-19),(l-19+5+1):dv,(l-19+1):(l-19+5))
            V<-Vtemp[reorder,reorder]
            
            
            # iters and objects
            pitr <- dims(b[[1]])$iter
            vitr <- dim(vcov(object)[[1]])[3]
            mitr <- max(c(nsim, pitr, vitr))
            it <- seq(mitr)
            
            # sanity checks - must be 1 or n
            if(!((nsim==1 | nsim==mitr) & (pitr==1 | pitr==mitr) & (vitr==1 | vitr==mitr)))
              stop("The number of iters and simulations must be 1 or n.")
            
            # if there are iters in pars or vcov simulation must be done by iter
            if(pitr!=1 | vitr != 1){
              stop("Ther2 should be not iterations")
             
            } else {
              if(!is.null(seed)) set.seed(seed)
              parsim <- mvrEmpT(mitr, as.numeric(unlist(b)), V, empirical=empirical)
            }
            
            # add to object and return
            if(mitr > pitr){
            l<-length(object@stkmodel@params)
            
            out@stkmodel@params <- propagate(object@stkmodel@params, mitr)
            out @ stkmodel@params[] <- parsim[1:l,]
            
            for (i in 1:length(object@qmodel)){
            li<-length(object@qmodel[[i]]@params)
            out@qmodel[[i]]@params <- propagate(object@qmodel[[i]]@params, mitr)
            out @ qmodel[[i]]@params[] <- parsim[(l+1):(l+li),]
            l<-l+li
            }
            
            for (i in 1:length(object@vmodel)){
              li<-length(object@vmodel[[i]]@params)
              out@vmodel[[i]]@params <- propagate(object@vmodel[[i]]@params, mitr)
              out @ vmodel[[i]]@params[] <- parsim[(l+1):(l+li),]
              l<-l+li
            }}
            
            out
          }


###corrected residuals computing
residuals2<-function(object, stock, indices, ...) { 
  # object holder 
  lst <- list() 
  length(lst) <- length(indices) + 2	 
  # catch 
  lst[[1]] <- stdlogres2(catch.n(stock), catch.n(object),predict(object)$vmodel$catch) 
  # indices 
  idx <- index(object) 
  for(i in 1:length(indices)){ 
    lst[[i+1]] <- stdlogres2(index(indices[[i]]), idx[[i]],predict(object)$vmodel[[i+1]]) 
  } 
  lst[[length(lst)]] <- stdlogres(catch(stock), computeCatch(stock + object)) 
  # out 
  names(lst) <- c("catch.n", names(indices), "catch") 
  new("a4aFitResiduals", FLQuants(lst)) 
} 


#################residuals#################################
stdlogres2<- function(obs, fit, vpred, ...){ 
  flq <- log(obs/fit)   
  #res <- apply(flq, c(1,3:6), scale, center=FALSE) 
  res <- flq/(vpred) 
  dimnames(res) <- dimnames(flq) 
  as(res, "FLQuant") 
} 



load("C:\\Leire\\Sardina\\Ispra Whorkshop\\to send\\BB_6plus.RData")
cv.ln<-function(x){
  round(sqrt(exp(1/x)-1),10)
}
prec.ln<-function(cv){
  round( 1/log(cv^2 +1),10)
}

a.gam<-function(cv){
  round(1/cv^2,10)
}
cv.gam<-function(x){
  round( 1/sqrt(x),10)
}

save_datalist<-function(stk,idx,changetau=F,cv=1000){
  dat<-list()
  
  # number(year,age)
  
  # read from FLr objects
  
  #bioman index
  #dat$Idepm<-c(7.8,3.3,7.8,11,3.8,2.3,11,6.1,10,4.3,5.6,5.5,8.1)*1000000 #(Ptot/DF,DF~10^6)
  #catch at age
  dat$C<-t(stk@catch.n[drop=T])
  #weigth at age
  dat$WC<-t(stk@stock.wt[drop=T])
  #maturity
  dat$mat<-as.numeric(stk@mat[,1])
  #natural mortality
  dat$M<-as.numeric(stk@m[,1])
  #acgas index
  dat$Iac<-t(index(idx[[1]])[drop=T])
  dat$Idepm<-(index(idx[[2]])[drop=T])
  #proportion m and harvest before spowning  
  dat$m.spwn<-stk@m.spwn[drop=T][,1]
  dat$h.spwn<-t(stk@harvest.spwn[drop=T])
  
  
  #for piors
  dat$logmu.sr<-14
  dat$tau.sr<-prec.ln(2)
  if(changetau){dat$tau.sr<-prec.ln(cv)}
  dat$logmu.n1<-14
  dat$tau.n1<-prec.ln(2)
  if(changetau){dat$tau.n1<-prec.ln(cv)}
  
  dat$logmu.f<--0.7
  dat$tau.f<-prec.ln(2)
  if(changetau){dat$tau.f<-prec.ln(cv)}
  
  dat$logmu.s<--0.7
  dat$tau.s<-prec.ln(2)
  if(changetau){dat$tau.s<-prec.ln(cv)}
  
  #dat$a.s<-0
  #dat$b.s<-8
  
  dat$a.Iac<-2.5
  dat$b.Iac<-0.025
  if(changetau){dat$a.Iac<-a.gam(cv)}
  
  dat$a.Idepm<-2.5
  dat$b.Idepm<-0.025
  if(changetau){dat$a.Idepm<-a.gam(cv)}
  
  dat$a.C<-2.5
  dat$b.C<-0.0125
  if(changetau){dat$a.C<-a.gam(cv)}
  
  
  dat$logmu.qac<-0
  dat$tau.qac<-prec.ln(2)
  if(changetau){dat$tau.qac<-prec.ln(cv)}
  
  dat$logmu.qdepm<-4
  dat$tau.qdepm<-prec.ln(2)
  if(changetau){dat$tau.qdepm<-prec.ln(cv)}
  
  #additional
  dat$refage<-2 # age 1
  
  dat$Y<-dim(dat$C)[1]
  dat$A<-dim(dat$C)[2]
  
  dat$t_depm<-mean(range(idx[[2]])[c('startf','endf')])  #bioman year time
  dat$t_ac<-mean(range(idx[[1]])[c('startf','endf')]) #acgas year time
  
  return(dat)
}

dat<-save_datalist(BB.stk,BB.idx)

Bayes2a4a<-function(a4a.MC,jag.fit){
  p.to.save<-c('F','N','muC')
  for (i in 1:length(p.to.save)){
    array_list<-jag.fit$BUGSoutput$sims.list
    assign(p.to.save[i],get(p.to.save[i],array_list))
  }
  
  ####all outputs in a data frame#####################
  b_output<-jag.fit$BUGSoutput$sims.array[,1,]
  
  Bayes.fit<-a4a.MC
  ###Stock.n -N
  namesdim<-dimnames(Bayes.fit@stock.n)
  N<-N[,-14,]
  N<-aperm(N,c(3,2,1))
  dimnames(N)<-list(namesdim$age,namesdim$year,1:dim(b_output)[1])
  mN<-melt(N)
  mN<-cbind(mN,'unique','all','unique')
  colnames(mN)<-c('age','year','iter','data','unit','season','area')
  Bayes.fit@stock.n<-as.FLQuant(mN)
  
  ###harvest -F
  namesdim<-dimnames(Bayes.fit@harvest)
  F<-aperm(F,c(3,2,1))
  dimnames(F)<-list(namesdim$age,namesdim$year,1:dim(b_output)[1])
  mF<-melt(F)
  mF<-cbind(mF,'unique','all','unique')
  colnames(mF)<-c('age','year','iter','data','unit','season','area')
  Bayes.fit@harvest<-as.FLQuant(mF)
  
  
  ###Catch.n - muC
  namesdim<-dimnames(Bayes.fit@catch.n)
  muC<-aperm(muC,c(3,2,1))
  dimnames(muC)<-list(namesdim$age,namesdim$year,1:dim(b_output)[1])
  mmuC<-melt(muC)
  mmuC<-cbind(mmuC,'unique','all','unique')
  colnames(mmuC)<-c('age','year','iter','data','unit','season','area')
  Bayes.fit@catch.n<-as.FLQuant(mmuC)
  
  
  #units(Bayes.fit)<-units(a4a.MC)
  Bayes.fit@harvest@units<-'f'
  return(Bayes.fit)
}


##for prior tables

#rec
fcv<-function(x){sd(x)/mean(x)}

ptable<-function(dat,a4a.sto,a4a.mc,jagsout){
  priors1<-data.frame(Param=NA,pior=NA,priorm=NA,cv=NA,baym=NA,baycv=NA,mcm=NA,mccv=NA,stom=NA,stocv=NA,mcint=NA)
  for (i in 1:(dat$Y)){
    qnt<-round(qlnorm(c(0.5,0.05,0.95),dat$logmu.sr,1/sqrt(dat$tau.sr))/simp,2)
    qntj<-round(quantile(jagsout$N[,i,1],c(0.5,0.05,0.95))/simp,2)
    qntmc<-round(quantile(stock.n(a4a.mc)[1,i,drop=T],c(0.5,0.05,0.95))/simp,2)
    fila<-c(
      paste0('rec_',i+2001),
      paste0('ln(',dat$logmu.sr,',',round(1/sqrt(dat$tau.sr),2),')'),
      paste0(qnt[1],'(',qnt[2],',',qnt[3],')'),
      round(cv.ln(dat$tau.sr),0),
      paste0(qntj[1],'(',qntj[2],',',qntj[3],')'),
      round(median(stock.n(a4a.mc)[1,i,])/simp,2),
      round(median(stock.n(a4a.sto)[1,i,])/simp,2),
      
      
      round(fcv(jagsout$N[,i,1]),3),
      round(fcv(stock.n(a4a.mc)[1,i,]),3),
      round(fcv(stock.n(a4a.sto)[1,i,]),3),
      paste0(qntmc[1],'(',qntmc[2],',',qntmc[3],')')
    )
    priors1[i,]<-fila
  }
  
  #n1
  
  priors2<-data.frame(Param=NA,pior=NA,priorm=NA,cv=NA,baym=NA,baycv=NA,mcm=NA,mccv=NA,stom=NA,stocv=NA,mcint=NA)
  for (i in 1:(dat$A)){
    qnt<-round(qlnorm(c(0.5,0.05,0.95),dat$logmu.n1,1/sqrt(dat$tau.n1))/simp,2)
    qntj<-round(quantile(jagsout$N[,1,i],c(0.5,0.05,0.95))/simp,2)
    qntmc<-round(quantile(stock.n(a4a.mc)[i,1,drop=T],c(0.5,0.05,0.95))/simp,2)
    fila<-c(
      paste0('N1_',i-1),
      paste0('ln(',dat$logmu.n1,',',round(1/sqrt(dat$tau.n1),2),')'),
      paste0(qnt[1],'(',qnt[2],',',qnt[3],')'),
      round(cv.ln(dat$tau.n1),0),
      paste0(qntj[1],'(',qntj[2],',',qntj[3],')'),
      round(median(stock.n(a4a.mc)[i,1,])/simp,2),
      round(median(stock.n(a4a.sto)[i,1,])/simp,2),
      
      round(fcv(jagsout$N[,1,i]),3),
      round(fcv(stock.n(a4a.mc)[i,1,]),3),
      round(fcv(stock.n(a4a.sto)[i,1,]),3),
      paste0(qntmc[1],'(',qntmc[2],',',qntmc[3],')')
    )
    priors2[i,]<-fila
  }
  
  #fy
  priors4<-data.frame(Param=NA,pior=NA,priorm=NA,cv=NA,baym=NA,baycv=NA,mcm=NA,mccv=NA,stom=NA,stocv=NA,mcint=NA)
  for (i in 1:(dim(jagsout$Fparams)[2])){
    qnt<-round(qnorm(c(0.5,0.05,0.95),dat$zero[1],sqrt(dat$Smatrix[1,1])),2)
    qntj<-round(quantile(jagsout$Fparams[,i],c(0.5,0.05,0.95)),2)
    qntmc<-round(quantile(pars(a4a.mc)@stkmodel@params[i,drop=T],c(0.5,0.05,0.95)),2)
    fila<-c(
      paste0('F_',i),
      paste0('N(',dat$zero[1],',',round(sqrt(dat$Smatrix[1,1]),2),')'),
      paste0(qnt[1],'(',qnt[2],',',qnt[3],')'),
      round(sqrt(dat$Smatrix[1,1])/0.1,0),
      paste0(qntj[1],'(',qntj[2],',',qntj[3],')'),
      round(median(pars(a4a.mc)@stkmodel@params[i,drop=T]),2),
      round(median(pars(a4a.sto)@stkmodel@params[i,drop=T]),2),
      
      
      round(fcv(jagsout$Fparams[,i]),3),
      round(fcv(pars(a4a.mc)@stkmodel@params[i,drop=T]),2),
      round(fcv(pars(a4a.sto)@stkmodel@params[i,drop=T]),2),
      paste0(qntmc[1],'(',qntmc[2],',',qntmc[3],')')
    )
    priors4[i,]<-fila
  }
  
  
  
  
  
  #q
  
  priors3<-data.frame(Param=NA,pior=NA,priorm=NA,cv=NA,baym=NA,baycv=NA,mcm=NA,mccv=NA,stom=NA,stocv=NA,mcint=NA)
  #q1
  qnt<-round(qlnorm(c(0.5,0.05,0.95),dat$logmu.qac,1/sqrt(dat$tau.qac)),0)
  qntj<-round(quantile(jagsout$qac,c(0.5,0.05,0.95)),2)
  qntmc<-round(quantile(predict(a4a.mc)$qmodel$pelgas[1,1,drop=T],c(0.5,0.05,0.95)),2)
  fila<-c(
    paste0('q1'),
    paste0('ln(',dat$logmu.qac,',',round(1/sqrt(dat$tau.qac),2),')'),
    paste0(qnt[1],'(',qnt[2],',',qnt[3],')'),
    round(cv.ln(dat$tau.n1),0),
    
    paste0(qntj[1],'(',qntj[2],',',qntj[3],')'),
    round(median(predict(a4a.mc)$qmodel$pelgas[1,1,]),2),
    round(median(predict(a4a.sto)$qmodel$pelgas[1,1,]),2),
    round(fcv(jagsout$qac),3),
    
    
    round(fcv(predict(a4a.mc)$qmodel$pelgas[1,1,]),3),
    
    
    round(fcv(predict(a4a.sto)$qmodel$pelgas[1,1,]),3),
    paste0(qntmc[1],'(',qntmc[2],',',qntmc[3],')')
  )
  priors3[1,]<-fila
  #q2
  qnt<-round(qlnorm(c(0.5,0.05,0.95),dat$logmu.qdepm,1/sqrt(dat$tau.qdepm)),0)
  qntj<-round(quantile(jagsout$qdepm,c(0.5,0.05,0.95)),2)
  qntmc<-round(quantile(predict(a4a.mc)$qmodel$dpm[1,1,drop=T],c(0.5,0.05,0.95)),2)
  fila<-c(
    paste0('q2'),
    paste0('ln(',dat$logmu.qdepm,',',round(1/sqrt(dat$tau.qdepm),2),')'),
    paste0(qnt[1],'(',qnt[2],',',qnt[3],')'),
    round(cv.ln(dat$tau.n1),0),
    
    paste0(qntj[1],'(',qntj[2],',',qntj[3],')'),
    round(median(predict(a4a.mc)$qmodel$dpm[1,1,]),2),
    round(median(predict(a4a.sto)$qmodel$dpm[1,1,]),2),
    round(fcv(jagsout$qdepm),3),
    
    
    round(fcv(predict(a4a.mc)$qmodel$dpm[1,1,]),3),
    
    
    round(fcv(predict(a4a.sto)$qmodel$dpm[1,1,]),3),
    paste0(qntmc[1],'(',qntmc[2],',',qntmc[3],')')
  )
  priors3[2,]<-fila
  
  #sd I1
  qnt<-round(1/sqrt(qgamma(c(0.5,0.95,0.05),dat$a.Iac,dat$b.Iac)),2)
  qntj<-round(quantile(1/sqrt(jagsout$tau.Iac),c(0.5,0.05,0.95)),2)
  qntmc<-round(quantile(predict(a4a.mc)$vmodel$pelgas[1,1,drop=T],c(0.5,0.05,0.95)),2)
  fila<-c(
    paste0('sdI1'),
    paste0('invgamma(',dat$a.Iac,',',dat$b.Iac,')'),
    paste0(qnt[1],'(',qnt[2],',',qnt[3],')'),
    round(cv.gam(dat$a.Iac),0),
    
    paste0(qntj[1],'(',qntj[2],',',qntj[3],')'),
    round(median(predict(a4a.mc)$vmodel$pelgas[1,1,]),2),
    round(median(predict(a4a.sto)$vmodel$pelgas[1,1,]),2),
    round(fcv(1/sqrt(jagsout$tau.Iac)),3),
    
    
    round(fcv(predict(a4a.mc)$vmodel$pelgas[1,1,]),3),
    
    
    round(fcv(predict(a4a.sto)$vmodel$pelgas[1,1,]),3),
    paste0(qntmc[1],'(',qntmc[2],',',qntmc[3],')')
  )
  priors3[3,]<-fila
  
  #sd I2
  qnt<-round(1/sqrt(qgamma(c(0.5,0.95,0.05),dat$a.Idepm,dat$b.Idepm)),2)
  qntj<-round(quantile(1/sqrt(jagsout$tau.Idepm),c(0.5,0.05,0.95)),2)
  qntmc<-round(quantile(predict(a4a.mc)$vmodel$dpm[1,1,drop=T],c(0.5,0.05,0.95)),2)
  fila<-c(
    paste0('sdI2'),
    paste0('invgamma(',dat$a.Idepm,',',dat$b.Idepm,')'),
    paste0(qnt[1],'(',qnt[2],',',qnt[3],')'),
    round(cv.gam(dat$a.Idepm),0),
    
    paste0(qntj[1],'(',qntj[2],',',qntj[3],')'),
    round(median(predict(a4a.mc)$vmodel$dpm[1,1,]),2),
    round(median(predict(a4a.sto)$vmodel$dpm[1,1,]),2),
    round(fcv(1/sqrt(jagsout$tau.Idepm)),3),
    
    
    round(fcv(predict(a4a.mc)$vmodel$dpm[1,1,]),3),
    
    
    round(fcv(predict(a4a.sto)$vmodel$dpm[1,1,]),3),
    paste0(qntmc[1],'(',qntmc[2],',',qntmc[3],')')
  )
  priors3[4,]<-fila
  
  #sd IC
  qnt<-round(1/sqrt(qgamma(c(0.5,0.95,0.05),dat$a.C,dat$b.C)),2)
  qntj<-round(quantile(1/sqrt(jagsout$tau.C),c(0.5,0.05,0.95)),2)
  qntmc<-round(quantile(predict(a4a.mc)$vmodel$catch[1,1,drop=T],c(0.5,0.05,0.95)),2)
  fila<-c(
    paste0('sdCatch'),
    paste0('invgamma(',dat$a.C,',',dat$b.C,')'),
    paste0(qnt[1],'(',qnt[2],',',qnt[3],')'),
    round(cv.gam(dat$a.Idepm),0),
    
    paste0(qntj[1],'(',qntj[2],',',qntj[3],')'),
    round(median(predict(a4a.mc)$vmodel$catch[1,1,]),2),
    round(median(predict(a4a.sto)$vmodel$catch[1,1,]),2),
    round(fcv(1/sqrt(jagsout$tau.C)),3),
    
    
    round(fcv(predict(a4a.mc)$vmodel$catch[1,1,]),3),
    
    
    round(fcv(predict(a4a.sto)$vmodel$catch[1,1,]),3),
    paste0(qntmc[1],'(',qntmc[2],',',qntmc[3],')')
  )
  priors3[5,]<-fila
  
  priors<-rbind(priors3,priors4,priors1,priors2)
  colnames(priors)<-c('Param','pior','prior_med','cv','bay_med','mc_med','sto_med','bay_cv','mc_cv','sto_cv','a4a_mc')
  
  return(priors)
}


# 
# 
# require(stats)
# formula(PlantGrowth)         # check the default formula
# pg <- unstack(PlantGrowth)   # unstack according to this formula
# pg
# stack(pg)                    # now put it back together
# stack(pg, select = -ctrl)    # omitting one vector



plot_quants<-function(quant_vec){
  obj<-as.data.frame(quant_vec)
  ggplot(obj, aes(y = data, x = year,colour=qname,fill=qname)) + 
    stat_summary(fun.y = median,
                 fun.ymin = function(x) quantile(x,0.05), 
                 fun.ymax = function(x) quantile(x,0.95), 
                 geom = "ribbon",alpha=0.4) +
    stat_summary(fun.y = median,
                 geom = "line") +
    facet_grid( age~., scales = "free")}

####for correlations plots#############
varmatrix<-function(object){
  forcor<-cbind(
    t(object$ssb[drop=T]),
    
    t(object$F[1,,drop=T]),
    t(object$F[2,,drop=T]),
    t(object$F[3,,drop=T]),
    t(object$F[4,,drop=T]),
    t(object$F[5,,drop=T]),
    t(object$F[6,,drop=T]),
    t(object$F[7,,drop=T]),
    
    as.matrix((object$q1[1,1,drop=T])),
    as.matrix((object$q2[1,1,drop=T])),
    as.matrix((object$vc[1,1,drop=T])),
    as.matrix((object$vq1[1,1,drop=T])),
    as.matrix((object$vq2[1,1,drop=T]))
  )
  colnames(forcor)<-paste0(c(rep('ssb',13),rep(paste0('F',1:7),each=13),'q1','q2','vc','vI1','vI2'),colnames(forcor))
  return(forcor)
}

alliters<-function(object){
  alli<-varmatrix(object[[1]])
  for(n in 2:nsim){
    alli<-rbind(alli,varmatrix(object[[n]]))}
  return(alli)
}

#levelplot(cor(varmatrix(sim.boots[[3]])),col.regions=rainbow(1000)[600:1],cuts=100,at=seq(-1,1, length.out=600),main='Bootstrap')

          