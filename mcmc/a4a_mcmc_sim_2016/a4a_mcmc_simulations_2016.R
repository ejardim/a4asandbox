######################################################################
### parts of code (from 2016) used for paper:
# Uncertainty estimation and model selection in stock assessment
# models with non-parametric effects on ﬁshing mortality
### Leire Citores                         
### 2024 note: this code is not working now 
### needs updates on libraries and functions...
### This code runs the simulation for the case where f is
### modeled using factors fmodel<-~factor(age)+factor(year).
### The different shapes presented in the paper are listed
### at the end of the script
######################################################################

###########load sardine8abd stock data###########
load("BB_6plus.RData")
#################################################
source("sim_functions.R")

############SIMULATIONS###########################
## a4a fit using multivariate normal for uncertainty
## a4a fit using mcmc option
## Bayesian fit with Jags

index.var(BB.idx[[1]])<-1
qmodel<-list(~1,~1)
vmodel<-list(~1,~1,~1)

fmodel<-~factor(age)+factor(year)
modelname<-'Ffac'
a4a.te <- a4aSCA(BB.stk, BB.idx,fmodel=fmodel,qmodel=qmodel,vmodel=vmodel)

nsim<-100

sim.mc<-list()
sim.sto<-list()
catch.sim<-list()
idx1.sim<-list()
idx2.sim<-list()


sdC<-1/sqrt(prec.ln(0.25))
sdI<-1/sqrt(prec.ln(0.2))

for (i in 1:nsim){
  catch.sim[[i]]<-catch.n(a4a.te)*rlnorm(dat$A*dat$Y,0,sdC)
  idx1.sim[[i]]<-index(a4a.te)[[1]]*rlnorm(dat$Y*(dat$A-1),0,sdI)
  idx2.sim[[i]]<-index(a4a.te)[[2]]*rlnorm(dat$Y,0,sdI)
}


####################################
###### a4a fit MN   ################
####################################
library(foreach)
library(doParallel)

#start time
strt<-Sys.time()

#loop
i=1

#sim.sto<-foreach(i=1:nsim,.export=c("a4aSCA","simulate","simulate_L0","mvrEmpT","pars","predict",'ssb','SCAMCMC'), .verbose=T) %dopar% {
i=1
for(i in 1:nsim){
  BBs.sim<-BB.stk
  BBs.sim@catch.n<-catch.sim[[i]]
  #sim.stk[[i]]<-BBs.sim
  
  BBi.sim<-BB.idx
  BBi.sim[[1]]@index<-idx1.sim[[i]]
  BBi.sim[[2]]@index<-idx2.sim[[i]]
  #sim.idx[[i]]<-BBi.sim
  
  unlink('C:\\probak\\sim_sto',recursive = T, force = T)
  a4a.normal <- a4aSCA(BBs.sim,BBi.sim,fmodel=fmodel,qmodel=qmodel,vmodel=vmodel,wkdir='C:\\probak\\sim_sto')
  a4a.sto <-simulate_L(a4a.normal,10000)
  refpts(brp(FLBRP(BBs.sim+a4a.normal)))
  plot( refpts(brp(FLBRP(BBs.sim+a4a.normal))))
  #assign(paste('a4a.sto',i,sep=''),simulate(a4a.normal,2000))
  #a4a.sto<-simulate(a4a.normal,2000)
  pre<-predict(a4a.sto)
  ssb1<-ssb(BBs.sim+a4a.sto)
  
  sim.sto[[i]]<-list(F=pre$stkmodel$harvest,rec=pre$stkmodel$rec,q1=pre$qmodel$pelgas,q2=pre$qmodel$dpm,ssb1=ssb1,vc=pre$vmodel$catch,vq1=pre$vmodel$pelgas,vq2=pre$vmodel$dpm)
  
}


####################################
#######MCMC in a4a  ################
####################################
cl<-makeCluster(4)
registerDoParallel(cl)

#start time
strt<-Sys.time()

sim.mc<-foreach(i=1:nsim,.export=c("a4aSCA","predict",'ssb','SCAMCMC'), .verbose=T) %dopar% {
  
  BBs.sim<-BB.stk
  BBs.sim@catch.n<-catch.sim[[i]]
  #sim.stk[[i]]<-BBs.sim
  
  BBi.sim<-BB.idx
  BBi.sim[[1]]@index<-idx1.sim[[i]]
  BBi.sim[[2]]@index<-idx2.sim[[i]]
  #sim.idx[[i]]<-BBi.sim
  
  
  a4a.mc <- a4aSCA(BBs.sim,BBi.sim,fmodel=fmodel,qmodel=qmodel,vmodel=vmodel,fit='MCMC',mcmc=SCAMCMC(mcmc=100000,mcsave=100,mcprobe=0.3,mcrb=7))
  
  #assign(paste('a4a.mc',i,sep=''),a4a.mc)
  pre<-predict(a4a.mc)
  ssb1<-ssb(BBs.sim+a4a.mc)
  list(F=pre$stkmodel$harvest,rec=pre$stkmodel$rec,q1=pre$qmodel$pelgas,q2=pre$qmodel$dpm,ssb1=ssb1,vc=pre$vmodel$catch,vq1=pre$vmodel$pelgas,vq2=pre$vmodel$dpm)
  
}

print(Sys.time()-strt)
stopCluster(cl)

####################################
#######Bayes in jags################
####################################
modelname<-'bayfac_hightau'
##run jags##########
nsim<-100
load('C:\\Leire\\Paper1\\simu\\simulationsFfac100.RData')
source("C:\\Leire\\Paper1\\simu\\sim_functions.R")
setwd('C:\\Leire\\Sardina\\Ispra Whorkshop\\Ernesto\\Comp')
index.var(BB.idx[[1]])<-1
qmodel<-list(~1,~1)
vmodel<-list(~1,~1,~1)
fmodel<-~factor(age)+factor(year)
a4a.MC <- a4aSCA(BB.stk, BB.idx,fmodel=fmodel,qmodel=qmodel,vmodel=vmodel,fit='MCMC',mcmc=SCAMCMC(mcmc=1000,mcsave=1))

p.to.save=c('ssb','qac','qdepm','tau.C','tau.Iac','tau.Idepm','F','N','muC','muIac','muIdepm')
p.to.save=c('qac','qdepm','tau.C','tau.Iac','tau.Idepm','fy','sa','rec','N0')

sim.bay<-list()

cl<-makeCluster(4)
registerDoParallel(cl)

strt<-Sys.time()

sim.bay<-foreach(j=1:nsim,.export=c("jags","units","as.FLQuant","melt","index","predict",'harvest','ssb','rec','save_datalist','Bayes2a4a'), .verbose=T) %dopar% {
  #for (j in 1:3){
  BBs.sim<-BB.stk
  BBs.sim@catch.n<-catch.sim[[j]]
  #sim.stk[[i]]<-BBs.sim
  
  BBi.sim<-BB.idx
  BBi.sim[[1]]@index<-idx1.sim[[j]]
  BBi.sim[[2]]@index<-idx2.sim[[j]]
  #sim.idx[[i]]<-BBi.sim
  dat<-save_datalist(BBs.sim,BBi.sim,changetau=T,cv=1000)
  
  jag.fit<-jags(model.file='model0_fsep_q1_v111.jags',data=dat,parameters.to.save=p.to.save,n.chains=1, n.iter=150000,n.burnin=50000, n.thin=100, working.directory=NULL,DIC=F, progress.bar = "text")
  
  q1ba<-predict(a4a.MC)$qmodel$pelgas
  q1ba[,,]<-NA
  for( y in 1:dat$Y){
    for (i in 1:(dat$A-1)){
      q1ba[i,y,]<-(jag.fit$BUGSoutput$sims.list$qac)
    }}
  
  q2ba<-predict(a4a.MC)$qmodel$dpm
  q2ba[,,]<-NA
  for( y in 1:dat$Y){
    q2ba[1,y,]<-(jag.fit$BUGSoutput$sims.list$qdepm)}
  
  
  
  vcba<-predict(a4a.MC)$vmodel$catch
  vcba[,,]<-NA
  for( y in 1:dat$Y){
    for (i in 1:(dat$A)){
      vcba[i,y,]<-1/sqrt(jag.fit$BUGSoutput$sims.list$tau.C)
    }}
  
  vI1ba<-predict(a4a.MC)$vmodel$pelgas
  vI1ba[,,]<-NA
  for( y in 1:dat$Y){
    for (i in 1:(dat$A-1)){
      vI1ba[i,y,]<-1/sqrt(jag.fit$BUGSoutput$sims.list$tau.Iac)
    }}
  
  vI2ba<-predict(a4a.MC)$vmodel$dpm
  vI2ba[,,]<-NA
  for( y in 1:dat$Y){
    vI2ba[1,y,]<-1/sqrt(jag.fit$BUGSoutput$sims.list$tau.Idepm)}
  
  Bayes.fit<-Bayes2a4a(a4a.MC,jag.fit)
  
  #sim.bay[[j]]<-
  list(F=harvest(Bayes.fit),rec1=rec(BBs.sim+Bayes.fit),q1=q1ba,q2=q2ba,ssb1=ssb(BBs.sim+Bayes.fit),vc=vcba,vq1=vI1ba,vq2=vI2ba)
}


print(Sys.time()-strt)

############################################# save

filename0<-paste0('C:\\Leire\\Paper1\\simu\\simulationsF',modelname,'100.RData')
save(list=c('sim.mc','sim.sto','sim.bay','catch.sim','idx1.sim','idx2.sim','a4a.te','BB.idx','BB.stk'),file=filename0)

########################
####BOOTSTRAP###########
########################

modelname<-'bootfac'
############
nsim<-100
load('C:\\Leire\\Paper1\\simu\\simulationsFfac100.RData')
rm(sim.mc)
rm(sim.sto)
source("C:\\Leire\\Paper1\\simu\\sim_functions.R")


###Conditioned Parametric bootstraping (punt_etal_1993)

index.var(BB.idx[[1]])<-1
qmodel<-list(~1,~1)
vmodel<-list(~1,~1,~1)

fmodel<-~factor(age)+factor(year)


start<-Sys.time()

cl<-makeCluster(4,outfile="foreach_fac.txt")
registerDoParallel(cl)

set.seed(1)

boot_size<-1000
#clusterExport(cl, c("nsim")) 

sim.boots<-foreach(i=1:nsim,.packages=c("FLa4a","tcltk"), .verbose=T) %dopar% {
  errors<-0
  
  #for(i in 1:1){
  BBs.sim<-BB.stk
  BBs.sim@catch.n<-catch.sim[[i]]
  #sim.stk[[i]]<-BBs.sim
  
  BBi.sim<-BB.idx
  BBi.sim[[1]]@index<-idx1.sim[[i]]
  BBi.sim[[2]]@index<-idx2.sim[[i]]
  #sim.idx[[i]]<-BBi.sim
  
  a4a.normal <- a4aSCA(BBs.sim,BBi.sim,fmodel=fmodel,qmodel=qmodel,vmodel=vmodel)
  
  boot_fit<-boot_samp(a4a.normal,BBs.sim,BBi.sim)[[1]]
  
  pre<-predict(boot_fit)
  ssb1<-ssb(BBs.sim+boot_fit)
  
  sim.boot<-list(F=pre$stkmodel$harvest,rec=pre$stkmodel$rec,q1=pre$qmodel$pelgas,q2=pre$qmodel$dpm,ssb1=ssb1,vc=pre$vmodel$catch,vq1=pre$vmodel$pelgas,vq2=pre$vmodel$dpm)
  sim.boot<-sapply(sim.boot,expand,iter=1:boot_size)
  
  
  for (j in 2:boot_size){
    
    output<-boot_samp(a4a.normal,BBs.sim,BBi.sim)
    boot_fit<-output[[1]]
    
    errors<-errors+output[[2]]
    
    pre<-predict(boot_fit)
    ssb1<-ssb(BBs.sim+boot_fit)
    
    sim.boot$F[,,,,,j]<-pre$stkmodel$harvest
    sim.boot$rec[,,,,,j]<-pre$stkmodel$rec
    sim.boot$q1[,,,,,j]<-pre$qmodel$pelgas
    sim.boot$q2[,,,,,j]<-pre$qmodel$dpm
    sim.boot$ssb1[,,,,,j]<-ssb1
    sim.boot$vc[,,,,,j]<-pre$vmodel$catch
    sim.boot$vq1[,,,,,j]<-pre$vmodel$pelgas
    sim.boot$vq2[,,,,,j]<-pre$vmodel$dpm
  }
  
  print(paste('nsim=',i,Sys.time(),', errors:',errors))
  sim.boot
}


time<-Sys.time()-start
time


filename0<-paste0('C:\\Leire\\Paper1\\simu\\simulationsF',modelname,'100.RData')
save(list=c('sim.boots','catch.sim','idx1.sim','idx2.sim','a4a.te','BB.idx','BB.stk'),file=filename0)

rm(sim.boots)

stopCluster(cl)

####################### end simulations ###############################



### F shapes used for the paper

shape1 <- a4aSCA(BB.stk, BB.idx,fmodel=~factor(age)+factor(year),qmodel=qmodel,vmodel=vmodel)
shape2 <- a4aSCA(BB.stk, BB.idx,fmodel=~s(age,k=6)+s(year,k=8),qmodel=qmodel,vmodel=vmodel)
shape3 <- a4aSCA(BB.stk, BB.idx,fmodel=~s(age,k=3)+s(year,k=4),qmodel=qmodel,vmodel=vmodel)
shape4 <- a4aSCA(BB.stk, BB.idx,fmodel=~te(age,year,k=c(6,5)),qmodel=qmodel,vmodel=vmodel)
shape5 <- a4aSCA(BB.stk, BB.idx,fmodel=~te(age,year,k=c(3,3)),qmodel=qmodel,vmodel=vmodel)
