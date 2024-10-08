
install.packages("msy", repos="http://flr-project.org/R")
install.packages("ggplotFL", repos="http://flr-project.org/R")
install.packages("FLa4a", repos="http://flr-project.org/R")
install.packages("FLasher", repos="http://flr-project.org/R")
install.packages("FLCore", repos="http://flr-project.org/R")
install.packages("FLa4a", repos="http://flr-project.org/R")
install.packages("FLEDA", repos="http://flr-project.org/R")

library(FLCore)
library(FLa4a)
library(FLEDA)
library(ggplotFL)
library(gridExtra)
library(icesAdvice)
library(plotly)
library(icesTAF)
library(FLasher)
library(msy)
library(ggplotFL)

setwd("D:/WorkingGroups/2024_WG/2024WGBIE/a4a/meg8c9a/test")

load("meg8c9ainputs.RData")
load("meg8c9aIndices.RData")


SavePlot<-function(plotname,width=6,height=4){
  file <- file.path(paste0('../test/meg8c9a_a4a_',plotname,'.png'))
  dev.print(png,file,width=width,height=height,units='in',res=600)
}


# specify submodels defined as follows:


#fmod, a formula object depicting the model for log fishing mortality at age
fmod <- ~ te(age, year, k = c(4,20))
#fmod <- ~factor(replace(age, age > 6, 6)) + factor(year)  

#srmod a formula object depicting the model for log recruitment
srmod <- ~factor(year) #this stock-recruitment model (srmod) is 'free'; i.e. there is no restriction on the estimated recruitment, based on the SSB. 

#qmod a list of formula objects depicting the models for log survey catchability at age
qmod <- list(~I(1/(1 + exp(-age))))


fit1 <- sca(stock,tun.sel[c(1)],fmodel=fmod,qmodel=qmod,srmodel=srmod)


save(stock,tun.sel,fit1,file='../test/meg8c9afitDef.Rdata')

submodels(fit1)



stk1 <- stock + fit1

plot(stk1)
save(stk1,file="../test/meg8c9a_stk1.RData")

name(stk1) <- 'Run1'

AIC(fit1)
# 843.0017
BIC(fit1)
#  1368.767

# hacked function from a4aFitresiduals-class
plotr<- function(x, y=missing, auxline="smooth",...){
  args <- list()
  args$data <- as.data.frame(x)
  args$x <- as.formula("data~year|factor(age)*qname")
  args$type=c("p", auxline)
  args$groups <- quote(qname)
  args$cex=0.3
  args$lwd=0
  args$ylab="standardized residuals"
  args$xlab=""
  args$panel=function(x,y,...){
    panel.abline(h=0, col.line="gray80")
    panel.xyplot(x,y,...)
  }
  args$par.settings=list(
    superpose.symbol=list(col=1, pch=19, cex=0.2), 
    superpose.line=list(col="gray75", lty=1, lwd=0, col=NA), 
    strip.background=list(col="gray90"), 
    strip.border=list(col="black"), 
    box.rectangle=list(col="gray90"))
  args$main="log residuals of catch and abundance indices by age"
  if(is(latticeExtra::useOuterStrips, "function")) latticeExtra::useOuterStrips(do.call("xyplot", args)) else do.call("xyplot", args)
}

res <- residuals(fit1, stock, tun.sel[c(1)])
plotr(res,auxline='none')
SavePlot('Residuals1',10,6)

res$catch.n[is.na(res$catch.n)] <- 0 #hack

res$SPGFS[is.na(res$SPGFS)] <- 0 #hack


bubbles(res)
SavePlot('Residuals2',10,6)

res <- residuals(fit1, stock, tun.sel[c(1)])
qqmath(res)
SavePlot('Residuals3')

wireframe(fit1@harvest)
a <- xyplot(data~age,groups=year,stk1@harvest,type='b',ylim=c(0,2),ylab='F',main='Fishing mortality')
a
SavePlot('F')

#


fitted <- predict(pars(fit1))
a1 <- xyplot(data~age,groups=year,data=fitted$qmodel[1],type='b',ylab='Catchabilty',main='Q SPGFS',ylim=c(0.000,0.003))
a1

grid.arrange(a,a1,ncol=2)
SavePlot('F and Q',6,3)


# simulate

fits <- simulate(fit1, 1000)
flqs <- FLQuants(sim=iterMedians(stock.n(fits)), det=stock.n(fit1))
keylst <- list(points=FALSE, lines=TRUE, space="right")
xyplot(data~year|factor(age), groups=qname, data=flqs, type="l", main="Median simulations VS fit", scales=list(y=list(relation="free")), auto.key=keylst)
stks <- stock + fits
plot(stks) + facet_wrap(~qname,scales='free_y')
SavePlot('summary0')


#summary table for the report
#F (Average F(last 3 years)), R geometric mean (1998:2020))

years <- stk1@range[4]:stk1@range[5]
nyears <- length(years)
GM <- exp(mean(log(window(stock.n(stk1)["1",], start=1998, end=2021)))) #since 1998
GM 
#[1] 4149.428
fsq <- fsq <- mean(fbar(stk1)[,nyears-2:0]) #F not scaled
#fsq <- fsq <- fbar(stk1)[,nyears]  #F SCALED TO THE LAST YEAR
fsq
#[1] 0.08991529


# 08/06/2022 in the ADG the 2sd is requested and it is included in all xxse variables for example: tsbse <- 2*apply(tsb(stks),2,sd).


#Summary table using coeficient interval of 1sd

tsb <- tsb(stk1)
tsbse <- apply(tsb(stks),2,sd)
ssb <- ssb(stk1)
ssbse <- apply(ssb(stks),2,sd)
ssbInt <- NA #get this from stf
catchobs <- stock@catch
discardsobs <- stock@discards
catch <- stk1@catch
recr <- stk1@stock.n[1,]
recrse <- apply(stks@stock.n[1,],2,sd)
fbar <- fbar(stk1)
fbarse <- apply(fbar(stks),2,sd)

#Summary table using coeficient interval of 2sd

tsb <- tsb(stk1)
tsb2sd <- 2*apply(tsb(stks),2,sd)
ssb <- ssb(stk1)
ssb2sd <- 2*apply(ssb(stks),2,sd)
ssbInt <- NA #get this from stf
catchobs <- stock@catch
discardsobs <- stock@discards
catch <- stk1@catch
recr <- stk1@stock.n[1,]
recr2sd <- 2*apply(stks@stock.n[1,],2,sd)
fbar <- fbar(stk1)
fbar2sd <- 2*apply(fbar(stks),2,sd)

sum1 <- data.frame(Year=c(years,paste0(max(years)+1,'*'))
                   ,Lan=c(catchobs-discardsobs,NA)
                   ,Dis=c(discardsobs,NA)
                   ,Cat=c(catchobs,NA)
                   ,CatEst=c(catch,NA)
                   ,Tsb=c(tsb,NA)
                   ,Ssb=c(ssb,ssbInt)
                   ,SsbCv=c(ssbse/ssb,NA)
                   ,Recr=c(recr,GM)
                   ,RecrCv=c(recrse/recr,NA)
                   ,Fbar=c(fbar,fsq)
                   ,FbarCv=c(fbarse/fbar,NA)
)

knitr::kable(sum1,row.names=F,digits=c(rep(3,12)))
write.csv(sum1,'../test/Summary_final.csv',row.names=F)

plot(years,tsb,type='l',ylim=c(0,120))
lines(years,ssb,type='l',col=2)
plot(years,ssb/tsb,type='l',ylim=c(0,1))

#summary for standard graphs template with coeficient interval of 1sd
sum2 <- data.frame(year=c(years,max(years+1))
                   ,recrlo=c(recr-recrse,NA)
                   ,recr=c(recr,GM)
                   ,recrhi=c(recr+recrse,NA)
                   ,tsblo=c(tsb-tsbse,NA)
                   ,tsb=c(tsb,NA)
                   ,tsbhi=c(tsb+tsbse,NA)
                   ,ssblo=c(ssb-ssbse,NA)
                   ,ssb=c(ssb,ssbInt)
                   ,ssbhi=c(ssb+ssbse,NA)
                   ,catch=c(catchobs,NA)
                   ,lan=c(catchobs-discardsobs,NA)
                   ,dis=c(discardsobs,NA)
                   ,ibc=NA
                   ,ur=NA
                   ,yssb=NA
                   ,flo=c(fbar-fbarse,NA)
                   ,f=c(fbar,fsq)
                   ,fhi=c(fbar+fbarse,NA)
)
knitr::kable(sum2,row.names=F,digits=c(rep(0,16),rep(3,3)))
write.csv(sum2,'../test/Summary_SAG_1sd.csv', row.names=F)

#summary for standard graphs template with coeficient interval of 2sd
sum3 <- data.frame(year=c(years,max(years+1))
                   ,recrlo=c(recr-recr2sd,NA)
                   ,recr=c(recr,GM)
                   ,recrhi=c(recr+recr2sd,NA)
                   ,tsblo=c(tsb-tsb2sd,NA)
                   ,tsb=c(tsb,NA)
                   ,tsbhi=c(tsb+tsb2sd,NA)
                   ,ssblo=c(ssb-ssb2sd,NA)
                   ,ssb=c(ssb,ssbInt)
                   ,ssbhi=c(ssb+ssb2sd,NA)
                   ,catch=c(catchobs,NA)
                   ,lan=c(catchobs-discardsobs,NA)
                   ,dis=c(discardsobs,NA)
                   ,ibc=NA
                   ,ur=NA
                   ,yssb=NA
                   ,flo=c(fbar-fbar2sd,NA)
                   ,f=c(fbar,fsq)
                   ,fhi=c(fbar+fbar2sd,NA)
)
knitr::kable(sum3,row.names=F,digits=c(rep(0,16),rep(3,3)))
write.csv(sum3,'../test/Summary_SAG_2sd.csv', row.names=F)

#sensitivity for template - not needed anymore?
ages <- stk1@range[1]:stk1@range[2]
stk3y <- window(stk1,start=max(years)-2)
p <- apply(landings.n(stk3y)/catch.n(stk3y),1,mean,na.rm=T)
p <- c(ifelse(is.na(p),0,p))
sen <- data.frame(Age=ages
                  ,M=c(apply(m(stk3y),1,mean))
                  ,Mat=c(apply(mat(stk3y),1,mean))
                  ,PF=c(apply(harvest.spwn(stk3y),1,mean))
                  ,PM=c(apply(m.spwn(stk3y),1,mean))
                  ,Sel=c(apply(harvest(stk3y),1,mean))*p
                  ,WeCa=c(apply(landings.wt(stk3y),1,mean))
                  ,Fd=c(apply(harvest(stk3y),1,mean))*(1-p)
                  ,WeCad=c(apply(discards.wt(stk3y),1,mean))
                  ,Fi=0
                  ,WeCai=0
)
knitr::kable(sen,row.names=F,digits=c(0,2,2,2,2,3,3,3,3,3,3))
write.csv(sen,'../test/Sen.csv',row.names=F)



plot(fit1, stock)
SavePlot('Fit1')
plot(fit1, tun.sel[1])
SavePlot('fit2')

fitSumm(fit1)


#retro
#fit1 <- sca(stock, tun.sel[c(1)], fmodel=fmod, qmodel=qmod, srmodel=srmod)
#stk1 <- stock + fit1
back <- 5
retro <- split(1:back, 1:back)
retro <- lapply(retro, function(x){
  yr <- range(stock)["maxyear"] - x
  stk <- window(stock, end=yr)
  tun <- window(tun.sel[c(1)], end=yr)
  stk + sca(stk, tun, fmodel=fmod, qmodel=qmod, srmodel=srmod)
})
retro$"0" <- stock + fit1
plot(FLStocks(retro),col=1,lwd=1)+ facet_wrap(~qname,scales='free_y') 

SavePlot('Retro',10,6.5)



Retro_F <- data.frame(Y0=c(fbar(retro$`0`))
                      ,Y1=c(fbar(retro$`1`),NA)
                      ,Y2=c(fbar(retro$`2`),NA,NA)
                      ,Y3=c(fbar(retro$`3`),NA,NA,NA)
                      ,Y4=c(fbar(retro$`4`),NA,NA,NA,NA)               
                      ,Y5=c(fbar(retro$`5`),NA,NA,NA,NA,NA)
                      
)

mohn(Retro_F,details=T) 
mohn(Retro_F,plot=T) 

#[1]  -0.4974443

error <- function(years,x,se){
  #  polygon(c(years,rev(years)),c(x-se,rev(c(x+se))),col='#0000FF20',border=NA)
  polygon(c(years,rev(years)),c(x-1.96*se,rev(c(x+1.96*se))),col='#0000FF20',border=NA)
}

pal <- c("#7570B3","#1B9E77","#E6AB02","#D95F02","#E7298A")
years <- stk1@range[4]:stk1@range[5]
nyears <-length(years)
fbar <- fbar(stk1)
fbarse <- apply(fbar(stks),2,sd)
ylim <- c(0,max(fbar+1.96*fbarse))
xlim <-range(years)
plot(NA,xlim=xlim,ylim=ylim,xlab='Year',ylab='F',main=paste0('Fbar ',stk1@range[6],'-',stk1@range[7]))
error(years,fbar,fbarse)
for(i in 1:5) lines(years[1:(nyears-i)],fbar(retro[[i]]),lwd=1,col=pal[i])
lines(years,fbar,lwd=2,col=4)
#abline(h=0.28,lty=2)
#abline(h=0.181,lty=3)
#abline(h=0.39,lty=3)
#lines(1986:2005,c(0.3486,0.3176,0.3423,0.3854,0.3835,0.3433,0.271,0.2041,0.2105,0.2664,0.3087,0.3065,0.2512,0.155,0.1261,0.182,0.1981,0.2328,0.2376,0.1622))
#legend('topleft',c('Fbar','Fmsy','Frange'),lty=c(1,2,3),col=c(4,1,1),lwd=c(2,1,1),bty='n',y.intersp=0.8)
SavePlot('Retro_F_boundaries',10,6.5)

Retro_SSB <- data.frame(Y0=c(ssb(retro$`0`))
                        ,Y1=c(ssb(retro$`1`),NA)
                        ,Y2=c(ssb(retro$`2`),NA,NA)
                        ,Y3=c(ssb(retro$`3`),NA,NA,NA)
                        ,Y4=c(ssb(retro$`4`),NA,NA,NA,NA)               
                        ,Y5=c(ssb(retro$`5`),NA,NA,NA,NA,NA)
)
mohn(Retro_SSB) 
mohn(Retro_SSB,details=T) 
mohn(Retro_SSB,plot=T) 


#[1] 1.030311

ssb <- ssb(stk1)/1000
ssbse <- apply(ssb(stks),2,sd)/1000
ylim <- c(0,max(ssb+1.96*ssbse))
xlim <-range(years)
plot(NA,xlim=xlim,ylim=ylim,xlab='Year',ylab='F',main=paste0('SSB ',stk1@range[6],'-',stk1@range[7]))
error(years,ssb,ssbse)
for(i in 1:5) lines(years[1:(nyears-i)],ssb(retro[[i]])/1000,lwd=1,col=pal[i])
lines(years,ssb,lwd=2,col=4)
#abline(h=0.28,lty=2)
#abline(h=0.181,lty=3)
#abline(h=0.39,lty=3)
#lines(1986:2005,c(0.3486,0.3176,0.3423,0.3854,0.3835,0.3433,0.271,0.2041,0.2105,0.2664,0.3087,0.3065,0.2512,0.155,0.1261,0.182,0.1981,0.2328,0.2376,0.1622))
#legend('topleft',c('Fbar','Fmsy','Frange'),lty=c(1,2,3),col=c(4,1,1),lwd=c(2,1,1),bty='n',y.intersp=0.8)
SavePlot('Retro_SSB_boundaries',10,6.5)

recr <- function(x) x@stock.n[1,]
Retro_R <- data.frame(Y0=c(recr(retro$`0`))
                      ,Y1=c(recr(retro$`1`),NA)
                      ,Y2=c(recr(retro$`2`),NA,NA)
                      ,Y3=c(recr(retro$`3`),NA,NA,NA)
                      ,Y4=c(recr(retro$`4`),NA,NA,NA,NA)               
                      ,Y5=c(recr(retro$`5`),NA,NA,NA,NA,NA)
)
mohn(Retro_R) 
mohn(Retro_R,details=T) 
mohn(Retro_R,plot=T) 



#$rho
#[1]  0.1287669


############### 


# do own summary plot with ref pts etc
error <- function(years,x,se){
  #  polygon(c(years,rev(years)),c(x-se,rev(c(x+se))),col='#0000FF20',border=NA)
  polygon(c(years,rev(years)),c(x-1.96*se,rev(c(x+1.96*se))),col='#0000FF20',border=NA)
}

pal <- c("#7570B3","#1B9E77","#E6AB02","#D95F02","#E7298A")

catch <- catch(stk1)/1000
#catchse <- apply(catch(stks),2,sd)/1000 
catchse <- 2*apply(catch(stks),2,sd)/1000 #2sd multiply apply *2
catchobs <- catch(stock)/1000
discardsobs<- discards(stock)/1000
landingsobs<- landings(stock)/1000
#dw <- discards.wt(stk1)
#dw[,as.character(1986:2002)] <- apply(window(dw,start=2003),1,mean)
discards <- discards(stk1)/1000
discardsse <- 2*apply(discards(stks),2,sd)/1000
#discards <- apply(catch.n(stk1)*dw*pdis[1:8,]/1000,2,sum)
plot(years,discardsobs,type='l')
lines(years,discards,col=2)

windows(10,8,10)
par(mfrow=c(2,2),mar=c(4.5,4,2,1))
xlim <- range(years)
ylim <- c(0,max(catch+1.96*catchse))
plot(NA,xlim=xlim,ylim=ylim,xlab='Year',ylab='Kt',main='Catch')
error(years,catch,catchse)
points(years,catchobs,cex=0.75)
for(i in 1:5) lines(years[1:(nyears-i)],retro[[i]]@catch/1000,lwd=1,col=pal[i])
lines(years,catch,lwd=2,col=4)
lines(years,discardsobs,lty=3,col=4)
lines(years,landingsobs,lty=3,col=7)
legend('topleft',c('Obs catch','Est catch','Discards','Landings'),lwd=c(NA,2,1,1),col=c(1,4,4,7),lty=c(NA,1,3,3),pch=c(1,NA,NA,NA),bty='n',y.intersp=0.8)

recr <- stock.n(stk1)[1,]/1000
recrse <- 2*apply(stock.n(stks)[1,],2,sd)/1000 #2sd multiply apply *2
ylim <- c(0,max(recr+1.96*recrse))
plot(NA,xlim=xlim,ylim=ylim,xlab='Year',ylab='Millions',main='Recruits age 0')
error(years,recr,recrse)
for(i in 1:5) lines(years[1:(nyears-i)],retro[[i]]@stock.n[1,]/1000,lwd=1,col=pal[i])
lines(years,recr,lwd=2,col=4)
abline(h=GM/1000,lty=2)
#lines(1986:2006,c(17137,11234,11069,13357,17553,23674,23224,21945,19216,17802,19940,24444,31931,37577,38408,37148,27821,25333,22081,22086,22086)/1000) # last xsa assessment 2006
legend('topleft',c('Recruits','GM'),lty=c(1,2),col=c(4,1),lwd=c(2,1),bty='n',y.intersp=0.8)


fbar <- fbar(stk1)
fbarse <- 2*apply(fbar(stks),2,sd)  #2sd
ylim <- c(0,max(fbar+1.96*fbarse))
plot(NA,xlim=xlim,ylim=ylim,xlab='Year',ylab='F',main=paste0('Fbar ',stk1@range[6],'-',stk1@range[7]))
error(years,fbar,fbarse)
for(i in 1:5) lines(years[1:(nyears-i)],fbar(retro[[i]]),lwd=1,col=pal[i])
lines(years,fbar,lwd=2,col=4)
abline(h=0.173,lty=2,col=4)
abline(h=0.45,lty=2)
abline(h=0.619,lty=3)
#lines(1986:2005,c(0.3486,0.3176,0.3423,0.3854,0.3835,0.3433,0.271,0.2041,0.2105,0.2664,0.3087,0.3065,0.2512,0.155,0.1261,0.182,0.1981,0.2328,0.2376,0.1622))
legend('bottomleft',c('Fmsy','Fpa','Flim'),lty=c(2,2,3),col=c(4,1,1),lwd=c(2,1,1),bty='n',y.intersp=0.8)

ssb <- ssb(stk1)/1000
ssbse <- 2*apply(ssb(stks),2,sd)/1000 #2sd
ylim <- c(0,max(ssb+1.96*ssbse))
plot(NA,xlim=xlim,ylim=ylim,xlab='Year',ylab='Kt',main='SSB')
error(years,ssb,ssbse)
for(i in 1:5) lines(years[1:(nyears-i)],ssb(retro[[i]])/1000,lwd=1,col=pal[i])
lines(years,ssb,col=4,lwd=2)
abline(h=725/1000,lty=2)
#lines(1986:2006,c(54219,48522,41531,38148,35287,38234,32447,30477,35398,45010,48558,45368,46144,45676,48715,52042,54714,56597,67360,78989,82348)/1000)
legend('topleft',c('SSB','Btrigger'),lty=c(1,2),col=c(4,1),lwd=c(2,1),bty='n',y.intersp=0.8)

par(mfrow=c(1,1))
plot.window(0:1,0:1)
box(col='#00000000')
legend('center',paste(-1:-5,'years'),col=pal,lty=1,bg='white',title='Retrospective',cex=0.8)

SavePlot('summary_updated',10,8)
dev.off()

write.csv(fit1@stock.n@.Data,'../test/stocknumbers.csv',row.names=T) 
write.csv(fit1@harvest@.Data,'../test/harvest.csv',row.names=T) 

#### 
par(mfrow=c(1,1))
kobe <- data.frame(Btrig=725,Fmsy=0.173,B=rev(c(ssb(stk1))),F=rev(c(fbar(stk1))))
with(kobe[1,],plot(B/Btrig,F/Fmsy,xlim=c(0,10),ylim=c(0,5.5),main='meg8c9a'))
rect(-1,-1,10,10,col='yellow',border=NA)
rect(-1,1,1,10,col='red',border=NA)
rect(1,-1,10,1,col='green',border=NA)
with(kobe[1,],text(B/Btrig,F/Fmsy,max(years),pos=3))
with(kobe[1,],points(B/Btrig,F/Fmsy,pch=16))
with(kobe,lines(B/Btrig,F/Fmsy))
SavePlot('Kobe',6,6)

#########################################################################################
