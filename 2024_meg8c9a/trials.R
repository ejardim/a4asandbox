# trial 2
fmod <- ~s(age, k=4) + s(year, k=10) + te(age, year, k=c(3,5))  
srmod <- ~factor(year) #this stock-recruitment model (srmod) is 'free'; i.e. there is no restriction on the estimated recruitment, based on the SSB. 
qmod <- list(~I(1/(1 + exp(-age))))
fit02 <- sca(stock,tun.sel[c(1)],fmodel=fmod,qmodel=qmod,srmodel=srmod)
fitSumm(fit02)
res02 <- residuals(fit02, stock, tun.sel[c(1)])
plotr(res02,auxline='none')
stk02 <- stock + fit02

xyplot(data~age,groups=year,data=predict(pars(fit02))$qmodel[1],type='b',ylab='Catchabilty',main='Q SPGFS',ylim=c(0.000,0.003))

sim02 <- simulate(fit02, 1000)
flqs <- FLQuants(sim=iterMedians(stock.n(sim02)), det=stock.n(fit02))
keylst <- list(points=FALSE, lines=TRUE, space="right")
xyplot(data~year|factor(age), groups=qname, data=flqs, type="l", main="Median simulations VS fit", scales=list(y=list(relation="free")), auto.key=keylst)
stks02 <- stock + sim02
plot(stks02) + facet_wrap(~qname,scales='free_y')

xyplot(data~year, groups=qname, data=FLQuants(f02=computeCatch(stk02)/catch(stock), f01=computeCatch(stk1)/catch(stock)), type="l")

fmod <- ~s(age, k=4) + s(year, k=10) + te(age, year, k=c(3,5))  
qmod <- list(~I(1/(1 + exp(-age))) + s(year, k=3))

back <- 5
retro02 <- split(1:back, 1:back)
retro02 <- lapply(retro02, function(x){
  yr <- range(stock)["maxyear"] - x
  stk <- window(stock, end=yr)
  tun <- window(tun.sel[c(1)], end=yr)
  stk + sca(stk, tun, fmodel=fmod, qmodel=qmod, srmodel=srmod)
})
retro02$"0" <- stock + sca(stock,tun.sel[c(1)],fmodel=fmod,qmodel=qmod,srmodel=srmod)
plot(FLStocks(retro02),col=1,lwd=1)+ facet_wrap(~qname,scales='free_y') 

yr <- range(stock)["maxyear"] - 2
stk. <- window(stock, end=yr)
tun. <- window(tun.sel[c(1)], end=yr)
fmod <- ~s(age, k=4) + s(year, k=9) + te(age, year, k=c(3,5))  
retro02$"2" <- stk. + sca(stk., tun., fmodel=fmod, qmodel=qmod, srmodel=srmod)

yr <- range(stock)["maxyear"] - 3
stk. <- window(stock, end=yr)
tun. <- window(tun.sel[c(1)], end=yr)
fmod <- ~s(age, k=4) + s(year, k=9) + te(age, year, k=c(3,5))  
retro02$"3" <- stk. + sca(stk., tun., fmodel=fmod, qmodel=qmod, srmodel=srmod)

yr <- range(stock)["maxyear"] - 4
stk. <- window(stock, end=yr)
tun. <- window(tun.sel[c(1)], end=yr)
fmod <- ~s(age, k=4) + s(year, k=8) + te(age, year, k=c(3,5))  
retro02$"4" <- stk. + sca(stk., tun., fmodel=fmod, qmodel=qmod, srmodel=srmod)

yr <- range(stock)["maxyear"] - 5
stk. <- window(stock, end=yr)
tun. <- window(tun.sel[c(1)], end=yr)
fmod <- ~s(age, k=4) + s(year, k=8) + te(age, year, k=c(3,5))  
retro02$"5" <- stk. + sca(stk., tun., fmodel=fmod, qmodel=qmod, srmodel=srmod)


m1 <- list(~I(1/(1 + exp(-age))))
m2 <- list(~age)
m3 <- list(~s(age, k = 4))
m4 <- list(~factor(age))
m5 <- list(~sqrt(age))


fitSumm(sca(stock,tun.sel[c(1)],fmodel=fmod,qmodel=m1,srmodel=srmod))
fitSumm(sca(stock,tun.sel[c(1)],fmodel=fmod,qmodel=m2,srmodel=srmod))
fitSumm(sca(stock,tun.sel[c(1)],fmodel=fmod,qmodel=m3,srmodel=srmod))
fitSumm(sca(stock,tun.sel[c(1)],fmodel=fmod,qmodel=m4,srmodel=srmod))
fitSumm(sca(stock,tun.sel[c(1)],fmodel=fmod,qmodel=m5,srmodel=srmod))
fitSumm(sca(stock,tun.sel[c(1)],fmodel=fmod,qmodel=m6,srmodel=srmod))


# trial 3
stk00 <- stock
idx00 <- tun.sel[1]
# cv of observed catches
varslt <- catch.n(stk00)
varslt[] <- 0.1
catch.n(stk00) <- FLQuantDistr(catch.n(stk00), varslt)
# cv of observed indices
varslt <- index(idx00[[1]])
varslt[] <- 0.8
index.var(idx00[[1]]) <- varslt

fmod <- ~s(age, k=4) + s(year, k=10) + te(age, year, k=c(3,5))  
srmod <- ~factor(year) #this stock-recruitment model (srmod) is 'free'; i.e. there is no restriction on the estimated recruitment, based on the SSB. 
qmod <- list(~I(1/(1 + exp(-age))))
fit03 <- sca(stk00, idx00,fmodel=fmod,qmodel=qmod,srmodel=srmod)
fitSumm(fit03)
res03 <- residuals(fit03, stk00, idx00[c(1)])
plotr(res03,auxline='none')
stk03 <- stock + fit03

xyplot(data~age,groups=year,data=predict(pars(fit03))$qmodel[1],type='b',ylab='Catchabilty',main='Q SPGFS',ylim=c(0.000,0.003))

sim03 <- simulate(fit03, 1000)
flqs <- FLQuants(sim=iterMedians(stock.n(sim03)), det=stock.n(fit03))
keylst <- list(points=FALSE, lines=TRUE, space="right")
xyplot(data~year|factor(age), groups=qname, data=flqs, type="l", main="Median simulations VS fit", scales=list(y=list(relation="free")), auto.key=keylst)
stks03 <- stock + sim03
plot(stks03) + facet_wrap(~qname,scales='free_y')

xyplot(data~year, groups=qname, data=FLQuants(f03=computeCatch(stk03)/catch(stock), f02=computeCatch(stk02)/catch(stock), f01=computeCatch(stk1)/catch(stock)), type="l", auto.key=TRUE)

# trial 4
fmod <- ~s(age, k=5) + s(year, k=10) + te(age, year, k=c(3,5))  
srmod <- ~factor(year) #this stock-recruitment model (srmod) is 'free'; i.e. there is no restriction on the estimated recruitment, based on the SSB. 
qmod <- list(~s(age, k=4, by=breakpts(year, 2013)))
fit04 <- sca(stock,tun.sel[c(1)],fmodel=fmod,qmodel=qmod,srmodel=srmod)
fitSumm(fit04)
res04 <- residuals(fit04, stock, tun.sel[c(1)])
plotr(res04,auxline='none')
stk04 <- stock + fit04

xyplot(data~age,groups=year,data=predict(pars(fit04))$qmodel[1],type='b',ylab='Catchabilty',main='Q SPGFS')

sim04 <- simulate(fit04, 1000)
flqs <- FLQuants(sim=iterMedians(stock.n(sim04)), det=stock.n(fit04))
keylst <- list(points=FALSE, lines=TRUE, space="right")
xyplot(data~year|factor(age), groups=qname, data=flqs, type="l", main="Median simulations VS fit", scales=list(y=list(relation="free")), auto.key=keylst)
stks04 <- stock + sim04
plot(stks04) + facet_wrap(~qname,scales='free_y')

xyplot(data~year, groups=qname, data=FLQuants(f04=computeCatch(stk04)/catch(stock), f03=computeCatch(stk03)/catch(stock), f02=computeCatch(stk02)/catch(stock), f01=computeCatch(stk1)/catch(stock)), type="l", auto.key=TRUE)
coef(pars(fit04)@qmodel)

