 
data(mut09)
 
mut09 <- dps567
mut09.idx <- dps567.idx


# mean model
fit01 <- sca(mut09, mut09.idx, fmod=~1, qmod=list(~1), srmod=~1, vmod=list(~1, ~1),  n1mod=~1)
res01 <- residuals(fit01, mut09, mut09.idx)
plot(res01)
plot(res01, auxline="spline", by="age")

# f age model
fit02 <- sca(mut09, mut09.idx, fmod=~factor(age), qmod=list(~1), srmod=~1, vmod=list(~1, ~1),  n1mod=~1)
res02 <- residuals(fit02, mut09, mut09.idx)
plot(res02)
plot(res02, auxline="spline", by="age")


# f age * year model
fit04 <- sca(mut09, mut09.idx, fmod=~factor(age) + factor(year), qmod=list(~1), srmod=~1, vmod=list(~1, ~1),  n1mod=~1)
res04 <- residuals(fit04, mut09, mut09.idx)
plot(res04)
plot(res04, auxline="spline", by="year")

# q age
fit05 <- sca(mut09, mut09.idx, fmod=~1, qmod=list(~factor(age)), srmod=~1, vmod=list(~1, ~1),  n1mod=~1)
res05 <- residuals(fit05, mut09, mut09.idx)
plot(res05)
plot(res05, auxline="spline", by="year")

# f and q age
fit06 <- sca(mut09, mut09.idx, fmod=~factor(age), qmod=list(~factor(age)), srmod=~1, vmod=list(~1, ~1),  n1mod=~1)
res06 <- residuals(fit06, mut09, mut09.idx)
plot(res06)
plot(res06, auxline="spline", by="year")


# f year model
fit03 <- sca(mut09, mut09.idx, fmod=~factor(age) + factor(year), qmod=list(~factor(age)), srmod=~1, vmod=list(~1, ~1),  n1mod=~1)
res03 <- residuals(fit03, mut09, mut09.idx)
plot(res03)
plot(res03, auxline="spline", by="age")




# q age
fit06 <- sca(mut09, mut09.idx, fmod=~factor(age) + s(year, k=8), qmod=list(~s(age, k=3)), srmod=~1, vmod=list(~1, ~1),  n1mod=~1)
res06 <- residuals(fit06, mut09, mut09.idx)
plot(res06)
plot(res06, auxline="spline", by="year")


# var age model

fit08 <- sca(dps567, dps567.idx, fmod=~factor(age) + factor(year), qmod=list(~factor(age)), srmod=~factor(year), vmod=list(~1, ~1), n1mod=~factor(age))

flqs <- FLQuants(Mod08=predict(fit08)$vmodel$catch[,"2021"], Mod07=predict(fit07)$vmodel$catch[,"2021"])
xyplot(data~age, data=flqs, group=qname, type="l", auto.key=T)

flqs <- FLQuants(Mod08=catch.n(dps567+simulate(fit08, nsim=500))[,"2021"], Mod07=catch.n(dps567+simulate(fit07, nsim=500))[,"2021"])
bwplot(data~qname|factor(age),  data=as.data.frame(flqs), scales="free", auto.key=T)

res07 <- residuals(fit07, mut09, mut09.idx)
plot(res07)

# f  + q age model
fit04 <- sca(mut09, mut09.idx, fmod=~factor(age), qmod=list(~factor(age)), srmod=~1, vmod=list(~1, ~1),  n1mod=~1)
plot(mut09+fit04)
wireframe(harvest(fit04))
wireframe(stock.n(fit04))
res04 <- residuals(fit04, mut09, mut09.idx)
plot(res04)


