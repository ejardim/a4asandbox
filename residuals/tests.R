library(FLa4a)
data(ple4)
ple4.fit <- sca(ple4, ple4.indices)
 
data(mut09)
 
# mean model
fit01 <- sca(mut09, mut09.idx, fmod=~1, qmod=list(~1), srmod=~1, vmod=list(~1, ~1),  n1mod=~1)
plot(mut09+fit01)
wireframe(harvest(fit01))
wireframe(stock.n(fit01))
res01 <- residuals(fit01, mut09, mut09.idx)
plot(res01)

# f age model
fit02 <- sca(mut09, mut09.idx, fmod=~factor(age), qmod=list(~1), srmod=~1, vmod=list(~1, ~1),  n1mod=~1)
plot(mut09+fit02)
wireframe(harvest(fit02))
wireframe(stock.n(fit02))
res02 <- residuals(fit02, mut09, mut09.idx)
plot(res02)

# q age model
fit03 <- sca(mut09, mut09.idx, fmod=~factor(age), qmod=list(~factor(age)), srmod=~1, vmod=list(~1, ~1),  n1mod=~1)
plot(mut09+fit03)
wireframe(harvest(fit03))
wireframe(stock.n(fit03))
res03 <- residuals(fit03, mut09, mut09.idx)
plot(res03)

# sr model
fit04 <- sca(mut09, mut09.idx, fmod=~factor(age), qmod=list(~factor(age)), srmod=~factor(year), vmod=list(~1, ~1),  n1mod=~1)
plot(mut09+fit04)
wireframe(harvest(fit04))
wireframe(stock.n(fit04))
res04 <- residuals(fit04, mut09, mut09.idx)
plot(res04)

# f year model
fit05 <- sca(mut09, mut09.idx, fmod=~factor(year), qmod=list(~1), srmod=~1, vmod=list(~1, ~1),  n1mod=~1)
plot(mut09+fit05)
wireframe(harvest(fit05))
wireframe(stock.n(fit05))
res05 <- residuals(fit05, mut09, mut09.idx)
plot(res05)

# f age * year model
fit06 <- sca(mut09, mut09.idx, fmod=~factor(age) + factor(year), qmod=list(~1), srmod=~1, vmod=list(~1, ~1),  n1mod=~1)
plot(mut09+fit06)
wireframe(harvest(fit06))
wireframe(stock.n(fit06))
res06 <- residuals(fit06, mut09, mut09.idx)
plot(res06)

# var age model
fit07 <- sca(mut09, mut09.idx, fmod=~factor(age) + factor(year), qmod=list(~factor(age)), srmod=~1, vmod=list(~factor(age), ~factor(age)),  n1mod=~1)

flqs <- FLQuants(Mod06=predict(fit06)$vmodel$catch[,"2021"], Mod07=predict(fit07)$vmodel$catch[,"2021"])
xyplot(data~age, data=flqs, group=qname, type="l", auto.key=T)

flqs <- FLQuants(Mod06=catch.n(mut09+simulate(fit06, nsim=500))[,"2021"], Mod07=catch.n(mut09+simulate(fit07, nsim=500))[,"2021"])
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


