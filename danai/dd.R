rm(list= ls())
getwd()
stkfit <- readRDS("github/a4asandbox/danai/STECF 23-09 - Annex_I/STECF 23-09 - Annex I/Stock_Objs_2309/HKE_1_5_6_7_final.rds")

idx <- readRDS("github/a4asandbox/danai/STECF 23-09 - Annex_I/STECF 23-09 - Annex I/Stock_Objs_2309/HKE_1_5_6_7_idx.rds")

stk <- readRDS("github/a4asandbox/danai/STECF 23-09 - Annex_I/STECF 23-09 - Annex I/Stock_Objs_2309/HKE_1_5_6_7_stk.rds")

fit <- readRDS("github/a4asandbox/danai/STECF 23-09 - Annex_I/STECF 23-09 - Annex I/Stock_Objs_2309/HKE_1_5_6_7_fit.rds")


stkf <- stk + fit


fmodel <-~1
qmodel <- list(~1)
srmodel <- ~1

fit01 <- sca(stk, idx, fmodel, qmodel, srmodel)
res01 <- residuals(fit01, stk, idx)

flq <- res01$catch.n

plot(res01, by = 'year', auxline = '')

library(randtests)

# Assume you have fitted a GAM model called `model`
# Extract Pearson residuals
residuals <- residuals(model, type = "pearson")

# Run the runs test
runs_test_result <- runs.test(res01)
print(runs_test_result)

library(ggplot2)
ggplot(as.data.frame(flq), aes(x = age, y = data)) + geom_point() + geom_smooth()+facet_grid(~year)

fit@pars@stkmodel@coefficients
coef(fit)
idx



stk <- readRDS("github/a4asandbox/danai/STECF 23-09 - Annex_I/STECF 23-09 - Annex I/Stock_Objs_2309/HKE_1_5_6_7_stk.rds")

idx <- readRDS("github/a4asandbox/danai/STECF 23-09 - Annex_I/STECF 23-09 - Annex I/Stock_Objs_2309/HKE_1_5_6_7_idx.rds")

fit <- readRDS("github/a4asandbox/danai/STECF 23-09 - Annex_I/STECF 23-09 - Annex I/Stock_Objs_2309/HKE_1_5_6_7_fit.rds")


fit01 <- sca(stk, idx, 
             fmodel = ~s(year, k =10) + s(age, k= 3),
             qmodel = list(~s(age, k =3)),
             srmodel = ~1)

res01 <- residuals(fit01, stk, idx)
plot(res01, by = 'age')

flq01 <- res01$catch.n

ggplot(as.data.frame(flq01), aes(x = age, y = data)) + geom_point() + 
  geom_smooth(method = 'lm')+facet_grid(~year)+ geom_hline(yintercept = 0, linetype = 'dotted')

fit02 <- sca(stk, idx, 
             fmodel = ~te(year,age, k = c(3,6)),
             qmodel = list(~s(age, k =3)),
             srmodel = ~1)

res02 <- residuals(fit02, stk, idx)
plot(res02, by = 'age', auxline = "")

flq02 <- res02$catch.n
flq03 <- res02$`1`
ggplot(as.data.frame(flq02), aes(x = age, y = data)) + geom_point() + 
  geom_smooth(method = 'lm')+facet_grid(~year) + geom_hline(yintercept = 0, linetype = 'dotted')
ggplot(as.data.frame(flq03), aes(x = age, y = data)) + geom_point() + 
  geom_smooth(method = 'lm')+facet_grid(~year) + geom_hline(yintercept = 0, linetype = 'dotted')

fitSumm(fit02)
fitSumm
# save(stk,fit,idx, file = "github/a4asandbox/danai/DPS_5_6_7.rda")


