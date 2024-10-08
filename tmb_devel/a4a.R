library(TMB)

rundir <- "../colin_local/runs/run"
# get data
Xs <- readRDS(file.path(rundir, "modelmatrices.rds"))
fit <- readRDS(file.path(rundir, "fit.rds"))
data_a4a <- readRDS(file.path(rundir, "data.rds"))

# get a4a results
n.out <- read.table(file.path(rundir, "n.out"))
logn.out <- log(n.out)

data_old <- list(Y_old = rnorm(10) + 1:10, x_old = 1:10)
parameters_old <- list(a_old = 0, b_old = 0, logSigma_old = 0)



aux <- data_a4a$obs[c("year", "fleet", "age")]
aux <- as.matrix(aux)


data_new <- list(
  obs = data_a4a$obs$observation,
  aux = aux,
  minYear = data_a4a$years[1],
  minAge = data_a4a$ages[1],
  M = exp(t(matrix(data_a4a$aux$m, nrow = diff(data_a4a$ages) + 1, ncol = diff(data_a4a$years) + 1))),
  designF = Xs$fmodel,
  designN1 = Xs$n1model,
  designR = Xs$rmodel
)

# set up pars

pars <- coef(fit)$stkmodel[drop = TRUE]
fpars <- pars[grepl("fMod", names(pars))]
n1pars <- pars[grepl("n1Mod", names(pars))]
rpars <- pars[grepl("rMod", names(pars))]

#n_centering <- pars(fit)@stkmodel@centering[drop = TRUE]

parameters_new <- list(
  # Fpar = rep(0, ncol(Xs$fmodel))
  Fpar = unname(fpars),
  N1par = unname(n1pars),
  Rpar = unname(rpars)
)


data <- c(data_old, data_new)
parameters <- c(parameters_new, parameters_old)

# compile and load
compile("a4a.cpp")
dyn.load(dynlib("a4a"))
obj <- MakeADFun(data, parameters, DLL = "a4a")

opt <- do.call("optim", obj)

obj$report()$logF
obj$report()$logN1
obj$report()$logR
obj$report()$logN

# check F
range(
  t(obj$report()$logF) - log(harvest(fit))[drop = TRUE]
)

# check N1
range(
  obj$report()$logN1 - logn.out[1, -1]
)

# check R
range(
  obj$report()$logR - logn.out[, 1]
)

# check N
range(
  obj$report()$logN - logn.out
)


saveRDS(opt, "a4a_opt.rds")
