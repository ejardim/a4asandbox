
library(FLa4a)
library(glue)

# read a stock assessment
stkdir <- "colin_local/stecf/stks"

stock <- "MUT_6"

files <- dir(stkdir, stock)

type <- gsub(glue("{stock}_|[.]rds"), "", files)

stk <- lapply(files, function(x) readRDS(file.path(stkdir, x)))
names(stk) <- type

if (TRUE) {
  data(ple4)
  data(ple4.indices)

  fmodel <- ~ factor(age) + factor(year)
  qmodel <- list(~ factor(age), ~ factor(age))
  ple4.fit <- sca(fmodel = fmodel, qmodel = qmodel, stock = ple4, indices = ple4.indices[2:3])

  stk <- list(
    stk = ple4,
    idx = ple4.indices[2:3],
    fit = ple4.fit
  )
}

wkdir <- "colin_local/runs/run"
unlink(wkdir, recursive = TRUE)

fit <-
  sca(
    stk$stk, stk$idx,
    fmodel = fmodel(stk$fit),
    qmodel = formula(qmodel(stk$fit)),
    n1model = n1model(stk$fit),
    vmodel = formula(vmodel(stk$fit)),
    wkdir = "colin_local/runs/run"
  )


get.line <- function(file, line, what = numeric()) {
  scan(file, what = what, nlines = 1, skip = line-1, quiet = TRUE)
}

# read model matrices
read.cfg <- function(model) {
  file <- file.path(wkdir, glue("{model}.cfg"))
  npars <- get.line(file, 3)
  ndata <- get.line(file, 5)
  data <- read.table(file, skip = 6, nrow = ndata)
  as.matrix(data)
}

modelmatrices <-
  list(
    fmodel = read.cfg("fmodel"),
    qmodel = read.cfg("qmodel"),
    n1model = read.cfg("ny1model"),
    vmodel = read.cfg("vmodel"),
    rmodel = read.cfg("rmodel")
  )

# read data
read.data <- function() {
  file <- file.path(wkdir, glue("a4a.dat"))

  ages <- get.line(file, 3)
  years <- get.line(file, 5)
  nsurveys <- get.line(file, 7)
  survey_minages <- get.line(file, 9)
  survey_maxages <- get.line(file, 11)
  survey_times <- get.line(file, 13)


  noobs <- get.line(file, 19)
  obs <- read.table(file, skip = 21, nrow = noobs)
  names(obs) <- c("fleet", "year", "age", "observation", "weights")

  naux <- get.line(file, noobs + 24)
  aux <- read.table(file, skip = noobs + 25, nrow = naux)
  names(aux) <- c("year", "age", "m", "m.spwn", "harvest.spwn", "mat.wt", "wt")

  list(
    ages = ages,
    years = years,
    nsurveys = nsurveys,
    survey_minages = survey_minages,
    survey_maxages = survey_maxages,
    survey_times = survey_times,
    obs = obs,
    aux = aux
  )
}

data <- read.data()

saveRDS(modelmatrices, file.path(wkdir, "modelmatrices.rds"))
saveRDS(fit, file.path(wkdir, "fit.rds"))
saveRDS(stk, file.path(wkdir, "stk.rds"))
saveRDS(data, file.path(wkdir, "data.rds"))
