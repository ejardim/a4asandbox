
library(FLa4a)
library(glue)

# read a stock assessment
stkdir <- "colin_local/stecf/stks"

stock <- "MUT_6"

files <- dir(stkdir, stock)

type <- gsub(glue("{stock}_|[.]rds"), "", files)

stk <- lapply(files, function(x) readRDS(file.path(stkdir, x)))
names(stk) <- type

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

#
read.cfg <- function(model) {
  file <- file.path(wkdir, glue("{model}.cfg"))
  npars <- scan(file, what = integer(), nlines = 1, skip = 2, quiet = TRUE)
  ndata <- scan(file, what = integer(), nlines = 1, skip = 4, quiet = TRUE)
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

saveRDS(modelmatrices, file.path(wkdir, "modelmatrices.rds"))
saveRDS(fit, file.path(wkdir, "fit.rds"))
saveRDS(stk, file.path(wkdir, "stk.rds"))
