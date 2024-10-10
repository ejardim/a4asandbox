
# get data
x <- system2("Rscript", c("tmb_devel/data_create_dataset.R"), stdout = TRUE, stderr = TRUE)
cat(x, sep = "\n")

setwd("tmb_devel")

# run linreg in fresh R session
x <- system2("Rscript", c("a4a.R"), stdout = TRUE, stderr = TRUE)
cat(x, sep = "\n")

setwd("..")
