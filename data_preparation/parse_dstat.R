library(boot)
library(pbapply)

calc_dstat <- function(d) {
  dstat = (sum(d$abba) - sum(d$baba)) / (sum(d$abba) + sum(d$baba))
  return(dstat)
}

process_d <- function(x) {
  d = read.csv(x, stringsAsFactors = F)
  # do ancestral allele setting
  d = d[d[, 6] %in% c(0, 1), ]
  # if outgroup allele == 0, something wrong
  # need to address
  d$abba = (1 - d[, 3]) * d[, 4] * d[, 5] * (1 - d[, 6])
  d$baba = d[, 3] * (1 - d[, 4]) * d[, 5] * (1 - d[, 6])
  
  D = calc_dstat(d)
  num_inf = nrow(d[d$abba > 0 | d$baba > 0, ])

  boots = lapply(1:100, function(x) {d[sample(1:nrow(d), nrow(d), replace =  T), ]})
  dboots = unlist(lapply(boots, calc_dstat))
  
  D = calc_dstat(d)
  d_sd = sd(dboots)
  d_z = abs(D) / d_sd
  vals = c(names(d)[3:6], num_inf, D, d_sd, d_z)
  
  return(vals)
}

setwd("~/Dropbox/Encelia/analysis/pop_gen/dstat/")
files = list.files(".", pattern = "csv")
dvals = pblapply(files, process_d)


xx = as.data.frame(do.call(rbind, dvals), stringsAsFactors = F)
xx1$p = 2*pnorm(-abs(as.numeric(xx1$V8)))
