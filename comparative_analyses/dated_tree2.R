# prep tree
t = read.tree("~/Dropbox/Encelia/analysis/phylogeny/astral/astral.miss0.6.brlen.tre")
t = root(t, c("Enceliopsis_covillei", "Xylorhiza_tortifolia"),
         resolve.root = T)
# our divergence dating says PAL - VEN = 1 myr
# smith et al tree says root height is 1.36 myr
rootenc1 = getMRCA(t, t$tip.label[grep("Encelia_", t$tip.label)])
calib = makeChronosCalib(t, node = rootenc1, age.min = 1.36)

lambda = c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0, 1, 10)
model = c("correlated", "relaxed", "discrete")

lnl = data.frame(lnl = rep(NA, length(lambda) * length(model)),
                 phiic = rep(NA, length(lambda) * length(model)),
                 lambda = rep(NA, length(lambda) * length(model)),
                 model = rep(NA, length(lambda) * length(model)))

for (j in 1:length(model)) {
  for (i in 1:length(lambda)) {
    
    if (j == 3) {
      ctrl <- chronos.control(nb.rate.cat = 1)
      tt = chronos(t, lambda[i], calibration = calib, model = model[j], control = ctrl)
    } else {
      tt = chronos(t, lambda[i], calibration = calib, model = model[j])
    }
    n = i + length(lambda) * (j - 1)
    
    lnl[n, "lnl"] = attr(tt, "ploglik")
    lnl[n, "lambda"] = lambda[i]
    lnl[n, "model"] = model[j]
    lnl[n, "phiic"] = attr(tt, "PHIIC")$PHIIC
  }
}

lnl2 = lnl %>% filter(lambda < 0.99)
aa = ggplot(lnl2, aes(lambda, phiic, color = model)) + 
  geom_point() + xlab(expression(lambda)) +
  ylab(expression(phi ~ "IC")) + ylim(0, 250)
save_plot("~/Dropbox/Encelia/manuscripts/Encelia_Phylogeny/figures/dating.png", aa)

ctrl <- chronos.control(nb.rate.cat = 1)
tt = chronos(t, 0, calibration = calib, model = model[j], control = ctrl)
write.tree(tt, "~/Dropbox/Encelia/analysis/phylogeny/astral/astral.miss0.6.brlen.dated.tre")
