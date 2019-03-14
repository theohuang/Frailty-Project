## Frailty Simulation Analysis
## Last updated: March 14, 2019

library(dplyr)
library(ggplot2)
library(data.table)
library(xtable)
library(pROC)
dir.res <- "/Users/Theo/Dropbox (Partners HealthCare)/Frailty Github/Simulations"


### Variances of 0.3
### Family-specific frailty distribution
for(i in 1:100){
  load(paste(dir.res, "/Family Frailty Distribution/Var03/simfamdistv03_", i, ".RData", sep = ""))
  if(i == 1){
    sim.fam.v03 <- res.post
  } else{
    sim.fam.v03 <- rbind(sim.fam.v03, res.post)
  }
}

### Getting the risk prediction results
for(i in 1:100){
  load(paste(dir.res, "/Risk Predictions/Var03/simresv03_", i, ".RData", sep = ""))
  if(i == 1){
    sim.res.v03 <- res
  } else{
    if(ncol(res) == 12){
      res <- res[, 1:11]
    }
    sim.res.v03 <- rbind(sim.res.v03, res)
  }
}


w.list.b <- seq(-1.5, 1.5, 0.5)
w.list.o <- seq(-1.5, 1.5, 0.5)
supp.w <- expand.grid(w.list.b, w.list.o)
names(supp.w) <- c("W.BC", "W.OC")
f.w <- rep(1, nrow(supp.w)) / nrow(supp.w) ## uniform


## Heatmaps of means and medians of posterior frailty probabilities
ggplot(cbind(supp.w, Mean = base::colMeans(sim.fam.v03[, 1:49], na.rm = TRUE)),
       aes(x = W.BC, y = W.OC, fill = Mean)) +
  geom_tile() + labs(x = "Breast Cancer Frailty",
                     y = "Ovarian Cancer Frailty") +
  scale_x_continuous(breaks = w.list.b) +
  scale_y_continuous(breaks = w.list.o)
ggplot(cbind(supp.w, Median = apply(sim.fam.v03[, 1:49], 2, median, na.rm = TRUE)),
       aes(x = W.BC, y = W.OC, fill = Median)) +
  geom_tile() + labs(x = "Breast Cancer Frailty",
                     y = "Ovarian Cancer Frailty") +
  scale_x_continuous(breaks = w.list.b) +
  scale_y_continuous(breaks = w.list.o)



load(paste(getwd(), "/Simulations/simdat_v03.RData", sep = ""))
sim.res.v03$BC.5 <- filter(fam.sim, FamID %in% sim.res.v03$FamID, isProband == 1)$AffectedBreast.fu
sim.res.v03.noc <- filter(sim.res.v03, FamID %in% filter(fam.sim, FamID %in% sim.res.v03$FamID, isProband == 1, AffectedOvary == 0)$FamID)


### Variances of 2
### Family-specific frailty distribution
for(i in 1:100){
  load(paste(dir.res, "/Family Frailty Distribution/Var2/simfamdistv2_", i, ".RData", sep = ""))
  if(i == 1){
    sim.fam.v2 <- res.post
  } else{
    sim.fam.v2 <- rbind(sim.fam.v2, res.post)
  }
}

### Getting the risk prediction results
for(i in 1:100){
  load(paste(dir.res, "/Risk Predictions/Var2/simresv2_", i, ".RData", sep = ""))
  if(i == 1){
    sim.res.v2 <- res
  } else{
    if(ncol(res) == 12){
      res <- res[, 1:11]
    }
    sim.res.v2 <- rbind(sim.res.v2, res)
  }
}

## Heatmaps of means and medians of posterior frailty probabilities
ggplot(cbind(supp.w, Mean = base::colMeans(sim.fam.v2[, 1:49], na.rm = TRUE)),
       aes(x = W.BC, y = W.OC, fill = Mean)) +
  geom_tile() + labs(x = "Breast Cancer Frailty",
                     y = "Ovarian Cancer Frailty") +
  scale_x_continuous(breaks = w.list.b) +
  scale_y_continuous(breaks = w.list.o)
ggplot(cbind(supp.w, Median = apply(sim.fam.v2[, 1:49], 2, median, na.rm = TRUE)),
       aes(x = W.BC, y = W.OC, fill = Median)) +
  geom_tile() + labs(x = "Breast Cancer Frailty",
                     y = "Ovarian Cancer Frailty") +
  scale_x_continuous(breaks = w.list.b) +
  scale_y_continuous(breaks = w.list.o)

load(paste(getwd(), "/Simulations/simdat_v2.RData", sep = ""))
sim.res.v2$BC.5 <- filter(fam.sim, FamID %in% sim.res.v2$FamID, isProband == 1)$AffectedBreast.fu
sim.res.v2.noc <- filter(sim.res.v2, FamID %in% filter(fam.sim, FamID %in% sim.res.v2$FamID, isProband == 1, AffectedOvary == 0)$FamID)

# 
# load(paste(getwd(), "/Simulations/simdat_v03.RData", sep = ""))
# nrow(filter(fam.sim, AffectedBreast >= 1, isProband != 1)) / length(unique(fam.sim$FamID))
# load(paste(getwd(), "/Simulations/simdat_v2.RData", sep = ""))
# nrow(filter(fam.sim, AffectedBreast >= 1, isProband != 1)) / length(unique(fam.sim$FamID))
# load(paste(getwd(), "/Simulations/simdat_dis.RData", sep = ""))
# nrow(filter(fam.sim, AffectedBreast >= 1, isProband != 1)) / length(unique(fam.sim$FamID))



### Discrete distribution
### Family-specific frailty distribution
for(i in 1:100){
  load(paste(dir.res, "/Family Frailty Distribution/Dis/simfamdistdis_", i, ".RData", sep = ""))
  if(i == 1){
    sim.fam.dis <- res.post
  } else{
    sim.fam.dis <- rbind(sim.fam.dis, res.post)
  }
}

### Getting the risk prediction results
for(i in 1:100){
  load(paste(dir.res, "/Risk Predictions/Dis/simresdis_", i, ".RData", sep = ""))
  if(i == 1){
    sim.res.dis <- res
  } else{
    if(ncol(res) == 12){
      res <- res[, 1:11]
    }
    sim.res.dis <- rbind(sim.res.dis, res)
  }
}


## Heatmaps of means and medians of posterior frailty probabilities
ggplot(cbind(supp.w, Mean = base::colMeans(sim.fam.dis[, 1:49], na.rm = TRUE)),
       aes(x = W.BC, y = W.OC, fill = Mean)) +
  geom_tile() + labs(x = "Breast Cancer Frailty",
                     y = "Ovarian Cancer Frailty") +
  scale_x_continuous(breaks = w.list.b) +
  scale_y_continuous(breaks = w.list.o)
ggplot(cbind(supp.w, Median = apply(sim.fam.dis[, 1:49], 2, median, na.rm = TRUE)),
       aes(x = W.BC, y = W.OC, fill = Median)) +
  geom_tile() + labs(x = "Breast Cancer Frailty",
                     y = "Ovarian Cancer Frailty") +
  scale_x_continuous(breaks = w.list.b) +
  scale_y_continuous(breaks = w.list.o)

load(paste(getwd(), "/Simulations/simdat_dis.RData", sep = ""))
sim.res.dis$BC.5 <- filter(fam.sim, FamID %in% sim.res.dis$FamID, isProband == 1)$AffectedBreast.fu
sim.res.dis.noc <- filter(sim.res.dis, FamID %in% filter(fam.sim, FamID %in% sim.res.dis$FamID, isProband == 1, AffectedOvary == 0)$FamID)


# mean.post.b <- function(post, w.list.b, supp.w){
#   res <- 0
#   for(j in 1:length(w.list.b)){
#     res <- res + sum(post[which(supp.w$W.BC == w.list.b[j])]) * w.list.b[j]
#   }
#   return(res)
# }
# mean.post.o <- function(post, w.list.o, supp.w){
#   res <- 0
#   for(j in 1:length(w.list.o)){
#     res <- res + sum(post[which(supp.w$W.OC == w.list.o[j])]) * w.list.o[j]
#   }
#   return(res)
# }
# 
# sim.res.dis.noc$wbc.mean <- apply(filter(sim.fam.dis, FamID %in% sim.res.dis.noc$FamID)[, 1:(ncol(sim.fam.dis) - 1)], 1,
#                                         mean.post.b, w.list.b = w.list.b, supp.w = supp.w)
# sim.res.dis.noc$woc.mean <- apply(filter(sim.fam.dis, FamID %in% sim.res.dis.noc$FamID)[, 1:(ncol(sim.fam.dis) - 1)], 1,
#                                         mean.post.o, w.list.o = w.list.o, supp.w = supp.w)
# sim.res.dis.noc$wbc.mean.dis <- sapply(sim.res.dis.noc$wbc.mean, function(x) w.list.b[which.min(abs(w.list.b - x))])
# sim.res.dis.noc$woc.mean.dis <- sapply(sim.res.dis.noc$woc.mean, function(x) w.list.o[which.min(abs(w.list.o - x))])
# ggplot(sim.res.dis.noc, aes(wbc.mean.dis, woc.mean.dis)) +
#   stat_bin2d(aes(fill = stat(count)), binwidth = c(0.5, 0.5))

# par(mfrow = c(3, 3))
# for(i in 1:length(w.list.b15)){
#   hist(post.prior.dis[, which(supp.w$W.BC == w.list.b[i] & supp.w$W.OC == 0)],
#        xlim = c(0, 0.1), ylim = c(0, 4500), xlab = "Posterior",
#        main = paste("BC Frailty = ", w.list.b[i], sep = ""))
#   abline(v = mean(post.all.dis[, which(supp.w$W.BC == w.list.b[i] & supp.w$W.OC == 0)]),
#          col = "blue", lwd = 3)
#   abline(v = median(post.all.dis[, which(supp.w$W.BC == w.list.b[i] & supp.w$W.OC == 0)]),
#          col = "red", lwd = 3)
# }
# par(mfrow = c(1, 1))



# # average number of relatives with BC
# dt <- data.table(fam.sim)
# dt <- dt[, length(which(.SD$isProband != 1 & .SD$AffectedBreast == 1)), by = FamID]
# mean(dt$V1)
# 
# dt <- data.table(fam.sim)
# dt <- dt[, nrow(.SD), by = FamID]
# mean(dt$V1)


getOE <- function(dat, out, pred){
  sum(select(dat, out)[, 1]) / sum(select(dat, pred)[, 1])
}

getAUC <- function(dat, out, pred){
  pROC::auc(select(dat, out)[, 1], select(dat, pred)[, 1], data = dat)
}

getRBS <- function(dat, out, pred){
  sqrt(mean((select(dat, out)[, 1] - select(dat, pred)[, 1])^2))
}


perf.boot.sim <- function(dat, out, pred, nboot){
  res.boot <- setNames(data.frame(matrix(NA, nboot, 3)), c("OE", "AUC", "rBS"))
  for(i in 1:nboot){
    dat.boot <- dat[sample(1:nrow(dat), nrow(dat), replace = TRUE), ]
    res.boot$OE[i] <- getOE(dat.boot, out, pred)
    res.boot$AUC[i] <- getAUC(dat.boot, out, pred)
    res.boot$rBS[i] <- getRBS(dat.boot, out, pred)
  }
  res <- c(getOE(dat, out, pred), quantile(res.boot$OE, 0.025), quantile(res.boot$OE, 0.975),
           getAUC(dat, out, pred), quantile(res.boot$AUC, 0.025), quantile(res.boot$AUC, 0.975),
           getRBS(dat, out, pred), quantile(res.boot$rBS, 0.025), quantile(res.boot$rBS, 0.975))
  return(res)
}


### overall results
# res.all <- setNames(data.frame(matrix(NA, 6, 3)), c("OE", "AUC", "rBS"))
# res.all$OE <- c(getOE(res.prior03.noc, "BC.5", "Prob.BC.5"),
#                 getOE(res.prior03.noc, "BC.5", "Prob.BC.5.nf"),
#                 getOE(res.prior2.noc, "BC.5", "Prob.BC.5"),
#                 getOE(res.prior2.noc, "BC.5", "Prob.BC.5.nf"),
#                 getOE(res.prior.dis.noc, "BC.5", "Prob.BC.5"),
#                 getOE(res.prior.dis.noc, "BC.5", "Prob.BC.5.nf"))
# res.all$AUC <- c(getAUC(res.prior03.noc, "BC.5", "Prob.BC.5"),
#                 getAUC(res.prior03.noc, "BC.5", "Prob.BC.5.nf"),
#                 getAUC(res.prior2.noc, "BC.5", "Prob.BC.5"),
#                 getAUC(res.prior2.noc, "BC.5", "Prob.BC.5.nf"),
#                 getAUC(res.prior.dis.noc, "BC.5", "Prob.BC.5"),
#                 getAUC(res.prior.dis.noc, "BC.5", "Prob.BC.5.nf"))
# res.all$rBS <- c(getRBS(res.prior03.noc, "BC.5", "Prob.BC.5"),
#                 getRBS(res.prior03.noc, "BC.5", "Prob.BC.5.nf"),
#                 getRBS(res.prior2.noc, "BC.5", "Prob.BC.5"),
#                 getRBS(res.prior2.noc, "BC.5", "Prob.BC.5.nf"),
#                 getRBS(res.prior.dis.noc, "BC.5", "Prob.BC.5"),
#                 getRBS(res.prior.dis.noc, "BC.5", "Prob.BC.5.nf"))

nboot <- 1000
res.all <- setNames(data.frame(matrix(NA, 6, 9)), c("OE", "OE_lo", "OE_hi",
                                                   "AUC", "AUC_lo", "AUC_hi",
                                                   "rBS", "rBS_lo", "rBS_hi"))

start <- Sys.time()
res.all[1, ] <- perf.boot.sim(sim.res.v03.noc, "BC.5", "Prob.BC.5", nboot)
res.all[2, ] <- perf.boot.sim(sim.res.v03.noc, "BC.5", "Prob.BC.5.nf", nboot)
res.all[3, ] <- perf.boot.sim(sim.res.v2.noc, "BC.5", "Prob.BC.5", nboot)
res.all[4, ] <- perf.boot.sim(sim.res.v2.noc, "BC.5", "Prob.BC.5.nf", nboot)
res.all[5, ] <- perf.boot.sim(sim.res.dis.noc, "BC.5", "Prob.BC.5", nboot)
res.all[6, ] <- perf.boot.sim(sim.res.dis.noc, "BC.5", "Prob.BC.5.nf", nboot)
print(difftime(Sys.time(), start, units = "secs"))

res.all.tab <- setNames(data.frame(matrix(NA, 6, 3)), c("OE", "AUC", "rBS"))
for(i in 1:6){
  res.all.tab$OE[i] <- paste(round(res.all$OE[i], 3), " (", round(res.all$OE_lo[i], 3),
                             ", ", round(res.all$OE_hi[i], 3), ")", sep = "")
  res.all.tab$AUC[i] <- paste(round(res.all$AUC[i], 3), " (", round(res.all$AUC_lo[i], 3),
                             ", ", round(res.all$AUC_hi[i], 3), ")", sep = "")
  res.all.tab$rBS[i] <- paste(round(res.all$rBS[i], 3), " (", round(res.all$rBS_lo[i], 3),
                             ", ", round(res.all$rBS_hi[i], 3), ")", sep = "")
}

xtable(res.all.tab)


save(sim.res.v03.noc, sim.res.v2.noc, sim.res.dis.noc,
     sim.fam.v03, sim.fam.v2, sim.fam.dis,
     res.all, res.all.tab,
     file = "Frailty_Sim_Analysis_Results.RData")





########## High risk families ##########

### Variances of 0.3
### Family-specific frailty distribution
for(i in 1:100){
  load(paste(dir.res, "/Family Frailty Distribution/High Risk/Var03/simfamdistv03hr_", i, ".RData", sep = ""))
  if(i == 1){
    sim.fam.v03.hr <- res.post
  } else{
    sim.fam.v03.hr <- rbind(sim.fam.v03.hr, res.post)
  }
}

### Getting the risk prediction results
for(i in 1:100){
  load(paste(dir.res, "/Risk Predictions/High Risk/Var03/simresv03hr_", i, ".RData", sep = ""))
  if(i == 1){
    sim.res.v03.hr <- res
  } else{
    if(ncol(res) == 12){
      res <- res[, 1:11]
    }
    sim.res.v03.hr <- rbind(sim.res.v03.hr, res)
  }
}

## Heatmaps of means and medians of posterior frailty probabilities
ggplot(cbind(supp.w, Mean = base::colMeans(sim.fam.v03.hr[, 1:49], na.rm = TRUE)),
       aes(x = W.BC, y = W.OC, fill = Mean)) +
  geom_tile() + labs(x = "Breast Cancer Frailty",
                     y = "Ovarian Cancer Frailty") +
  scale_x_continuous(breaks = w.list.b) +
  scale_y_continuous(breaks = w.list.o)
ggplot(cbind(supp.w, Median = apply(sim.fam.v03.hr[, 1:49], 2, median, na.rm = TRUE)),
       aes(x = W.BC, y = W.OC, fill = Median)) +
  geom_tile() + labs(x = "Breast Cancer Frailty",
                     y = "Ovarian Cancer Frailty") +
  scale_x_continuous(breaks = w.list.b) +
  scale_y_continuous(breaks = w.list.o)



load(paste(getwd(), "/Simulations/simdat_v03_hr.RData", sep = ""))
sim.res.v03.hr$BC.5 <- filter(fam.sim, FamID %in% sim.res.v03.hr$FamID, isProband == 1)$AffectedBreast.fu
sim.res.v03.noc.hr <- filter(sim.res.v03.hr, FamID %in% filter(fam.sim, FamID %in% sim.res.v03.hr$FamID, isProband == 1, AffectedOvary == 0)$FamID)



### Variances of 2
### Family-specific frailty distribution
for(i in 1:99){
  load(paste(dir.res, "/Family Frailty Distribution/High Risk/Var2/simfamdistv2hr_", i, ".RData", sep = ""))
  if(i == 1){
    sim.fam.v2.hr <- res.post
  } else{
    sim.fam.v2.hr <- rbind(sim.fam.v2.hr, res.post)
  }
}

### Getting the risk prediction results
for(i in 1:99){
  load(paste(dir.res, "/Risk Predictions/High Risk/Var2/simresv2hr_", i, ".RData", sep = ""))
  if(i == 1){
    sim.res.v2.hr <- res
  } else{
    if(ncol(res) == 12){
      res <- res[, 1:11]
    }
    sim.res.v2.hr <- rbind(sim.res.v2.hr, res)
  }
}


## Heatmaps of means and medians of posterior frailty probabilities
ggplot(cbind(supp.w, Mean = base::colMeans(sim.fam.v2.hr[, 1:49], na.rm = TRUE)),
       aes(x = W.BC, y = W.OC, fill = Mean)) +
  geom_tile() + labs(x = "Breast Cancer Frailty",
                     y = "Ovarian Cancer Frailty") +
  scale_x_continuous(breaks = w.list.b) +
  scale_y_continuous(breaks = w.list.o)
ggplot(cbind(supp.w, Median = apply(sim.fam.v2.hr[, 1:49], 2, median, na.rm = TRUE)),
       aes(x = W.BC, y = W.OC, fill = Median)) +
  geom_tile() + labs(x = "Breast Cancer Frailty",
                     y = "Ovarian Cancer Frailty") +
  scale_x_continuous(breaks = w.list.b) +
  scale_y_continuous(breaks = w.list.o)

load(paste(getwd(), "/Simulations/simdat_v2_hr.RData", sep = ""))
sim.res.v2.hr$BC.5 <- filter(fam.sim, FamID %in% sim.res.v2.hr$FamID, isProband == 1)$AffectedBreast.fu
sim.res.v2.noc.hr <- filter(sim.res.v2.hr, FamID %in% filter(fam.sim, FamID %in% sim.res.v2.hr$FamID, isProband == 1, AffectedOvary == 0)$FamID)


### Discrete distribution
### Family-specific frailty distribution
for(i in 1:99){
  load(paste(dir.res, "/Family Frailty Distribution/High Risk/Dis/simfamdistdishr_", i, ".RData", sep = ""))
  if(i == 1){
    sim.fam.dis.hr <- res.post
  } else{
    sim.fam.dis.hr <- rbind(sim.fam.dis.hr, res.post)
  }
}

### Getting the risk prediction results
for(i in 1:99){
  load(paste(dir.res, "/Risk Predictions/High Risk/Dis/simresdishr_", i, ".RData", sep = ""))
  if(i == 1){
    sim.res.dis.hr <- res
  } else{
    if(ncol(res) == 12){
      res <- res[, 1:11]
    }
    sim.res.dis.hr <- rbind(sim.res.dis.hr, res)
  }
}


## Heatmaps of means and medians of posterior frailty probabilities
ggplot(cbind(supp.w, Mean = base::colMeans(sim.fam.dis.hr[, 1:49], na.rm = TRUE)),
       aes(x = W.BC, y = W.OC, fill = Mean)) +
  geom_tile() + labs(x = "Breast Cancer Frailty",
                     y = "Ovarian Cancer Frailty") +
  scale_x_continuous(breaks = w.list.b) +
  scale_y_continuous(breaks = w.list.o)
ggplot(cbind(supp.w, Median = apply(sim.fam.dis.hr[, 1:49], 2, median, na.rm = TRUE)),
       aes(x = W.BC, y = W.OC, fill = Median)) +
  geom_tile() + labs(x = "Breast Cancer Frailty",
                     y = "Ovarian Cancer Frailty") +
  scale_x_continuous(breaks = w.list.b) +
  scale_y_continuous(breaks = w.list.o)

load(paste(getwd(), "/Simulations/simdat_dis_hr.RData", sep = ""))
sim.res.dis.hr$BC.5 <- filter(fam.sim, FamID %in% sim.res.dis.hr$FamID, isProband == 1)$AffectedBreast.fu
sim.res.dis.noc.hr <- filter(sim.res.dis.hr, FamID %in% filter(fam.sim, FamID %in% sim.res.dis.hr$FamID, isProband == 1, AffectedOvary == 0)$FamID)


nboot <- 1000
res.all.hr <- setNames(data.frame(matrix(NA, 6, 9)), c("OE", "OE_lo", "OE_hi",
                                                    "AUC", "AUC_lo", "AUC_hi",
                                                    "rBS", "rBS_lo", "rBS_hi"))

start <- Sys.time()
res.all.hr[1, ] <- perf.boot.sim(sim.res.v03.noc.hr, "BC.5", "Prob.BC.5", nboot)
res.all.hr[2, ] <- perf.boot.sim(sim.res.v03.noc.hr, "BC.5", "Prob.BC.5.nf", nboot)
res.all.hr[3, ] <- perf.boot.sim(sim.res.v2.noc.hr, "BC.5", "Prob.BC.5", nboot)
res.all.hr[4, ] <- perf.boot.sim(sim.res.v2.noc.hr, "BC.5", "Prob.BC.5.nf", nboot)
res.all.hr[5, ] <- perf.boot.sim(sim.res.dis.noc.hr, "BC.5", "Prob.BC.5", nboot)
res.all.hr[6, ] <- perf.boot.sim(sim.res.dis.noc.hr, "BC.5", "Prob.BC.5.nf", nboot)
print(difftime(Sys.time(), start, units = "secs"))

res.all.tab.hr <- setNames(data.frame(matrix(NA, 6, 3)), c("OE", "AUC", "rBS"))
for(i in 1:6){
  res.all.tab.hr$OE[i] <- paste(round(res.all.hr$OE[i], 3), " (", round(res.all.hr$OE_lo[i], 3),
                             ", ", round(res.all.hr$OE_hi[i], 3), ")", sep = "")
  res.all.tab.hr$AUC[i] <- paste(round(res.all.hr$AUC[i], 3), " (", round(res.all.hr$AUC_lo[i], 3),
                              ", ", round(res.all.hr$AUC_hi[i], 3), ")", sep = "")
  res.all.tab.hr$rBS[i] <- paste(round(res.all.hr$rBS[i], 3), " (", round(res.all.hr$rBS_lo[i], 3),
                              ", ", round(res.all.hr$rBS_hi[i], 3), ")", sep = "")
}

xtable(res.all.tab.hr)


save(sim.res.v03.noc.hr, sim.res.v2.noc.hr, sim.res.dis.noc.hr,
     sim.fam.v03.hr, sim.fam.v2.hr, sim.fam.dis.hr,
     res.all.hr, res.all.tab.hr,
     file = "Frailty_Sim_Analysis_Results_hr.RData")
