## Frailty Simulation Analysis
## Last updated: March 10, 2020

library(dplyr)
library(ggplot2)
library(data.table)
library(xtable)
library(plyr)
library(pROC)
dir.sim <- "/Users/thuang/Dropbox (Partners HealthCare)/Frailty Project/Simulations"

theme_update(plot.title = element_text(hjust = 0.5))



###### Loading simulated data ######

### Variances of 0.3
### Family-specific frailty distribution
for(i in 1:200){
  load(paste(dir.sim, "/Family Frailty Distribution/Population/Var03/simfamdistv03_", i, ".RData", sep = ""))
  if(i == 1){
    sim.fam.v03 <- res.post
  } else{
    sim.fam.v03 <- rbind(sim.fam.v03, res.post)
  }
}

### Getting the risk prediction results
for(i in 1:200){
  load(paste(dir.sim, "/Risk Predictions/Population/Var03/simresv03_", i, ".RData", sep = ""))
  if(i == 1){
    sim.res.v03 <- res
  } else{
    if(ncol(res) == 12){
      res <- res[, 1:11]
    }
    sim.res.v03 <- rbind(sim.res.v03, res)
  }
}


### High risk
### Family-specific frailty distribution
for(i in 1:198){
  load(paste(dir.sim, "/Family Frailty Distribution/High Risk/Var03/simfamdistv03hr_", i, ".RData", sep = ""))
  if(i == 1){
    sim.fam.v03.hr <- res.post
  } else{
    sim.fam.v03.hr <- rbind(sim.fam.v03.hr, res.post)
  }
}

### Getting the risk prediction results
for(i in 1:198){
  load(paste(dir.sim, "/Risk Predictions/High Risk/Var03/simresv03hr_", i, ".RData", sep = ""))
  if(i == 1){
    sim.res.v03.hr <- res
  } else{
    if(ncol(res) == 12){
      res <- res[, 1:11]
    }
    sim.res.v03.hr <- rbind(sim.res.v03.hr, res)
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
                     y = "Ovarian Cancer Frailty",
                     title = "Bivariate Normal, Var = 0.3") +
  scale_x_continuous(breaks = w.list.b) +
  scale_y_continuous(breaks = w.list.o)
ggplot(cbind(supp.w, Median = apply(sim.fam.v03[, 1:49], 2, median, na.rm = TRUE)),
       aes(x = W.BC, y = W.OC, fill = Median)) +
  geom_tile() + labs(x = "Breast Cancer Frailty",
                     y = "Ovarian Cancer Frailty",
                     title = "Bivariate Normal, Var = 0.3") +
  scale_x_continuous(breaks = w.list.b) +
  scale_y_continuous(breaks = w.list.o)

## high risk
ggplot(cbind(supp.w, Mean = base::colMeans(sim.fam.v03.hr[, 1:49], na.rm = TRUE)),
       aes(x = W.BC, y = W.OC, fill = Mean)) +
  geom_tile() + labs(x = "Breast Cancer Frailty",
                     y = "Ovarian Cancer Frailty",
                     title = "Bivariate Normal, Var = 0.3, High-risk allele frequency") +
  scale_x_continuous(breaks = w.list.b) +
  scale_y_continuous(breaks = w.list.o)
ggplot(cbind(supp.w, Median = apply(sim.fam.v03.hr[, 1:49], 2, median, na.rm = TRUE)),
       aes(x = W.BC, y = W.OC, fill = Median)) +
  geom_tile() + labs(x = "Breast Cancer Frailty",
                     y = "Ovarian Cancer Frailty",
                     title = "Bivariate Normal, Var = 0.3, High-risk allele frequency") +
  scale_x_continuous(breaks = w.list.b) +
  scale_y_continuous(breaks = w.list.o)



load(paste(getwd(), "/Simulations/simdat_v03.RData", sep = ""))
sim.res.v03$BC.5 <- filter(fam.sim, FamID %in% sim.res.v03$FamID, isProband == 1)$AffectedBreast.fu
sim.res.v03.noc <- filter(sim.res.v03, FamID %in% filter(fam.sim, FamID %in% sim.res.v03$FamID, isProband == 1, AffectedOvary == 0)$FamID)

## family size
dt <- data.table(fam.sim)
dt <- dt[, nrow(.SD), by = FamID]
sim.res.v03.noc$FamSize <- merge(sim.res.v03.noc, dt, by = "FamID")$V1

load(paste(getwd(), "/Simulations/simdat_v03_hr.RData", sep = ""))
sim.res.v03.hr$BC.5 <- filter(fam.sim, FamID %in% sim.res.v03.hr$FamID, isProband == 1)$AffectedBreast.fu
sim.res.v03.hr.noc <- filter(sim.res.v03.hr, FamID %in% filter(fam.sim, FamID %in% sim.res.v03.hr$FamID, isProband == 1, AffectedOvary == 0)$FamID)

## family size
dt <- data.table(fam.sim)
dt <- dt[, nrow(.SD), by = FamID]
sim.res.v03.hr.noc$FamSize <- merge(sim.res.v03.hr.noc, dt, by = "FamID")$V1


### Variances of 2
### Family-specific frailty distribution
for(i in 1:200){
  load(paste(dir.sim, "/Family Frailty Distribution/Population/Var2/simfamdistv2_", i, ".RData", sep = ""))
  if(i == 1){
    sim.fam.v2 <- res.post
  } else{
    sim.fam.v2 <- rbind(sim.fam.v2, res.post)
  }
}

### Getting the risk prediction results
for(i in 1:200){
  load(paste(dir.sim, "/Risk Predictions/Population/Var2/simresv2_", i, ".RData", sep = ""))
  if(i == 1){
    sim.res.v2 <- res
  } else{
    if(ncol(res) == 12){
      res <- res[, 1:11]
    }
    sim.res.v2 <- rbind(sim.res.v2, res)
  }
}

### high risk
### using discrete uniform population-level frailty distribution
for(i in 1:198){
  load(paste(dir.sim, "/Family Frailty Distribution/High Risk/Var2/simfamdistv2hr_", i, ".RData", sep = ""))
  if(i == 1){
    sim.fam.v2.hr <- res.post
  } else{
    sim.fam.v2.hr <- rbind(sim.fam.v2.hr, res.post)
  }
}

### Getting the risk prediction results
for(i in 1:198){
  load(paste(dir.sim, "/Risk Predictions/High Risk/Var2/simresv2hr_", i, ".RData", sep = ""))
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
ggplot(cbind(supp.w, Mean = base::colMeans(sim.fam.v2[, 1:49], na.rm = TRUE)),
       aes(x = W.BC, y = W.OC, fill = Mean)) +
  geom_tile() + labs(x = "Breast Cancer Frailty",
                     y = "Ovarian Cancer Frailty",
                     title = "Bivariate normal, Var = 2") +
  scale_x_continuous(breaks = w.list.b) +
  scale_y_continuous(breaks = w.list.o)
ggplot(cbind(supp.w, Median = apply(sim.fam.v2[, 1:49], 2, median, na.rm = TRUE)),
       aes(x = W.BC, y = W.OC, fill = Median)) +
  geom_tile() + labs(x = "Breast Cancer Frailty",
                     y = "Ovarian Cancer Frailty",
                     title = "Bivariate normal, Var = 2") +
  scale_x_continuous(breaks = w.list.b) +
  scale_y_continuous(breaks = w.list.o)

ggplot(cbind(supp.w, Mean = base::colMeans(sim.fam.v2.hr[, 1:49], na.rm = TRUE)),
       aes(x = W.BC, y = W.OC, fill = Mean)) +
  geom_tile() + labs(x = "Breast Cancer Frailty",
                     y = "Ovarian Cancer Frailty",
                     title = "Bivariate normal, Var = 2, High-risk allele frequency") +
  scale_x_continuous(breaks = w.list.b) +
  scale_y_continuous(breaks = w.list.o)
ggplot(cbind(supp.w, Median = apply(sim.fam.v2.hr[, 1:49], 2, median, na.rm = TRUE)),
       aes(x = W.BC, y = W.OC, fill = Median)) +
  geom_tile() + labs(x = "Breast Cancer Frailty",
                     y = "Ovarian Cancer Frailty",
                     title = "Bivariate normal, Var = 2, High-risk allele frequency") +
  scale_x_continuous(breaks = w.list.b) +
  scale_y_continuous(breaks = w.list.o)

load(paste(getwd(), "/Simulations/simdat_v2.RData", sep = ""))
sim.res.v2$BC.5 <- filter(fam.sim, FamID %in% sim.res.v2$FamID, isProband == 1)$AffectedBreast.fu
sim.res.v2.noc <- filter(sim.res.v2, FamID %in% filter(fam.sim, FamID %in% sim.res.v2$FamID, isProband == 1, AffectedOvary == 0)$FamID)

## family size
dt <- data.table(fam.sim)
dt <- dt[, nrow(.SD), by = FamID]
sim.res.v2.noc$FamSize <- merge(sim.res.v2.noc, dt, by = "FamID")$V1


load(paste(getwd(), "/Simulations/simdat_v2_hr.RData", sep = ""))
sim.res.v2.hr$BC.5 <- filter(fam.sim, FamID %in% sim.res.v2.hr$FamID, isProband == 1)$AffectedBreast.fu
sim.res.v2.hr.noc <- filter(sim.res.v2.hr, FamID %in% filter(fam.sim, FamID %in% sim.res.v2.hr$FamID, isProband == 1, AffectedOvary == 0)$FamID)

## family size
dt <- data.table(fam.sim)
dt <- dt[, nrow(.SD), by = FamID]
sim.res.v2.hr.noc$FamSize <- merge(sim.res.v2.hr.noc, dt, by = "FamID")$V1




### Discrete distribution
### Family-specific frailty distribution
for(i in 1:200){
  load(paste(dir.sim, "/Family Frailty Distribution/Population/Dis/simfamdistdis_", i, ".RData", sep = ""))
  if(i == 1){
    sim.fam.dis <- res.post
  } else{
    sim.fam.dis <- rbind(sim.fam.dis, res.post)
  }
}

### Getting the risk prediction results
for(i in 1:200){
  load(paste(dir.sim, "/Risk Predictions/Population/Dis/simresdis_", i, ".RData", sep = ""))
  if(i == 1){
    sim.res.dis <- res
  } else{
    if(ncol(res) == 12){
      res <- res[, 1:11]
    }
    sim.res.dis <- rbind(sim.res.dis, res)
  }
}

## Correlated frailty distribution
### Family-specific frailty distribution
for(i in 1:200){
  load(paste(dir.sim, "/Family Frailty Distribution/Population/Dis/Cor/simfamdistdiscor_", i, ".RData", sep = ""))
  if(i == 1){
    sim.fam.dis.cor <- res.post
  } else{
    sim.fam.dis.cor <- rbind(sim.fam.dis.cor, res.post)
  }
}

### Getting the risk prediction results
for(i in 1:200){
  load(paste(dir.sim, "/Risk Predictions/Population/Dis/Cor/simresdiscor_", i, ".RData", sep = ""))
  if(i == 1){
    sim.res.dis.cor <- res
  } else{
    if(ncol(res) == 12){
      res <- res[, 1:11]
    }
    sim.res.dis.cor <- rbind(sim.res.dis.cor, res)
  }
}

ggplot(cbind(supp.w, f.w),
       aes(x = W.BC, y = W.OC, fill = f.w)) +
  geom_tile() + geom_text(aes(label = round(f.w, 3)), col = "white") +
  labs(x = "Breast Cancer Frailty",
       y = "Ovarian Cancer Frailty",
       title = "Correlated Discrete Distribution") +
  scale_x_continuous(breaks = w.list.b) +
  scale_y_continuous(breaks = w.list.o) +
  scale_fill_continuous(name = expression(f[W](w)))

## High risk
### Family-specific frailty distribution
for(i in 1:198){
  load(paste(dir.sim, "/Family Frailty Distribution/High Risk/Dis/simfamdistdishr_", i, ".RData", sep = ""))
  if(i == 1){
    sim.fam.dis.hr <- res.post
  } else{
    sim.fam.dis.hr <- rbind(sim.fam.dis.hr, res.post)
  }
}

### Getting the risk prediction results
for(i in 1:198){
  load(paste(dir.sim, "/Risk Predictions/High Risk/Dis/simresdishr_", i, ".RData", sep = ""))
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
ggplot(cbind(supp.w, Mean = base::colMeans(sim.fam.dis[, 1:49], na.rm = TRUE)),
       aes(x = W.BC, y = W.OC, fill = Mean)) +
  geom_tile() + labs(x = "Breast Cancer Frailty",
                     y = "Ovarian Cancer Frailty",
                     title = "Discrete uniform") +
  scale_x_continuous(breaks = w.list.b) +
  scale_y_continuous(breaks = w.list.o)
ggplot(cbind(supp.w, Median = apply(sim.fam.dis[, 1:49], 2, median, na.rm = TRUE)),
       aes(x = W.BC, y = W.OC, fill = Median)) +
  geom_tile() + labs(x = "Breast Cancer Frailty",
                     y = "Ovarian Cancer Frailty",
                     title = "Discrete uniform") +
  scale_x_continuous(breaks = w.list.b) +
  scale_y_continuous(breaks = w.list.o)

hist(sim.fam.dis[, which(supp.w$W.BC == 1.5 & supp.w$W.OC == 0)])



## Heatmaps of means and medians of posterior frailty probabilities
ggplot(cbind(supp.w, Mean = base::colMeans(sim.fam.dis.hr[, 1:49], na.rm = TRUE)),
       aes(x = W.BC, y = W.OC, fill = Mean)) +
  geom_tile() + labs(x = "Breast Cancer Frailty",
                     y = "Ovarian Cancer Frailty",
                     title = "Discrete uniform, High-risk allele frequency") +
  scale_x_continuous(breaks = w.list.b) +
  scale_y_continuous(breaks = w.list.o)
ggplot(cbind(supp.w, Median = apply(sim.fam.dis.hr[, 1:49], 2, median, na.rm = TRUE)),
       aes(x = W.BC, y = W.OC, fill = Median)) +
  geom_tile() + labs(x = "Breast Cancer Frailty",
                     y = "Ovarian Cancer Frailty",
                     title = "Discrete uniform, High-risk allele frequency") +
  scale_x_continuous(breaks = w.list.b) +
  scale_y_continuous(breaks = w.list.o)

load(paste(getwd(), "/Simulations/simdat_dis.RData", sep = ""))
sim.res.dis$BC.5 <- filter(fam.sim, FamID %in% sim.res.dis$FamID, isProband == 1)$AffectedBreast.fu
sim.res.dis.noc <- filter(sim.res.dis, FamID %in% filter(fam.sim, FamID %in% sim.res.dis$FamID, isProband == 1, AffectedOvary == 0)$FamID)

## Correlated frailty distribution
load(paste(dir.sim, "/simdat_dis_cor.RData", sep = ""))
sim.res.dis.cor$BC.5 <- filter(fam.sim, FamID %in% sim.res.dis.cor$FamID, isProband == 1)$AffectedBreast.fu
sim.res.dis.cor.noc <- filter(sim.res.dis.cor, FamID %in% filter(fam.sim, FamID %in% sim.res.dis.cor$FamID, isProband == 1, AffectedOvary == 0)$FamID)


## family size
dt <- data.table(fam.sim)
dt <- dt[, nrow(.SD), by = FamID]
sim.res.dis.noc$FamSize <- merge(sim.res.dis.noc, dt, by = "FamID")$V1

load(paste(getwd(), "/Simulations/simdat_dis_hr.RData", sep = ""))
sim.res.dis.hr$BC.5 <- filter(fam.sim, FamID %in% sim.res.dis.hr$FamID, isProband == 1)$AffectedBreast.fu
sim.res.dis.hr.noc <- filter(sim.res.dis.hr, FamID %in% filter(fam.sim, FamID %in% sim.res.dis.hr$FamID, isProband == 1, AffectedOvary == 0)$FamID)

## family size
dt <- data.table(fam.sim)
dt <- dt[, nrow(.SD), by = FamID]
sim.res.dis.hr.noc$FamSize <- merge(sim.res.dis.hr.noc, dt, by = "FamID")$V1

## plotting true frailties vs frailty means for the discrete uniform
mean.post.b <- function(post, w.list.b, supp.w){
  res <- 0
  for(j in 1:length(w.list.b)){
    res <- res + sum(post[which(supp.w$W.BC == w.list.b[j])]) * w.list.b[j]
  }
  return(res)
}
mean.post.o <- function(post, w.list.o, supp.w){
  res <- 0
  for(j in 1:length(w.list.o)){
    res <- res + sum(post[which(supp.w$W.OC == w.list.o[j])]) * w.list.o[j]
  }
  return(res)
}

sim.res.dis.noc$wbc.mean <- apply(filter(sim.fam.dis, FamID %in% sim.res.dis.noc$FamID)[, 1:(ncol(sim.fam.dis) - 1)], 1,
                                  mean.post.b, w.list.b = w.list.b, supp.w = supp.w)
sim.res.dis.noc$woc.mean <- apply(filter(sim.fam.dis, FamID %in% sim.res.dis.noc$FamID)[, 1:(ncol(sim.fam.dis) - 1)], 1,
                                  mean.post.o, w.list.o = w.list.o, supp.w = supp.w)
sim.res.dis.noc$wbc.mean.dis <- sapply(sim.res.dis.noc$wbc.mean, function(x) w.list.b[which.min(abs(w.list.b - x))])
sim.res.dis.noc$woc.mean.dis <- sapply(sim.res.dis.noc$woc.mean, function(x) w.list.o[which.min(abs(w.list.o - x))])


sim.res.dis.noc$wbc.mode <- apply(filter(sim.fam.dis, FamID %in% sim.res.dis.noc$FamID)[, 1:(ncol(sim.fam.dis) - 1)], 1,
                                  function(x) supp.w$W.BC[which.max(x)])
sim.res.dis.noc$woc.mode <- apply(filter(sim.fam.dis, FamID %in% sim.res.dis.noc$FamID)[, 1:(ncol(sim.fam.dis) - 1)], 1,
                                  function(x) supp.w$W.OC[which.max(x)])

load(paste(getwd(), "/Simulations/simdat_bc_dis.RData", sep = ""))
sim.res.dis.noc <- merge(sim.res.dis.noc,
                         select(filter(fam.sim.bc, FamID %in% sim.res.dis.noc$FamID, isProband == 1),
                                W.BC, W.OC, FamID),
                         by = "FamID")

dt.b <- data.table(sim.res.dis.noc)
dt.b <- dt.b[, mean(.SD$wbc.mean), by = list(W.BC, W.OC)]
dt.b <- join(supp.w, dt.b)

dt.o <- data.table(sim.res.dis.noc)
dt.o <- dt.o[, mean(.SD$woc.mean), by = list(W.BC, W.OC)]
dt.o <- join(supp.w, dt.o)

ggplot(cbind(dt.b, V2 = dt.o[, 3], V3 = rowMeans(cbind(dt.b[, 3], dt.o[, 3]))),
       aes(x = W.BC, y = W.OC, fill = V3)) +
  geom_tile() + labs(x = "True Breast Cancer Frailty",
                     y = "True Ovarian Cancer Frailty") +
  geom_text(aes(label = paste(round(V1, 2), ", ", round(V2, 2))), col = "white") +
  scale_x_continuous(breaks = w.list.b) +
  scale_y_continuous(breaks = w.list.o) +
  scale_fill_continuous(name = "Mean of Mean BC and \n Mean OC Frailties")
ggplot(dt.o, aes(x = W.BC, y = W.OC, fill = V1)) +
  geom_tile() + labs(x = "True Breast Cancer Frailty",
                     y = "True Ovarian Cancer Frailty") +
  geom_text(aes(label = round(V1, 3)), col = "white") +
  scale_x_continuous(breaks = w.list.b) +
  scale_y_continuous(breaks = w.list.o) +
  scale_fill_continuous(name = "Mean of Mean OC Frailty")

ggplot(sim.res.dis.noc, aes(factor(W.BC), wbc.mean)) +
  geom_boxplot()

ggplot(sim.res.dis.noc, aes(factor(W.BC), wbc.mode)) +
  geom_boxplot()



getOE <- function(dat, out, pred){
  sum(select(dat, out)[, 1]) / sum(select(dat, pred)[, 1])
}

getAUC <- function(dat, out, pred){
  pROC::auc(select(dat, out)[, 1], select(dat, pred)[, 1], data = dat)
}

getRBS <- function(dat, out, pred){
  sqrt(mean((select(dat, out)[, 1] - select(dat, pred)[, 1])^2))
}


# perf.boot.sim <- function(dat, out, pred, nboot){
#   res.boot <- setNames(data.frame(matrix(NA, nboot, 3)), c("OE", "AUC", "rBS"))
#   for(i in 1:nboot){
#     dat.boot <- dat[sample(1:nrow(dat), nrow(dat), replace = TRUE), ]
#     res.boot$OE[i] <- getOE(dat.boot, out, pred)
#     res.boot$AUC[i] <- getAUC(dat.boot, out, pred)
#     res.boot$rBS[i] <- getRBS(dat.boot, out, pred)
#   }
#   res <- c(getOE(dat, out, pred), quantile(res.boot$OE, 0.025), quantile(res.boot$OE, 0.975),
#            getAUC(dat, out, pred), quantile(res.boot$AUC, 0.025), quantile(res.boot$AUC, 0.975),
#            getRBS(dat, out, pred), quantile(res.boot$rBS, 0.025), quantile(res.boot$rBS, 0.975))
#   return(res)
# }

perf.boot.sim.comp <- function(dat, out, pred.f, pred.nf, nboot, seed = 12345){
  set.seed(seed)
  res.boot <- setNames(data.frame(matrix(NA, nboot, 3)), c("OE", "AUC", "rBS"))
  oe <- auc <- rbs <- setNames(data.frame(matrix(NA, nboot, 2)), c("Frailty", "NoFrailty"))
  for(i in 1:nboot){
    dat.boot <- dat[sample(1:nrow(dat), nrow(dat), replace = TRUE), ]
    oe$Frailty[i] <- getOE(dat.boot, out, pred.f)
    oe$NoFrailty[i] <- getOE(dat.boot, out, pred.nf)
    auc$Frailty[i] <- getAUC(dat.boot, out, pred.f)
    auc$NoFrailty[i] <- getAUC(dat.boot, out, pred.nf)
    rbs$Frailty[i] <- getRBS(dat.boot, out, pred.f)
    rbs$NoFrailty[i] <- getRBS(dat.boot, out, pred.nf)
  }
  
  res <- rbind(c(getOE(dat, out, pred.f), quantile(oe$Frailty, 0.025), quantile(oe$Frailty, 0.975),
           getAUC(dat, out, pred.f), quantile(auc$Frailty, 0.025), quantile(auc$Frailty, 0.975),
           getRBS(dat, out, pred.f), quantile(rbs$Frailty, 0.025), quantile(rbs$Frailty, 0.975)),
           c(getOE(dat, out, pred.nf), quantile(oe$NoFrailty, 0.025), quantile(oe$NoFrailty, 0.975),
             getAUC(dat, out, pred.nf), quantile(auc$NoFrailty, 0.025), quantile(auc$NoFrailty, 0.975),
             getRBS(dat, out, pred.nf), quantile(rbs$NoFrailty, 0.025), quantile(rbs$NoFrailty, 0.975)))
  
  return(list(res = res, oe = oe, auc = auc, rbs = rbs))
}

perf.fs <- function(dat, out, pred.f, pred.nf, famsize){
  res <- setNames(data.frame(matrix(NA, 2, 6)),
                  c("OE.F", "AUC.F", "rBS.F", "OE.NF", "AUC.NF", "rBS.NF"))
  res$OE.F[1] <- getOE(filter(dat, FamSize <= famsize), out, pred.f)
  res$AUC.F[1] <- getAUC(filter(dat, FamSize <= famsize), out, pred.f)
  res$rBS.F[1] <- getRBS(filter(dat, FamSize <= famsize), out, pred.f)
  res$OE.NF[1] <- getOE(filter(dat, FamSize <= famsize), out, pred.nf)
  res$AUC.NF[1] <- getAUC(filter(dat, FamSize <= famsize), out, pred.nf)
  res$rBS.NF[1] <- getRBS(filter(dat, FamSize <= famsize), out, pred.nf)
  res$OE.F[2] <- getOE(filter(dat, FamSize > famsize), out, pred.f)
  res$AUC.F[2] <- getAUC(filter(dat, FamSize > famsize), out, pred.f)
  res$rBS.F[2] <- getRBS(filter(dat, FamSize > famsize), out, pred.f)
  res$OE.NF[2] <- getOE(filter(dat, FamSize > famsize), out, pred.nf)
  res$AUC.NF[2] <- getAUC(filter(dat, FamSize > famsize), out, pred.nf)
  res$rBS.NF[2] <- getRBS(filter(dat, FamSize > famsize), out, pred.nf)
  rownames(res) <- c(paste("FamSize <= ", famsize, sep = ""), paste("FamSize > ", famsize, sep = ""))
  return(res)
}

### performance results by family size

cutoff <- 30

xtable(rbind(perf.fs(sim.res.v03.noc, "BC.5", "Prob.BC.5", "Prob.BC.5.nf", cutoff),
             perf.fs(sim.res.v2.noc, "BC.5", "Prob.BC.5", "Prob.BC.5.nf", cutoff),
             perf.fs(sim.res.dis.noc, "BC.5", "Prob.BC.5", "Prob.BC.5.nf", cutoff),
             perf.fs(sim.res.v03.hr.noc, "BC.5", "Prob.BC.5", "Prob.BC.5.nf", cutoff),
             perf.fs(sim.res.v2.hr.noc, "BC.5", "Prob.BC.5", "Prob.BC.5.nf", cutoff),
             perf.fs(sim.res.dis.hr.noc, "BC.5", "Prob.BC.5", "Prob.BC.5.nf", cutoff))[, 1:3],
       digits = 3)



#### Performance Metrics ######

## BRCAPRO allele frequencies

nboot <- 1000
sim.v03.boot <- perf.boot.sim.comp(sim.res.v03.noc, "BC.5", "Prob.BC.5", "Prob.BC.5.nf", nboot)
sim.v2.boot <- perf.boot.sim.comp(sim.res.v2.noc, "BC.5", "Prob.BC.5", "Prob.BC.5.nf", nboot)
sim.dis.boot <- perf.boot.sim.comp(sim.res.dis.noc, "BC.5", "Prob.BC.5", "Prob.BC.5.nf", nboot)

mean(abs(sim.v03.boot$oe$Frailty - 1) <= abs(sim.v03.boot$oe$NoFrailty - 1))
mean(sim.v03.boot$auc$Frailty >= sim.v03.boot$auc$NoFrailty)
mean(sim.v03.boot$rbs$Frailty <= sim.v03.boot$rbs$NoFrailty)

mean(abs(sim.v2.boot$oe$Frailty - 1) <= abs(sim.v2.boot$oe$NoFrailty - 1))
mean(sim.v2.boot$auc$Frailty >= sim.v2.boot$auc$NoFrailty)
mean(sim.v2.boot$rbs$Frailty <= sim.v2.boot$rbs$NoFrailty)

mean(abs(sim.dis.boot$oe$Frailty - 1) <= abs(sim.dis.boot$oe$NoFrailty - 1))
mean(sim.dis.boot$auc$Frailty >= sim.dis.boot$auc$NoFrailty)
mean(sim.dis.boot$rbs$Frailty <= sim.dis.boot$rbs$NoFrailty)

res.all[1:2, ] <- sim.v03.boot$res
res.all[3:4, ] <- sim.v2.boot$res
res.all[5:6, ] <- sim.dis.boot$res

# Correlated frailty distribution
sim.cor.boot <- perf.boot.sim.comp(sim.res.dis.cor.noc, "BC.5", "Prob.BC.5", "Prob.BC.5.nf", nboot)

sim.cor.boot$res

mean(abs(sim.cor.boot$oe$Frailty - 1) <= abs(sim.cor.boot$oe$NoFrailty - 1))
mean(sim.cor.boot$auc$Frailty >= sim.cor.boot$auc$NoFrailty)
mean(sim.cor.boot$rbs$Frailty <= sim.cor.boot$rbs$NoFrailty)

save(sim.res.v03.noc, sim.res.v2.noc, sim.res.dis.noc, sim.res.dis.cor.noc,
     sim.fam.v03, sim.fam.v2, sim.fam.dis, sim.fam.dis.cor,
     res.all, sim.v03.boot, sim.v2.boot, sim.dis.boot, sim.cor.boot,
     file = paste0(dir.sim, "/Frailty_Sim_Analysis_Results.RData"))


## High-risk allele frequencies

sim.v03.hr.boot <- perf.boot.sim.comp(sim.res.v03.hr.noc, "BC.5", "Prob.BC.5", "Prob.BC.5.nf", nboot)
sim.v2.hr.boot <- perf.boot.sim.comp(sim.res.v2.hr.noc, "BC.5", "Prob.BC.5", "Prob.BC.5.nf", nboot)
sim.dis.hr.boot <- perf.boot.sim.comp(sim.res.dis.hr.noc, "BC.5", "Prob.BC.5", "Prob.BC.5.nf", nboot)

mean(abs(sim.v03.hr.boot$oe$Frailty - 1) <= abs(sim.v03.hr.boot$oe$NoFrailty - 1))
mean(sim.v03.hr.boot$auc$Frailty >= sim.v03.hr.boot$auc$NoFrailty)
mean(sim.v03.hr.boot$rbs$Frailty <= sim.v03.hr.boot$rbs$NoFrailty)

mean(abs(sim.v2.hr.boot$oe$Frailty - 1) <= abs(sim.v2.hr.boot$oe$NoFrailty - 1))
mean(sim.v2.hr.boot$auc$Frailty >= sim.v2.hr.boot$auc$NoFrailty)
mean(sim.v2.hr.boot$rbs$Frailty <= sim.v2.hr.boot$rbs$NoFrailty)

mean(abs(sim.dis.hr.boot$oe$Frailty - 1) <= abs(sim.dis.hr.boot$oe$NoFrailty - 1))
mean(sim.dis.hr.boot$auc$Frailty >= sim.dis.hr.boot$auc$NoFrailty)
mean(sim.dis.hr.boot$rbs$Frailty <= sim.dis.hr.boot$rbs$NoFrailty)

res.all.hr[1:2, ] <- sim.v03.hr.boot$res
res.all.hr[3:4, ] <- sim.v2.hr.boot$res
res.all.hr[5:6, ] <- sim.dis.hr.boot$res

save(sim.res.v03.hr.noc, sim.res.v2.hr.noc, sim.res.dis.hr.noc,
     sim.fam.v03.hr, sim.fam.v2.hr, sim.fam.dis.hr,
     res.all.hr, sim.v03.hr.boot, sim.v2.hr.boot, sim.dis.hr.boot,
     file = paste0(dir.sim, "/Frailty_Sim_Analysis_Results_hr.RData"))


