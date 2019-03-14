## Obtaining Posterior Distribution for Frailty Simulation
## Last Updated: March 13, 2019

library(dplyr)
library(ggplot2)
library(data.table)
dir.res <- "/Users/Theo/Dropbox (Partners HealthCare)/Frailty Github/Simulations"


### Results for all data (not just breast cancer subset)
for(i in 1:200){
  load(paste(dir.res, "/Population Frailty Distribution/Var03/simpopv03_", i, ".RData", sep = ""))
  if(i == 1){
    pop.v03 <- res.post
  } else{
    pop.v03 <- rbind(pop.v03, res.post)
  }
}

### Getting the new prior frailty distribution
f.w.pop.v03 <- apply(pop.v03[, 1:49], 2, median, na.rm = TRUE)
f.w.pop.v03 <- f.w.pop.v03 / sum(f.w.pop.v03)
save(f.w.pop.v03, file = paste(dir.res, "/Population Frailty Distribution/fwpop.v03.RData", sep = ""))

# load(paste(getwd(), "/Simulations/simdat_v03.RData", sep = ""))
# 
# pop.v03[, c("W.BC", "W.OC")] <- select(filter(fam.sim, isProband == 1), W.BC, W.OC)
# pop.v03[, c("W.BC.D", "W.OC.D")] <- cbind(apply(select(filter(fam.sim, isProband == 1), W.BC), 1, function(x) w.list.b[which.min(abs(x - w.list.b))]),
#                                                apply(select(filter(fam.sim, isProband == 1), W.OC), 1, function(x) w.list.o[which.min(abs(x - w.list.o))]))


### Heatmap of new prior distribution
w.list.b <- seq(-1.5, 1.5, 0.5)
w.list.o <- seq(-1.5, 1.5, 0.5)
supp.w <- expand.grid(w.list.b, w.list.o)
names(supp.w) <- c("W.BC", "W.OC")
ggplot(cbind(supp.w, Median = f.w.pop.v03),
       aes(x = W.BC, y = W.OC, fill = Median)) +
  geom_tile() + labs(x = "Breast Cancer Frailty",
                     y = "Ovarian Cancer Frailty") +
  scale_x_continuous(breaks = w.list.b) +
  scale_y_continuous(breaks = w.list.o)

# ## heatmap of mixture distribution with equal weights (mean instead of median)
# ggplot(cbind(supp.w, Mean = colMeans(pop.v03[, 1:49])),
#        aes(x = W.BC, y = W.OC, fill = Mean)) +
#   geom_tile() + labs(x = "Breast Cancer Frailty",
#                      y = "Ovarian Cancer Frailty") +
#   scale_x_continuous(breaks = w.list.b) +
#   scale_y_continuous(breaks = w.list.o)
# 
# 
# ## Heatmap of true (discretized) frailty distribution
# dist.w <- rep(0, nrow(supp.w))
# for(i in 1:nrow(supp.w)){
#   dist.w[i] <- length(which(pop.v03$W.BC.D == supp.w$W.BC[i] &
#                               pop.v03$W.OC.D == supp.w$W.OC[i]))
# }
# ggplot(cbind(supp.w, Sum = dist.w),
#        aes(x = W.BC, y = W.OC, fill = Sum)) +
#   geom_tile() + labs(x = "Breast Cancer Frailty",
#                      y = "Ovarian Cancer Frailty") +
#   scale_x_continuous(breaks = w.list.b) +
#   scale_y_continuous(breaks = w.list.o)
# 
# ## Scatterplot of true frailties
# ggplot(filter(fam.sim, isProband == 1), aes(W.BC, W.OC)) + geom_point(alpha = 0.3)
# 

# get.mean <- function(post, w.list, can){
#   dist.w <- rep(0, length(w.list))
#   for(i in 1:length(w.list)){
#     dist.w[i] <- sum(post[which(supp.w[, can] == w.list[i])])
#   }
#   return(sum(w.list * dist.w))
# }

# res <- pop.v03
# res$W.BC.mean <- res$W.OC.mean <- 0
# res$W.BC.mean <- apply(pop.v03, 1, get.mean, w.list = w.list.b, can = 1)
# res$W.OC.mean <- apply(pop.v03, 1, get.mean, w.list = w.list.o, can = 2)
# 
# pop.v03$W.BC.mean <- res$W.BC.mean
# pop.v03$W.OC.mean <- res$W.OC.mean
# pop.v03[, c("W.BC.mean.D", "W.OC.mean.D")] <- cbind(apply(matrix(pop.v03$W.BC.mean, ncol = 1), 1, function(x) w.list.b[which.min(abs(x - w.list.b))]),
#                                                          apply(matrix(pop.v03$W.OC.mean, ncol = 1), 1, function(x) w.list.o[which.min(abs(x - w.list.o))]))
# 
# 
# dist.w.mean <- rep(0, nrow(supp.w))
# for(i in 1:nrow(supp.w)){
#   dist.w.mean[i] <- length(which(pop.v03$W.BC.mean.D == supp.w$W.BC[i] &
#                                   pop.v03$W.OC.mean.D == supp.w$W.OC[i]))
# }
# ggplot(cbind(supp.w, Sum = dist.w.mean),
#        aes(x = W.BC, y = W.OC, fill = Sum)) +
#   geom_tile() + labs(x = "Breast Cancer Frailty",
#                      y = "Ovarian Cancer Frailty") +
#   scale_x_continuous(breaks = w.list.b) +
#   scale_y_continuous(breaks = w.list.o)
# 




### Variances of 2

### Results for all of cgn data (not just breast cancer subset)
for(i in 1:200){
  load(paste(dir.res, "/Population Frailty Distribution/Var2/simpopv2_", i, ".RData", sep = ""))
  if(i == 1){
    pop.v2 <- res.post
  } else{
    pop.v2 <- rbind(pop.v2, res.post)
  }
}

### Getting the new prior frailty distribution
f.w.pop.v2 <- apply(pop.v2[, 1:49], 2, median, na.rm = TRUE)
f.w.pop.v2 <- f.w.pop.v2 / sum(f.w.pop.v2)
save(f.w.pop.v2, file = paste(dir.res, "/Population Frailty Distribution/fwpop.v2.RData", sep = ""))

### Heatmap of new prior distribution
ggplot(cbind(supp.w, Median = f.w.pop.v2),
       aes(x = W.BC, y = W.OC, fill = Median)) +
  geom_tile() + labs(x = "Breast Cancer Frailty",
                     y = "Ovarian Cancer Frailty") +
  scale_x_continuous(breaks = w.list.b) +
  scale_y_continuous(breaks = w.list.o)

ggplot(cbind(supp.w, Median = apply(pop.v03[, 1:49], 2, mean, na.rm = TRUE)),
       aes(x = W.BC, y = W.OC, fill = Median)) +
  geom_tile() + labs(x = "Breast Cancer Frailty",
                     y = "Ovarian Cancer Frailty") +
  scale_x_continuous(breaks = w.list.b) +
  scale_y_continuous(breaks = w.list.o)


# res <- pop.v2
# res$W.BC.mean <- res$W.OC.mean <- 0
# res$W.BC.mean <- apply(pop.v2, 1, get.mean, w.list = w.list.b, can = 1)
# res$W.OC.mean <- apply(pop.v2, 1, get.mean, w.list = w.list.o, can = 2)
# 
# pop.v2$W.BC.mean <- res$W.BC.mean
# pop.v2$W.OC.mean <- res$W.OC.mean
# pop.v2[, c("W.BC", "W.OC")] <- select(filter(fam.sim, isProband == 1, FamID %in% pop.v2$FamID), W.BC, W.OC)



### Discrete distribution

### Results for all of cgn data (not just breast cancer subset)
for(i in 1:200){
  load(paste(dir.res, "/Population Frailty Distribution/Dis/simpopdis_", i, ".RData", sep = ""))
  if(i == 1){
    pop.dis <- res.post
  } else{
    pop.dis <- rbind(pop.dis, res.post)
  }
}

### Getting the new prior frailty distribution
f.w.pop.dis <- apply(pop.dis[, 1:49], 2, median, na.rm = TRUE)
f.w.pop.dis <- f.w.pop.dis / sum(f.w.pop.dis)
save(f.w.pop.dis, file = paste(dir.res, "/Population Frailty Distribution/fwpop.dis.RData", sep = ""))


### Heatmap of new prior distribution
ggplot(cbind(supp.w, Median = f.w.pop.dis),
       aes(x = W.BC, y = W.OC, fill = Median)) +
  geom_tile() + labs(x = "Breast Cancer Frailty",
                     y = "Ovarian Cancer Frailty") +
  scale_x_continuous(breaks = w.list.b) +
  scale_y_continuous(breaks = w.list.o)

# ggplot(cbind(supp.w, Median = apply(filter(pop.dis, W.BC == 1.5)[, 1:49], 2, median)),
#        aes(x = W.BC, y = W.OC, fill = Median)) +
#   geom_tile() + labs(x = "Breast Cancer Frailty",
#                      y = "Ovarian Cancer Frailty") +
#   scale_x_continuous(breaks = w.list.b) +
#   scale_y_continuous(breaks = w.list.o)



# get.mean <- function(post, w.list, can){
#   dist.w <- rep(0, length(w.list))
#   for(i in 1:length(w.list)){
#     dist.w[i] <- sum(post[which(supp.w[, can] == w.list[i])])
#   }
#   return(sum(w.list * dist.w))
# }
# 
# res <- pop.dis
# res$W.BC.mean <- res$W.OC.mean <- 0
# res$W.BC.mean <- apply(pop.dis, 1, get.mean, w.list = w.list.b, can = 1)
# res$W.OC.mean <- apply(pop.dis, 1, get.mean, w.list = w.list.o, can = 2)
# 
# pop.dis$W.BC.mean <- res$W.BC.mean
# pop.dis$W.OC.mean <- res$W.OC.mean
# pop.dis[, c("W.BC", "W.OC")] <- select(filter(fam.sim, isProband == 1, FamID %in% pop.dis$FamID), W.BC, W.OC)


save(pop.v03, pop.v2, pop.dis,
     file = paste(dir.res, "/Frailty Simulation Population Frailty Distribution.RData", sep = ""))



#### High risk families #####

### Results for all data (not just breast cancer subset)
for(i in 1:200){
  load(paste(dir.res, "/Population Frailty Distribution/High Risk/Var03/simpopv03hr_", i, ".RData", sep = ""))
  if(i == 1){
    pop.v03.hr <- res.post
  } else{
    pop.v03.hr <- rbind(pop.v03.hr, res.post)
  }
}
for(i in 1:200){
  load(paste(dir.res, "/Population Frailty Distribution/High Risk/Var2/simpopv2hr_", i, ".RData", sep = ""))
  if(i == 1){
    pop.v2.hr <- res.post
  } else{
    pop.v2.hr <- rbind(pop.v2.hr, res.post)
  }
}
for(i in 1:200){
  load(paste(dir.res, "/Population Frailty Distribution/High Risk/Dis/simpopdishr_", i, ".RData", sep = ""))
  if(i == 1){
    pop.dis.hr <- res.post
  } else{
    pop.dis.hr <- rbind(pop.dis.hr, res.post)
  }
}

### Getting the new prior frailty distribution
f.w.pop.v03.hr <- apply(pop.v03.hr[, 1:49], 2, median, na.rm = TRUE)
f.w.pop.v03.hr <- f.w.pop.v03.hr / sum(f.w.pop.v03.hr)
f.w.pop.v2.hr <- apply(pop.v2.hr[, 1:49], 2, median, na.rm = TRUE)
f.w.pop.v2.hr <- f.w.pop.v2.hr / sum(f.w.pop.v2.hr)
f.w.pop.dis.hr <- apply(pop.dis.hr[, 1:49], 2, median, na.rm = TRUE)
f.w.pop.dis.hr <- f.w.pop.dis.hr / sum(f.w.pop.dis.hr)
save(f.w.pop.v03.hr, file = paste(dir.res, "/Population Frailty Distribution/High Risk/fwpop.v03.hr.RData", sep = ""))
save(f.w.pop.v2.hr, file = paste(dir.res, "/Population Frailty Distribution/High Risk/fwpop.v2.hr.RData", sep = ""))
save(f.w.pop.dis.hr, file = paste(dir.res, "/Population Frailty Distribution/High Risk/fwpop.dis.hr.RData", sep = ""))


# load(paste(getwd(), "/Simulations/simdat_v03.RData", sep = ""))
# 
# pop.v03[, c("W.BC", "W.OC")] <- select(filter(fam.sim, isProband == 1), W.BC, W.OC)
# pop.v03[, c("W.BC.D", "W.OC.D")] <- cbind(apply(select(filter(fam.sim, isProband == 1), W.BC), 1, function(x) w.list.b[which.min(abs(x - w.list.b))]),
#                                                apply(select(filter(fam.sim, isProband == 1), W.OC), 1, function(x) w.list.o[which.min(abs(x - w.list.o))]))


### Heatmap of new prior distribution
w.list.b <- seq(-1.5, 1.5, 0.5)
w.list.o <- seq(-1.5, 1.5, 0.5)
supp.w <- expand.grid(w.list.b, w.list.o)
names(supp.w) <- c("W.BC", "W.OC")
ggplot(cbind(supp.w, Median = f.w.pop.v03.hr),
       aes(x = W.BC, y = W.OC, fill = Median)) +
  geom_tile() + labs(x = "Breast Cancer Frailty",
                     y = "Ovarian Cancer Frailty") +
  scale_x_continuous(breaks = w.list.b) +
  scale_y_continuous(breaks = w.list.o)
ggplot(cbind(supp.w, Median = f.w.pop.v2.hr),
       aes(x = W.BC, y = W.OC, fill = Median)) +
  geom_tile() + labs(x = "Breast Cancer Frailty",
                     y = "Ovarian Cancer Frailty") +
  scale_x_continuous(breaks = w.list.b) +
  scale_y_continuous(breaks = w.list.o)
ggplot(cbind(supp.w, Median = f.w.pop.dis.hr),
       aes(x = W.BC, y = W.OC, fill = Median)) +
  geom_tile() + labs(x = "Breast Cancer Frailty",
                     y = "Ovarian Cancer Frailty") +
  scale_x_continuous(breaks = w.list.b) +
  scale_y_continuous(breaks = w.list.o)


dt <- data.table()

