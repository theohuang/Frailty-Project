## Frailty Simulations
## Simulating with different frailty distributions
## Getting population-level frailty distribution
## Variance parameters of 2
## Higher risk families
## Last updated: March 14, 2019

rm(list = ls())
a1 <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

net <- TRUE

library(BayesMendel)
library(BMmultigene)
library(dplyr)
library(data.table)
library(mvtnorm)
library(doParallel)
registerDoParallel(cores = 4)


## loading simulated data
load(paste(getwd(), "/Frailty/Simulations/simdat_v2_hr.RData", sep = ""))
famid.list <- unique(fam.sim$FamID)

n.sim <- 500
a1.max <- ceiling(length(famid.list) / n.sim); n.max <- length(famid.list)
b1 <- (a1 - 1) * n.sim + 1
b2 <- ifelse(a1 == a1.max, n.max, a1 * n.sim)
int <- b1:b2

dat <- filter(fam.sim, FamID %in% famid.list[int])

## frailty distribution
w.list.b <- seq(-1.5, 1.5, 0.5)
w.list.o <- seq(-1.5, 1.5, 0.5)
supp.w <- expand.grid(w.list.b, w.list.o)
names(supp.w) <- c("W.BC", "W.OC")
f.w <- rep(1, nrow(supp.w)) / nrow(supp.w) ## uniform
f.w.b <- rep(0, length(w.list.b))
for(i in 1:length(w.list.b)){
  f.w.b[i] <- sum(f.w[supp.w$W.BC == w.list.b[i]])
}
f.w.o <- rep(0, length(w.list.o))
for(i in 1:length(w.list.o)){
  f.w.o[i] <- sum(f.w[supp.w$W.OC == w.list.o[i]])
}

source(paste(getwd(), "/Estimating Functions Discrete.R", sep = ""))
source(paste(getwd(), "/Parameters Discrete.R", sep = ""))



start <- Sys.time()
res.post <- foreach(i = 1:nrow(supp.w), .combine = cbind) %dopar% {
  post.fam(dat, i, int, net = net)
}
difftime(Sys.time(), start)

## normalizing
for(i in 1:nrow(res.post)){
  res.post[i, ] <- res.post[i, ] / sum(res.post[i, ])
}

res.post <- data.frame(res.post)
res.post$FamID <- unique(dat$FamID)

save(res.post, file = paste(getwd(), "/Frailty/Simulations/Population Frailty Distribution/High Risk/Var2/simpopv2hr_", a1, ".RData", sep = ""))






