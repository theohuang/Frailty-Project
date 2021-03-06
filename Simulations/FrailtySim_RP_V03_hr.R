## Running the frailty model on the simulated data
## Using variance parameters of 0.3
## Net Risk
## Higher risk
## Last updated: May 1, 2019

rm(list = ls())
a1 <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

net <- TRUE

library(BayesMendel)
library(dplyr)
library(data.table)
library(abind)
library(doParallel)
registerDoParallel(cores = 4)

## loading simulated data (BC subset)
load(paste(getwd(), "/Frailty/Simulations/simdat_bc_v03_hr.RData", sep = ""))

## only using subset of families
famid.list <- unique(fam.sim.bc$FamID)

## high-risk allele frequencies
af.aj <- c(0.05, 0.05) ## AJ allele frequency
af.naj <- c(0.02, 0.02) ## non-AJ allele frequency


n.sim <- 500
a1.max <- ceiling(length(famid.list) / n.sim); n.max <- length(famid.list)
b1 <- (a1 - 1) * n.sim + 1
b2 <- ifelse(a1 == a1.max, n.max, a1 * n.sim)
int <- b1:b2

dat <- filter(fam.sim.bc, FamID %in% famid.list[int])


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

save(res.post, file = paste(getwd(), "/Frailty/Simulations/Family Frailty Distribution/High Risk/Var03/simfamdistv03hr_", a1, ".RData", sep = ""))

## now getting risk predictions
start <- Sys.time()
res <- foreach(i = 1:length(int), .combine = rbind) %dopar% {
  rp.fam(dat, i, int, res.post, net = net, gt = FALSE,
         af.aj = af.aj, af.naj = af.naj)
}
difftime(Sys.time(), start, units = "secs")

res <- data.frame(res)
names(res) <- c("FamID", "Prob.BC.5", "Prob.OC.5", "Prob.BRCA",
                "Prob.BRCA1", "Prob.BRCA2", "Prob.BC.5.nf", "Prob.OC.5.nf",
                "Prob.BRCA.nf", "Prob.BRCA1.nf", "Prob.BRCA2.nf")

save(res, file = paste(getwd(), "/Frailty/Simulations/Risk Predictions/High Risk/Var03/simresv03hr_", a1, ".RData", sep = ""))




