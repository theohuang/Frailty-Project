## Frailty simulation
## Generating the families
## Variance parameters of 0.3
## Higher risk families
## Last udpated: March 11, 2019

source("Frailty Simulation Functions.R")


library(BayesMendel)
library(BMmultigene)
library(mvtnorm)
library(dplyr)
library(doParallel)
registerDoParallel(cores = 4)



cancers <- c("BC", "OC")
genes <- c("BRCA1", "BRCA2")
mu <- c(0, 0)
var <- c(0.3, 0.3)
rho <- 0
## BRCAPRO penetrances (using the BMmultigene assumptions)
penCancersF <- list(BC = penet.brca.net$fFX[, c("B00", "B10", "B01")],
                    OC = penet.brca.net$fFY[, c("B00", "B10", "B01")])
penCancersM <- list(BC = penet.brca.net$fMX[, c("B00", "B10", "B01")],
                    OC = penet.brca.net$fMY[, c("B00", "B10", "B01")])

## baseline hazards
hzd0.b <- cbhf.norm(mu[1], var[1], penCancersF$BC)
hzd0.o <- cbhf.norm(mu[2], var[2], penCancersF$OC)
hzd0.b.m <- cbhf.norm(mu[1], var[1], penCancersM$BC)
hzd0.o.m <- cbhf.norm(mu[2], var[2], penCancersM$OC) 


af.aj <- setNames(c(0.05, 0.05), genes) ## AJ allele frequency
af.naj <- setNames(c(0.02, 0.02), genes) ## non-AJ allele frequency

n.parallel <- 100
n.fam <- 1e3


## Generating the families
start <- Sys.time()
fam.sim <- foreach(i = 1:n.parallel, .combine = append) %dopar% {
  gen.fam.norm(n.fam, genes, cancers, penCancersF, penCancersM, af.aj, af.naj,
               mu = mu, var = var, rho = rho, hzd0.b = hzd0.b, hzd0.o = hzd0.o,
               seed = 999 + i)
}
print(difftime(Sys.time(), start, units = "secs"))


## adding FamIDs
for(i in 1:(n.parallel * n.fam)){
  fam.sim[[i]]$FamID <- i
}
## combining into a single data frame
fam.sim <- bind_rows(fam.sim)


## subset for BC risk prediction
fam.sim.bc <- filter(fam.sim, FamID %in% filter(fam.sim, isProband == 1, Gender == 0, AffectedBreast == 0)$FamID)


save(fam.sim, file = paste(getwd(), "/Frailty/Simulations/simdat_v03_hr.RData", sep = ""))
save(fam.sim.bc, file = paste(getwd(), "/Frailty/Simulations/simdat_bc_v03_hr.RData", sep = ""))




