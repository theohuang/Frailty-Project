## Frailty simulation
## Generating the families
## Generating from discrete distribution
## Last udpated: April 11, 2019

source("Frailty Simulation Functions.R")


library(BayesMendel)
library(dplyr)
library(abind)
library(doParallel)
registerDoParallel(cores = 4)

source(paste(getwd(), "/Simulations/Generating Families Functions/sim.simFam.R", sep = ""))
source(paste(getwd(), "/Simulations/Generating Families Functions/genCancerPen.R", sep = ""))
source(paste(getwd(), "/Simulations/Generating Families Functions/sim.buildGenoMat.R", sep = ""))
source(paste(getwd(), "/Simulations/Generating Families Functions/sim.linkParents.R", sep = ""))
source(paste(getwd(), "/Simulations/Generating Families Functions/sim.simCurAgeVar.R", sep = ""))
source(paste(getwd(), "/Simulations/Generating Families Functions/sim.simCancerVars.R", sep = ""))
source(paste(getwd(), "/Simulations/Generating Families Functions/sim.buildBranchOfAlleleMats.R", sep = ""))
source(paste(getwd(), "/Simulations/Generating Families Functions/helpers.R", sep = ""))


cancers <- c("BC", "OC")
genes <- c("BRCA1", "BRCA2")
## BRCAPRO penetrances (using the BMmultigene assumptions)
penCancersF <- list(BC = penet.brca.net$fFX[, c("B00", "B10", "B01")],
                    OC = penet.brca.net$fFY[, c("B00", "B10", "B01")])
penCancersM <- list(BC = penet.brca.net$fMX[, c("B00", "B10", "B01")],
                    OC = penet.brca.net$fMY[, c("B00", "B10", "B01")])

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

## baseline hazards
hzd0.b <- cbhf.dis(w.list.b, f.w.b, penCancersF$BC)
hzd0.o <- cbhf.dis(w.list.o, f.w.o, penCancersF$OC)
hzd0.b.m <- cbhf.dis(w.list.b, f.w.b, penCancersM$BC)
hzd0.o.m <- cbhf.dis(w.list.o, f.w.o, penCancersM$OC)


af.aj <- setNames(c(0.01366243, 0.01168798), genes) ## AJ allele frequency
af.naj <- setNames(c(0.0005829, 0.0006760), genes) ## non-AJ allele frequency

n.parallel <- 100
n.fam <- 1e3


## Generating the families
start <- Sys.time()
fam.sim <- foreach(i = 1:n.parallel, .combine = append) %dopar% {
  gen.fam.dis(n.fam, genes, cancers, penCancersF, penCancersM, af.aj, af.naj,
              supp.w = supp.w, f.w = f.w, hzd0.b = hzd0.b, hzd0.o = hzd0.o,
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


save(fam.sim, file = paste(getwd(), "/Frailty/Simulations/simdat_dis.RData", sep = ""))
save(fam.sim.bc, file = paste(getwd(), "/Frailty/Simulations/simdat_bc_dis.RData", sep = ""))




