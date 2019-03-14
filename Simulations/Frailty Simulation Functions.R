## Functions for the Frailty Simulations
## Last updated: November 28, 2018

## obtains the covariance matrix for a bivariate normal
mvn.sigma <- function(var, rho){
  matrix(c(var[1], rho * sqrt(var[1] * var[2]),
           rho * sqrt(var[1] * var[2]), var[2]), 2, 2)
}

## Obtaining the conditional penetrance from the
## conditional baseline hazard function
mpfa <- function(w, hzd0, nr = 94){
  ## mpfa = make penetrance frailty-adjusted
  ## hzd0 = conditional baseline hazard function
  hzd.cond <- 1 - (1 - hzd0)^exp(w) ## C-log-log link
  surv.cond <- apply(1 - hzd.cond, 2, cumprod)
  pen.cond <- surv.cond[1:(nr - 1), ] * hzd.cond[2:nr, ]
  pen.cond <- rbind(hzd.cond[1, ], pen.cond) ## S(0) = P(T > 0) = 1
  return(pen.cond)
}

## obtains the baseline hazards
cbhf.norm <- function(mu, var, pen, nr = 94){
  ## cbhf = conditional baseline hazard function
  ## pen is the marginal penetrance function
  ## mu is a 2-vector for the frailty means for BC and OC
  ## var is a 2-vector for the frailty variances for BC and OC
  hzd0.cond <- matrix(0, nr, 3)
  h.prev <- matrix(0, nr, 3)
  for(i in 1:nr){
    for(j in 1:3){
      if(i == 1){
        hzd0.cond[i, j] <- uniroot(function(h){
          integrate(f = function(w) (1 - (1 - h)^exp(w)) * dnorm(w, mean = mu[1], sd = sqrt(var[1])),
                    lower = -Inf, upper = Inf)$value - pen[i, j]
        }, interval = c(0, 1), tol = 1e-20)$root
        h.prev[i, j] <- hzd0.cond[i, j]
        # hfunc <- nlminb(0.5, heq, lower = 0, upper = 1,
        #                 surv = surv[j, i], hprev = 0, supp.w = supp.w, f.w = f.w,
        #                 dis = dis)
      } else {
        # hfunc <- nlminb(0.5, heq, lower = 0, upper = 1,
        #                 surv = surv[j, i], hprev = h.prev[1:(j - 1), i],
        #                 supp.w = supp.w, f.w = f.w, dis = dis)
        hzd0.cond[i, j] <- uniroot(function(h){
          integrate(f = function(w) (1 - (1 - h)^exp(w)) * prod(1 - h.prev[1:(i - 1), j])^exp(w) *
                      dnorm(w, mean = mu[1], sd = sqrt(var[1])),
                    lower = -Inf, upper = Inf)$value - pen[i, j]
        }, interval = c(0, 1), tol = 1e-20)$root
        h.prev[i, j] <- hzd0.cond[i, j]
      }
      # hzd0.cond[j, i] <- h.prev[j, i] <- hfunc$par
    }
  }
  return(hzd0.cond)
}

cbhf.dis <- function(w.list, f.w, pen, nr = 94){
  ## cbhf = conditional baseline hazard function
  ## pen is the marginal penetrance function
  hzd0.cond <- matrix(0, nr, 3)
  h.prev <- matrix(0, nr, 3)
  for(i in 1:nr){
    for(j in 1:3){
      if(i == 1){
        hzd0.cond[i, j] <- uniroot(function(h){
          sum((1 - (1 - h)^exp(w.list)) * f.w) - pen[i, j]
        }, interval = c(0, 1), tol = 1e-20)$root
        h.prev[i, j] <- hzd0.cond[i, j]
        # hfunc <- nlminb(0.5, heq, lower = 0, upper = 1,
        #                 surv = surv[j, i], hprev = 0, supp.w = supp.w, f.w = f.w,
        #                 dis = dis)
      } else {
        # hfunc <- nlminb(0.5, heq, lower = 0, upper = 1,
        #                 surv = surv[j, i], hprev = h.prev[1:(j - 1), i],
        #                 supp.w = supp.w, f.w = f.w, dis = dis)
        hzd0.cond[i, j] <- uniroot(function(h){
          sum((1 - (1 - h)^exp(w.list)) * prod(1 - h.prev[1:(i - 1), j])^exp(w.list) *
                f.w) - pen[i, j]
        }, interval = c(0, 1), tol = 1e-20)$root
        h.prev[i, j] <- hzd0.cond[i, j]
      }
      # hzd0.cond[j, i] <- h.prev[j, i] <- hfunc$par
    }
  }
  return(hzd0.cond)
}





## getting the Cancer Penetrance from the frailty-adjusted penetrance
genCP <- function(genes, cancers, penCancersF, penCancersM, w, hzd0.b, hzd0.o){
  ## mu is a 2-vector for the frailty means for BC and OC
  ## var is a 2-vector for the frailty variances for BC and OC
  penCancersF.cond <- penCancersF
  penCancersF.cond$BC <- mpfa(w[1], hzd0.b)
  penCancersF.cond$OC <- mpfa(w[2], hzd0.o)
  penCancersM.cond <- penCancersM
  penCancersM.cond$BC <- mpfa(w[1], hzd0.b.m)
  penCancersM.cond$OC <- mpfa(w[2], hzd0.o.m)
  CP.cond <- genCancerPen(genes, cancers, penCancersF.cond, penCancersM.cond, maxK = length(genes), age.last = 95)
  return(CP.cond)
}

## generating the families
gen.fam.norm <- function(n.sim, genes, cancers, penCancersF, penCancersM, af.aj, af.naj, p.aj = 5.5 / 285,
                        mu, var, rho, hzd0.b, hzd0.o, seed = NULL, age.max = 94, age.min = 7,
                        censoring = TRUE, t.fu = 5){
  ## t.fu = number of years of follow-up for risk prediction
  if(!is.null(seed)) set.seed(seed)
  
  w <- rmvnorm(n.sim, mu, mvn.sigma(var, rho))
  aj <- rbinom(n.sim, 1, p.aj)
  nSibsPatern <- matrix(sample(0:3, 2 * n.sim, replace = TRUE), n.sim, byrow = TRUE)
  nSibsMatern <- matrix(sample(0:3, 2 * n.sim, replace = TRUE), n.sim, byrow = TRUE)
  nSibs <- matrix(sample(0:3, 2 * n.sim, replace = TRUE), n.sim, byrow = TRUE)
  
  fam.sim <- list()
  for(i in 1:n.sim){
    ## AJ status
    if(aj[i] == 1){
      af <- af.aj
    } else{
      af <- af.naj
    }
    
    ## generating family from frailty-adjusted penetrance and family-specific allele frequency
    nGrandchild <- matrix(sample(0:3, 2 * (sum(nSibs[i, ]) + 1), replace = TRUE), sum(nSibs[i, ]) + 1, 2)
    fam.sim[[i]] <- tryCatch(sim.simFam(nSibsPatern[i, ], nSibsMatern[i, ], nSibs[i, ],
                                        nGrandchild, af,
                                        genCP(genes, cancers, penCancersF, penCancersM, w[i, ], hzd0.b, hzd0.o),
                                        includeGeno = TRUE, age.max = age.max,
                                        age.min = age.min, censoring = censoring),
                             error = function(e) NULL)
    fam.sim[[i]]$W.BC <- w[i, 1]
    fam.sim[[i]]$W.OC <- w[i, 2]
    fam.sim[[i]]$ethnic <- ifelse(aj[i] == 1, "AJ", "nonAJ")
    names(fam.sim[[i]])[names(fam.sim[[i]]) %in% c("isAffBC", "isAffOC", "AgeBC", "AgeOC")] <-
      c("AffectedBreast", "AffectedOvary", "AgeBreast", "AgeOvary")
    fam.sim[[i]]$AgeBreastContralateral <- 0
    fam.sim[[i]]$Twins <- 0
    fam.sim[[i]]$RELATION <- ifelse(fam.sim[[i]]$isProband == 1, 0, 1)
    fam.sim[[i]]$race <- "Unknown"
    
    ## for 5-year risk prediction, using the ages and affection statuses as the
    ## follow-up information and going back 5 years
    fam.sim[[i]]$AgeBreast.fu <- fam.sim[[i]]$AgeBreast
    fam.sim[[i]]$AgeOvary.fu <- fam.sim[[i]]$AgeOvary
    fam.sim[[i]]$AgeBreastContralateral.fu <- 0
    fam.sim[[i]]$AffectedBreast.fu <- fam.sim[[i]]$AffectedBreast
    fam.sim[[i]]$AffectedOvary.fu <- fam.sim[[i]]$AffectedOvary
    fam.sim[[i]]$AgeBreast <- ifelse(fam.sim[[i]]$AffectedBreast.fu == 0, fam.sim[[i]]$AgeBreast.fu - t.fu,
                                     ifelse(fam.sim[[i]]$AgeBreast.fu > fam.sim[[i]]$CurAge - t.fu, fam.sim[[i]]$CurAge - t.fu,
                                            fam.sim[[i]]$AgeBreast.fu))
    fam.sim[[i]]$AgeOvary <- ifelse(fam.sim[[i]]$AffectedOvary.fu == 0, fam.sim[[i]]$AgeOvary.fu - t.fu,
                                    ifelse(fam.sim[[i]]$AgeOvary.fu > fam.sim[[i]]$CurAge - t.fu, fam.sim[[i]]$CurAge - t.fu,
                                           fam.sim[[i]]$AgeOvary.fu))
    fam.sim[[i]]$AffectedBreast <- ifelse(fam.sim[[i]]$AffectedBreast.fu == 0, 0,
                                          ifelse(fam.sim[[i]]$AgeBreast.fu > fam.sim[[i]]$CurAge - t.fu, 0, 1))
    fam.sim[[i]]$AffectedOvary <- ifelse(fam.sim[[i]]$AffectedOvary.fu == 0, 0,
                                          ifelse(fam.sim[[i]]$AgeOvary.fu > fam.sim[[i]]$CurAge - t.fu, 0, 1))
  }
  return(fam.sim)
}

## using discrete distribution
gen.fam.dis <- function(n.sim, genes, cancers, penCancersF, penCancersM, af.aj, af.naj, p.aj = 5.5 / 285,
                         supp.w, f.w, hzd0.b, hzd0.o, seed = NULL, age.max = 94, age.min = 7,
                         censoring = TRUE, t.fu = 5){
  ## t.fu = number of years of follow-up for risk prediction
  if(!is.null(seed)) set.seed(seed)
  
  w <- as.matrix(supp.w[base::sample(1:49, n.sim, replace = TRUE, prob = f.w), ])
  aj <- rbinom(n.sim, 1, p.aj)
  nSibsPatern <- matrix(sample(0:3, 2 * n.sim, replace = TRUE), n.sim, byrow = TRUE)
  nSibsMatern <- matrix(sample(0:3, 2 * n.sim, replace = TRUE), n.sim, byrow = TRUE)
  nSibs <- matrix(sample(0:3, 2 * n.sim, replace = TRUE), n.sim, byrow = TRUE)
  
  fam.sim <- list()
  for(i in 1:n.sim){
    ## AJ status
    if(aj[i] == 1){
      af <- af.aj
    } else{
      af <- af.naj
    }
    
    ## generating family from frailty-adjusted penetrance and family-specific allele frequency
    nGrandchild <- matrix(sample(0:3, 2 * (sum(nSibs[i, ]) + 1), replace = TRUE), sum(nSibs[i, ]) + 1, 2)
    fam.sim[[i]] <- tryCatch(sim.simFam(nSibsPatern[i, ], nSibsMatern[i, ], nSibs[i, ],
                                        nGrandchild, af,
                                        genCP(genes, cancers, penCancersF, penCancersM, w[i, ], hzd0.b, hzd0.o),
                                        includeGeno = TRUE, age.max = age.max,
                                        age.min = age.min, censoring = censoring),
                             error = function(e) NULL)
    fam.sim[[i]]$W.BC <- w[i, 1]
    fam.sim[[i]]$W.OC <- w[i, 2]
    fam.sim[[i]]$ethnic <- ifelse(aj[i] == 1, "AJ", "nonAJ")
    names(fam.sim[[i]])[names(fam.sim[[i]]) %in% c("isAffBC", "isAffOC", "AgeBC", "AgeOC")] <-
      c("AffectedBreast", "AffectedOvary", "AgeBreast", "AgeOvary")
    fam.sim[[i]]$AgeBreastContralateral <- 0
    fam.sim[[i]]$Twins <- 0
    fam.sim[[i]]$RELATION <- ifelse(fam.sim[[i]]$isProband == 1, 0, 1)
    fam.sim[[i]]$race <- "Unknown"
    
    ## for 5-year risk prediction, using the ages and affection statuses as the
    ## follow-up information and going back 5 years
    fam.sim[[i]]$AgeBreast.fu <- fam.sim[[i]]$AgeBreast
    fam.sim[[i]]$AgeOvary.fu <- fam.sim[[i]]$AgeOvary
    fam.sim[[i]]$AgeBreastContralateral.fu <- 0
    fam.sim[[i]]$AffectedBreast.fu <- fam.sim[[i]]$AffectedBreast
    fam.sim[[i]]$AffectedOvary.fu <- fam.sim[[i]]$AffectedOvary
    fam.sim[[i]]$AgeBreast <- ifelse(fam.sim[[i]]$AffectedBreast.fu == 0, fam.sim[[i]]$AgeBreast.fu - t.fu,
                                     ifelse(fam.sim[[i]]$AgeBreast.fu > fam.sim[[i]]$CurAge - t.fu, fam.sim[[i]]$CurAge - t.fu,
                                            fam.sim[[i]]$AgeBreast.fu))
    fam.sim[[i]]$AgeOvary <- ifelse(fam.sim[[i]]$AffectedOvary.fu == 0, fam.sim[[i]]$AgeOvary.fu - t.fu,
                                    ifelse(fam.sim[[i]]$AgeOvary.fu > fam.sim[[i]]$CurAge - t.fu, fam.sim[[i]]$CurAge - t.fu,
                                           fam.sim[[i]]$AgeOvary.fu))
    fam.sim[[i]]$AffectedBreast <- ifelse(fam.sim[[i]]$AffectedBreast.fu == 0, 0,
                                          ifelse(fam.sim[[i]]$AgeBreast.fu > fam.sim[[i]]$CurAge - t.fu, 0, 1))
    fam.sim[[i]]$AffectedOvary <- ifelse(fam.sim[[i]]$AffectedOvary.fu == 0, 0,
                                         ifelse(fam.sim[[i]]$AgeOvary.fu > fam.sim[[i]]$CurAge - t.fu, 0, 1))
  }
  return(fam.sim)
}




