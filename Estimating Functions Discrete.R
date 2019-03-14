## Estimating Functions for Discrete Frailty, Bayesian Analysis
## Using BayesMendel package v2.1-5 -- crude/net penetrances
## Last updated: March 13, 2019

## Obtaining the conditional baseline hazard function
cbhf <- function(w.list, f.w, pen, nr = 94){
  ## cbhf = conditional baseline hazard function
  ## pen is the marginal penetrance function
  hzd0.cond <- matrix(0, nr, 5)
  h.prev <- matrix(0, nr, 5)
  for(i in 1:nr){
    for(j in c(1, 2, 4, 5)){
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

## This runs BRCAPRO, incorporating germline testing and interventions
## It doesn't include race since this would overwrite the penetrance
brcapro.fam <- function(fam, fut, id0, newpenet, race = NA, seed = 1,
                        net = TRUE, gt = TRUE){
  set.seed(seed)
  ## fut = follow-up time
  ## id0 = ID of proband
  ## newpenet = modified penetrance
  ## net = net risk
  ## gt = use genetic testing results
  brca.var <- c("ID", "Gender", "FatherID", "MotherID", "AffectedBreast",
                "AffectedOvary", "AgeBreast", "AgeOvary", "AgeBreastContralateral",
                "Twins", "ethnic", "Death", "AgeDeath")
  if(gt == TRUE){
    if(all(c("BRCA1", "BRCA2") %in% names(fam))){
      if(any(fam$BRCA1 %in% 1:2 | fam$BRCA2 %in% 1:2)){
        germline.testing <- data.frame(BRCA1 = fam$BRCA1, BRCA2 = fam$BRCA2, TestOrder = rep(0, nrow(fam)))
      } else{
        germline.testing <- NULL
      } 
    } else{
      germline.testing <- NULL
    }
  } else{
    germline.testing <- NULL
  }
  if(all(c("Mastectomy", "AgeMastectomy") %in% names(fam))){
    mastectomy <- data.frame(Mastectomy = fam$Mastectomy, AgeMastectomy = fam$AgeMastectomy) 
  } else{
    mastectomy <- NULL
  }
  if(all(c("Oophorectomy", "AgeOophorectomy") %in% names(fam))){
    oophorectomy <- data.frame(Oophorectomy = fam$Oophorectomy, AgeOophorectomy = fam$AgeOophorectomy)
  } else {
    oophorectomy <- NULL
  }
  
  if(net == TRUE){
    params <- brcaparams(penetrance.net = newpenet, age.by = fut)
  } else{
    params <- brcaparams(penetrance.crude = newpenet, age.by = fut)
  }
  
  if(is.na(race)){
    result <- brcapro(family = fam[, which(names(fam) %in% brca.var)], counselee.id = id0,
                      germline.testing = germline.testing,
                      oophorectomy = oophorectomy, mastectomy = mastectomy,
                      params = params,
                      print = FALSE, net = net)
  } else{
    result <- brcapro(family = fam[, which(names(fam) %in% brca.var)], counselee.id = id0,
                      germline.testing = germline.testing,
                      oophorectomy = oophorectomy, mastectomy = mastectomy,
                      params = params,
                      print = FALSE, race = race, net = net)
  }
  return(result)
}

race.suf <- function(dat, gen, cancer, type, cond = TRUE){
  ## returns the conditional penetrance or survival for the corresponding race,
  ## given the gender and the cancer
  ## dat = family data
  ## gen = gender ("Female" or "Male")
  ## cancer is either "Breast" or "Ovarian"
  ## Cond is TRUE if we want conditional (on the frailty). If FALSE,
  ## we obtain the marginal (over the frailty).
  if(dat$race[1] == "Asian"){
    rc <- ".as"
  } else if(dat$race[1] == "Black"){
    rc <- ".bl"
  } else if(dat$race[1] == "Hispanic"){
    rc <- ".his"
  } else if(dat$race[1] == "NativeAmerican"){
    rc <- ".na"
  } else if(dat$race[1] == "White"){
    rc <- ".wh"
  } else {
    rc <- ""
  }
  if(cancer == "Breast"){
    if(gen == "Female"){
      if(type == "pen"){
        if(cond == TRUE){
          return(eval(parse(text = paste("pen.b.cond", rc, sep = ""))))
        } else{
          return(eval(parse(text = paste("pen.b.mar", rc, sep = ""))))
        }
      } else if(type == "surv"){
        if(cond == TRUE){
          return(eval(parse(text = paste("surv.b.cond", rc, sep = ""))))
        } else{
          return(eval(parse(text = paste("surv.b.mar", rc, sep = ""))))
        }
      }
    } else if(gen == "Male"){
      if(type == "pen"){
        if(cond == TRUE){
          return(eval(parse(text = paste("pen.b.cond.m", rc, sep = ""))))
        } else{
          return(eval(parse(text = paste("pen.b.mar.m", rc, sep = ""))))
        }
      } else if(type == "surv"){
        if(cond == TRUE){
          return(eval(parse(text = paste("surv.b.cond.m", rc, sep = ""))))
        } else{
          return(eval(parse(text = paste("surv.b.mar.m", rc, sep = ""))))
        }
      }
    }
  } else if(cancer == "Ovarian"){
    if(gen == "Female"){
      if(type == "pen"){
        if(cond == TRUE){
          return(eval(parse(text = paste("pen.o.cond", rc, sep = ""))))
        } else{
          return(eval(parse(text = paste("pen.o.mar", rc, sep = ""))))
        }
      } else if(type == "surv"){
        if(cond == TRUE){
          return(eval(parse(text = paste("surv.o.cond", rc, sep = ""))))
        } else{
          return(eval(parse(text = paste("surv.o.mar", rc, sep = ""))))
        }
      }
    }
  }
}

## obtaining the family with the age imputations
fam.new <- function(fam){
  dat.new <- fam
  rc <- dat.new$race[1]
  rel <- dat.new$RELATION
  id.rel <- dat.new$ID
  id.pro <- dat.new$ID[dat.new$RELATION == 0]
  if(!is.null(dat.new$Mastectomy)){
    mast <- dat.new$Mastectomy
  }
  if(!is.null(dat.new$Oophorectomy)){
    ooph <- dat.new$Oophorectomy
  }
  if(!is.null(dat.new$AgeMastectomy)){
    agemast <- dat.new$AgeMastectomy
  }
  if(!is.null(dat.new$AgeOophorectomy)){
    ageooph <- dat.new$AgeOophorectomy
  }
  if(!is.null(dat.new$BRCA1)){
    brca1 <- dat.new$BRCA1
  }
  if(!is.null(dat.new$BRCA2)){
    brca2 <- dat.new$BRCA2
  }
  ## imputing missing ages
  dat.new <- ImputeAge(
    fff = CheckFamStructure(model = "brcapro", fff = as.data.frame(dat.new), counselee.id = id.pro,
                            germline.testing = NULL, marker.testing = NULL,
                            oophorectomy = NULL, mastectomy = NULL,
                            imputeAges = TRUE, imputeRelatives = TRUE,
                            params = brcaparams()),
    params = brcaparams(), model = "brcapro"
  )$fff
  ## CheckFamStructure will remove the some columns if it needs to add relatives
  ## Here we add them back in because they are used in the estimating equation function
  ## and in brcapro
  dat.new$race <- rc
  dat.new$RELATION <- 999; dat.new$RELATION[dat.new$ID %in% id.rel] <- rel
  if(!is.null(fam$Mastectomy)){
    dat.new$Mastectomy <- 0; dat.new$Mastectomy[dat.new$ID %in% id.rel] <- mast
  }
  if(!is.null(fam$Oophorectomy)){
    dat.new$Oophorectomy <- 0; dat.new$Oophorectomy[dat.new$ID %in% id.rel] <- ooph
  }
  if(!is.null(fam$AgeMastectomy)){
    dat.new$AgeMastectomy <- 1; dat.new$AgeMastectomy[dat.new$ID %in% id.rel] <- agemast
  }
  if(!is.null(fam$AgeOophorectomy)){
    dat.new$AgeOophorectomy <- 1; dat.new$AgeOophorectomy[dat.new$ID %in% id.rel] <- ageooph
  }
  if(!is.null(dat.new$BRCA1)){
    dat.new$BRCA1 <- 0; dat.new$BRCA1[dat.new$ID %in% id.rel] <- brca1
  } 
  if(!is.null(dat.new$BRCA2)){
    dat.new$BRCA2 <- 0; dat.new$BRCA2[dat.new$ID %in% id.rel] <- brca2
  }
  dat.new$Death <- rep(NA, nrow(dat.new))
  dat.new$AgeDeath <- rep(NA, nrow(dat.new))
  return(dat.new)
}

## changing the family IDs so that they correspond to the row numbers
## also adding fathers/mothers if family members have fathers listed but not mothers,
## or vice versa
changeID <- function(fam){
  ## first re-ordering IDs
  fam$ID2 <- 1:nrow(fam)
  mid <- fid <- rep(0, nrow(fam))
  mid[fam$MotherID != 0] <- match(fam$MotherID, fam$ID)[fam$MotherID != 0]
  fid[fam$FatherID != 0] <- match(fam$FatherID, fam$ID)[fam$FatherID != 0]
  fam$MotherID <- mid
  fam$FatherID <- fid
  fam$ID <- fam$ID2
  fam$ID2 <- NULL
  
  ## now adding fathers/mothers
  n.fam <- nrow(fam)
  ct <- 0
  fid.miss <- mid.miss <- matrix(nrow = 2, ncol = 0)
  for(i in 1:n.fam){
    if(fam$MotherID[i] == 0 & fam$FatherID[i] != 0){
      if(fam$FatherID[i] %in% fid.miss[1, ]){
        ## if we already created the missing mother, use that
        fam$MotherID[i] <- fid.miss[2, which(fid.miss[1, ] == fam$FatherID[i])]
      } else{
        ## if not, create the missing mother
        ct <- ct + 1
        fid.miss <- cbind(fid.miss, c(fam$FatherID[i], n.fam + ct))
        fam <- rbind(fam, fam[1, ])
        fam$ID[n.fam + ct] <- n.fam + ct
        fam$RELATION[n.fam + ct] <- 99
        fam$Gender[n.fam + ct] <- 0
        fam$Twins[n.fam + ct] <- 0
        fam$MotherID[n.fam + ct] <- 0
        fam$FatherID[n.fam + ct] <- 0
        fam$AffectedBreast[n.fam + ct] <- 0
        fam$AffectedOvary[n.fam + ct] <- 0
        fam$AgeBreast[n.fam + ct] <- NA
        fam$AgeOvary[n.fam + ct] <- NA
        fam$AgeBreastContralateral[n.fam + ct] <- NA
        fam$Death[n.fam + ct] <- NA
        fam$AgeDeath[n.fam + ct] <- NA
        fam$MotherID[i] <- n.fam + ct
      }
    } else if(fam$MotherID[i] != 0 & fam$FatherID[i] == 0){
      if(fam$MotherID[i] %in% mid.miss[1, ]){
        ## if we already created the missing father, use that
        fam$FatherID[i] <- mid.miss[2, which(mid.miss[1, ] == fam$MotherID[i])]
      } else{
        ## if not, create the missing father
        ct <- ct + 1
        mid.miss <- cbind(mid.miss, c(fam$MotherID[i], n.fam + ct))
        fam <- rbind(fam, fam[1, ])
        fam$ID[n.fam + ct] <- n.fam + ct
        fam$RELATION[n.fam + ct] <- 99
        fam$Gender[n.fam + ct] <- 1
        fam$Twins[n.fam + ct] <- 0
        fam$MotherID[n.fam + ct] <- 0
        fam$FatherID[n.fam + ct] <- 0
        fam$AffectedBreast[n.fam + ct] <- 0
        fam$AffectedOvary[n.fam + ct] <- 0
        fam$AgeBreast[n.fam + ct] <- NA
        fam$AgeOvary[n.fam + ct] <- NA
        fam$AgeBreastContralateral[n.fam + ct] <- NA
        fam$Death[n.fam + ct] <- NA
        fam$AgeDeath[n.fam + ct] <- NA
        fam$FatherID[i] <- n.fam + ct
      }
    }
  }
  return(fam)
}

## obtaining the likelihood matrix, given the frailty-adjusted penetrance
getLIK <- function(fam, pen, surv){
  LIK <- matrix(0, nrow(fam), 4)
  for(i in 1:nrow(fam)){
    for(j in c(1, 2, 4, 5)){
      if(fam$Gender[i] == 0){
        LIK[i, which(c(1, 2, 4, 5) == j)] <- pen$fFX[fam$AgeBreast[i], j]^(fam$AffectedBreast[i]) *
          surv$fFX[fam$AgeBreast[i], j]^(1 - fam$AffectedBreast[i]) *
          pen$fFY[fam$AgeOvary[i], j]^(fam$AffectedOvary[i]) *
          surv$fFY[fam$AgeOvary[i], j]^(1 - fam$AffectedOvary[i])
      } else{
        LIK[i, which(c(1, 2, 4, 5) == j)] <- pen$fMX[fam$AgeBreast[i], j]^(fam$AffectedBreast[i]) *
          surv$fMX[fam$AgeBreast[i], j]^(1 - fam$AffectedBreast[i])
      }
    }
  }
  return(LIK)
}


## posterior distribution
post <- function(w, f.w, dat, seed = 1, net = TRUE, nIter = 100){
  ## w = frailty vector
  ## dat = family data
  ## f.w = discrete frailty pmf
  
  ## setting seed because of possible age imputation
  set.seed(seed)
  
  age.max <- 94
  dat$AgeBreast[dat$AgeBreast > age.max] <- age.max
  dat$AgeBreastContralateral[dat$AgeBreastContralateral > age.max] <- age.max
  dat$AgeOvary[dat$AgeOvary > age.max] <- age.max
  if(!is.null(dat$AgeMasteectomy)){
    dat$AgeMastectomy[dat$AgeMastectomy > age.max] <- age.max
  }
  if(!is.null(dat$AgeOophorectomy)){
    dat$AgeOophorectomy[dat$AgeOophorectomy > age.max] <- age.max
  }
  
  ## getting the frailty-adjusted penetrance
  ## these are net penetrances since we are getting carrier probabilities
  pen.b.dat <- race.suf(dat, "Female", "Breast", "pen")
  pen.b.m.dat <- race.suf(dat, "Male", "Breast", "pen")
  pen.o.dat <- race.suf(dat, "Female", "Ovarian", "pen")
  
  newpenet <- penet.brca.net
  newpenet$fFX[, 1:5] <- pen.b.dat[, , which(w.list.b == w[1])]
  newpenet$fMX[, 1:5] <- pen.b.m.dat[, , which(w.list.b == w[1])]
  newpenet$fFY[, 1:5] <- pen.o.dat[, , which(w.list.o == w[2])]
  
  ## if penetrance has a 0, change it to 1e-5 (to avoid errors in peeling)
  for(k in c(1, 3)){
    newpenet[[k]][, c(1, 2, 4, 5)][which(newpenet[[k]][, c(1, 2, 4, 5)] == 0)] <- 1e-5
  }
  ## if penetrance sums to 1, normalizing so that it sums to 0.99999
  ## (to avoid errors in peeling)
  for(j in 1:3){
    if(all(!is.na(newpenet[[j]]))){
      if(any(colSums(newpenet[[j]]) >= 1)){
        for(k in c(1, 2, 4, 5)){
          if(sum(newpenet[[j]][, k]) >= 1){
            newpenet[[j]][, k] <- newpenet[[j]][, k] * 0.99999 / sum(newpenet[[j]][, k])
          }
        }
      }
    }
  }
  
  newsurv <- list(fFX = 1 - apply(newpenet$fFX, 2, cumsum),
                  fMX = 1 - apply(newpenet$fMX, 2, cumsum),
                  fFY = 1 - apply(newpenet$fFY, 2, cumsum))
  
  
  ## re-ordering the IDs and imputing ages
  ## first use fam.new to avoid issues with Mother/FatherIDs that don't correspond to any ID
  ped <- fam.new(changeID(fam.new(dat)))
  
  ## multiple imputation for unknown affected ages
  ped$isProband <- 0; ped$isProband[ped$RELATION == 0] <- 1
  if(ped$ethnic[1] == "nonAJ"){
    prevs <- c(1 - (1 - 0.0005829)^2, 1 - (1 - 0.0006760)^2)
  } else if(ped$ethnic[1] == "Italian"){
    prevs <- c(1 - (1 - 0.001673779)^2, 1 - (1 - 0.00132622)^2)
  } else{
    prevs <- c(1 - (1 - 0.01366243)^2, 1 - (1 - 0.01168798)^2)
  }
  
  age.miss <- any((ped$AffectedBreast %in% 1:2 & ped$AgeBreast == 1) |
                    (ped$AffectedOvary == 1 & ped$AgeOvary == 1))
  
  if(age.miss == TRUE){
    p.hw.impute <- rep(0, nIter)
    ifamily = ImputeAge(fff = CheckFamStructure(model = "brcapro", fff = ped, counselee.id = filter(ped, isProband == 1)$ID,
                                                germline.testing = NULL, marker.testing = NULL,
                                                oophorectomy = NULL, mastectomy = NULL,
                                                imputeAges = TRUE, imputeRelatives = TRUE,
                                                params = brcaparams(penetrance.net = newpenet)),
                        params = brcaparams(penetrance.net = newpenet),
                        model = "brcapro")
    for(i in 1:nIter){
      
      # Run each imputed age through runPeeling if there are affected relatives
      # who had age imputed
      nf <- ifamily$fff
      nf$AgeBreast[ifamily$fff$AffectedBreast == 1 & ifamily$fff$AgeBreast == 1] <-
        ifamily$age_bc[i, ifamily$fff$AffectedBreast == 1 & ifamily$fff$AgeBreast == 1]
      nf$AgeBreast[ifamily$fff$AffectedBreast == 2 & ifamily$fff$AgeBreast == 1] <-
        ifamily$age_bcc[i, ifamily$fff$AffectedBreast == 2 & ifamily$fff$AgeBreast == 1]
      nf$AgeOvary[ifamily$fff$AffectedOvary == 1 & ifamily$fff$AgeOvary == 1] <-
        ifamily$age_oc[i, ifamily$fff$AffectedOvary == 1 & ifamily$fff$AgeOvary == 1]
      
      LIK <- getLIK(nf, newpenet, newsurv)
      p.hw.impute[i] <- sum(pp.peelingParing(nf, prevs, LIK, T = 2, normalize = FALSE))
    }
    p.hw <- mean(p.hw.impute)
  } else{
    LIK <- getLIK(ped, newpenet, newsurv)
    p.hw <- sum(pp.peelingParing(ped, prevs, LIK, T = 2, normalize = FALSE))
  }
  ind.w <- which(supp.w$W.BC == w[1] & supp.w$W.OC == w[2])
  return(f.w[ind.w] * p.hw)
}

## posterior distribution for each family
post.fam <- function(dat, i, int, net = TRUE){
  temp <- vector()
  for(j in 1:length(int)){
    temp <- c(temp,
              tryCatch(post(as.numeric(supp.w[i, ]), f.w,
                            filter(dat, FamID == unique(dat$FamID)[j]),
                            net = net), error = function(e) NA))
  }
  return(temp)
}

## risk predictions for each family
rp.fam <- function(dat, i, int, res.post, net = TRUE, gt = TRUE){
  fam <- tryCatch(fam.new(dat[dat$FamID == famid.list[int][i], ]), error = function(e) NULL)
  if(is.null(fam)){
    return(c(famid.list[int][i], rep(NA, 11)))
  }
  fam2 <- dat[dat$FamID == famid.list[int][i], ]
  id.pro <- fam2$ID[fam2$RELATION == 0]
  
  ## return NAs if proband is older than 93
  if(fam2$AgeBreast[fam2$ID == id.pro] > 93){
    return(c(fam2$FamID[1], rep(NA, 10)))
  }
  
  res.brca <- data.frame(matrix(0, nrow(supp.w), 5))
  names(res.brca) <- c("P.BC.5", "P.OC.5", "P.BRCA", "P.BRCA1", "P.BRCA2")
  for(ww in 1:nrow(supp.w)){
    ind.wb <- which(w.list.b == supp.w[ww, 1])
    ind.wo <- which(w.list.o == supp.w[ww, 2])
    
    ## getting the frailty-adjusted penetrance
    pen.b.dat <- race.suf(dat, "Female", "Breast", "pen")
    pen.b.m.dat <- race.suf(dat, "Male", "Breast", "pen")
    pen.o.dat <- race.suf(dat, "Female", "Ovarian", "pen")
    
    newpenet <- penet.brca.net
    newpenet$fFX[, 1:5] <- pen.b.cond[, , ind.wb]
    newpenet$fMX[, 1:5] <- pen.b.cond.m[, , ind.wb]
    newpenet$fFY[, 1:5] <- pen.o.cond[, , ind.wo]
    
    ## if penetrance has a 0, change it to 1e-5 (to avoid errors in peeling)
    for(j in 1:3){
      if(all(!is.na(newpenet[[j]]))){
        newpenet[[j]][, c(1, 2, 4, 5)][which(newpenet[[j]][, c(1, 2, 4, 5)] == 0)] <- 1e-5
      }
    }
    ## if penetrance sums to 1, normalizing so that it sums to 0.99999
    ## (to avoid errors in peeling)
    for(j in 1:3){
      if(all(!is.na(newpenet[[j]]))){
        if(any(colSums(newpenet[[j]]) >= 1)){
          for(k in c(1, 2, 4, 5)){
            if(sum(newpenet[[j]][, k]) >= 1){
              newpenet[[j]][, k] <- newpenet[[j]][, k] * 0.99999 / sum(newpenet[[j]][, k])
            }
          }
        }
      }
    }
    
    result <- brcapro.fam(fam, fut = 5, id0 = id.pro, newpenet = newpenet, net = net, gt = gt)
    res.brca$P.BC.5[ww] <- tryCatch(result@predictions[1, 2], error = function(e) NA)
    res.brca$P.OC.5[ww] <- tryCatch(result@predictions[1, 3], error = function(e) NA) 
    res.brca$P.BRCA[ww] <- tryCatch(result@probs[1], error = function(e) NA)
    res.brca$P.BRCA1[ww] <- tryCatch(result@probs[2], error = function(e) NA)
    res.brca$P.BRCA2[ww] <- tryCatch(result@probs[3], error = function(e) NA)
  }
  
  ## normalizing the posterior distribution so it sums to 1 for each family
  res.post.norm <- res.post[, 1:nrow(supp.w)]
  for(k in 1:length(int)){
    res.post.norm[k, ] <-
      res.post[k, 1:nrow(supp.w)] / sum(res.post[k, 1:nrow(supp.w)])
  }
  
  res.fam <- rep(0, 11)
  res.fam[1] <- fam2$FamID[1]
  
  ## results with frailty
  res.fam[2:6] <- as.numeric(c(sum(res.brca$P.BC.5 * res.post.norm[i, ]),
                               sum(res.brca$P.OC.5 * res.post.norm[i, ]),
                               sum(res.brca$P.BRCA * res.post.norm[i, ]),
                               sum(res.brca$P.BRCA1 * res.post.norm[i, ]),
                               sum(res.brca$P.BRCA2 * res.post.norm[i, ])))
  
  ## results without frailty
  pen.b.dat.mar <- race.suf(dat, "Female", "Breast", "pen", cond = FALSE)
  pen.b.m.dat.mar <- race.suf(dat, "Male", "Breast", "pen", cond = FALSE)
  pen.o.dat.mar <- race.suf(dat, "Female", "Ovarian", "pen", cond = FALSE)
  
  defaultpenet <- penet.brca.net
  defaultpenet$fFX[, 1:5] <- pen.b.dat.mar
  defaultpenet$fMX[, 1:5] <- pen.b.m.dat.mar
  defaultpenet$fFY[, 1:5] <- pen.o.dat.mar
  
  res.nf <- brcapro.fam(fam, fut = 5, id0 = id.pro,
                        newpenet = defaultpenet, net = net, gt = gt)
  res.fam[7:11] <- as.numeric(c(tryCatch(res.nf@predictions[1, 2], error = function(e) NA),
                                tryCatch(res.nf@predictions[1, 3], error = function(e) NA),
                                tryCatch(res.nf@probs[1], error = function(e) NA),
                                tryCatch(res.nf@probs[2], error = function(e) NA),
                                tryCatch(res.nf@probs[3], error = function(e) NA)))
  
  return(res.fam)
}

