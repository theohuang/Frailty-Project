## Obtaining the baseline hazard functions
## Discrete frailty distribution
## Last Updated: November 16, 2018
## Theo Huang


## Loading packages
library(BayesMendel)

if(net == TRUE){
  penet <- penet.brca.net
  baseline.race <- baseline.race.brca.net
} else{
  penet <- penet.brca.crude
  baseline.race <- baseline.race.brca.crude
}


## Marginal (with respect to frailty variate) penetrance and survival function

pen.b.mar <- penet$fFX[, 1:5]
pen.b.mar.as <- pen.b.mar; pen.b.mar.as[, 1] <- baseline.race$fFX[, 1]
pen.b.mar.bl <- pen.b.mar; pen.b.mar.bl[, 1] <- baseline.race$fFX[, 2]
pen.b.mar.his <- pen.b.mar; pen.b.mar.his[, 1] <- baseline.race$fFX[, 3]
pen.b.mar.na <- pen.b.mar; pen.b.mar.na[, 1] <- baseline.race$fFX[, 4]
pen.b.mar.wh <- pen.b.mar; pen.b.mar.wh[, 1] <- baseline.race$fFX[, 5]

pen.b.mar.m <- penet$fMX[, 1:5]
pen.b.mar.m.as <- pen.b.mar.m; pen.b.mar.m.as[, 1] <- baseline.race$fMX[, 1]
pen.b.mar.m.bl <- pen.b.mar.m; pen.b.mar.m.bl[, 1] <- baseline.race$fMX[, 2]
pen.b.mar.m.his <- pen.b.mar.m; pen.b.mar.m.his[, 1] <- baseline.race$fMX[, 3]
pen.b.mar.m.na <- pen.b.mar.m; pen.b.mar.m.na[, 1] <- baseline.race$fMX[, 4]
pen.b.mar.m.wh <- pen.b.mar.m; pen.b.mar.m.wh[, 1] <- baseline.race$fMX[, 5]

pen.o.mar <- penet$fFY[, 1:5]
pen.o.mar.as <- pen.o.mar; pen.o.mar.as[, 1] <- baseline.race$fFY[, 1]
pen.o.mar.bl <- pen.o.mar; pen.o.mar.bl[, 1] <- baseline.race$fFY[, 2]
pen.o.mar.his <- pen.o.mar; pen.o.mar.his[, 1] <- baseline.race$fFY[, 3]
pen.o.mar.na <- pen.o.mar; pen.o.mar.na[, 1] <- baseline.race$fFY[, 4]
pen.o.mar.wh <- pen.o.mar; pen.o.mar.wh[, 1] <- baseline.race$fFY[, 5]

## Extracting conditional baseline hazard functions

## Female
hzd0.b.cond <- cbhf(w.list.b, f.w.b, pen.b.mar)
hzd0.b.cond.as <- cbhf(w.list.b, f.w.b, pen.b.mar.as)
hzd0.b.cond.bl <- cbhf(w.list.b, f.w.b, pen.b.mar.bl)
hzd0.b.cond.his <- cbhf(w.list.b, f.w.b, pen.b.mar.his)
hzd0.b.cond.na <- cbhf(w.list.b, f.w.b, pen.b.mar.na)
hzd0.b.cond.wh <- cbhf(w.list.b, f.w.b, pen.b.mar.wh)

hzd0.b.cond.m <- cbhf(w.list.b, f.w.b, pen.b.mar.m)
hzd0.b.cond.m.as <- cbhf(w.list.b, f.w.b, pen.b.mar.m.as)
hzd0.b.cond.m.bl <- cbhf(w.list.b, f.w.b, pen.b.mar.m.bl)
hzd0.b.cond.m.his <- cbhf(w.list.b, f.w.b, pen.b.mar.m.his)
hzd0.b.cond.m.na <- cbhf(w.list.b, f.w.b, pen.b.mar.m.na)
hzd0.b.cond.m.wh <- cbhf(w.list.b, f.w.b, pen.b.mar.m.wh)

hzd0.o.cond <- cbhf(w.list.o, f.w.o, pen.o.mar)
hzd0.o.cond.as <- cbhf(w.list.b, f.w.o, pen.o.mar.as)
hzd0.o.cond.bl <- cbhf(w.list.b, f.w.o, pen.o.mar.bl)
hzd0.o.cond.his <- cbhf(w.list.b, f.w.o, pen.o.mar.his)
hzd0.o.cond.na <- cbhf(w.list.b, f.w.o, pen.o.mar.na)
hzd0.o.cond.wh <- cbhf(w.list.b, f.w.o, pen.o.mar.wh)


pen.b.cond <- pen.b.cond.as <- pen.b.cond.bl <- pen.b.cond.his <-
  pen.b.cond.na <- pen.b.cond.wh <- pen.b.cond.m <- pen.b.cond.m.as <-
  pen.b.cond.m.bl <- pen.b.cond.m.his <- pen.b.cond.m.na <- pen.b.cond.m.wh <-
  surv.b.cond <- surv.b.cond.as <- surv.b.cond.bl <- surv.b.cond.his <-
  surv.b.cond.na <- surv.b.cond.wh <- surv.b.cond.m <- surv.b.cond.m.as <-
  surv.b.cond.m.bl <- surv.b.cond.m.his <- surv.b.cond.m.na <- surv.b.cond.m.wh <-
  array(0, dim = c(94, 5, length(w.list.b)))
for(i in 1:length(w.list.b)){
  pen.b.cond[, , i] <- mpfa(w.list.b[i], hzd0.b.cond)
  pen.b.cond.as[, , i] <- mpfa(w.list.b[i], hzd0.b.cond.as)
  pen.b.cond.bl[, , i] <- mpfa(w.list.b[i], hzd0.b.cond.bl)
  pen.b.cond.his[, , i] <- mpfa(w.list.b[i], hzd0.b.cond.his)
  pen.b.cond.na[, , i] <- mpfa(w.list.b[i], hzd0.b.cond.na)
  pen.b.cond.wh[, , i] <- mpfa(w.list.b[i], hzd0.b.cond.wh)
  pen.b.cond.m[, , i] <- mpfa(w.list.b[i], hzd0.b.cond.m)
  pen.b.cond.m.as[, , i] <- mpfa(w.list.b[i], hzd0.b.cond.m.as)
  pen.b.cond.m.bl[, , i] <- mpfa(w.list.b[i], hzd0.b.cond.m.bl)
  pen.b.cond.m.his[, , i] <- mpfa(w.list.b[i], hzd0.b.cond.m.his)
  pen.b.cond.m.na[, , i] <- mpfa(w.list.b[i], hzd0.b.cond.m.na)
  pen.b.cond.m.wh[, , i] <- mpfa(w.list.b[i], hzd0.b.cond.m.wh)
  surv.b.cond[, , i] <- 1 - apply(pen.b.cond[, , i], 2, cumsum)
  surv.b.cond.as[, , i] <- 1 - apply(pen.b.cond.as[, , i], 2, cumsum)
  surv.b.cond.bl[, , i] <- 1 - apply(pen.b.cond.bl[, , i], 2, cumsum)
  surv.b.cond.his[, , i] <- 1 - apply(pen.b.cond.his[, , i], 2, cumsum)
  surv.b.cond.na[, , i] <- 1 - apply(pen.b.cond.na[, , i], 2, cumsum)
  surv.b.cond.wh[, , i] <- 1 - apply(pen.b.cond.wh[, , i], 2, cumsum)
  surv.b.cond.m[, , i] <- 1 - apply(pen.b.cond.m[, , i], 2, cumsum)
  surv.b.cond.m.as[, , i] <- 1 - apply(pen.b.cond.m.as[, , i], 2, cumsum)
  surv.b.cond.m.bl[, , i] <- 1 - apply(pen.b.cond.m.bl[, , i], 2, cumsum)
  surv.b.cond.m.his[, , i] <- 1 - apply(pen.b.cond.m.his[, , i], 2, cumsum)
  surv.b.cond.m.na[, , i] <- 1 - apply(pen.b.cond.m.na[, , i], 2, cumsum)
  surv.b.cond.m.wh[, , i] <- 1 - apply(pen.b.cond.m.wh[, , i], 2, cumsum)
}
pen.o.cond <- pen.o.cond.as <- pen.o.cond.bl <- pen.o.cond.his <-
  pen.o.cond.na <- pen.o.cond.wh <- surv.o.cond <- surv.o.cond.as <-
  surv.o.cond.bl <- surv.o.cond.his <- surv.o.cond.na <- surv.o.cond.wh <-
  array(0, dim = c(94, 5, length(w.list.o)))
for(i in 1:length(w.list.o)){
  pen.o.cond[, , i] <- mpfa(w.list.o[i], hzd0.o.cond)
  pen.o.cond.as[, , i] <- mpfa(w.list.b[i], hzd0.o.cond.as)
  pen.o.cond.bl[, , i] <- mpfa(w.list.b[i], hzd0.o.cond.bl)
  pen.o.cond.his[, , i] <- mpfa(w.list.b[i], hzd0.o.cond.his)
  pen.o.cond.na[, , i] <- mpfa(w.list.b[i], hzd0.o.cond.na)
  pen.o.cond.wh[, , i] <- mpfa(w.list.b[i], hzd0.o.cond.wh)
  surv.o.cond[, , i] <- 1 - apply(pen.o.cond[, , i], 2, cumsum)
  surv.o.cond.as[, , i] <- 1 - apply(pen.o.cond.as[, , i], 2, cumsum)
  surv.o.cond.bl[, , i] <- 1 - apply(pen.o.cond.bl[, , i], 2, cumsum)
  surv.o.cond.his[, , i] <- 1 - apply(pen.o.cond.his[, , i], 2, cumsum)
  surv.o.cond.na[, , i] <- 1 - apply(pen.o.cond.na[, , i], 2, cumsum)
  surv.o.cond.wh[, , i] <- 1 - apply(pen.o.cond.wh[, , i], 2, cumsum)
}