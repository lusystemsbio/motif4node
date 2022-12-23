rm(list = ls())
print("REPLACEME")
Sys.time()
load("/work/lulab/Ben/BenStuff/research/projects/oldhome/scripts/networks/rset/rsetlists/output/rsetlistREPLACEME")
library(parallel)
library(sRACIPE)
library(combinat)
detectCores()
ks.test.self.imp <- function(a,b){
  x <- a
  y <- b
  n.x <- length(x)
  n.y <- length(y)
  n <- n.x * n.y/(n.x + n.y)
  w <- c(x, y)
  z <- cumsum(ifelse(order(w) <= n.x, 1/n.x, -1/n.y))
  return(max(abs(z)))
}
best.dist.final.test <- function(nums){
  x <- log2(t(assay(rsetlistREPLACEME[[nums[1]]])))
  x <- x[rowSums(x) != -Inf,]
  x <- scale(x)
  y <- log2(t(assay(rsetlistREPLACEME[[nums[2]]])))
  y <- y[rowSums(y) != -Inf,]
  y <- scale(y)
  scores <- vector(length = 24)
  a.1 <- x[,"A"]
  b.1 <- x[,"B"]
  c.1 <- x[,"C"]
  d.1 <- x[,"D"]
  AB.1 <- a.1*b.1
  AC.1 <- a.1*c.1
  AD.1 <- a.1*d.1
  BC.1 <- b.1*c.1
  BD.1 <- b.1*d.1
  CD.1 <- c.1*d.1
  i.a <- y[,"A"]
  i.b <- y[,"B"]
  i.c <- y[,"C"]
  i.d <- y[,"D"]
  i.AB <- i.a *i.b
  i.AC <- i.a*i.c
  i.AD <- i.a*i.d
  i.BC <- i.b*i.c
  i.BD <- i.b*i.d
  i.CD <- i.c*i.d
  #1 ABCD
  ks.AB.1 <- ks.test.self.imp(AB.1, i.AB)
  ks.AC.1 <- ks.test.self.imp(AC.1, i.AC)
  ks.AD.1 <- ks.test.self.imp(AD.1, i.AD)
  ks.BC.1 <- ks.test.self.imp(BC.1, i.BC)
  ks.BD.1 <- ks.test.self.imp(BD.1, i.BD)
  ks.CD.1 <- ks.test.self.imp(CD.1, i.CD)
  ks.A <- ks.test.self.imp(a.1, i.a)
  ks.B <- ks.test.self.imp(b.1, i.b)
  ks.C <- ks.test.self.imp(c.1, i.c)
  ks.D <- ks.test.self.imp(d.1, i.d)
  scores[[1]] <- mean(c(ks.AB.1, ks.AC.1, ks.AD.1, ks.BC.1, ks.BD.1, ks.CD.1, ks.A, ks.B, ks.C, ks.D))
  #2 ABDC
  ks.AC.2 <- ks.test.self.imp(AC.1, i.AD)
  ks.AD.2 <- ks.test.self.imp(AD.1, i.AC)
  ks.BC.2 <- ks.test.self.imp(BC.1, i.BD)
  ks.BD.2 <- ks.test.self.imp(BD.1, i.BC)
  ks.A <- ks.test.self.imp(a.1, i.a)
  ks.B <- ks.test.self.imp(b.1, i.b)
  ks.C <- ks.test.self.imp(c.1, i.d)
  ks.D <- ks.test.self.imp(d.1, i.c)
  scores[[2]] <- mean(c(ks.AB.1, ks.AC.2, ks.AD.2, ks.BC.2, ks.BD.2, ks.CD.1, ks.A, ks.B, ks.C, ks.D))
  #3 ADBC
  ks.AB.3 <- ks.test.self.imp(AB.1, i.AD)
  ks.AC.3 <- ks.test.self.imp(AC.1, i.AB)
  ks.BD.3 <- ks.test.self.imp(BD.1, i.CD)
  ks.CD.3 <- ks.test.self.imp(CD.1, i.BC)
  ks.A <- ks.test.self.imp(a.1, i.a)
  ks.B <- ks.test.self.imp(b.1, i.d)
  ks.C <- ks.test.self.imp(c.1, i.b)
  ks.D <- ks.test.self.imp(d.1, i.c)
  scores[[3]] <- mean(c(ks.AB.3, ks.AC.3, ks.AD.2, ks.BC.2, ks.BD.3, ks.CD.3, ks.A, ks.B, ks.C, ks.D))
  #4 DABC
  ks.AC.4 <- ks.test.self.imp(AC.1, i.BD)
  ks.AD.4 <- ks.test.self.imp(AD.1, i.CD)
  ks.BC.4 <- ks.test.self.imp(BC.1, i.AB)
  ks.BD.4 <- ks.test.self.imp(BD.1, i.AC)
  ks.A <- ks.test.self.imp(a.1, i.d)
  ks.B <- ks.test.self.imp(b.1, i.a)
  ks.C <- ks.test.self.imp(c.1, i.b)
  ks.D <- ks.test.self.imp(d.1, i.c)
  scores[[4]] <- mean(c(ks.AB.3, ks.AC.4, ks.AD.4, ks.BC.4, ks.BD.4, ks.CD.3, ks.A, ks.B, ks.C, ks.D))
  #5 DACB
  ks.AC.5 <- ks.test.self.imp(AC.1, i.CD)
  ks.AD.5 <- ks.test.self.imp(AD.1, i.BD)
  ks.BC.5 <- ks.test.self.imp(BC.1, i.AC)
  ks.BD.5 <- ks.test.self.imp(BD.1, i.AB)
  ks.A <- ks.test.self.imp(a.1, i.d)
  ks.B <- ks.test.self.imp(b.1, i.a)
  ks.C <- ks.test.self.imp(c.1, i.c)
  ks.D <- ks.test.self.imp(d.1, i.b)
  scores[[5]] <- mean(c(ks.AB.3, ks.AC.5, ks.AD.5, ks.BC.5, ks.BD.5, ks.CD.3, ks.A, ks.B, ks.C, ks.D))
  #6 ADCB
  ks.AD.6 <- ks.test.self.imp(AD.1, i.AB)
  ks.BC.6 <- ks.test.self.imp(BC.1, i.CD)
  ks.A <- ks.test.self.imp(a.1, i.a)
  ks.B <- ks.test.self.imp(b.1, i.d)
  ks.C <- ks.test.self.imp(c.1, i.c)
  ks.D <- ks.test.self.imp(d.1, i.b)
  scores[[6]] <- mean(c(ks.AB.3, ks.AC.1, ks.AD.6, ks.BC.6, ks.BD.1, ks.CD.3, ks.A, ks.B, ks.C, ks.D))
  #7 ACDB
  ks.AB.7 <- ks.test.self.imp(AB.1, i.AC)
  ks.CD.7 <- ks.test.self.imp(CD.1, i.BD)
  ks.A <- ks.test.self.imp(a.1, i.a)
  ks.B <- ks.test.self.imp(b.1, i.c)
  ks.C <- ks.test.self.imp(c.1, i.d)
  ks.D <- ks.test.self.imp(d.1, i.b)
  scores[[7]] <- mean(c(ks.AB.7, ks.AC.2, ks.AD.6, ks.BC.6, ks.BD.2, ks.CD.7, ks.A, ks.B, ks.C, ks.D))
  #8 ACBD
  ks.A <- ks.test.self.imp(a.1, i.a)
  ks.B <- ks.test.self.imp(b.1, i.c)
  ks.C <- ks.test.self.imp(c.1, i.b)
  ks.D <- ks.test.self.imp(d.1, i.d)
  scores[[8]] <- mean(c(ks.AB.7, ks.AC.3, ks.AD.1, ks.BC.1, ks.BD.3, ks.CD.7, ks.A, ks.B, ks.C, ks.D))
  #9 CABD
  ks.AC.9 <- ks.test.self.imp(AC.1, i.BC)
  ks.BD.9 <- ks.test.self.imp(BD.1, i.AD)
  ks.A <- ks.test.self.imp(a.1, i.c)
  ks.B <- ks.test.self.imp(b.1, i.a)
  ks.C <- ks.test.self.imp(c.1, i.b)
  ks.D <- ks.test.self.imp(d.1, i.d)
  scores[[9]] <- mean(c(ks.AB.7, ks.AC.9, ks.AD.4, ks.BC.4, ks.BD.9, ks.CD.7, ks.A, ks.B, ks.C, ks.D))
  #10 CADB
  ks.AD.10 <- ks.test.self.imp(AD.1, i.BC)
  ks.BC.10 <- ks.test.self.imp(BC.1, i.AD)
  ks.A <- ks.test.self.imp(a.1, i.c)
  ks.B <- ks.test.self.imp(b.1, i.a)
  ks.C <- ks.test.self.imp(c.1, i.d)
  ks.D <- ks.test.self.imp(d.1, i.b)
  scores[[10]] <- mean(c(ks.AB.7, ks.AC.5, ks.AD.10, ks.BC.10, ks.BD.5, ks.CD.7, ks.A, ks.B, ks.C, ks.D))
  #11 CDAB
  ks.AB.11 <- ks.test.self.imp(AB.1, i.CD)
  ks.CD.11 <- ks.test.self.imp(CD.1, i.AB)
  ks.A <- ks.test.self.imp(a.1, i.c)
  ks.B <- ks.test.self.imp(b.1, i.d)
  ks.C <- ks.test.self.imp(c.1, i.a)
  ks.D <- ks.test.self.imp(d.1, i.b)
  scores[[11]] <- mean(c(ks.AB.11, ks.AC.1, ks.AD.10, ks.BC.10, ks.BD.1, ks.CD.11, ks.A, ks.B, ks.C, ks.D))
  #12 DCAB
  ks.A <- ks.test.self.imp(a.1, i.d)
  ks.B <- ks.test.self.imp(b.1, i.c)
  ks.C <- ks.test.self.imp(c.1, i.a)
  ks.D <- ks.test.self.imp(d.1, i.b)
  scores[[12]] <- mean(c(ks.AB.11, ks.AC.2, ks.AD.5, ks.BC.5, ks.BD.2, ks.CD.11, ks.A, ks.B, ks.C, ks.D))
  #13 DCBA
  ks.A <- ks.test.self.imp(a.1, i.d)
  ks.B <- ks.test.self.imp(b.1, i.c)
  ks.C <- ks.test.self.imp(c.1, i.b)
  ks.D <- ks.test.self.imp(d.1, i.a)
  scores[[13]] <- mean(c(ks.AB.11, ks.AC.4, ks.AD.1, ks.BC.1, ks.BD.4, ks.CD.11, ks.A, ks.B, ks.C, ks.D))
  #14 CDBA
  ks.A <- ks.test.self.imp(a.1, i.c)
  ks.B <- ks.test.self.imp(b.1, i.d)
  ks.C <- ks.test.self.imp(c.1, i.b)
  ks.D <- ks.test.self.imp(d.1, i.a)
  scores[[14]] <- mean(c(ks.AB.11, ks.AC.9, ks.AD.2, ks.BC.2, ks.BD.9, ks.CD.11, ks.A, ks.B, ks.C, ks.D))
  #15 CBDA
  ks.AB.15 <- ks.test.self.imp(AB.1, i.BC)
  ks.CD.15 <- ks.test.self.imp(CD.1, i.AD)
  ks.A <- ks.test.self.imp(a.1, i.c)
  ks.B <- ks.test.self.imp(b.1, i.b)
  ks.C <- ks.test.self.imp(c.1, i.d)
  ks.D <- ks.test.self.imp(d.1, i.a)
  scores[[15]] <- mean(c(ks.AB.15, ks.AC.5, ks.AD.2, ks.BC.2, ks.BD.5, ks.CD.15, ks.A, ks.B, ks.C, ks.D))
  #16 CBAD
  ks.A <- ks.test.self.imp(a.1, i.c)
  ks.B <- ks.test.self.imp(b.1, i.b)
  ks.C <- ks.test.self.imp(c.1, i.a)
  ks.D <- ks.test.self.imp(d.1, i.d)
  scores[[16]] <- mean(c(ks.AB.15, ks.AC.1, ks.AD.4, ks.BC.4, ks.BD.1, ks.CD.15, ks.A, ks.B, ks.C, ks.D))
  #17 BCAD
  ks.A <- ks.test.self.imp(a.1, i.b)
  ks.B <- ks.test.self.imp(b.1, i.c)
  ks.C <- ks.test.self.imp(c.1, i.a)
  ks.D <- ks.test.self.imp(d.1, i.d)
  scores[[17]] <- mean(c(ks.AB.15, ks.AC.3, ks.AD.5, ks.BC.5, ks.BD.3, ks.CD.15, ks.A, ks.B, ks.C, ks.D))
  #18 BCDA
  ks.A <- ks.test.self.imp(a.1, i.b)
  ks.B <- ks.test.self.imp(b.1, i.c)
  ks.C <- ks.test.self.imp(c.1, i.d)
  ks.D <- ks.test.self.imp(d.1, i.a)
  scores[[18]] <- mean(c(ks.AB.15, ks.AC.4, ks.AD.6, ks.BC.6, ks.BD.4, ks.CD.15, ks.A, ks.B, ks.C, ks.D))
  #19 BDCA
  ks.AB.19 <- ks.test.self.imp(AB.1, i.BD)
  ks.CD.19 <- ks.test.self.imp(CD.1, i.AC)
  ks.A <- ks.test.self.imp(a.1, i.b)
  ks.B <- ks.test.self.imp(b.1, i.d)
  ks.C <- ks.test.self.imp(c.1, i.c)
  ks.D <- ks.test.self.imp(d.1, i.a)
  scores[[19]] <- mean(c(ks.AB.19, ks.AC.9, ks.AD.6, ks.BC.6, ks.BD.9, ks.CD.19, ks.A, ks.B, ks.C, ks.D))
  #20 DBCA
  ks.A <- ks.test.self.imp(a.1, i.d)
  ks.B <- ks.test.self.imp(b.1, i.b)
  ks.C <- ks.test.self.imp(c.1, i.c )
  ks.D <- ks.test.self.imp(d.1, i.a)
  scores[[20]] <- mean(c(ks.AB.19, ks.AC.5, ks.AD.1, ks.BC.1, ks.BD.5, ks.CD.19, ks.A, ks.B, ks.C, ks.D))
  #21 DBAC
  ks.A <- ks.test.self.imp(a.1, i.d)
  ks.B <- ks.test.self.imp(b.1, i.b)
  ks.C <- ks.test.self.imp(c.1, i.a)
  ks.D <- ks.test.self.imp(d.1, i.c)
  scores[[21]] <- mean(c(ks.AB.19, ks.AC.2, ks.AD.4, ks.BC.4, ks.BD.2, ks.CD.19, ks.A, ks.B, ks.C, ks.D))
  #22 BDAC
  ks.A <- ks.test.self.imp(a.1, i.b)
  ks.B <- ks.test.self.imp(b.1, i.d)
  ks.C <- ks.test.self.imp(c.1, i.a)
  ks.D <- ks.test.self.imp(d.1, i.c)
  scores[[22]] <- mean(c(ks.AB.19, ks.AC.3, ks.AD.10, ks.BC.10, ks.BD.3, ks.CD.19, ks.A, ks.B, ks.C, ks.D))
  #23 BADC
  ks.A <- ks.test.self.imp(a.1, i.b)
  ks.B <- ks.test.self.imp(b.1, i.a)
  ks.C <- ks.test.self.imp(c.1, i.d)
  ks.D <- ks.test.self.imp(d.1, i.c)
  scores[[23]] <- mean(c(ks.AB.1, ks.AC.4, ks.AD.10, ks.BC.10, ks.BD.4, ks.CD.1, ks.A, ks.B, ks.C, ks.D))
  #24 BACD
  ks.A <- ks.test.self.imp(a.1, i.b)
  ks.B <- ks.test.self.imp(b.1, i.a)
  ks.C <- ks.test.self.imp(c.1, i.c)
  ks.D <- ks.test.self.imp(d.1, i.d)
  scores[[24]] <- mean(c(ks.AB.1, ks.AC.9, ks.AD.5, ks.BC.5, ks.BD.9, ks.CD.1, ks.A, ks.B, ks.C, ks.D))
  scores.1 <<- scores
  return(scores[which.min(scores)])
}

dist.combos.REPLACEME <- combn2(1:1000)
dist.combo.list.REPLACEME <- list(length = length(dist.combos.REPLACEME))
for (i in 1:nrow(dist.combos.REPLACEME)){
  dist.combo.list.REPLACEME[[i]] <- dist.combos.REPLACEME[i,]
}
dist.REPLACEME <- mclapply(dist.combo.list.REPLACEME, best.dist.final.test, mc.preschedule = TRUE, mc.cores = 28)
save(dist.REPLACEME, file ="/work/lulab/Ben/distmat/NEWMAN_SUBSAMPLING/output/dist.REPLACEME.R")
save(dist.combos.REPLACEME, file ="/work/lulab/Ben/distmat/NEWMAN_SUBSAMPLING/output/dist.combos.REPLACEME.R")
Sys.time()

