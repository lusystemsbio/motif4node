rm(list = ls())
library(sRACIPE)
library(combinat)
mk.ecdf <- function(rset){
  ecdfs <- vector(mode = "list", length = 10)
  x <- scale(log2(t(assay(rset)))[rowSums(log2(t(assay(rset)))) != -Inf,])
  a <- x[,"A"]
  b <- x[,"B"]
  c <- x[,"C"]
  d <- x[,"D"]
  ecdfs[[1]] <- ecdf(a)
  ecdfs[[2]] <- ecdf(b)
  ecdfs[[3]] <- ecdf(c)
  ecdfs[[4]] <- ecdf(d)
  ecdfs[[5]] <- ecdf(a*b)
  ecdfs[[6]] <- ecdf(a*c)
  ecdfs[[7]] <- ecdf(a*d)
  ecdfs[[8]] <- ecdf(b*c)
  ecdfs[[9]] <- ecdf(b*d)
  ecdfs[[10]] <- ecdf(c*d)
  return(ecdfs)
}

for (i in 1:61){
  load(paste0("/work/lulab/Ben/BenStuff/research/projects/oldhome/scripts/networks/rset/rsetlists/output/","rsetlist", i))
  nme <- paste0("rsetlist", i)
  nme1 <- paste0("rsetlist.ecdf.", i)
  nme2 <- paste0("/work/lulab/Ben/distmat/output/ecdfs/", nme1)
  assign(nme1, lapply(get(nme), mk.ecdf))
  save(list = nme1, file = nme2)
  rm(list = c(nme, nme1))
}
