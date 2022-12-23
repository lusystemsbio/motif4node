rm(list = ls())
library(igraph)
library(parallel)
load("/work/lulab/Ben/distmat/NEWMAN_SUBSAMPLING/inception/input/combos.list")
load("/work/lulab/Ben/distmat/NEWMAN_SUBSAMPLING/inception/input/ks.dist.mat.all.complete")
get.average.distance <- function(combos){
  tmp1 <- names.list[[combos[1]]] #this is id of batch and community
  tmp2 <- names.list[[combos[2]]]  
  membership.1 <- get(tmp1)#this is vecotr of memebrship
  membership.2 <- get(tmp2)
  batch.1 <- strsplit(tmp1, split = "[.]")[[1]][2]#subset(batch)
  batch.2 <- strsplit(tmp2, split = "[.]")[[1]][2]
  index.1 <- membership.1 +((as.numeric(batch.1) - 1) *1000) 
  index.2 <- membership.2 +((as.numeric(batch.2) - 1) *1000)
  index.3 <- c(index.1, index.2)
  sub <- ks.dist.mat.all[index.3,index.3]#this makes smaller distance matrix of everything in both communities
  colnames(sub) <- index.3
  rownames(sub) <- index.3
  a <- length(index.1)
  b <- length(index.2)
  c <- expand.grid(1:a, 1:b) 
  comparison.vec <- vector(length = nrow(c))
  for (i in 1:nrow(c)){
    comparison.vec[[i]] <- sub[as.character(index.1[[c[i,1]]]),as.character(index.2[[c[i,2]]])]
  }
  return(mean(comparison.vec))
}
files <- list.files("/work/lulab/Ben/distmat/NEWMAN_SUBSAMPLING/inception/input/clusters/")
names.list <- vector(mode = "list", length = 1)
k <- 0
l <- 0
for (i in files){
  load(paste0("/work/lulab/Ben/distmat/NEWMAN_SUBSAMPLING/inception/input/clusters/", i))
  nme <- i
  for(j in 1:length(get(nme))){
    l <- l + 1
    nme1 <- paste0(nme, ".", j)
    assign(nme1, which(membership(get(nme)) == names(sizes(get(nme))[j])))
    names.list[[l]] <- nme1
  }
  rm(list = nme)
  k <- k+1
}
distances <- mclapply(combos.list, get.average.distance, mc.preschedule = TRUE, mc.cores = 64)
save(distances, file = "/work/lulab/Ben/distmat/NEWMAN_SUBSAMPLING/inception/output/community.distances")

