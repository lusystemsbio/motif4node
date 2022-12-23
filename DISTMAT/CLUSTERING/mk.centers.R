load("/work/lulab/Ben/distmat/NEWMAN_SUBSAMPLING/inception/input/ks.dist.mat.all.complete")
load("/work/lulab/Ben/distmat/NEWMAN_SUBSAMPLING/inception/split/Iterative_Subsampling/ASSOCIATION/clean/input/clusters.g10.100k")


for(i in 1:length(clusters.g10.100k)){
  tmp.cluster <- clusters.g10.100k[[i]]
  tmp <- ks.dist.mat.all[tmp.cluster,tmp.cluster]
  tmp[is.na(tmp)] <- 0
  get.average.dist <- function(ind){
    return(mean(tmp[,ind]))
  }
  average.dist <- lapply(1:ncol(tmp), get.average.dist)
  a.d <- unlist(average.dist)
  print(paste0("center of ", names(clusters.g10.100k)[[i]], " is ", tmp.cluster[[which(a.d == min(a.d))]]))
}



