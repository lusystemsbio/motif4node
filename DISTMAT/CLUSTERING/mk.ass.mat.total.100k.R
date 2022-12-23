library(sRACIPE)
library(igraph)
library(combinat)


print("ass")
for(i in 1:10){
  if (i == 1){
    load("/work/lulab/Ben/distmat/NEWMAN_SUBSAMPLING/inception/split/Iterative_Subsampling/ASSOCIATION/output/ass.mat")
    ass.mat.total <- ass.mat
    print(i - 1)
  } else {
    load(paste0("/work/lulab/Ben/distmat/NEWMAN_SUBSAMPLING/inception/split/Iterative_Subsampling/ASSOCIATION/output/ass.mat.", i-1))
    ass.mat.total <- ass.mat.total + ass.mat
    print(i-1)
  }
}
save(ass.mat.total, file = "/work/lulab/Ben/distmat/NEWMAN_SUBSAMPLING/inception/split/Iterative_Subsampling/ASSOCIATION/output/ass.mat.total.100k")

