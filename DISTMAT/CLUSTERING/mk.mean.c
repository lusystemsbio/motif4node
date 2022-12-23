#mean.c
library(sRACIPE)
get.info <- function(find.me){
  find.me <- as.numeric(find.me)
  if(length(unlist(strsplit(as.character(find.me), split = ""))) <= 3){
    return(c(0, find.me))
  } else if((length(unlist(strsplit(as.character(find.me), split = ""))) <= 4) & find.me %% 1000 != 0){
    return(c(as.numeric(unlist(strsplit(as.character(find.me), split = ""))[1]), as.numeric(paste(as.numeric(unlist(strsplit(as.character(find.me), split = ""))[2:4]), collapse = ""))))
  } else if((length(unlist(strsplit(as.character(find.me), split = ""))) <= 5)& find.me %% 1000 != 0){
    return(c(as.numeric(paste(as.numeric(unlist(strsplit(as.character(find.me), split = ""))[1:2]), collapse = "")), as.numeric(paste(as.numeric(unlist(strsplit(as.character(find.me), split = ""))[3:5]), collapse = ""))))
  } else if(as.numeric(find.me) %% 1000 == 0){
    return(c(as.numeric(find.me)/1000 - 1, 1000))
  }
}
for(i in 1:61){
  load(paste0("/work/lulab/Ben/BenStuff/research/projects/oldhome/scripts/networks/rset/rsetlists/output/rsetlist", i))
}
ind.to.mean.c <- function(ind){
  info <- get.info(ind)
  tmp <- get(paste0("rsetlist", info[[1]]+1))[[info[[2]]]]
  tmp <- as.data.frame(log2(t(assay(tmp))))
  return(mean(tmp$C))
}
ind <- 1:60212
mean.c <- lapply(ind, ind.to.mean.c)
save(mean.c, file = "mean.c")

