load("/home/clausb/scripts/objects/tmat.Rdata")


nameslist <- list()
length(nameslist) <- length(unlist(as.list(unique(tmat['Trace']))))
names(nameslist) <- unlist(as.list(unique(tmat['Trace'])))
for (i in 0:8){
  nameslist[[i+1]] <- tmat[tmat['Trace'] == i,]
}
save(nameslist, file = "nameslist.Rdata")
