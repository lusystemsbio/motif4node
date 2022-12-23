bensface <- list()
variabless <- 1:4
danya <- as.matrix(expand.grid(lapply(numeric(length(variabless)), function(x) c(0,1,2))), ncol=length(variabless))
tricep <- as.matrix(expand.grid(lapply(numeric(length(variabless)), function(x) c(1:81))), ncol=length(variabless))

bensface <- list()
for (i in 1:nrow(tricep)){
  temp <- matrix(nrow = 4, ncol = 4)
  temp[1,] <-danya[tricep[i,1],] 
  temp[2,] <-danya[tricep[i,2],]
  temp[3,] <-danya[tricep[i,3],]
  temp[4,] <-danya[tricep[i,4],]
  bensface[[i]] <- temp
}
save(bensface, file="allnetworks.RData")
