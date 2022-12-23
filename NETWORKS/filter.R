load("allnetworks.RData", verbose = FALSE)
done <- "models removed"
doner <- "models remaining"
trash <- list()
trashcount <- 0
pb <- txtProgressBar(min = 0, max = length(bensface), style = 3)
for (i in (1:length(bensface))){
  if (colSums(bensface[[i]])[1] == 0){
    trashcount <- trashcount+1
    trash[trashcount] <- i
  } else if (colSums(bensface[[i]])[1] == 1 & bensface[[i]][1,1] == 1){
    trashcount <- trashcount+1
    trash[trashcount] <- i
  } else if (colSums(bensface[[i]])[1] == 2 & bensface[[i]][1,1] == 2){
    trashcount <- trashcount+1
    trash[trashcount] <- i
  } else if (colSums(bensface[[i]])[2] == 0){
    trashcount <- trashcount+1
    trash[trashcount] <- i
  } else if (colSums(bensface[[i]])[2] == 1 & bensface[[i]][2,2] == 1){
    trashcount <- trashcount+1
    trash[trashcount] <- i
    setTxtProgressBar(pb, i) 
  } else if (colSums(bensface[[i]])[2] == 2 & bensface[[i]][2,2] == 2){
    trashcount <- trashcount+1
    trash[trashcount] <- i
  } else if (colSums(bensface[[i]])[3] == 0){
    trashcount <- trashcount+1
    trash[trashcount] <- i
  } else if (colSums(bensface[[i]])[3] == 1 & bensface[[i]][3,3] == 1){
    trashcount <- trashcount+1
    trash[trashcount] <- i
    setTxtProgressBar(pb, i) 
  } else if (colSums(bensface[[i]])[3] == 2 & bensface[[i]][3,3] == 2){
    trashcount <- trashcount+1
    trash[trashcount] <- i
  } else if (colSums(bensface[[i]])[4] == 0){
    trashcount <- trashcount+1
    trash[trashcount] <- i
  } else if (colSums(bensface[[i]])[4] == 1 & bensface[[i]][4,4] == 1){
    trashcount <- trashcount+1
    trash[trashcount] <- i
  } else if (colSums(bensface[[i]])[4] == 2 & bensface[[i]][4,4] == 2){
    trashcount <- trashcount+1
    trash[trashcount] <- i
  } else if ((rowSums(bensface[[i]])[1] == 0) | (rowSums(bensface[[i]])[2] == 0) | (rowSums(bensface[[i]])[3] == 0) | (rowSums(bensface[[i]])[4] == 0)){
    trashcount <- trashcount+1
    trash[trashcount] <- i
  }
  setTxtProgressBar(pb, i)
}

  
beep()
  #length(trash)
  #trashcount
  filter2 <- trash
  remaining <- length(bensface) - trashcount
  paste(trashcount, done, remaining, doner, sep = " ")
  save(trash, file = "badmodels.Rdata" )
  
