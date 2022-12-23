load("/home/clausb/scripts/objects/strace.Rdata")
for (i in 1:355983){
  j <- i+1
  if (strace[[i,1]] == strace[[j,1]]) {
    if (strace[[i,2]] == strace[[j,2]]) {
      if (strace[[i,3]] == strace[[j,3]]) {
        if (strace[[i,4]] == strace[[j,4]]) {
          if (strace[[i,5]] == strace[[j,5]]) {
            if (strace[[i,6]] == strace[[j,6]]) {
              strace[[i,1]] <- "A"
            }
          }
        }
      }
    }
  }
}

strace$Eigenvalue_1 == "A"
strace <- strace[tmp$Eigenvalue_1 != "A",]
save(strace, file = "removed.Rdata")
