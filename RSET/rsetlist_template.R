library(sRACIPE)

indexlist_to_rset <- function(index){
  tmp <- ALLUNIQUEMODELS[[index]]
  colnames(tmp) <- c("A", "B", "C", "D")
  rownames(tmp) <- c("A", "B", "C", "D")
  tmp2 <- as.data.frame(matrix(ncol = 5))
  colnames(tmp2) <- c("from", "to", "arrows.to.type", "color", "type")
  k <- 0
  for (i in 1:4){
    for (j in 1:4){
      if (tmp[[i,j]] == 1){
        k <- k+1
        tmp2[k,1] <- colnames(tmp)[j]
        tmp2[k,2] <- rownames(tmp)[i]
        tmp2[k,3] <- "arrow"
        tmp2[k,4] <- "blue"
        tmp2[k,5] <- 1
        
      } else {
        if (tmp[[i,j]] == 2){
          k <- k+1
          tmp2[k,1] <- colnames(tmp)[j]
          tmp2[k,2] <- rownames(tmp)[i]
          tmp2[k,3] <- "circle"
          tmp2[k,4] <- "red"
          tmp2[k,5] <- 2
        }
      }
    }
  }
  edges <- tmp2
  Edges <- data.frame(edges$from, edges$to, edges$type)
  colnames(Edges) <- c("Source", "Target","Type")
  nodes <- data.frame(id = c("A", "B", "C", "D"), label = c("A", "B", "C", "D"))
  rSet <- sRACIPE::sracipeSimulate(circuit = Edges, numModels = 10000, plots = FALSE)
  return(rSet)
}


if (!file.exists("/projects/lu-lab/BenStuff/research/projects/oldhome/scripts/networks/rset/rsetlists/output/rsetlistNUMBER1")){
  print("creating rsetlistNUMBER1")
  load("/projects/lu-lab/BenStuff/research/projects/oldhome/scripts/networks/rset/rsetlists/input/ALLUNIQUEMODELS_final.Rdata")
  load("/projects/lu-lab/BenStuff/research/projects/oldhome/scripts/networks/rset/indexlists/output/indexlistNUMBER1.Rdata")
  rsetlistNUMBER1 <- lapply(indexlistNUMBER1, indexlist_to_rset)
  save(rsetlistNUMBER1, file ="/projects/lu-lab/BenStuff/research/projects/oldhome/scripts/networks/rset/rsetlists/output/rsetlistNUMBER1" )
}

