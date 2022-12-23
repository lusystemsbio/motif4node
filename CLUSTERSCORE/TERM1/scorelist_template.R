library(sRACIPE)

if (!file.exists("/projects/lu-lab/BenStuff/research/projects/oldhome/scripts/networks/clustscore/output/scorelistNUMBER1")){
  load("/projects/lu-lab/BenStuff/research/projects/oldhome/scripts/networks/rset/rsetlists/output/rsetlistNUMBER1")
  scorelistNUMBER1 <- list()
  length(scorelistNUMBER1) <- length(rsetlistNUMBER1)
  failedlistNUMBER1 <- list()
  clusterlistNUMBER1 <- list()
  length(clusterlistNUMBER1) <- length(scorelistNUMBER1)
  n <- 1
  for (i in 1:length(rsetlistNUMBER1)){
    tlog_gex <- t(log2(assay(rsetlistNUMBER1[[i]])))
    tlog_gex <- tlog_gex[rowSums(tlog_gex) != -Inf,]
    clusters <- kmeans(tlog_gex, centers = 3, iter.max = 50)
    #individual distances
    if (clusters$ifault != 4){
      dist_AB <- dist(clusters$centers, method = "euclidean")[1]
      dist_BC <- dist(clusters$centers, method = "euclidean")[3]
      dist_AC <- dist(clusters$centers, method = "euclidean")[2]
      #individual r terms
      rA <- sqrt(clusters$withinss[[1]]/clusters$size[[1]])
      rB <- sqrt(clusters$withinss[[2]]/clusters$size[[2]])
      rC <- sqrt(clusters$withinss[[3]]/clusters$size[[3]])
      #indvidual scoring terms
      score_AB <- dist_AB/(rA*rB)
      score_BC <- dist_BC/(rB*rC)
      score_AC <- dist_AC/(rA*rC)
      score_vector <-c(score_AB, score_BC, score_AC)
      score_vector <- score_vector[order(score_vector)]
      scorelistNUMBER1[[i]] <- score_vector
      clusterlistNUMBER1[[i]] <- clusters
    } else if (clusters$ifault == 4){
      failedlistNUMBER1[[n]] <- i
      n <- n +1
    }
  }
  save(scorelistNUMBER1, file = "/projects/lu-lab/BenStuff/research/projects/oldhome/scripts/networks/clustscore/output/scorelistNUMBER1")
  save(failedlistNUMBER1, file = "/projects/lu-lab/BenStuff/research/projects/oldhome/scripts/networks/clustscore/output/failedlistNUMBER1")
  save(clusterlistNUMBER1, file = "/projects/lu-lab/BenStuff/research/projects/oldhome/scripts/networks/clustscore/output/clusterlistNUMBER1")
}
