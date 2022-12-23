load("/projects/lu-lab/BenStuff/research/projects/oldhome/scripts/networks/clustscore/term2/input/clusterlistNUMBER1")
scorelist_term2_NUMBER1 <- list()
length(scorelist_term2_NUMBER1) <- length(clusterlistNUMBER1)
for (i in 1:length(clusterlistNUMBER1)){
  if (!is.null(clusterlistNUMBER1[[i]])){
    dist_AB <- dist(clusterlistNUMBER1[[i]]$centers, method = "euclidean")[1]
    dist_BC <- dist(clusterlistNUMBER1[[i]]$centers, method = "euclidean")[3]
    dist_AC <- dist(clusterlistNUMBER1[[i]]$centers, method = "euclidean")[2]
    rA <- sqrt(clusterlistNUMBER1[[i]]$withinss[[1]]/clusterlistNUMBER1[[i]]$size[[1]])
    rB <- sqrt(clusterlistNUMBER1[[i]]$withinss[[2]]/clusterlistNUMBER1[[i]]$size[[2]])
    rC <- sqrt(clusterlistNUMBER1[[i]]$withinss[[3]]/clusterlistNUMBER1[[i]]$size[[3]])
    score_one <- abs((dist_AB *(rA+rB)) - ((dist_AC *(rA+rC))+(dist_BC*(rB+rC))))
    score_two <- abs((dist_AC *(rA+rC)) - ((dist_AB *(rA+rB))+(dist_BC*(rB+rC))))
    score_three <- abs((dist_BC *(rB+rC)) - ((dist_AB *(rA+rB))+(dist_AC*(rA+rC))))
    score_vector <-c(score_one, score_two, score_three)
    score_vector <- score_vector[order(score_vector)]
    score_matrix <- t(as.matrix(score_vector))
    colnames(score_matrix) <- c('score_one', 'score_two', 'score_three')
    score_matrix <- as.data.frame(score_matrix)
    scorelist_term2_NUMBER1[[i]] <- score_vector
  }
  save(scorelist_term2_NUMBER1, file = "/projects/lu-lab/BenStuff/research/projects/oldhome/scripts/networks/clustscore/term2/output/scorelistNUMBER1")
}

