if (!file.exists("/projects/lu-lab/BenStuff/research/projects/oldhome/scripts/networks/clustscore/slap_index/scoredf_term1_NUMBER1")){
  load("/projects/lu-lab/BenStuff/research/projects/oldhome/scripts/networks/clustscore/output/scorelistNUMBER1")
  num <- NUMBER1
  nums <- (((num - 1)*1000)+1):(num * 1000)
  scorematNUMBER1 <- matrix(nrow = length(scorelistNUMBER1), ncol = 4)
  for (i in (1:length(scorelistNUMBER1))){
    tmp <- scorelistNUMBER1[[i]]
    scorematNUMBER1[i,4] <- nums[[i]]
    if (!is.null(tmp)){
      scorematNUMBER1[i,1] <- tmp[[1]]
      scorematNUMBER1[i,2] <- tmp[[2]]
      scorematNUMBER1[i,3] <- tmp[[3]]
    }
  }
  colnames(scorematNUMBER1) <- c("score1","score2","score3", "index")
  scoredf_term1_NUMBER1 <- as.data.frame(scorematNUMBER1)
  save(scoredf_term1_NUMBER1, file = "/projects/lu-lab/BenStuff/research/projects/oldhome/scripts/networks/clustscore/slap_index/output/scoredf_term1_NUMBER1")
  rm(scoredfNUMBER1, scorematNUMBER1, scorelistNUMBER1)
}
