# check m.l.
#z.score = function(mean, sd, score){
#  return((score - mean)/sd)
#}

#' m.l. never used!
get_from_index <- function(num){
  if(num%%1000 == 0){list <- paste0("rsetlist.ecdf.", (num%/%1000))
  ind <- 1000
  } else {list <- paste0("rsetlist.ecdf.", (num%/%1000) + 1)
  ind <- num%%1000
  }
  return(get(list)[[ind]])
}

#' m.l. an old function, is it the same as get_motif but one less argument? 
get_motif_rename <- function(number){
  A_1 <- matrix(c(0,1,0,0), nrow = 2, ncol = 2)
  A_2 <- matrix(c(1,1,0,0), nrow = 2, ncol = 2)
  A_3 <- matrix(c(1,1,0,1), nrow = 2, ncol = 2)
  A_4 <- matrix(c(0,1,0,1), nrow = 2, ncol = 2)
  A_5 <- matrix(c(2,1,0,0), nrow = 2, ncol = 2)
  A_6 <- matrix(c(0,1,0,2), nrow = 2, ncol = 2)
  A_7 <- matrix(c(2,1,0,2), nrow = 2, ncol = 2)
  A_8 <- matrix(c(1,1,0,2), nrow = 2, ncol = 2)
  A_9 <- matrix(c(2,1,0,1), nrow = 2, ncol = 2)
  
  ##B 10 : 18
  B_1 <- matrix(c(0,2,0,0), nrow = 2, ncol = 2)
  B_2 <- matrix(c(1,2,0,0), nrow = 2, ncol = 2)
  B_3 <- matrix(c(1,2,0,1), nrow = 2, ncol = 2)
  B_4 <- matrix(c(0,2,0,1), nrow = 2, ncol = 2)
  B_5 <- matrix(c(2,2,0,0), nrow = 2, ncol = 2)
  B_6 <- matrix(c(0,2,0,2), nrow = 2, ncol = 2)
  B_7 <- matrix(c(2,2,0,2), nrow = 2, ncol = 2)
  B_8 <- matrix(c(1,2,0,2), nrow = 2, ncol = 2)
  B_9 <- matrix(c(2,2,0,1), nrow = 2, ncol = 2)
  
  #C 19 : 27
  C_1 <- matrix(c(0,1,1,0), nrow = 2, ncol = 2)
  C_2 <- matrix(c(1,1,1,0), nrow = 2, ncol = 2)
  C_3 <- matrix(c(1,1,1,1), nrow = 2, ncol = 2)
  C_4 <- matrix(c(0,1,1,1), nrow = 2, ncol = 2)
  C_5 <- matrix(c(2,1,1,0), nrow = 2, ncol = 2)
  C_6 <- matrix(c(0,1,1,2), nrow = 2, ncol = 2)
  C_7 <- matrix(c(2,1,1,2), nrow = 2, ncol = 2)
  C_8 <- matrix(c(1,1,1,2), nrow = 2, ncol = 2)
  C_9 <- matrix(c(2,1,1,1), nrow = 2, ncol = 2)
  
  #D 28 : 36
  D_1 <- matrix(c(0,1,2,0), nrow = 2, ncol = 2)
  D_2 <- matrix(c(1,1,2,0), nrow = 2, ncol = 2)
  D_3 <- matrix(c(1,1,2,1), nrow = 2, ncol = 2)
  D_4 <- matrix(c(0,1,2,1), nrow = 2, ncol = 2)
  D_5 <- matrix(c(2,1,2,0), nrow = 2, ncol = 2)
  D_6 <- matrix(c(0,1,2,2), nrow = 2, ncol = 2)
  D_7 <- matrix(c(2,1,2,2), nrow = 2, ncol = 2)
  D_8 <- matrix(c(1,1,2,2), nrow = 2, ncol = 2)
  D_9 <- matrix(c(2,1,2,1), nrow = 2, ncol = 2)
  
  #E 37 : 45
  E_1 <- matrix(c(0,2,2,0), nrow = 2, ncol = 2)
  E_2 <- matrix(c(1,2,2,0), nrow = 2, ncol = 2)
  E_3 <- matrix(c(1,2,2,1), nrow = 2, ncol = 2)
  E_4 <- matrix(c(0,2,2,1), nrow = 2, ncol = 2)
  E_5 <- matrix(c(2,2,2,0), nrow = 2, ncol = 2)
  E_6 <- matrix(c(0,2,2,2), nrow = 2, ncol = 2)
  E_7 <- matrix(c(2,2,2,2), nrow = 2, ncol = 2)
  E_8 <- matrix(c(1,2,2,2), nrow = 2, ncol = 2)
  E_9 <- matrix(c(2,2,2,1), nrow = 2, ncol = 2)
  
  #F 46 : 54
  F_1 <- matrix(c(0,2,1,0), nrow = 2, ncol = 2)
  F_2 <- matrix(c(1,2,1,0), nrow = 2, ncol = 2)
  F_3 <- matrix(c(1,2,1,1), nrow = 2, ncol = 2)
  F_4 <- matrix(c(0,2,1,1), nrow = 2, ncol = 2)
  F_5 <- matrix(c(2,2,1,0), nrow = 2, ncol = 2)
  F_6 <- matrix(c(0,2,1,2), nrow = 2, ncol = 2)
  F_7 <- matrix(c(2,2,1,2), nrow = 2, ncol = 2)
  F_8 <- matrix(c(1,2,1,2), nrow = 2, ncol = 2)
  F_9 <- matrix(c(2,2,1,1), nrow = 2, ncol = 2)
  
  # G 55:63
  G_1 <- matrix(c(0,0,1,0), nrow = 2, ncol = 2)
  G_2 <- matrix(c(1,0,1,0), nrow = 2, ncol = 2)
  G_3 <- matrix(c(1,0,1,1), nrow = 2, ncol = 2)
  G_4 <- matrix(c(0,0,1,1), nrow = 2, ncol = 2)
  G_5 <- matrix(c(2,0,1,0), nrow = 2, ncol = 2)
  G_6 <- matrix(c(0,0,1,2), nrow = 2, ncol = 2)
  G_7 <- matrix(c(2,0,1,2), nrow = 2, ncol = 2)
  G_8 <- matrix(c(1,0,1,2), nrow = 2, ncol = 2)
  G_9 <- matrix(c(2,0,1,1), nrow = 2, ncol = 2)
  
  #H 64:72
  H_1 <- matrix(c(0,0,2,0), nrow = 2, ncol = 2)
  H_2 <- matrix(c(1,0,2,0), nrow = 2, ncol = 2)
  H_3 <- matrix(c(1,0,2,1), nrow = 2, ncol = 2)
  H_4 <- matrix(c(0,0,2,1), nrow = 2, ncol = 2)
  H_5 <- matrix(c(2,0,2,0), nrow = 2, ncol = 2)
  H_6 <- matrix(c(0,0,2,2), nrow = 2, ncol = 2)
  H_7 <- matrix(c(2,0,2,2), nrow = 2, ncol = 2)
  H_8 <- matrix(c(1,0,2,2), nrow = 2, ncol = 2)
  H_9 <- matrix(c(2,0,2,1), nrow = 2, ncol = 2)
  
  k <- 0
  motif_list <- list()
  for(i in 1:8){
    for (j in 1:9){
      k <- k + 1
      motif_list[k] <- paste0(LETTERS[i], "_", j)
    }
  }
  
  #new naming convention
  motif_list[[22]] <- motif_list[[41]]
  motif_list[[24]] <- motif_list[[43]]
  motif_list[[27]] <- motif_list[[44]]
  
  tmp <- get(motif_list[[number]])
  colnames(tmp) <- c("A", "B")
  rownames(tmp) <- c("A", "B")
  tmp2 <- as.data.frame(matrix(ncol = 4))
  colnames(tmp2) <- c("from", "to", "arrows.to.type", "color")
  k <- 0
  for (i in 1:2){
    for (j in 1:2){
      if (tmp[[i,j]] == 1){
        k <- k+1
        tmp2[k,1] <- colnames(tmp)[j]
        tmp2[k,2] <- rownames(tmp)[i]
        tmp2[k,3] <- "arrow"
        tmp2[k,4] <- "blue"
      } else {
        if (tmp[[i,j]] == 2){
          k <- k+1
          tmp2[k,1] <- colnames(tmp)[j]
          tmp2[k,2] <- rownames(tmp)[i]
          tmp2[k,3] <- "circle"
          tmp2[k,4] <- "red"
        }
      }
    }
  }
  #tmp2
  edges <- tmp2
  nodes <- data.frame(id = c("A", "B"), label = c("A", "B"))
  return(c(nodes,edges))
}

#' m.l. what does this function do? any use?
get_info <- function(find.me){
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

#'m.l. a function that has never been used. undefined EMat_pca and pcaRotation
gen_perm <- function(meh){
  meh = meh
  tmp.emat = apply(EMat_pca, 2, sample)
  ecdf.1  <- cal_ecdf_pca(scale(as.data.frame(scale(tmp.emat, rep(0, 1448), FALSE) %*% -pcaRotation)))
  return(ecdf.1)
}
