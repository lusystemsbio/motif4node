############ stats for 2d enrichments ############
#########
rm(list = ls())
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
##### generate a list of all 72 2-node circuit motifs (including redundant ones)
generate_motif_list <- function() {
  motif_list = list()
  motif_list[[1]] <- matrix(c(0,1,0,0), nrow = 2, ncol = 2)
  motif_list[[2]] <- matrix(c(1,1,0,0), nrow = 2, ncol = 2)
  motif_list[[3]] <- matrix(c(1,1,0,1), nrow = 2, ncol = 2)
  motif_list[[4]] <- matrix(c(0,1,0,1), nrow = 2, ncol = 2)
  motif_list[[5]] <- matrix(c(2,1,0,0), nrow = 2, ncol = 2)
  motif_list[[6]] <- matrix(c(0,1,0,2), nrow = 2, ncol = 2)
  motif_list[[7]] <- matrix(c(2,1,0,2), nrow = 2, ncol = 2)
  motif_list[[8]] <- matrix(c(1,1,0,2), nrow = 2, ncol = 2)
  motif_list[[9]] <- matrix(c(2,1,0,1), nrow = 2, ncol = 2)
  
  ##B 10 : 18
  motif_list[[10]] <- matrix(c(0,2,0,0), nrow = 2, ncol = 2)
  motif_list[[11]] <- matrix(c(1,2,0,0), nrow = 2, ncol = 2)
  motif_list[[12]] <- matrix(c(1,2,0,1), nrow = 2, ncol = 2)
  motif_list[[13]] <- matrix(c(0,2,0,1), nrow = 2, ncol = 2)
  motif_list[[14]] <- matrix(c(2,2,0,0), nrow = 2, ncol = 2)
  motif_list[[15]] <- matrix(c(0,2,0,2), nrow = 2, ncol = 2)
  motif_list[[16]] <- matrix(c(2,2,0,2), nrow = 2, ncol = 2)
  motif_list[[17]] <- matrix(c(1,2,0,2), nrow = 2, ncol = 2)
  motif_list[[18]] <- matrix(c(2,2,0,1), nrow = 2, ncol = 2)
  
  #C 19 : 27
  motif_list[[19]] <- matrix(c(0,1,1,0), nrow = 2, ncol = 2)
  motif_list[[20]] <- matrix(c(1,1,1,0), nrow = 2, ncol = 2)
  motif_list[[21]] <- matrix(c(1,1,1,1), nrow = 2, ncol = 2)
  motif_list[[22]] <- matrix(c(0,1,1,1), nrow = 2, ncol = 2)
  motif_list[[23]] <- matrix(c(2,1,1,0), nrow = 2, ncol = 2)
  motif_list[[24]] <- matrix(c(0,1,1,2), nrow = 2, ncol = 2)
  motif_list[[25]] <- matrix(c(2,1,1,2), nrow = 2, ncol = 2)
  motif_list[[26]] <- matrix(c(1,1,1,2), nrow = 2, ncol = 2)
  motif_list[[27]] <- matrix(c(2,1,1,1), nrow = 2, ncol = 2)
  
  #D 28 : 36
  motif_list[[28]] <- matrix(c(0,1,2,0), nrow = 2, ncol = 2)
  motif_list[[29]] <- matrix(c(1,1,2,0), nrow = 2, ncol = 2)
  motif_list[[30]] <- matrix(c(1,1,2,1), nrow = 2, ncol = 2)
  motif_list[[31]] <- matrix(c(0,1,2,1), nrow = 2, ncol = 2)
  motif_list[[32]] <- matrix(c(2,1,2,0), nrow = 2, ncol = 2)
  motif_list[[33]] <- matrix(c(0,1,2,2), nrow = 2, ncol = 2)
  motif_list[[34]] <- matrix(c(2,1,2,2), nrow = 2, ncol = 2)
  motif_list[[35]] <- matrix(c(1,1,2,2), nrow = 2, ncol = 2)
  motif_list[[36]] <- matrix(c(2,1,2,1), nrow = 2, ncol = 2)
  
  #E 37 : 45
  motif_list[[37]] <- matrix(c(0,2,2,0), nrow = 2, ncol = 2)
  motif_list[[38]] <- matrix(c(1,2,2,0), nrow = 2, ncol = 2)
  motif_list[[39]] <- matrix(c(1,2,2,1), nrow = 2, ncol = 2)
  motif_list[[40]] <- matrix(c(0,2,2,1), nrow = 2, ncol = 2)
  motif_list[[41]] <- matrix(c(2,2,2,0), nrow = 2, ncol = 2)
  motif_list[[42]] <- matrix(c(0,2,2,2), nrow = 2, ncol = 2)
  motif_list[[43]] <- matrix(c(2,2,2,2), nrow = 2, ncol = 2)
  motif_list[[44]] <- matrix(c(1,2,2,2), nrow = 2, ncol = 2)
  motif_list[[45]] <- matrix(c(2,2,2,1), nrow = 2, ncol = 2)
  
  #F 46 : 54
  motif_list[[46]] <- matrix(c(0,2,1,0), nrow = 2, ncol = 2)
  motif_list[[47]] <- matrix(c(1,2,1,0), nrow = 2, ncol = 2)
  motif_list[[48]] <- matrix(c(1,2,1,1), nrow = 2, ncol = 2)
  motif_list[[49]] <- matrix(c(0,2,1,1), nrow = 2, ncol = 2)
  motif_list[[50]] <- matrix(c(2,2,1,0), nrow = 2, ncol = 2)
  motif_list[[51]] <- matrix(c(0,2,1,2), nrow = 2, ncol = 2)
  motif_list[[52]] <- matrix(c(2,2,1,2), nrow = 2, ncol = 2)
  motif_list[[53]] <- matrix(c(1,2,1,2), nrow = 2, ncol = 2)
  motif_list[[54]] <- matrix(c(2,2,1,1), nrow = 2, ncol = 2)
  
  # G 55:63
  motif_list[[55]] <- matrix(c(0,0,1,0), nrow = 2, ncol = 2)
  motif_list[[56]] <- matrix(c(1,0,1,0), nrow = 2, ncol = 2)
  motif_list[[57]] <- matrix(c(1,0,1,1), nrow = 2, ncol = 2)
  motif_list[[58]] <- matrix(c(0,0,1,1), nrow = 2, ncol = 2)
  motif_list[[59]] <- matrix(c(2,0,1,0), nrow = 2, ncol = 2)
  motif_list[[60]] <- matrix(c(0,0,1,2), nrow = 2, ncol = 2)
  motif_list[[61]] <- matrix(c(2,0,1,2), nrow = 2, ncol = 2)
  motif_list[[62]] <- matrix(c(1,0,1,2), nrow = 2, ncol = 2)
  motif_list[[63]] <- matrix(c(2,0,1,1), nrow = 2, ncol = 2)
  
  #H 64:72
  motif_list[[64]] <- matrix(c(0,0,2,0), nrow = 2, ncol = 2)
  motif_list[[65]] <- matrix(c(1,0,2,0), nrow = 2, ncol = 2)
  motif_list[[66]] <- matrix(c(1,0,2,1), nrow = 2, ncol = 2)
  motif_list[[67]] <- matrix(c(0,0,2,1), nrow = 2, ncol = 2)
  motif_list[[68]] <- matrix(c(2,0,2,0), nrow = 2, ncol = 2)
  motif_list[[69]] <- matrix(c(0,0,2,2), nrow = 2, ncol = 2)
  motif_list[[70]] <- matrix(c(2,0,2,2), nrow = 2, ncol = 2)
  motif_list[[71]] <- matrix(c(1,0,2,2), nrow = 2, ncol = 2)
  motif_list[[72]] <- matrix(c(2,0,2,1), nrow = 2, ncol = 2)
  
  return(motif_list)
}

#### convert redundant indices to the indices for all 39 non-redundant motifs 
generate_index_conversion  <- function(){
  
  #initial indices
  new_ind = seq(1,72)
  
  new_ind[46] <- 28
  new_ind[49] <- 29
  new_ind[48] <- 30
  new_ind[47] <- 31
  new_ind[51] <- 32
  new_ind[50] <- 33
  new_ind[52] <- 34
  new_ind[54] <- 35
  new_ind[53] <- 36
  
  new_ind[55] <- 1
  new_ind[56] <- 4
  new_ind[57] <- 3
  new_ind[58] <- 2
  new_ind[59] <- 6
  new_ind[60] <- 5
  new_ind[61] <- 7
  new_ind[62] <- 9
  new_ind[63] <- 8
  
  new_ind[64] <- 10
  new_ind[65] <- 13
  new_ind[66] <- 12
  new_ind[67] <- 11
  new_ind[68] <- 15
  new_ind[69] <- 14
  new_ind[70] <- 16
  new_ind[71] <- 18
  new_ind[72] <- 17
  
  new_ind[22] <- 20
  new_ind[24] <- 23
  new_ind[27] <- 26
  new_ind[40] <- 38
  
  # replace the indices for the following motifs to initially unused indices #22, #24, #27
  new_ind[42] <- 22
  new_ind[45] <- 27
  new_ind[41] <- 22
  new_ind[43] <- 24
  new_ind[44] <- 27
  
  return(new_ind)
}

### generate the data to check whether two motifs overlapping or not, based on the location of the motifs in the 4-node circuit
generate_overlap_data <- function() {
  
  ### location definition:
  ### 1: AB, 2: AC, 3: AD, 4: BC, 5: BD, 6: CD
  overlap_all = matrix(1, nrow = 6, ncol = 6)
  overlap_all[1,6] = 0
  overlap_all[2,5] = 0
  overlap_all[3,4] = 0
  overlap_all[4,3] = 0
  overlap_all[5,2] = 0
  overlap_all[6,1] = 0
  
  return(overlap_all)
}

##### find all 2-node circuit motifs from a 4-node circuit & the location of the motifs in the 4-node circuit
find_motif_2node_full <- function(mat, motif_list, new_ind){
  AA <- mat[1,1]
  AB <- mat[2,1]
  AC <- mat[3,1]
  AD <- mat[4,1]
  BA <- mat[1,2]
  BB <- mat[2,2]
  BC <- mat[3,2]
  BD <- mat[4,2]
  CA <- mat[1,3]
  CB <- mat[2,3]
  CC <- mat[3,3]
  CD <- mat[4,3]
  DA <- mat[1,4]
  DB <- mat[2,4]
  DC <- mat[3,4]
  DD <- mat[4,4]
  
  #turn 4x4 matrix into 6 2x2
  mat22 = list()
  mat22[[1]] <- matrix(c(AA, AB, BA, BB), ncol = 2, nrow = 2)
  mat22[[2]] <- matrix(c(AA, AC, CA, CC), ncol = 2, nrow = 2)
  mat22[[3]] <- matrix(c(AA, AD, DA, DD), ncol = 2, nrow = 2)
  mat22[[4]] <- matrix(c(BB, BC, CB, CC), ncol = 2, nrow = 2)
  mat22[[5]] <- matrix(c(BB, BD, DB, DD), ncol = 2, nrow = 2)
  mat22[[6]] <- matrix(c(CC, CD, DC, DD), ncol = 2, nrow = 2)
  
  #find the indices of all 2-node motifs (redundant version, a total of 72)
  # then convert the indices to the non-redundant version
  index <- integer(6)
  index2 <- integer(6)
  k <- 0 
  for (i in 1:length(mat22)){
    #    ref <- mat22[[i]]
    for(j in 1:length(motif_list)){
      #      test <- motif_list[[j]]
      if (all(mat22[[i]] == motif_list[[j]])){
        k <- k + 1
        index[[k]] <- new_ind[j]   ### assign non-redundant indices
        index2[[k]] <- i   ### which 2-node circuit this belongs to
        break
      }
    }
  }
  
  motifs <- index[1:k]
  location <- index2[1:k]
  return(list(motifs = motifs, location = location))
}

##### find all 2-node circuit motifs from a 4-node circuit (a short version of the previous function)
find_motif_2node <- function(mat, motif_list, new_ind){
  results_motifs = find_motif_2node_full(mat, motif_list, new_ind)
  return(results_motifs$motifs)
}

#### find the relationship between two 2-node circuit motifs 
find_interactions_2node_new <- function(mat, motif_list, new_ind, overlap_list){
  
  results_motifs = find_motif_2node_full(mat, motif_list, new_ind)
  
  num_motifs = length(results_motifs$motifs)
  
  combo_list <- matrix(nrow = num_motifs*(num_motifs-1)/2, ncol = 3)
  colnames(combo_list) <- c("motif1", "motif2", "interaction") 
  
  ind_combo = 0
  for(i in 1:(num_motifs-1)){
    motif_i = results_motifs$motifs[i]
    location_i = results_motifs$location[i]
    for(j in (i+1):num_motifs){
      motif_j = results_motifs$motifs[j]
      location_j = results_motifs$location[j]
      
      ind_combo = ind_combo + 1
      if_interact = overlap_list[location_i, location_j]
      
      if(motif_i > motif_j){
        combo_list[ind_combo,] = c(motif_j, motif_i, if_interact)
      } else {
        combo_list[ind_combo,] = c(motif_i, motif_j, if_interact)
      }
    }
  }
  
  return(combo_list)
}

#### network plotting function
plot.net <- function(tf_links = tf_links){
  require(visNetwork)
  tf_links <- tf_links[which(tf_links[,3] !=0),]
  topology=data.frame(as.matrix(tf_links), stringsAsFactors = F)
  node_list <- unique(c(topology[,1], topology[,2]))
  nodes <- data.frame(id = node_list,  font.size =30, value=c(rep(1,length(node_list))))
  label <- as.character(nodes$id)
  nodes <- cbind(nodes, label)
  edge_col <- data.frame(c(1,2),c("blue","darkred"))
  colnames(edge_col) <- c("relation", "color")
  arrow_type <- data.frame(c(1,2),c("arrow","circle"))
  colnames(arrow_type) <- c("type", "color")
  edges <- data.frame(from =c(topology[,1]), to = c(topology[,2])
                      , arrows.to.type	=arrow_type$color[c(as.numeric(topology[,3]))]
                      , width = 3
                      , color = edge_col$color[c(as.numeric(topology[,3]))]
  )
  visNetwork(nodes, edges, height = "1000px", width = "100%") %>%
    visEdges(arrows = "to") %>%
    visOptions(manipulation = F) %>%
    visLayout(randomSeed = 123) %>%
    visNodes(scaling = list(label = list(enabled = F))) %>%
    visPhysics(solver = "forceAtlas2Based", stabilization = FALSE)%>%
    visNodes(size = 10)
}

#### obtain the topology data for a specific circuit motif (number should be from 1 to 39)
get_motif <- function(number, motif_list){
  if (number > 39){
    print("Warning, the motif number is more than 39")
    return()
  } else if (number == 22){
    motif = motif_list[[41]]
  } else if (number == 24) {
    motif = motif_list[[43]]
  } else if (number == 27) {
    motif = motif_list[[44]]
  } else {
    motif = motif_list[[number]]
  }
  
  colnames(motif) <- c("A", "B")
  rownames(motif) <- c("A", "B")
  
  edges <- as.data.frame(matrix(ncol = 4))
  colnames(edges) <- c("from", "to", "arrows.to.type", "color")
  
  k <- 0
  for (i in 1:2){
    for (j in 1:2){
      if (motif[[i,j]] == 1){
        k <- k+1
        edges[k,1] <- colnames(motif)[j]
        edges[k,2] <- rownames(motif)[i]
        edges[k,3] <- "arrow"
        edges[k,4] <- "blue"
      } else {
        if (motif[[i,j]] == 2){
          k <- k+1
          edges[k,1] <- colnames(motif)[j]
          edges[k,2] <- rownames(motif)[i]
          edges[k,3] <- "circle"
          edges[k,4] <- "red"
        }
      }
    }
  }
  nodes <- data.frame(id = c("A", "B"), label = c("A", "B"))
  return(c(nodes,edges))
}
#### plot a specific circuit motif
plot.motif <- function(number, motif_list){
  tmp <- get_motif(number, motif_list)
  Source <- tmp$from
  Target <- tmp$to
  Type <- tmp$color
  tmp <- as.data.frame(cbind(Source,Target,Type))
  i.r <- which(tmp$Type == "red")
  i.b <- which(tmp$Type == "blue")
  Type <- vector(length = sum(length(i.r),length(i.b)))
  Type[i.r] <- 2
  Type[i.b] <- 1
  tmp$Type <- Type
  plot.net(tmp)
}

##### plot a network from the adjacency matrix 
plot.adj <- function(adj){
  seq.adj <- adj
  tmp2 <- as.data.frame(matrix(ncol = 3))
  colnames(tmp2) <- c("from", "to", "arrows.to.type")
  k <- 0
  for (i in 1:ncol(seq.adj)){
    for (j in 1:ncol(seq.adj)){
      if (seq.adj[[i,j]] == 1){
        k <- k+1
        tmp2[k,1] <- colnames(seq.adj)[j]
        tmp2[k,2] <- rownames(seq.adj)[i]
        tmp2[k,3] <- 1
      } else if(seq.adj[[i,j]] == 2){
        k <- k+1
        tmp2[k,1] <- colnames(seq.adj)[j]
        tmp2[k,2] <- rownames(seq.adj)[i]
        tmp2[k,3] <- 2
      }else if(seq.adj[[i,j]] == 0){
        k <- k+1
        tmp2[k,1] <- colnames(seq.adj)[j]
        tmp2[k,2] <- rownames(seq.adj)[i]
        tmp2[k,3] <- 0
      }
    }
  }
  
  plot.net(tmp2)
}

#### Hill function
hill <- function(x, x0, n) {
  y = (x/x0)**n
  return(1/(1+y))
}

#### making the full symmetric matrix from the upper diagonal matrix
makeSymm <- function(m) {
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  return(m)
}

#### enrichment, single motif
#### input: all.circuits -- a list of all 60212 non-redundant four-node circuits
####        all.scores -- an array of scores for the above circuits.
####        motif_list -- 2-node motif info generated by the function: generate_motif_list()
####        new_ind -- mapping to non-redundant 2-node motifs generated by the function: generate_index_conversion()
###         decreasing -- whether the scores are ordered from high to low (T) or from low to high (F)
###         topCircuits -- the number of top circuits the enrichment analysis emphasizing on
###         nhill -- Hill coefficient for the Hill function as the weighting factor for motif coutns, large n makes the Hill function more binary
enrichment_single <- function(all.circuits, all.scores, motif_list, new_ind, decreasing = T, topCircuits  = 600, nhill = 20) {
  
  #if(any(is.na(all.scores))) print("Warning: NA is found in the scores!")
  all.scores <- all.scores[-which(is.na(all.scores$score)),]
  ntot_circuit = nrow(all.scores)
  #if(length(all.circuits) != ntot_circuit) print("Warning: the size of the circuits and scores are not consistant!")
  
  ntot_motif = 39  # a total of 39 non-redundant 2-node motifs
  
  # weighted counts saved in freq_1 for the first 1% 
  freq_1 = integer(ntot_motif) 
  freq_99 = integer(ntot_motif) 
  
  if(decreasing){
    all.order = order(all.scores$score, decreasing = T)
    all.index = all.scores$index[all.order]
    all.scores_sorted = all.scores$score[all.order]
    score_threshold = all.scores_sorted[topCircuits]
    
    for (i in 1:ntot_circuit){
      circuit = all.circuits[[all.index[i]]]
      motif = find_motif_2node(circuit, motif_list, new_ind)
      for (j in 1:length(motif)) {
        freq_1[motif[j]] = freq_1[motif[j]] + 1 - hill(all.scores_sorted[i], score_threshold, nhill)
        freq_99[motif[j]] = freq_99[motif[j]] + hill(all.scores_sorted[i], score_threshold, nhill)
      }
    }  
  }else {
    all.order = order(all.scores$score, decreasing = F)
    all.index = all.scores$index[all.order]
    all.scores_sorted = all.scores$score[all.order]
    score_threshold = all.scores_sorted[topCircuits]
    
    for (i in 1:ntot_circuit){
      circuit = all.circuits[[all.index[i]]]
      motif = find_motif_2node(circuit, motif_list, new_ind)
      for (j in 1:length(motif)) {
        freq_1[motif[j]] = freq_1[motif[j]] + hill(all.scores_sorted[i], score_threshold, nhill)
        freq_99[motif[j]] = freq_99[motif[j]] + 1 - hill(all.scores_sorted[i], score_threshold, nhill)
      }
    }  
  }
  
  freq_1 = freq_1 / sum(freq_1)
  freq_99 = freq_99 / sum(freq_99)
  lfc = log2(freq_1 / freq_99)
  
  outcomes = data.frame(lfc)
  #  save(outcomes, file = "enrich.1.Rdata")
  
  return(outcomes)
}


#### enrichment, coupling between two motifs
#### input: all.circuits -- a list of all 60212 non-redundant four-node circuits
####        all.scores -- an array of scores for the above circuits.
####        motif_list -- 2-node motif info generated by the function: generate_motif_list()
####        new_ind -- mapping to non-redundant 2-node motifs generated by the function: generate_index_conversion()
####        overlap_list -- info of overlapping from the motif location in the 4-node circuit generated by the function: generate_overlap_data()
####        decreasing -- whether the scores are ordered from high to low (T) or from low to high (F)
####        if_overlap -- whether consider two motifs with overlapping (1), without overlapping (0), or both (2)
####        topCircuits -- the number of top circuits the enrichment analysis emphasizing on
####        nhill -- Hill coefficient for the Hill function as the weighting factor for motif coutns, large n makes the Hill function more binary
enrichment_coupling <- function(all.circuits, all.scores, motif_list, new_ind, overlap_list, decreasing = T, if_overlap = 2, 
                                topCircuits  = 600, nhill = 20) {
  
  #if(any(is.na(all.scores))) print("Warning: NA is found in the scores!")
  all.scores <- all.scores[-which(is.na(all.scores$score)),]
  ntot_circuit = nrow(all.scores)
  #if(length(all.circuits) != ntot_circuit) print("Warning: the size of the circuits and scores are not consistant!")
  
  ntot_motif = 39  # a total of 39 non-redundant 2-node motifs
  
  # weighted counts saved in freq_1 for the first 1% 
  freq_1_mat = matrix(0, ntot_motif, ntot_motif)
  freq_99_mat = matrix(0, ntot_motif, ntot_motif) 
  
  if(decreasing){
    all.order = order(all.scores$score, decreasing = T)
    all.index = all.scores$index[all.order]
    all.scores_sorted = all.scores$score[all.order]
    score_threshold = all.scores_sorted[topCircuits]
    
    for (i in 1:ntot_circuit){
      circuit = all.circuits[[all.index[i]]]
      motifs_coupling = find_interactions_2node_new(circuit, motif_list, new_ind, overlap_list)
      
      for(j in 1: nrow(motifs_coupling)){
        motif1 = motifs_coupling[j,1]
        motif2 = motifs_coupling[j,2]
        if( (if_overlap == 2) ||
            ((if_overlap == 0) && (motifs_coupling[j,3] == 0)) ||
            ((if_overlap == 1) && (motifs_coupling[j,3] == 1))){
          freq_1_mat[motif1, motif2] = freq_1_mat[motif1, motif2] + 1 - hill(all.scores_sorted[i], score_threshold, nhill)
          freq_99_mat[motif1, motif2] = freq_99_mat[motif1, motif2] + hill(all.scores_sorted[i], score_threshold, nhill)
        }
      }
    }  
  }else {
    all.order = order(all.scores$score, decreasing = F)
    all.index = all.scores$index[all.order]
    all.scores_sorted = all.scores$score[all.order]
    score_threshold = all.scores_sorted[topCircuits]
    
    for (i in 1:ntot_circuit){
      circuit = all.circuits[[all.index[i]]]
      motifs_coupling = find_interactions_2node_new(circuit, motif_list, new_ind, overlap_list)
      
      for(j in 1: nrow(motifs_coupling)){
        motif1 = motifs_coupling[j,1]
        motif2 = motifs_coupling[j,2]
        if( (if_overlap == 2) ||
            ((if_overlap == 0) && (motifs_coupling[j,3] == 0)) ||
            ((if_overlap == 1) && (motifs_coupling[j,3] == 1))){
          freq_1_mat[motif1, motif2] = freq_1_mat[motif1, motif2] + hill(all.scores_sorted[i], score_threshold, nhill)
          freq_99_mat[motif1, motif2] = freq_99_mat[motif1, motif2] + 1 - hill(all.scores_sorted[i], score_threshold, nhill)
        }
      }
    }
  }
  
  freq_1_mat = freq_1_mat / sum(freq_1_mat)
  freq_99_mat = freq_99_mat / sum(freq_99_mat)
  lfc = log2(freq_1_mat / freq_99_mat)
  lfc = makeSymm(lfc)
  
  outcomes = as.data.frame(lfc)
  colnames(outcomes) = seq(1,ntot_motif)
  
  return(outcomes)
}

#### a convenient function to perform different motif coupling analyses altogether
enrichment_coupling_all_cases <- function(all.circuits, all.scores, motif_list, new_ind, overlap_list, decreasing = T, 
                                          topCircuits  = 600, nhill = 20) {
  outcomes_0 = enrichment_coupling(all.circuits, all.scores, motif_list, new_ind, overlap_list, decreasing, if_overlap = 0, 
                                   topCircuits, nhill)
  outcomes_1 = enrichment_coupling(all.circuits, all.scores, motif_list, new_ind, overlap_list, decreasing, if_overlap = 1, 
                                   topCircuits, nhill)
  outcomes_2 = enrichment_coupling(all.circuits, all.scores, motif_list, new_ind, overlap_list, decreasing, if_overlap = 2, 
                                   topCircuits, nhill)
  
  #  save(outcomes_0, outcomes_1, outcomes_2, file = "enrich.2.Rdata")
  return(list(t0 = outcomes_0, t1 = outcomes_1, t2 = outcomes_2))
}

#### function to perform the whole enrichment analysis (1 motif, motif coupling)
#### all.circuits: a list of the 4x4 topology matrices for all 4-node circuits
#### all.scores: a vector the scores for the corresponding circuits
#### ylim: a vector of two elements specifying the y range for ploting the bar plot; e.g. c(-2,2)
#### colors_breaks: the color breaks used to plot heatmaps; e.g. seq(-4,4,by=0.2)

my_analysis <- function(all.circuits, all.scores, ylim = NULL, color_breaks = NULL, filename = NULL){
  motif_list = generate_motif_list()
  new_ind = generate_index_conversion()
  overlap_list = generate_overlap_data()
  
  outcome_single = enrichment_single(ALLUNIQUEMODELS, all.scores, motif_list, new_ind, decreasing = T, topCircuit = 600) 
  ggplot(outcome_single, aes(x= reorder(row.names(outcome_single), -lfc), lfc))  +
    geom_bar(stat="identity", fill = "skyblue")+
    theme_classic(base_size = 30) +
    xlab("Motif index") + ylab("Enrichment score") + coord_cartesian(ylim=ylim)
  
  outcome_coupling_all = enrichment_coupling_all_cases(ALLUNIQUEMODELS, all.scores, motif_list, new_ind, overlap_list, 
                                                       decreasing = T, topCircuit = 600)
  
  breaksList = color_breaks
  plot_t0 = pheatmap(outcome_coupling_all$t0,
                     color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
                     breaks = breaksList) 
  plot_t1 = pheatmap(outcome_coupling_all$t1,
                     color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
                     breaks = breaksList) 
  plot_t2 = pheatmap(outcome_coupling_all$t2,
                     color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
                     breaks = breaksList)
  
  if(!is.null(filename)){
    ggsave(paste0(filename, "_1D_enrcih.pdf"), plot = last_plot(), device = "pdf", width = 20, height = 10)
    save_pheatmap_pdf(plot_t0, paste0(filename, "_2D_enrcih_0.pdf"), width = 7, height = 7)
    save_pheatmap_pdf(plot_t1, paste0(filename, "_2D_enrcih_1.pdf"), width = 7, height = 7)
    save_pheatmap_pdf(plot_t2, paste0(filename, "_2D_enrcih_2.pdf"), width = 7, height = 7)
  } else {
    plot_t0
    plot_t1
    plot_t2
  }
  return()
}



#### save pheatmap plot to PDF
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
#### term1 stats - 10k runs ####
load("/work/lulab/Ben/BenStuff/research/projects/oldhome/scripts/networks/rset/rsetlists/input/ALLUNIQUEMODELS_final.Rdata")
load("/work/lulab/Ben/BenStuff/research/projects/oldhome/scripts/networks/clustscore/term1/slap_index/output/term1_single_sorted.Rdata")
motif_list = generate_motif_list()
new_ind = generate_index_conversion()

term1_single_sorted$index = sample(term1_single_sorted$index)
term1_stats_df = as.data.frame(matrix(nrow =39, ncol = 1000))
for(i in 1:1000){
  term1_single_sorted$index = sample(term1_single_sorted$index)
  term1_stats_df[,i] = enrichment_single(all.circuits= ALLUNIQUEMODELS, all.scores=term1_single_sorted, motif_list, new_ind=new_ind, decreasing = T, topCircuits  = 600, nhill = 20)[,1]
}

save(term1_stats_df, file = "term1_stats_REPLACEME")








