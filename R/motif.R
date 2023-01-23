#' A combined motif enrichment analysis for single two-node circuit motifs and motif coupling
#' @param all.circuits List of the topologies of all 60212 non-redundant 4-node circuits. Default "all.circuits" from the package data.
#' @param all.scores Data frame. 1st column: scores; 2nd column: circuit index. The data frame is ordered by the scores.
#' @param ylim Vector of numerics (2). Y axis limit for single motif enrichment. Default: NULL.
#' @param color_breaks Vector that defines color scaling for pheatmap. Default: NULL
#' @param filename Character. Prefix of filenames for plotting. Default: NULL. If provided, plots are also saved to files.
#' @param decreasing Logical. Whether circuits are ranked by the scores in a decreasing order (T) or not (F). Default T.
#' @param topCircuits Integer. Number of top circuits for the enrichment analysis. Default 600.
#' @return List of plotting objects for single motif and motif coupling enrichment analyses.
#' @import ggplot2
#' @import pheatmap
#' @import RColorBrewer
#' @import grDevices
#' @importFrom stats reorder
#' @export
motif_analysis <- function(all.circuits = all.circuits, all.scores, ylim = NULL, color_breaks = NULL, filename = NULL, 
                           decreasing = T, topCircuits = 600){
  motif_list = generate_motif_list()
  new_ind = generate_index_conversion()
  overlap_list = generate_overlap_data()
  grouping = circuit_grouping()
  
  ntot_circuit = nrow(all.scores)
  if(length(all.circuits) != ntot_circuit) print("Warning: the size of the circuits and scores are not consistant!")
  if(any(is.na(all.scores))) print("Warning: NA is found in the scores! Circuits with NA scores will be removed from the analysis")
  if(sum(is.na(all.scores$score)) >0){
    all.scores <- all.scores[-which(is.na(all.scores$score)),]
  }
  
  outcome_single = enrichment_single(all.circuits, all.scores, motif_list, new_ind, decreasing = decreasing, topCircuits = topCircuits)
  
  outcome_single$Groups = grouping[as.numeric(rownames(outcome_single))]
  cbPalette <- c("#56B4E9", "#CC79A7", "#D55E00", "#009E73", "#E69F00", "#999999", "#F0E442", "#0072B2")
  
  outcome_single$index = row.names(outcome_single)
  motif_order <- order(outcome_single$lfc, decreasing = T)
  outcome_single$index <- factor(outcome_single$index, levels = outcome_single$index[motif_order])
   
  p = ggplot(outcome_single, aes_string(x= "index",y = "lfc", fill = "Groups"))  +
    geom_bar(stat="identity") +
    xlab("Motif index") + ylab("Enrichment score") + coord_cartesian(ylim=ylim) +
    scale_fill_manual(values=cbPalette)
  
#  p = ggplot(outcome_single, aes(x= reorder(row.names(outcome_single), -lfc), y = lfc, fill = Groups))  +
#    geom_bar(stat="identity") +
#    xlab("Motif index") + ylab("Enrichment score") + coord_cartesian(ylim=ylim) +
#    scale_fill_manual(values=cbPalette)
  
  outcome_coupling_all = enrichment_coupling_all_cases(all.circuits, all.scores, motif_list, new_ind, overlap_list, 
                                                       decreasing = decreasing, topCircuits = topCircuits)

  breaksList = color_breaks
  colors = cbPalette[1:5];names(colors) <- c(1:5)
  ann_colors = list(groups = colors)
  plot_t0 = pheatmap(outcome_coupling_all$t0,
                     color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
                     breaks = breaksList, annotation_col = as.data.frame(grouping), 
                     annotation_colors = ann_colors)
  plot_t1 = pheatmap(outcome_coupling_all$t1,
                     color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
                     breaks = breaksList, annotation_col = as.data.frame(grouping), 
                     annotation_colors = ann_colors)
  plot_t2 = pheatmap(outcome_coupling_all$t2,
                     color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
                     breaks = breaksList, annotation_col = as.data.frame(grouping), 
                     annotation_colors = ann_colors)


  if(!is.null(filename)){
    ggsave(paste0(filename, "_1D_enrcih.pdf"), plot = last_plot(), device = "pdf", width = 20, height = 10)
    save_pheatmap_pdf(plot_t0, paste0(filename, "_2D_enrcih_0.pdf"), width = 7.5, height = 7)
    save_pheatmap_pdf(plot_t1, paste0(filename, "_2D_enrcih_1.pdf"), width = 7.5, height = 7)
    save_pheatmap_pdf(plot_t2, paste0(filename, "_2D_enrcih_2.pdf"), width = 7.5, height = 7)
  }
   
  return(list(single = p, coupling_0 = plot_t0, coupling_1 = plot_t1, coupling_2 = plot_t2))
}

#' Generate permutations for p-value calculations
#' @param all.circuits List of the topologies of all 60212 non-redundant 4-node circuits. Default "all.circuits" from the package data.
#' @param all.scores Data frame. 1st column: scores; 2nd column: circuit index. The data frame is ordered by the scores.
#' @param decreasing  Logical. Whether circuits are ranked by the scores in a decreasing order (T) or not (F). Default T.
#' @param topCircuits Integer. Number of top circuits for the enrichment analysis. Default 600.
#' @param no_perm Integer. Number of permutations.
#' @return List of permuted enrichment scores.
#' @export
single_motif_permute <- function(all.circuits = all.circuits, all.scores, decreasing = T, topCircuits = 600, 
                                 no_perm){
  new_ind = generate_index_conversion()
  motif_list = generate_motif_list()
  outcome_single = enrichment_single(all.circuits, all.scores, motif_list, new_ind, decreasing = decreasing, topCircuits = topCircuits)
  perm.values = lapply(1:no_perm, gen_ran_enrich, all.circuits = all.circuits, all.scores = all.scores, 
                       motif_list, new_ind, topCircuits = topCircuits)
  return(list(unlist(list(outcome_single, unlist(perm.values)))))
}

#' Motif enrichment analysis for single two-node circuit motifs
#' @param all.circuits List of the topologies of all 60212 non-redundant 4-node circuits. Default "all.circuits" from the package data.
#' @param all.scores Data frame. 1st column: scores; 2nd column: circuit index. The data frame is ordered by the scores.
#' @param motif_list List of 2 by 2 integer matrix. 2-node motif info generated by the function: generate_motif_list.
#' @param new_ind Vector of integer. mapping to non-redundant 2-node motifs generated by the function: generate_index_conversion.
#' @param decreasing Logic. Whether circuits are ranked by the scores in a decreasing order (T) or not (F). Default T.
#' @param topCircuits Integer. Number of top circuits for the enrichment analysis. Default 600.
#' @param nhill Integer/Numeric. Hill coefficient for the Hill function as the weighting factor for motif counts, 
#' large n makes the Hill function more binary. Default 20.
#' @return Data frame containing the enrichment scores of each 2-node circuit motif (39 by 1)
#' @export
enrichment_single <- function(all.circuits = all.circuits, all.scores, motif_list, new_ind, 
                              decreasing = T, topCircuits  = 600, nhill = 20) {
  
  if(any(is.na(all.scores))) print("Warning: NA is found in the scores!")
  
  ntot_circuit = nrow(all.scores)
  #if(length(all.circuits) != ntot_circuit) print("Warning: the size of the circuits and scores are not consistent!")
  
  ntot_motif = 39  # a total of 39 non-redundant 2-node motifs
  
  # weighted counts saved in freq_1 for the first 1%
  freq_1 = integer(ntot_motif)
  freq_99 = integer(ntot_motif)
  
  if(decreasing){
    all.order = order(all.scores$score, decreasing = decreasing)
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
    all.order = order(all.scores$score, decreasing = decreasing)
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
  
  return(outcomes)
}

#' Motif enrichment analysis for the coupling of two two-node circuit motifs
#' @param all.circuits List of the topologies of all 60212 non-redundant 4-node circuits. Default "all.circuits" from the package data.
#' @param all.scores Data frame. 1st column: scores; 2nd column: circuit index. The data frame is ordered by the scores.
#' @param motif_list List of 2 by 2 integer matrix. 2-node motif info generated by the function: generate_motif_list.
#' @param new_ind Vector of integer. mapping to non-redundant 2-node motifs generated by the function: generate_index_conversion.
#' @param overlap_list Matrix (6 by 6 integers). Info of overlapping from the motif location in the 4-node circuit
#' generated by the function: generate_overlap_data.
#' @param decreasing Logic. Whether circuits are ranked by the scores in a decreasing order (T) or not (F). Default T.
#' @param topCircuits Integer. Number of top circuits for the enrichment analysis. Default 600.
#' @param if_overlap Whether consider two motifs with overlapping (1), without overlapping (0), or both (2). Default 2.
#' @param nhill Integer/Numeric. Hill coefficient for the Hill function as the weighting factor for motif counts, 
#' large n makes the Hill function more binary. Default 20.
#' @return A data frame containing the enrichment scores of the coupling between 2-node circuit motifs (39 by 39)
#' @export
enrichment_coupling <- function(all.circuits = all.circuits, all.scores, motif_list, new_ind, overlap_list, 
                                decreasing = T, if_overlap = 2, topCircuits  = 600, nhill = 20) {
  
  if(any(is.na(all.scores))) print("Warning: NA is found in the scores!")
 
  ntot_circuit = nrow(all.scores)
  #if(length(all.circuits) != ntot_circuit) print("Warning: the size of the circuits and scores are not consistant!")
  
  ntot_motif = 39  # a total of 39 non-redundant 2-node motifs
  
  # weighted counts saved in freq_1 for the first 1%
  freq_1_mat = matrix(0, ntot_motif, ntot_motif)
  freq_99_mat = matrix(0, ntot_motif, ntot_motif)
  
  if(decreasing){
    all.order = order(all.scores$score, decreasing = decreasing)
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
    all.order = order(all.scores$score, decreasing = decreasing)
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

#' A convenient function to perform different motif coupling analyses altogether
#' @param all.circuits List of the topologies of all 60212 non-redundant 4-node circuits. Default "all.circuits" from the package data.
#' @param all.scores Data frame. 1st column: scores; 2nd column: circuit index. The data frame is ordered by the scores.
#' @param motif_list List of 2 by 2 integer matrix. 2-node motif info generated by the function: generate_motif_list.
#' @param new_ind Vector of integer. mapping to non-redundant 2-node motifs generated by the function: generate_index_conversion.
#' @param overlap_list Matrix (6 by 6 integers). Info of overlapping from the motif location in the 4-node circuit
#' generated by the function: generate_overlap_data.
#' @param decreasing Logic. Whether circuits are ranked by the scores in a decreasing order (T) or not (F). Default T.
#' @param topCircuits Integer. Number of top circuits for the enrichment analysis. Default 600.
#' @param nhill Integer/Numeric. Hill coefficient for the Hill function as the weighting factor for motif counts, 
#' large n makes the Hill function more binary. Default 20.
#' @return List containing the enrichment of coupling between 2-node motifs for all cases 
#'         (without overlapping, overlapping, and both)
#' @export
enrichment_coupling_all_cases <- function(all.circuits, all.scores, motif_list, new_ind, overlap_list, 
                                          decreasing = decreasing, topCircuits  = 600, nhill = 20){
  outcomes_0 = enrichment_coupling(all.circuits, all.scores, motif_list, new_ind, overlap_list, 
                                   decreasing, if_overlap = 0, topCircuits, nhill = nhill)
  outcomes_1 = enrichment_coupling(all.circuits, all.scores, motif_list, new_ind, overlap_list, 
                                   decreasing, if_overlap = 1, topCircuits, nhill = nhill)
  outcomes_2 = enrichment_coupling(all.circuits, all.scores, motif_list, new_ind, overlap_list, 
                                   decreasing, if_overlap = 2, topCircuits, nhill = nhill)
  return(list(t0 = outcomes_0, t1 = outcomes_1, t2 = outcomes_2))
}

#' Grouping two-node circuit motifs by their types
#' @return Vector of factor: types for all two-node motifs.
#' @export
circuit_grouping <- function(){
  grouping = integer(39)
  grouping[c(1:9)] = 1   #single activation, blue
  grouping[c(10:18)] = 2   # single inhibition, purple
  grouping[c(19,20,21,23,25,26)] = 3 # double activation, red
  grouping[c(22,24,27,37,38,39)] = 4  # double inhibition, green
  grouping[c(28:36)] = 5  # activation/inhibition, yellow
  return(as.factor(grouping))
}

# Generate a vector of enrichment scores by a random permutation
gen_ran_enrich <- function(all.circuits = all.circuits, all.scores, motif_list, new_ind, topCircuits = 600){
  scoremat[,2] = sample(all.scores[,2])
  return(enrichment_single(all.circuits, scoremat, motif_list, new_ind, decreasing = F, topCircuits = topCircuits))
}

# Find all 2-node circuit motifs from a 4-node circuit & the location of the motifs in the 4-node circuit
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

# Find all 2-node circuit motifs from a 4-node circuit
find_motif_2node <- function(mat, motif_list, new_ind){
  results_motifs = find_motif_2node_full(mat, motif_list, new_ind)
  return(results_motifs$motifs)
}

# Find the relationship between two 2-node circuit motifs
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

# Generate a list of all 72 2-node circuit motifs (including redundant ones)
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

#' Convert redundant indices to the indices for all 39 non-redundant motifs 
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

# Generate the data to check whether two motifs overlapping or not, based on the location of the motifs in the 4-node circuit
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

# Assign five types of circuit coupling (Double Cross, Cross, Line, Ship and Box) to a four-node circuit
#' @import igraph
classify_topo_coupling <- function(adj){
  tmp = adj
  dimnames(tmp) = list(LETTERS[1:4],LETTERS[1:4])
  #plot_net(adj_to_tpo(tmp))
  tmp[tmp != 0] <- 1
  diag(tmp) = 0
  g =graph_from_adjacency_matrix(tmp)
  u = as.undirected(g)
  e = as_edgelist(u)
  n = nrow(e)
  l = length(unique(e[,1])) + length(unique(e[,2]))
  if(n == 6){
    return("Double Cross")
  } else if(n == 5){
    return("Cross")
  } else if (n == 3){
    return("Line")
  } else if(l == 5){
    return("Ship")
  } else {
    return("Box")
  }
}

# Inhibitory Hill function
hill <- function(x, x0, n) {
  y = (x/x0)**n
  return(1/(1+y))
}

# Obtain the adjacency matrix for a specific circuit motif (number should be from 1 to 39)
get_motif_adj <- function(number, motif_list){
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
  return(motif)
}

# Obtain the topology data for a specific circuit motif (number should be from 1 to 39)
get_motif <- function(number, motif_list){
  
  motif = get_motif_adj(number, motif_list)
  
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

# Make the full symmetric matrix from the upper diagonal matrix.
makeSymm <- function(m) {
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  return(m)
}

#' Generate a random scale-free gene network consisting the two-node motifs of choice
#' @param num_nodes Integer. Number of nodes of the gene network to be generated.
#' @param motif_list List of 2 by 2 integer matrix. 2-node motif info generated by the function: generate_motif_list.
#' @param motif_choice Vector of integer. A vector of indices of the selected two-node circuit motifs
#' @return Matrix. Adjacency matrix of the generated network
#' @import igraph
#' @export
gen_network_scalefree <- function(num_nodes, motif_list, motif_choice){
  n.nodes <- num_nodes/2
  g <- sample_pa(n.nodes)
  g.adj <- as_adjacency_matrix(g, sparse = F)
  inter <- which(g.adj == 1, arr.ind = T)
  dimnames(g.adj) <- list(as.character(1:nrow(g.adj)), as.character(1:nrow(g.adj)))
  ###################expand matrix to double
  new.mat <- matrix(ncol = n.nodes*2, nrow = n.nodes*2, rep(0,n.nodes*4))
  for(i in 1: nrow(inter)){
    source = inter[i,2]
    target = inter[i,1]
    new.source = source + source
    new.target = (target - 1) +target
    new.mat[new.target, new.source] <-1
  }
  ##########################add in inhibition
  new.inter = which(new.mat == 1, arr.ind = T)
  for(i in 1:nrow(new.inter)){
    new.mat[new.inter[i,1], new.inter[i,2]] <- sample(1:2, 1, replace = T, prob = NULL)
  }
  ################# replace interactions with motifs
  for(i in 1:n.nodes){
    x = i
    y1 = (x-1)+x
    y2 = x+x
    if(length(motif_choice) == 1){
      choice <- sample(1:2, 1)
      if(choice == 1){
        new.mat[y1:y2, y1:y2] <- get_motif_adj(motif_choice)
      } else {
        new.mat[y2:y1,y2:y1] <- get_motif_adj(motif_choice)
      }
    } else {
      choice <- sample(1:2, 1)
      if(choice == 1){
        new.mat[y1:y2, y1:y2] <- get_motif_adj(sample(motif_choice,1))
      } else {
        new.mat[y2:y1,y2:y1] <- get_motif_adj(sample(motif_choice,1))
      }
    }
    
    
  }
  dimnames(new.mat) <- list(as.character(1:nrow(new.mat)), as.character(1:nrow(new.mat)))
  return(new.mat)
}

# Convert an adjacency matrix to the format of circuit topology
# A data frame of three columns: Source, Target, Interaction type
adj_to_tpo <- function(adj){
  seq.adj <- adj
  tmp2 <- as.data.frame(matrix(ncol = 3))
  colnames(tmp2) <- c("source", "target", "type")
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
      }
    }
  }
  #   else if(seq.adj[[i,j]] == 0){
  #     k <- k+1
  #     tmp2[k,1] <- colnames(seq.adj)[j]
  #     tmp2[k,2] <- rownames(seq.adj)[i]
  #     tmp2[k,3] <- 0
  #   }
  # }
  #  }
  #  tmp2 <- tmp2[tmp2$type != 0,]
  return(tmp2)
}
