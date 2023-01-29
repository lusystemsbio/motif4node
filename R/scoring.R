#' Simulate one of the 4-node circuits from the 60212 unique circuits 
#' @param index Numeric index number of 4-node circuit to be simulate. Takes values from 1:60212
#' @param Gaussian Logical. If T, kinetic parameters will be sampled from a gaussian distribution. If F, kinetic parameters will be sampled from a uniform distribution
#' @param numModels Numeric. Number of models to be simulated. Default: 10000
#' @param all.circuits List. The topology of all 60212 circuit motifs. Default "all.circuits" from the package data.
#' @return rset: sRACIPE object. RACIPE simulation results for a circuit
#' @import sRACIPE
#' @export
sim_4node <- function(index, Gaussian = F, numModels = 10000, all.circuits = all.circuits){
  adj = all.circuits[[index]]
  dimnames(adj) = list(LETTERS[1:4],LETTERS[1:4])
  if(Gaussian){
    rset = sracipeSimulate(adj_to_tpo(adj), numModels = numModels, plotToFile = F, plots = F, integrate  = F)
    rset.g = simu_rnorm(rset, numModels)
    return(rset.g)
  } else {
    rset = sracipeSimulate(adj_to_tpo(adj), numModels = numModels, plotToFile = F, plots = F, integrate  = T)
    return(rset)
  }
}

#' The scoring function for ranking circuits with a triangular state distribution 
#' @param rset The sRACIPE object of the simulated circuit.
#' @return return(min(score_vector)): returns the score for the triangular state distribution.
#' @import sRACIPE
#' @import SummarizedExperiment
#' @export
trig_score <- function(rset){
  tlog_gex <- t(log2(assay(rset)))
  tlog_gex <- tlog_gex[rowSums(tlog_gex) != -Inf,]
  clusters <- kmeans(tlog_gex, centers = 3, iter.max = 50)
  #individual distances
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
  return(min(score_vector))
}

#' The scoring function for ranking circuits with a linear state distribution 
#' @param rset The sRACIPE object of the simulated circuit.
#' @return return(min(score_vector)): returns the score for the linear state distribution.
#' @importFrom stats kmeans dist
#' @import sRACIPE
#' @import SummarizedExperiment
#' @export
lin_score <- function(rset){
  tlog_gex <- t(log2(assay(rset)))
  tlog_gex <- tlog_gex[rowSums(tlog_gex) != -Inf,]
  clusters <- kmeans(tlog_gex, centers = 3, iter.max = 50)
  dist_AB <- dist(clusters$centers, method = "euclidean")[1]
  dist_BC <- dist(clusters$centers, method = "euclidean")[3]
  dist_AC <- dist(clusters$centers, method = "euclidean")[2]
  rA <- sqrt(clusters$withinss[[1]]/clusters$size[[1]])
  rB <- sqrt(clusters$withinss[[2]]/clusters$size[[2]])
  rC <- sqrt(clusters$withinss[[3]]/clusters$size[[3]])
  score_one <- abs((dist_AB *(rA+rB)) - ((dist_AC *(rA+rC))+(dist_BC*(rB+rC))))
  score_two <- abs((dist_AC *(rA+rC)) - ((dist_AB *(rA+rB))+(dist_BC*(rB+rC))))
  score_three <- abs((dist_BC *(rB+rC)) - ((dist_AB *(rA+rB))+(dist_AC*(rA+rC))))
  score_vector <-c(score_one, score_two, score_three)
  return(min(score_vector))
}

#' Calculate th KS distance of two gene expression distributions 
#' @param query sRACIPE object or PCA matrix of the query data
#' @param reference sRACIPE object of the reference data
#' @param experimental Logical. T: query is the PCA matrix from an experimental dataset. F: query is an sRACIPE object
#' @return the distance between the gene expression distributions from the query and reference.
#' @export
dist_ks <- function(query, reference, experimental){
  tmp2 = cal_ecdf(reference)
  if(experimental){
    tmp1 = cal_ecdf_pca(query)
    return(dist_to_experimental_2node(tmp1, tmp2))
  } else {
    tmp1 = cal_ecdf(query)
    return(best_dist_final_from_list(tmp1, tmp2))
  }
}

# Compute the minimum distances of the state distributions between two four-node circuits.
# All possible permutations are considered.
best_dist_final_from_list <- function(x,y){
  scores <- vector(length = 24)
  a.1 <- x[[1]]
  b.1 <- x[[2]]
  c.1 <- x[[3]]
  d.1 <- x[[4]]
  AB.1 <- x[[5]]
  AC.1 <- x[[6]]
  AD.1 <- x[[7]]
  BC.1 <- x[[8]]
  BD.1 <- x[[9]]
  CD.1 <- x[[10]]
  i.a <- y[[1]]
  i.b <- y[[2]]
  i.c <- y[[3]]
  i.d <- y[[4]]
  i.AB <- y[[5]]
  i.AC <- y[[6]]
  i.AD <- y[[7]]
  i.BC <- y[[8]]
  i.BD <- y[[9]]
  i.CD <- y[[10]]
  #1 ABCD
  ks.AB.1 <- cal_ks_test(AB.1, i.AB, -5,5,.01)
  ks.AC.1 <- cal_ks_test(AC.1, i.AC, -5,5,.01)
  ks.AD.1 <- cal_ks_test(AD.1, i.AD, -5,5,.01)
  ks.BC.1 <- cal_ks_test(BC.1, i.BC, -5,5,.01)
  ks.BD.1 <- cal_ks_test(BD.1, i.BD, -5,5,.01)
  ks.CD.1 <- cal_ks_test(CD.1, i.CD, -5,5,.01)
  ks.A <- cal_ks_test(a.1, i.a, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.b, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.c, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.d, -5,5,.01)
  scores[[1]] <- mean(c(ks.AB.1, ks.AC.1, ks.AD.1, ks.BC.1, ks.BD.1, ks.CD.1, ks.A, ks.B, ks.C, ks.D))
  #2 ABDC
  ks.AC.2 <- cal_ks_test(AC.1, i.AD, -5,5,.01)
  ks.AD.2 <- cal_ks_test(AD.1, i.AC, -5,5,.01)
  ks.BC.2 <- cal_ks_test(BC.1, i.BD, -5,5,.01)
  ks.BD.2 <- cal_ks_test(BD.1, i.BC, -5,5,.01)
  ks.A <- cal_ks_test(a.1, i.a, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.b, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.d, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.c, -5,5,.01)
  scores[[2]] <- mean(c(ks.AB.1, ks.AC.2, ks.AD.2, ks.BC.2, ks.BD.2, ks.CD.1, ks.A, ks.B, ks.C, ks.D))
  #3 ADBC
  ks.AB.3 <- cal_ks_test(AB.1, i.AD, -5,5,.01)
  ks.AC.3 <- cal_ks_test(AC.1, i.AB, -5,5,.01)
  ks.BD.3 <- cal_ks_test(BD.1, i.CD, -5,5,.01)
  ks.CD.3 <- cal_ks_test(CD.1, i.BC, -5,5,.01)
  ks.A <- cal_ks_test(a.1, i.a, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.d, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.b, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.c, -5,5,.01)
  scores[[3]] <- mean(c(ks.AB.3, ks.AC.3, ks.AD.2, ks.BC.2, ks.BD.3, ks.CD.3, ks.A, ks.B, ks.C, ks.D))
  #4 DABC
  ks.AC.4 <- cal_ks_test(AC.1, i.BD, -5,5,.01)
  ks.AD.4 <- cal_ks_test(AD.1, i.CD, -5,5,.01)
  ks.BC.4 <- cal_ks_test(BC.1, i.AB, -5,5,.01)
  ks.BD.4 <- cal_ks_test(BD.1, i.AC, -5,5,.01)
  ks.A <- cal_ks_test(a.1, i.d, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.a, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.b, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.c, -5,5,.01)
  scores[[4]] <- mean(c(ks.AB.3, ks.AC.4, ks.AD.4, ks.BC.4, ks.BD.4, ks.CD.3, ks.A, ks.B, ks.C, ks.D))
  #5 DACB
  ks.AC.5 <- cal_ks_test(AC.1, i.CD, -5,5,.01)
  ks.AD.5 <- cal_ks_test(AD.1, i.BD, -5,5,.01)
  ks.BC.5 <- cal_ks_test(BC.1, i.AC, -5,5,.01)
  ks.BD.5 <- cal_ks_test(BD.1, i.AB, -5,5,.01)
  ks.A <- cal_ks_test(a.1, i.d, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.a, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.c, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.b, -5,5,.01)
  scores[[5]] <- mean(c(ks.AB.3, ks.AC.5, ks.AD.5, ks.BC.5, ks.BD.5, ks.CD.3, ks.A, ks.B, ks.C, ks.D))
  #6 ADCB
  ks.AD.6 <- cal_ks_test(AD.1, i.AB, -5,5,.01)
  ks.BC.6 <- cal_ks_test(BC.1, i.CD, -5,5,.01)
  ks.A <- cal_ks_test(a.1, i.a, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.d, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.c, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.b, -5,5,.01)
  scores[[6]] <- mean(c(ks.AB.3, ks.AC.1, ks.AD.6, ks.BC.6, ks.BD.1, ks.CD.3, ks.A, ks.B, ks.C, ks.D))
  #7 ACDB
  ks.AB.7 <- cal_ks_test(AB.1, i.AC, -5,5,.01)
  ks.CD.7 <- cal_ks_test(CD.1, i.BD, -5,5,.01)
  ks.A <- cal_ks_test(a.1, i.a, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.c, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.d, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.b, -5,5,.01)
  scores[[7]] <- mean(c(ks.AB.7, ks.AC.2, ks.AD.6, ks.BC.6, ks.BD.2, ks.CD.7, ks.A, ks.B, ks.C, ks.D))
  #8 ACBD
  ks.A <- cal_ks_test(a.1, i.a, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.c, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.b, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.d, -5,5,.01)
  scores[[8]] <- mean(c(ks.AB.7, ks.AC.3, ks.AD.1, ks.BC.1, ks.BD.3, ks.CD.7, ks.A, ks.B, ks.C, ks.D))
  #9 CABD
  ks.AC.9 <- cal_ks_test(AC.1, i.BC, -5,5,.01)
  ks.BD.9 <- cal_ks_test(BD.1, i.AD, -5,5,.01)
  ks.A <- cal_ks_test(a.1, i.c, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.a, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.b, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.d, -5,5,.01)
  scores[[9]] <- mean(c(ks.AB.7, ks.AC.9, ks.AD.4, ks.BC.4, ks.BD.9, ks.CD.7, ks.A, ks.B, ks.C, ks.D))
  #10 CADB
  ks.AD.10 <- cal_ks_test(AD.1, i.BC, -5,5,.01)
  ks.BC.10 <- cal_ks_test(BC.1, i.AD, -5,5,.01)
  ks.A <- cal_ks_test(a.1, i.c, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.a, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.d, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.b, -5,5,.01)
  scores[[10]] <- mean(c(ks.AB.7, ks.AC.5, ks.AD.10, ks.BC.10, ks.BD.5, ks.CD.7, ks.A, ks.B, ks.C, ks.D))
  #11 CDAB
  ks.AB.11 <- cal_ks_test(AB.1, i.CD, -5,5,.01)
  ks.CD.11 <- cal_ks_test(CD.1, i.AB, -5,5,.01)
  ks.A <- cal_ks_test(a.1, i.c, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.d, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.a, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.b, -5,5,.01)
  scores[[11]] <- mean(c(ks.AB.11, ks.AC.1, ks.AD.10, ks.BC.10, ks.BD.1, ks.CD.11, ks.A, ks.B, ks.C, ks.D))
  #12 DCAB
  ks.A <- cal_ks_test(a.1, i.d, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.c, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.a, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.b, -5,5,.01)
  scores[[12]] <- mean(c(ks.AB.11, ks.AC.2, ks.AD.5, ks.BC.5, ks.BD.2, ks.CD.11, ks.A, ks.B, ks.C, ks.D))
  #13 DCBA
  ks.A <- cal_ks_test(a.1, i.d, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.c, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.b, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.a, -5,5,.01)
  scores[[13]] <- mean(c(ks.AB.11, ks.AC.4, ks.AD.1, ks.BC.1, ks.BD.4, ks.CD.11, ks.A, ks.B, ks.C, ks.D))
  #14 CDBA
  ks.A <- cal_ks_test(a.1, i.c, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.d, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.b, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.a, -5,5,.01)
  scores[[14]] <- mean(c(ks.AB.11, ks.AC.9, ks.AD.2, ks.BC.2, ks.BD.9, ks.CD.11, ks.A, ks.B, ks.C, ks.D))
  #15 CBDA
  ks.AB.15 <- cal_ks_test(AB.1, i.BC, -5,5,.01)
  ks.CD.15 <- cal_ks_test(CD.1, i.AD, -5,5,.01)
  ks.A <- cal_ks_test(a.1, i.c, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.b, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.d, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.a, -5,5,.01)
  scores[[15]] <- mean(c(ks.AB.15, ks.AC.5, ks.AD.2, ks.BC.2, ks.BD.5, ks.CD.15, ks.A, ks.B, ks.C, ks.D))
  #16 CBAD
  ks.A <- cal_ks_test(a.1, i.c, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.b, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.a, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.d, -5,5,.01)
  scores[[16]] <- mean(c(ks.AB.15, ks.AC.1, ks.AD.4, ks.BC.4, ks.BD.1, ks.CD.15, ks.A, ks.B, ks.C, ks.D))
  #17 BCAD
  ks.A <- cal_ks_test(a.1, i.b, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.c, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.a, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.d, -5,5,.01)
  scores[[17]] <- mean(c(ks.AB.15, ks.AC.3, ks.AD.5, ks.BC.5, ks.BD.3, ks.CD.15, ks.A, ks.B, ks.C, ks.D))
  #18 BCDA
  ks.A <- cal_ks_test(a.1, i.b, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.c, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.d, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.a, -5,5,.01)
  scores[[18]] <- mean(c(ks.AB.15, ks.AC.4, ks.AD.6, ks.BC.6, ks.BD.4, ks.CD.15, ks.A, ks.B, ks.C, ks.D))
  #19 BDCA
  ks.AB.19 <- cal_ks_test(AB.1, i.BD, -5,5,.01)
  ks.CD.19 <- cal_ks_test(CD.1, i.AC, -5,5,.01)
  ks.A <- cal_ks_test(a.1, i.b, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.d, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.c, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.a, -5,5,.01)
  scores[[19]] <- mean(c(ks.AB.19, ks.AC.9, ks.AD.6, ks.BC.6, ks.BD.9, ks.CD.19, ks.A, ks.B, ks.C, ks.D))
  #20 DBCA
  ks.A <- cal_ks_test(a.1, i.d, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.b, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.c , -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.a, -5,5,.01)
  scores[[20]] <- mean(c(ks.AB.19, ks.AC.5, ks.AD.1, ks.BC.1, ks.BD.5, ks.CD.19, ks.A, ks.B, ks.C, ks.D))
  #21 DBAC
  ks.A <- cal_ks_test(a.1, i.d, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.b, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.a, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.c, -5,5,.01)
  scores[[21]] <- mean(c(ks.AB.19, ks.AC.2, ks.AD.4, ks.BC.4, ks.BD.2, ks.CD.19, ks.A, ks.B, ks.C, ks.D))
  #22 BDAC
  ks.A <- cal_ks_test(a.1, i.b ,-5,5,.01)
  ks.B <- cal_ks_test(b.1, i.d, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.a, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.c, -5,5,.01)
  scores[[22]] <- mean(c(ks.AB.19, ks.AC.3, ks.AD.10, ks.BC.10, ks.BD.3, ks.CD.19, ks.A, ks.B, ks.C, ks.D))
  #23 BADC
  ks.A <- cal_ks_test(a.1, i.b, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.a, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.d, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.c, -5,5,.01)
  scores[[23]] <- mean(c(ks.AB.1, ks.AC.4, ks.AD.10, ks.BC.10, ks.BD.4, ks.CD.1, ks.A, ks.B, ks.C, ks.D))
  #24 BACD
  ks.A <- cal_ks_test(a.1, i.b, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.a, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.c, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.d, -5,5,.01)
  scores[[24]] <- mean(c(ks.AB.1, ks.AC.9, ks.AD.5, ks.BC.5, ks.BD.9, ks.CD.1, ks.A, ks.B, ks.C, ks.D))
  return(scores[which.min(scores)])
}

#' Calculate z-score for a specific distance
#' @param distances list of distances
#' @param score Numeric. score to calculate the z-score
#' @return z-score
#' @export
z_score = function(distances, score){
  m = mean(unlist(distances))
  s = sd(unlist(distances))
  return((score - m)/s)
}

# Compute Empirical Cumpulative Distirbution Function (ECDF) of RACIPE-simulated gene expression data.
# Gene expression data is log transformed and standardized.
# The function computes 10 ecdfs:
# 1 - 4: the distributions of the gene expression for genes A, B, C, D, respectively. 
# 5 - 10: the distributions of the products of gene expression for two genes 
#         A & B, A & C, A & D, B & C, B & D, C & D, respectively. 
#' @import sRACIPE
#' @import SummarizedExperiment
#' @importFrom stats ecdf
cal_ecdf <- function(rset){
  ecdfs <- vector(mode = "list", length = 10)
  x <- scale(log2(t(assay(rset)))[rowSums(log2(t(assay(rset)))) != -Inf,])
  a <- x[,"A"]
  b <- x[,"B"]
  c <- x[,"C"]
  d <- x[,"D"]
  ecdfs[[1]] <- ecdf(a)
  ecdfs[[2]] <- ecdf(b)
  ecdfs[[3]] <- ecdf(c)
  ecdfs[[4]] <- ecdf(d)
  ecdfs[[5]] <- ecdf(a*b)
  ecdfs[[6]] <- ecdf(a*c)
  ecdfs[[7]] <- ecdf(a*d)
  ecdfs[[8]] <- ecdf(b*c)
  ecdfs[[9]] <- ecdf(b*d)
  ecdfs[[10]] <- ecdf(c*d)
  return(ecdfs)
}

# Compute Empirical Cumpulative Distirbution Function (ECDF) of data using the first four principal components.
# The function computes 10 ecdfs:
# 1 - 4: the distributions of the gene expression for PC1, PC2, PC3, PC4, respectively. 
# 5 - 10: the distributions of the products of gene expression for two PCs
#         PC1 & PC2, PC1 & PC3, PC1 & PC4, PC2 & PC3, PC2 & PC4, PC3 & PC4, respectively. 
#' @importFrom stats ecdf
cal_ecdf_pca <- function(PCA){
  ecdfs <- vector(mode = "list", length = 10)
  x <- PCA
  a <- x[,1]
  b <- x[,2]
  c <- x[,3]
  d <- x[,4]
  ecdfs[[1]] <- ecdf(a)
  ecdfs[[2]] <- ecdf(b)
  ecdfs[[3]] <- ecdf(c)
  ecdfs[[4]] <- ecdf(d)
  ecdfs[[5]] <- ecdf(a*b)
  ecdfs[[6]] <- ecdf(a*c)
  ecdfs[[7]] <- ecdf(a*d)
  ecdfs[[8]] <- ecdf(b*c)
  ecdfs[[9]] <- ecdf(b*d)
  ecdfs[[10]] <- ecdf(c*d)
  return(ecdfs)
}

#' Compare the state distributions of two four-node circuits, find the most matched genes, 
#' and project the simulated gene expression data of the 2nd circuit to the PCs of the 1st circuit
#' @param rset1 sRACIPE object. RACIPE simulation data for the first circuit
#' @param rset2 sRACIPE object. RACIPE simulation data for the second circuit
#' @return Numeric matrix. PC coordinates of RACIPE simulated gene expression of 
#' the second circuit projected onto the PCs of the RACIPE simulated gene expression of 
#' the first circuit
#' @import sRACIPE
#' @import SummarizedExperiment
#' @export
map_and_project <- function(rset1, rset2){
  tmp <- rset1
  tlog_gex <- log2(t(assay(tmp)))
  tlog_gex <- tlog_gex[rowSums(tlog_gex) != -Inf,]
  #tlog_gex_pca <- prcomp(tlog_gex, center = TRUE, scale = TRUE)
  #ggbiplot(tlog_gex_pca, choices = c(1,2), var.axes = FALSE, ellipse = TRUE)
  
  tmp.2 <- rset2
  tlog_gex.2 <- log2(t(assay(tmp.2)))
  tlog_gex.2 <- tlog_gex.2[rowSums(tlog_gex.2) != -Inf,]
  #tlog_gex_pca.2 <- prcomp(tlog_gex.2, center = TRUE, scale = TRUE)
  
  #pcs.switched <- scale(tlog_gex.2, tlog_gex_pca$center, tlog_gex_pca$scale) %*% tlog_gex_pca$rotation
  #plot(pcs.switched)
  
  tmp.1 <- cal_ecdf(tmp)
  tmp.2 <- cal_ecdf(tmp.2)
  mapping <- gene_mapping_ecdf(tmp.1, tmp.2)
  
  a.2 <- strsplit(mapping, split = "")[[1]][[1]]
  b.2 <- strsplit(mapping, split = "")[[1]][[2]]
  c.2 <- strsplit(mapping, split = "")[[1]][[3]]
  d.2 <- strsplit(mapping, split = "")[[1]][[4]]
  
  #colnames(tlog_gex)
  a.1 <- grep("A", colnames(tlog_gex))
  b.1 <- grep("B", colnames(tlog_gex))
  c.1 <- grep("C", colnames(tlog_gex))
  d.1 <- grep("D", colnames(tlog_gex))
  
  tlog_gex.corrected <- tlog_gex[,c(a.1,b.1,c.1,d.1)]
  tlog_gex.2.corrected <- tlog_gex.2[,c(a.2,b.2,c.2,d.2)]
  
  #colnames(tlog_gex.corrected)
  #colnames(tlog_gex.2.corrected)
  
  tlog_gex_pca <- prcomp(tlog_gex.corrected, center = TRUE, scale = TRUE)
  pcs.switched <- scale(tlog_gex.2.corrected, tlog_gex_pca$center, tlog_gex_pca$scale) %*% tlog_gex_pca$rotation
  return(pcs.switched)
  
}

# Generate Gaussian random numbers according to the statistics of random numbers of other types.
#' @importFrom stats sd rnorm
#' @importFrom scales rescale
gene_rnorm_col <- function(colno, tmp.params){
  tmp =  tmp.params[,colno]
  tmp.norm = rnorm(n=length(tmp), mean =mean(tmp), sd = sd(tmp) )
  tmp = rescale(tmp.norm, to = range(tmp))
  return(tmp)
  
}

#' RACIPE simulations with random kinetic parameters from Gaussian distributions
#' @param rset sRACIPE object. sRACIPE output from a standard RACIPE simulation 
#' @param numModels Numeric. Number of models to be simulated. Default: 10000.
#' @import sRACIPE
#' @export
simu_rnorm <- function(rset, numModels = 10000){
  tmp.params = sracipeParams(rset)
  g.me = setdiff(1:ncol(tmp.params), grep(x = colnames(tmp.params), pattern = "N"))
  mid = lapply(g.me, gene_rnorm_col, tmp.params = tmp.params)
  tmp.params[,g.me] = mid
  sracipeParams(rset) = tmp.params
  return(sracipeSimulate(rset, numModels = numModels, genParams = F, plots = F, integrate = T))
}

# Compute the KS distance between two ECDFs
cal_ks_test <- function(ECDF1, ECDF2, start, stop, by){
  diff <- abs(ECDF1(seq(from = start, to = stop, by = by)) - ECDF2(seq(from = start, to = stop, by = by)))
  return(diff[[which.max(diff)]])
}

# Map genes in two four-node circuits according to the KS distances of ECDFs
#' @importFrom combinat permn
gene_mapping_ecdf <- function(ECDF1, ECDF2){
  nmes <- lapply(permn(LETTERS[1:4]), paste, collapse = "" )
  x <- ECDF1
  y <- ECDF2
  scores <- vector(length = 24)
  a.1 <- x[[1]]
  b.1 <- x[[2]]
  c.1 <- x[[3]]
  d.1 <- x[[4]]
  AB.1 <- x[[5]]
  AC.1 <- x[[6]]
  AD.1 <- x[[7]]
  BC.1 <- x[[8]]
  BD.1 <- x[[9]]
  CD.1 <- x[[10]]
  i.a <- y[[1]]
  i.b <- y[[2]]
  i.c <- y[[3]]
  i.d <- y[[4]]
  i.AB <- y[[5]]
  i.AC <- y[[6]]
  i.AD <- y[[7]]
  i.BC <- y[[8]]
  i.BD <- y[[9]]
  i.CD <- y[[10]]
  #1 ABCD
  ks.AB.1 <- cal_ks_test(AB.1, i.AB, -5,5,.01)
  ks.AC.1 <- cal_ks_test(AC.1, i.AC, -5,5,.01)
  ks.AD.1 <- cal_ks_test(AD.1, i.AD, -5,5,.01)
  ks.BC.1 <- cal_ks_test(BC.1, i.BC, -5,5,.01)
  ks.BD.1 <- cal_ks_test(BD.1, i.BD, -5,5,.01)
  ks.CD.1 <- cal_ks_test(CD.1, i.CD, -5,5,.01)
  ks.A <- cal_ks_test(a.1, i.a, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.b, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.c, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.d, -5,5,.01)
  scores[[1]] <- mean(c(ks.AB.1, ks.AC.1, ks.AD.1, ks.BC.1, ks.BD.1, ks.CD.1, ks.A, ks.B, ks.C, ks.D))
  #2 ABDC
  ks.AC.2 <- cal_ks_test(AC.1, i.AD, -5,5,.01)
  ks.AD.2 <- cal_ks_test(AD.1, i.AC, -5,5,.01)
  ks.BC.2 <- cal_ks_test(BC.1, i.BD, -5,5,.01)
  ks.BD.2 <- cal_ks_test(BD.1, i.BC, -5,5,.01)
  ks.A <- cal_ks_test(a.1, i.a, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.b, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.d, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.c, -5,5,.01)
  scores[[2]] <- mean(c(ks.AB.1, ks.AC.2, ks.AD.2, ks.BC.2, ks.BD.2, ks.CD.1, ks.A, ks.B, ks.C, ks.D))
  #3 ADBC
  ks.AB.3 <- cal_ks_test(AB.1, i.AD, -5,5,.01)
  ks.AC.3 <- cal_ks_test(AC.1, i.AB, -5,5,.01)
  ks.BD.3 <- cal_ks_test(BD.1, i.CD, -5,5,.01)
  ks.CD.3 <- cal_ks_test(CD.1, i.BC, -5,5,.01)
  ks.A <- cal_ks_test(a.1, i.a, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.d, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.b, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.c, -5,5,.01)
  scores[[3]] <- mean(c(ks.AB.3, ks.AC.3, ks.AD.2, ks.BC.2, ks.BD.3, ks.CD.3, ks.A, ks.B, ks.C, ks.D))
  #4 DABC
  ks.AC.4 <- cal_ks_test(AC.1, i.BD, -5,5,.01)
  ks.AD.4 <- cal_ks_test(AD.1, i.CD, -5,5,.01)
  ks.BC.4 <- cal_ks_test(BC.1, i.AB, -5,5,.01)
  ks.BD.4 <- cal_ks_test(BD.1, i.AC, -5,5,.01)
  ks.A <- cal_ks_test(a.1, i.d, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.a, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.b, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.c, -5,5,.01)
  scores[[4]] <- mean(c(ks.AB.3, ks.AC.4, ks.AD.4, ks.BC.4, ks.BD.4, ks.CD.3, ks.A, ks.B, ks.C, ks.D))
  #5 DACB
  ks.AC.5 <- cal_ks_test(AC.1, i.CD, -5,5,.01)
  ks.AD.5 <- cal_ks_test(AD.1, i.BD, -5,5,.01)
  ks.BC.5 <- cal_ks_test(BC.1, i.AC, -5,5,.01)
  ks.BD.5 <- cal_ks_test(BD.1, i.AB, -5,5,.01)
  ks.A <- cal_ks_test(a.1, i.d, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.a, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.c, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.b, -5,5,.01)
  scores[[5]] <- mean(c(ks.AB.3, ks.AC.5, ks.AD.5, ks.BC.5, ks.BD.5, ks.CD.3, ks.A, ks.B, ks.C, ks.D))
  #6 ADCB
  ks.AD.6 <- cal_ks_test(AD.1, i.AB, -5,5,.01)
  ks.BC.6 <- cal_ks_test(BC.1, i.CD, -5,5,.01)
  ks.A <- cal_ks_test(a.1, i.a, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.d, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.c, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.b, -5,5,.01)
  scores[[6]] <- mean(c(ks.AB.3, ks.AC.1, ks.AD.6, ks.BC.6, ks.BD.1, ks.CD.3, ks.A, ks.B, ks.C, ks.D))
  #7 ACDB
  ks.AB.7 <- cal_ks_test(AB.1, i.AC, -5,5,.01)
  ks.CD.7 <- cal_ks_test(CD.1, i.BD, -5,5,.01)
  ks.A <- cal_ks_test(a.1, i.a, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.c, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.d, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.b, -5,5,.01)
  scores[[7]] <- mean(c(ks.AB.7, ks.AC.2, ks.AD.6, ks.BC.6, ks.BD.2, ks.CD.7, ks.A, ks.B, ks.C, ks.D))
  #8 ACBD
  ks.A <- cal_ks_test(a.1, i.a, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.c, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.b, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.d, -5,5,.01)
  scores[[8]] <- mean(c(ks.AB.7, ks.AC.3, ks.AD.1, ks.BC.1, ks.BD.3, ks.CD.7, ks.A, ks.B, ks.C, ks.D))
  #9 CABD
  ks.AC.9 <- cal_ks_test(AC.1, i.BC, -5,5,.01)
  ks.BD.9 <- cal_ks_test(BD.1, i.AD, -5,5,.01)
  ks.A <- cal_ks_test(a.1, i.c, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.a, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.b, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.d, -5,5,.01)
  scores[[9]] <- mean(c(ks.AB.7, ks.AC.9, ks.AD.4, ks.BC.4, ks.BD.9, ks.CD.7, ks.A, ks.B, ks.C, ks.D))
  #10 CADB
  ks.AD.10 <- cal_ks_test(AD.1, i.BC, -5,5,.01)
  ks.BC.10 <- cal_ks_test(BC.1, i.AD, -5,5,.01)
  ks.A <- cal_ks_test(a.1, i.c, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.a, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.d, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.b, -5,5,.01)
  scores[[10]] <- mean(c(ks.AB.7, ks.AC.5, ks.AD.10, ks.BC.10, ks.BD.5, ks.CD.7, ks.A, ks.B, ks.C, ks.D))
  #11 CDAB
  ks.AB.11 <- cal_ks_test(AB.1, i.CD, -5,5,.01)
  ks.CD.11 <- cal_ks_test(CD.1, i.AB, -5,5,.01)
  ks.A <- cal_ks_test(a.1, i.c, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.d, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.a, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.b, -5,5,.01)
  scores[[11]] <- mean(c(ks.AB.11, ks.AC.1, ks.AD.10, ks.BC.10, ks.BD.1, ks.CD.11, ks.A, ks.B, ks.C, ks.D))
  #12 DCAB
  ks.A <- cal_ks_test(a.1, i.d, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.c, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.a, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.b, -5,5,.01)
  scores[[12]] <- mean(c(ks.AB.11, ks.AC.2, ks.AD.5, ks.BC.5, ks.BD.2, ks.CD.11, ks.A, ks.B, ks.C, ks.D))
  #13 DCBA
  ks.A <- cal_ks_test(a.1, i.d, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.c, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.b, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.a, -5,5,.01)
  scores[[13]] <- mean(c(ks.AB.11, ks.AC.4, ks.AD.1, ks.BC.1, ks.BD.4, ks.CD.11, ks.A, ks.B, ks.C, ks.D))
  #14 CDBA
  ks.A <- cal_ks_test(a.1, i.c, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.d, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.b, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.a, -5,5,.01)
  scores[[14]] <- mean(c(ks.AB.11, ks.AC.9, ks.AD.2, ks.BC.2, ks.BD.9, ks.CD.11, ks.A, ks.B, ks.C, ks.D))
  #15 CBDA
  ks.AB.15 <- cal_ks_test(AB.1, i.BC, -5,5,.01)
  ks.CD.15 <- cal_ks_test(CD.1, i.AD, -5,5,.01)
  ks.A <- cal_ks_test(a.1, i.c, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.b, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.d, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.a, -5,5,.01)
  scores[[15]] <- mean(c(ks.AB.15, ks.AC.5, ks.AD.2, ks.BC.2, ks.BD.5, ks.CD.15, ks.A, ks.B, ks.C, ks.D))
  #16 CBAD
  ks.A <- cal_ks_test(a.1, i.c, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.b, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.a, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.d, -5,5,.01)
  scores[[16]] <- mean(c(ks.AB.15, ks.AC.1, ks.AD.4, ks.BC.4, ks.BD.1, ks.CD.15, ks.A, ks.B, ks.C, ks.D))
  #17 BCAD
  ks.A <- cal_ks_test(a.1, i.b, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.c, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.a, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.d, -5,5,.01)
  scores[[17]] <- mean(c(ks.AB.15, ks.AC.3, ks.AD.5, ks.BC.5, ks.BD.3, ks.CD.15, ks.A, ks.B, ks.C, ks.D))
  #18 BCDA
  ks.A <- cal_ks_test(a.1, i.b, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.c, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.d, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.a, -5,5,.01)
  scores[[18]] <- mean(c(ks.AB.15, ks.AC.4, ks.AD.6, ks.BC.6, ks.BD.4, ks.CD.15, ks.A, ks.B, ks.C, ks.D))
  #19 BDCA
  ks.AB.19 <- cal_ks_test(AB.1, i.BD, -5,5,.01)
  ks.CD.19 <- cal_ks_test(CD.1, i.AC, -5,5,.01)
  ks.A <- cal_ks_test(a.1, i.b, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.d, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.c, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.a, -5,5,.01)
  scores[[19]] <- mean(c(ks.AB.19, ks.AC.9, ks.AD.6, ks.BC.6, ks.BD.9, ks.CD.19, ks.A, ks.B, ks.C, ks.D))
  #20 DBCA
  ks.A <- cal_ks_test(a.1, i.d, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.b, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.c , -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.a, -5,5,.01)
  scores[[20]] <- mean(c(ks.AB.19, ks.AC.5, ks.AD.1, ks.BC.1, ks.BD.5, ks.CD.19, ks.A, ks.B, ks.C, ks.D))
  #21 DBAC
  ks.A <- cal_ks_test(a.1, i.d, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.b, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.a, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.c, -5,5,.01)
  scores[[21]] <- mean(c(ks.AB.19, ks.AC.2, ks.AD.4, ks.BC.4, ks.BD.2, ks.CD.19, ks.A, ks.B, ks.C, ks.D))
  #22 BDAC
  ks.A <- cal_ks_test(a.1, i.b ,-5,5,.01)
  ks.B <- cal_ks_test(b.1, i.d, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.a, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.c, -5,5,.01)
  scores[[22]] <- mean(c(ks.AB.19, ks.AC.3, ks.AD.10, ks.BC.10, ks.BD.3, ks.CD.19, ks.A, ks.B, ks.C, ks.D))
  #23 BADC
  ks.A <- cal_ks_test(a.1, i.b, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.a, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.d, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.c, -5,5,.01)
  scores[[23]] <- mean(c(ks.AB.1, ks.AC.4, ks.AD.10, ks.BC.10, ks.BD.4, ks.CD.1, ks.A, ks.B, ks.C, ks.D))
  #24 BACD
  ks.A <- cal_ks_test(a.1, i.b, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.a, -5,5,.01)
  ks.C <- cal_ks_test(c.1, i.c, -5,5,.01)
  ks.D <- cal_ks_test(d.1, i.d, -5,5,.01)
  scores[[24]] <- mean(c(ks.AB.1, ks.AC.9, ks.AD.5, ks.BC.5, ks.BD.9, ks.CD.1, ks.A, ks.B, ks.C, ks.D))
  return(nmes[[which.min(scores)]])
}

# Compute the minimum distances of the state distributions between a four-node circuit and experimental data.
# ECDF1: experimental; ECDF2: circuit simulations.
dist_to_experimental_2node <- function(ECDF1, ECDF2){
  x <- ECDF1
  y <- ECDF2
  scores <- vector(length = 12)
  a.1 <- x[[1]]
  b.1 <- x[[2]]
  AB.1 <- x[[5]]
  i.a <- y[[1]]
  i.b <- y[[2]]
  i.c <- y[[3]]
  i.d <- y[[4]]
  i.AB <- y[[5]]
  i.AC <- y[[6]]
  i.AD <- y[[7]]
  i.BC <- y[[8]]
  i.BD <- y[[9]]
  i.CD <- y[[10]]
  #1 AB
  ks.AB.1 <- cal_ks_test(AB.1, i.AB, -5,5,.01)
  ks.A <- cal_ks_test(a.1, i.a, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.b, -5,5,.01)
  scores[[1]] <- mean(c(ks.AB.1, ks.A, ks.B))
  #2 AD
  ks.AB.3 <- cal_ks_test(AB.1, i.AD, -5,5,.01)
  ks.A <- cal_ks_test(a.1, i.a, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.d, -5,5,.01)
  scores[[2]] <- mean(c(ks.AB.3,ks.A, ks.B))
  #3 DA
  ks.A <- cal_ks_test(a.1, i.d, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.a, -5,5,.01)
  scores[[3]] <- mean(c(ks.AB.3,ks.A, ks.B))
  #4 AC
  ks.AB.7 <- cal_ks_test(AB.1, i.AC, -5,5,.01)
  ks.A <- cal_ks_test(a.1, i.a, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.c, -5,5,.01)
  scores[[4]] <- mean(c(ks.AB.7, ks.A, ks.B))
  #5 CA
  ks.A <- cal_ks_test(a.1, i.c, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.a, -5,5,.01)
  scores[[5]] <- mean(c(ks.AB.7, ks.A, ks.B))
  #6 CD
  ks.AB.11 <- cal_ks_test(AB.1, i.CD, -5,5,.01)
  ks.A <- cal_ks_test(a.1, i.c, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.d, -5,5,.01)
  scores[[6]] <- mean(c(ks.AB.11, ks.A, ks.B))
  #7 DC
  ks.A <- cal_ks_test(a.1, i.d, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.c, -5,5,.01)
  scores[[7]] <- mean(c(ks.AB.11, ks.A, ks.B))
  #8 CB
  ks.AB.15 <- cal_ks_test(AB.1, i.BC, -5,5,.01)
  ks.A <- cal_ks_test(a.1, i.c, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.b, -5,5,.01)
  scores[[8]] <- mean(c(ks.AB.15, ks.A, ks.B))
  #9 BC
  ks.A <- cal_ks_test(a.1, i.b, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.c, -5,5,.01)
  scores[[9]] <- mean(c(ks.AB.15, ks.A, ks.B))
  #10 BD
  ks.AB.19 <- cal_ks_test(AB.1, i.BD, -5,5,.01)
  ks.A <- cal_ks_test(a.1, i.b, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.d, -5,5,.01)
  scores[[10]] <- mean(c(ks.AB.19, ks.A, ks.B))
  #11 DB
  ks.A <- cal_ks_test(a.1, i.d, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.b, -5,5,.01)
  scores[[11]] <- mean(c(ks.AB.19, ks.A, ks.B))
  #12 BA
  ks.A <- cal_ks_test(a.1, i.b, -5,5,.01)
  ks.B <- cal_ks_test(b.1, i.a, -5,5,.01)
  scores[[12]] <- mean(c(ks.AB.1, ks.A, ks.B))
  return(scores[[which.min(scores)]])
}
