#' Analysis script to evaluate the state distribution of a four-node gene circuit
#' @param rSet sRACIPE object. RACIPE simulation data.
#' @param numModels Numeric. Number of models to be simulated. Default: 10000
#' @return List. plot objects.
#' @importFrom webshot webshot
#' @import sRACIPE
#' @import SummarizedExperiment
#' @export
analysis_circuit_4node <- function(rSet, numModels = 10000){
  gex = log2(t(assay(rSet)))
  pca = prcomp(gex, center = T, scale = T)
  eigs = pca$sdev^2
  eigs = eigs / sum(eigs)
  ###
  
  tmp.g = simu_rnorm(rSet, numModels)
  
  gex.g = log2(t(assay(tmp.g)))
  pca.g = prcomp(gex.g, center = T, scale = T)
  eigs.g = pca.g$sdev^2
  eigs.g = eigs.g / sum(eigs.g)
  ###
  x = plot_net(sracipeCircuit(rSet))
  
  ###
  gex.1  = ggplot(data = as.data.frame(gex), aes_string(x = "A", y = "B")) +
    geom_point() +
    theme_classic()
  
  gex.2  = ggplot(data = as.data.frame(gex), aes_string(x = "A", y = "C")) +
    geom_point() +
    theme_classic()
  
  gex.3  = ggplot(data = as.data.frame(gex), aes_string(x = "A", y = "D")) +
    geom_point() +
    theme_classic()
  
  gex.4  = ggplot(data = as.data.frame(gex), aes_string(x = "B", y = "C")) +
    geom_point() +
    theme_classic()
  
  gex.5  = ggplot(data = as.data.frame(gex), aes_string(x = "B", y = "D")) +
    geom_point() +
    theme_classic()
  
  gex.6  = ggplot(data = as.data.frame(gex), aes_string(x = "C", y = "D")) +
    geom_point() +
    theme_classic()

  rslt =vector(mode = "list", length = 6) 
  rslt[[1]] = plot_grid(gex.1, gex.2, gex.3, gex.4, gex.5, gex.6)
  
  ###
  pca.1  = ggplot(data = as.data.frame(pca$x), aes_string(x = "PC1", y = "PC2")) +
    geom_point() +
    theme_classic() +
    labs( x = paste0("PC1 ",round(eigs[[1]],3)*100, "%" ), y =paste0("PC2 ",round(eigs[[2]],3)*100, "%" ) )
  
  pca.2 = ggplot(data = as.data.frame(pca$x), aes_string(x = "PC1", y = "PC3")) +
    geom_point() +
    theme_classic()+
    labs( x = paste0("PC1 ",round(eigs[[1]],3)*100, "%" ), y =paste0("PC3 ",round(eigs[[3]],3)*100, "%" ) )
  
  pca.3 = ggplot(data = as.data.frame(pca$x), aes_string(x = "PC1", y = "PC4")) +
    geom_point() +
    theme_classic()+
    labs( x = paste0("PC1 ",round(eigs[[1]],3)*100, "%" ), y =paste0("PC4 ",round(eigs[[4]],3)*100, "%" ) )
  
  pca.4 = ggplot(data = as.data.frame(pca$x), aes_string(x = "PC2", y = "PC3")) +
    geom_point() +
    theme_classic()+
    labs( x = paste0("PC2 ",round(eigs[[2]],3)*100, "%" ), y =paste0("PC3 ",round(eigs[[3]],3)*100, "%" ) )
 
  rslt[[2]] = ggdraw() +
    draw_plot(pca.1, 0, .5, .5, .5) +
    draw_plot(pca.2, .5, .5, .5, .5) +
    draw_plot(pca.3, 0, 0, .5, .5) +
    draw_plot(pca.4, .5, 0, .5, .5) +
    draw_plot_label(c("A", "B", "C", "D"), c(0, 0.5, 0, 0.5), c(1,1, 0.5, 0.5), size = 15)
  
  ####
  c.1 = ggplot(data = as.data.frame(pca$x), aes_string(x = "PC1", y  = "PC2")) + 
    stat_density_2d(aes_string (fill = "..level.."), geom = "polygon", bins = 30) +
    geom_density_2d(colour = "white", bins = 30) +
    theme_classic() +
    theme(legend.position="none")
  c.2 = ggplot(data = as.data.frame(pca$x), aes_string(x = "PC1", y  = "PC3")) + 
    stat_density_2d(aes_string (fill = "..level.."), geom = "polygon", bins = 30) +
    geom_density_2d(colour = "white", bins = 30) +
    theme_classic() +
    theme(legend.position="none")
  c.3 = ggplot(data = as.data.frame(pca$x), aes_string(x = "PC1", y  = "PC4")) + 
    stat_density_2d(aes_string (fill = "..level.."), geom = "polygon", bins = 30) +
    geom_density_2d(colour = "white", bins = 30) +
    theme_classic() +
    theme(legend.position="none")
  c.4 = ggplot(data = as.data.frame(pca$x), aes_string(x = "PC2", y  = "PC3")) + 
    stat_density_2d(aes_string (fill = "..level.."), geom = "polygon", bins = 30) +
    geom_density_2d(colour = "white", bins = 30) +
    theme_classic() +
    theme(legend.position="none")
  
  rslt[[3]] = ggdraw() +
    draw_plot(c.1, 0, .5, .5, .5) +
    draw_plot(c.2, .5, .5, .5, .5) +
    draw_plot(c.3, 0, 0, .5, .5) +
    draw_plot(c.4, .5, 0, .5, .5) +
    draw_plot_label(c("E", "F", "G", "H"), c(0, 0.5, 0, 0.5), c(1,1, 0.5, 0.5), size = 15)
  
  ####
  pca.1.g  = ggplot(data = as.data.frame(pca.g$x), aes_string(x = "PC1", y = "PC2")) +
    geom_point() +
    theme_classic() +
    labs( x = paste0("PC1 ",round(eigs.g[[1]],3)*100, "%" ), y =paste0("PC2 ",round(eigs.g[[2]],3)*100, "%" ) )
  
  pca.2.g = ggplot(data = as.data.frame(pca.g$x), aes_string(x = "PC1", y = "PC3")) +
    geom_point() +
    theme_classic()+
    labs( x = paste0("PC1 ",round(eigs.g[[1]],3)*100, "%" ), y =paste0("PC3 ",round(eigs.g[[3]],3)*100, "%" ) )
  
  pca.3.g = ggplot(data = as.data.frame(pca.g$x), aes_string(x = "PC1", y = "PC4")) +
    geom_point() +
    theme_classic()+
    labs( x = paste0("PC1 ",round(eigs.g[[1]],3)*100, "%" ), y =paste0("PC4 ",round(eigs.g[[4]],3)*100, "%" ) )
  
  pca.4.g= ggplot(data = as.data.frame(pca.g$x), aes_string(x = "PC2", y = "PC3")) +
    geom_point() +
    theme_classic()+
    labs( x = paste0("PC2 ",round(eigs.g[[2]],3)*100, "%" ), y =paste0("PC3 ",round(eigs.g[[3]],3)*100, "%" ) )
  
 rslt[[4]] =  ggdraw() +
    draw_plot(pca.1.g, 0, .5, .5, .5) +
    draw_plot(pca.2.g, .5, .5, .5, .5) +
    draw_plot(pca.3.g, 0, 0, .5, .5) +
    draw_plot(pca.4.g, .5, 0, .5, .5) +
    draw_plot_label(c("A", "B", "C", "D"), c(0, 0.5, 0, 0.5), c(1,1, 0.5, 0.5), size = 15)
  
  
  ####
  
  c.1.g = ggplot(data = as.data.frame(pca.g$x), aes_string(x = "PC1", y  = "PC2")) + 
    stat_density_2d(aes_string(fill = "..level.."), geom = "polygon", bins = 30) +
    geom_density_2d(colour = "white", bins = 30) +
    theme_classic() +
    theme(legend.position="none")
  c.2.g = ggplot(data = as.data.frame(pca.g$x), aes_string(x = "PC1", y  = "PC3")) + 
    stat_density_2d(aes_string(fill = "..level.."), geom = "polygon", bins = 30) +
    geom_density_2d(colour = "white", bins = 30) +
    theme_classic() +
    theme(legend.position="none")
  c.3.g = ggplot(data = as.data.frame(pca.g$x), aes_string(x = "PC1", y  = "PC4")) + 
    stat_density_2d(aes_string(fill = "..level.."), geom = "polygon", bins = 30) +
    geom_density_2d(colour = "white", bins = 30) +
    theme_classic() +
    theme(legend.position="none")
  c.4.g = ggplot(data = as.data.frame(pca.g$x), aes_string(x = "PC2", y  = "PC3")) + 
    stat_density_2d(aes_string(fill = "..level.."), geom = "polygon", bins = 30) +
    geom_density_2d(colour = "white", bins = 30) +
    theme_classic() +
    theme(legend.position="none")
  
  rslt[[5]] = ggdraw() +
    draw_plot(c.1.g, 0, .5, .5, .5) +
    draw_plot(c.2.g, .5, .5, .5, .5) +
    draw_plot(c.3.g, 0, 0, .5, .5) +
    draw_plot(c.4.g, .5, 0, .5, .5) +
    draw_plot_label(c("E", "F", "G", "H"), c(0, 0.5, 0, 0.5), c(1,1, 0.5, 0.5), size = 15)
  rslt[[6]] = x
return(rslt)
}

#' Analysis script to evaluate the state distribution of a two-node circuit motif
#' @param tpo Data frame. Topology data of a two-node circuit motif. 
#' The data frame contains three columns: Source, Target, Interaction type
#' @param numModels Numeric. Number of models to be simulated. Default: 10000
#' @import SummarizedExperiment
#' @import sRACIPE
#' @import htmlwidgets
#' @import ggplot2
#' @import cowplot
#' @importFrom webshot webshot
#' @importFrom stats prcomp
#' @export
analysis_circuit_2node <- function(tpo, numModels = 10000){
  rSet = sracipeSimulate(tpo, numModels = numModels)
  gex = log2(t(assay(rSet)))
  pca = prcomp(gex, center = T, scale = T)
  eigs = pca$sdev^2
  eigs = eigs / sum(eigs)
  ###
  tmp.g = simu_rnorm(rSet, numModels)
  
  gex.g = log2(t(assay(tmp.g)))
  pca.g = prcomp(gex.g, center = T, scale = T)
  eigs.g = pca.g$sdev^2
  eigs.g = eigs.g / sum(eigs.g)
  ###
  x = plot_net(sracipeCircuit(rSet))
  ###
  gex.1  = ggplot(data = as.data.frame(gex), aes_string(x = "A", y = "B")) +
    geom_point() +
    theme_classic()
  
  ####
  c.1 = ggplot(data = as.data.frame(gex), aes_string(x = "A", y = "B"))  +
    stat_density_2d(aes_string(fill = "..level.."), geom = "polygon", bins = 30) +
    geom_density_2d(colour = "white", bins = 30) +
    theme_classic() +
    theme(legend.position="none")
  
  ####
  gex.3  = ggplot(data = as.data.frame(gex.g), aes_string(x = "A", y = "B")) +
    geom_point() +
    theme_classic()
  
  ####
  
  c.1.g = ggplot(data = as.data.frame(gex.g), aes_string(x = "A", y = "B"))  +
    stat_density_2d(aes_string(fill = "..level.."), geom = "polygon", bins = 30) +
    geom_density_2d(colour = "white", bins = 30) +
    theme_classic()
  
  return(plot_grid(gex.1, c.1, gex.3, c.1.g))

}

#' Convert an adjacency matrix to a topology
#' @param adj the adjacency martix to be converted to a topology file
#' @export
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
      }else if(seq.adj[[i,j]] == 0){
        k <- k+1
        tmp2[k,1] <- colnames(seq.adj)[j]
        tmp2[k,2] <- rownames(seq.adj)[i]
        tmp2[k,3] <- 0
      }
    }
  }
  tmp2 <- tmp2[tmp2$type != 0,]
  return(tmp2)
}
