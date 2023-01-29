## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE, fig.width = 7, fig.height = 6, fig.align = "center")

## -----------------------------------------------------------------------------
library(motif4node)
library(ggplot2)
library(pheatmap)
library(sRACIPE)
set.seed(43)
data("all_circuits")

## -----------------------------------------------------------------------------
test = sim_4node(index = 11940, Gaussian = F, numModels = 200, all.circuits = all.circuits)

## ---- fig.width = 4.5, fig.height = 3-----------------------------------------
plt = plot_RACIPE(test)
plt[[1]] 

## ---- fig.width = 3, fig.height = 3-------------------------------------------
plt[[2]]

## ---- results = "hide"--------------------------------------------------------
circuit.list = sample(1:length(all.circuits), 10, replace = F)
test = lapply(X = circuit.list, FUN = sim_4node, Gaussian = F, numModels = 200, all.circuits = all.circuits)

## -----------------------------------------------------------------------------
scores_test = lapply(test, trig_score)
scoremat = data.frame(score = unlist(scores_test), index = 1:length(scores_test))
scoremat[order(scoremat$score, decreasing = T),]

## -----------------------------------------------------------------------------
load("distances.2node.rdata")
scoremat = data.frame(score = unlist(distances.2node), index = 1:length(distances.2node))

## -----------------------------------------------------------------------------
motif_results = motif_analysis(all.circuits = all.circuits, all.scores = scoremat, ylim = c(-2,2), 
                               color_breaks = seq(-4,4,by=0.2), filename = "Test", decreasing = F, topCircuits = 218)

## ---- fig.width = 9, fig.height = 5-------------------------------------------
motif_results$single # single motif enrichment

## -----------------------------------------------------------------------------
m = generate_motif_list()
#plot_motif(36, m)

## -----------------------------------------------------------------------------
net = gen_network_scalefree(num_nodes = 20, motif_list = m, motif_choice = c(36,1,31))
#plot_adj(net)

## -----------------------------------------------------------------------------
sessionInfo()

