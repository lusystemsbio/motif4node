


load("/home/clausb/R_image/allnetworks.Rdata", verbose = FALSE)
load("/home/clausb/R_image/remainingmodels.Rdata")

eigenvalues <- list()
for (i in 1:length(remainingmodels)) {
  e <- eigen(bensface[[remainingmodels[[i]]]])
  v <- e$values
  eigenvalues[[i]] <- v
}
save(eigenvalues, file = "eigenvalues.Rdata")
