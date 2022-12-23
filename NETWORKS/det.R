library(matlib)
load("/home/clausb/R_image/allnetworks.Rdata", verbose = FALSE)
load("/home/clausb/R_image/remainingmodels.Rdata")
detvalues <- list()
for (i in 1:length(remainingmodels)) {
	detvalues[[i]] <- det(bensface[[remainingmodels[[i]]]])
}
save(detvalues, file = "detvalues.Rdata")
