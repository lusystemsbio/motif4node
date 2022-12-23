
load("/home/clausb/R_image/allnetworks.Rdata", verbose = FALSE)
load("/home/clausb/R_image/remainingmodels.Rdata")

tvalues <- list()
for (i in 1:length(remainingmodels)) {
	tvalues[[i]] <- sum(diag(bensface[[remainingmodels[[i]]]]))
}
save(tvalues, file = "tvalues.Rdata")
