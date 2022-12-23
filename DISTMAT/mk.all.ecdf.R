rm(list = ls())
for(i in 1:61){
  load(paste0("/work/lulab/Ben/distmat/output/ecdfs/rsetlist.ecdf.", i))
print(paste0("loaded ", i))
}
all.ecdf <- c(rsetlist.ecdf.1,rsetlist.ecdf.2,rsetlist.ecdf.3,rsetlist.ecdf.4,rsetlist.ecdf.5,rsetlist.ecdf.6,rsetlist.ecdf.7,rsetlist.ecdf.8,rsetlist.ecdf.9,rsetlist.ecdf.10,rsetlist.ecdf.11,rsetlist.ecdf.12,rsetlist.ecdf.13,rsetlist.ecdf.14,rsetlist.ecdf.15,rsetlist.ecdf.16,rsetlist.ecdf.17,rsetlist.ecdf.18,rsetlist.ecdf.19,rsetlist.ecdf.20,rsetlist.ecdf.21,rsetlist.ecdf.22,rsetlist.ecdf.23,rsetlist.ecdf.24,rsetlist.ecdf.25,rsetlist.ecdf.26,rsetlist.ecdf.27,rsetlist.ecdf.28,rsetlist.ecdf.29,rsetlist.ecdf.30,rsetlist.ecdf.31,rsetlist.ecdf.32,rsetlist.ecdf.33,rsetlist.ecdf.34,rsetlist.ecdf.35,rsetlist.ecdf.36,rsetlist.ecdf.37,rsetlist.ecdf.38,rsetlist.ecdf.39,rsetlist.ecdf.40,rsetlist.ecdf.41,rsetlist.ecdf.42,rsetlist.ecdf.43,rsetlist.ecdf.44,rsetlist.ecdf.45,rsetlist.ecdf.46,rsetlist.ecdf.47,rsetlist.ecdf.48,rsetlist.ecdf.49,rsetlist.ecdf.50,rsetlist.ecdf.51,rsetlist.ecdf.52,rsetlist.ecdf.53,rsetlist.ecdf.54,rsetlist.ecdf.55,rsetlist.ecdf.56,rsetlist.ecdf.57,rsetlist.ecdf.58,rsetlist.ecdf.59,rsetlist.ecdf.60,rsetlist.ecdf.61)
save(all.ecdf, file = "/work/lulab/Ben/distmat/output/all.ecdf")
length(all.ecdf)
