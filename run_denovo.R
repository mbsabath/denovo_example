## Code to execute the Denovo method on a matched dataset
##
##

devtools::install_github("kwonsang/denovo")
library(denovo)

info <- read.csv("lalonde_matched.csv")

out <- denovo(info, split.ratio=0.25, total.significance=0.05, gamma=0.01))
saveRDS(out, "denovo_out.rds")

sensi <- denovo.sensi(info, Gamma.vec=c(1.01, 1.02, 1.05, 1.08), 
			total.significance=0.05, gamma=0.0001)
