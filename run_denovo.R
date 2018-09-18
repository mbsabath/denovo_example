## Code to execute the Denovo method on a matched dataset
##
##

devtools::install_github("kwonsang/denovo")
library(denovo)

info <- read.csv("lalonde_matched.csv")

out <- denovo(info)
saveRDS(out, "denovo_out.rds")