library(haplo.stats)

load("MB.rda")

indiv = nrow(Y)
snps = ncol(Y)
perc = (100 * length(which(is.na(Y))))/length(Y)

nPossibleHaplotypes = 2 ^ snps

Geno = cbind(substr(Y[,1],1,1),substr(Y[,1],2,2))

for (i in 2:snps) {
  Geno = cbind(Geno, substr(Y[,i],1,1) ,substr(Y[,i],2,2))
}
