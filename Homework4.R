library(haplo.stats)
library(xtable)

load("MB.rda")

#2
indiv = nrow(Y)
snps = ncol(Y)
perc = (100 * length(which(is.na(Y))))/length(Y)

tab = data.frame(indiv, snps, perc)

#3
nPossibleHaplotypes = 2 ^ snps

#4
Geno = cbind(substr(Y[,1],1,1),substr(Y[,1],2,2))

for (i in 2:snps) {
  Geno = cbind(Geno, substr(Y[,i],1,1) ,substr(Y[,i],2,2))
}

Snpnames <- paste("SNP",1:snps,sep="")
Snpnames
HaploRes <- haplo.em(Geno,locus.label=Snpnames,control=haplo.em.control(min.posterior=1e-3))
HaploRes

nHaplo = nrow(HaploRes$haplotype)
haploList = HaploRes$haplotype
indexMax = which(HaploRes$hap.prob == max(HaploRes$hap.prob))
mostCommonHaplotype = HaploRes$haplotype[indexMax,]
xtable(haploList[,20:28])

#5
names(table(Y[4,]))[2]
substr(names(table(Y[4,]))[2],1,1) != substr(names(table(Y[4,]))[2],2,2)

listIds=c()
for(i in 1:nrow(Y)) {
  cont=0
  for(j in 1:length(table(Y[i,]))) {
    if(substr(names(table(Y[i,]))[j],1,1) != substr(names(table(Y[i,]))[j],2,2)) {
      cont=cont+table(Y[i,])[names(table(Y[i,]))[j]]
    }
  }
  if(cont>1) {
    listIds[length(listIds)+1]=i;
  }
}
listIds


mostlikely = t(cbind(HaploRes$hap1code[listIds], HaploRes$hap2code[listIds]))
mostlikely[1,]=paste0(mostlikely[,1],"/",mostlikely[,2])
mostlikely = as.data.frame(t(mostlikely[1,]))
colnames(mostlikely)=listIds
xtable(mostlikely)

#6
Ydf = as.data.frame(Y)
drops <- c("rs5999890")
Ynone = Ydf[ , !(names(Ydf) %in% drops)]
Ynone = as.matrix(Ynone)
snpsnone = ncol(Ynone)

Genone = cbind(substr(Ynone[,1],1,1),substr(Ynone[,1],2,2))

for (i in 2:(snpsnone)) {
  Genone = cbind(Genone, substr(Ynone[,i],1,1) ,substr(Ynone[,i],2,2))
}
Snpnamesnone <- paste("SNP",1:snpsnone,sep="")
Snpnamesnone
HaploResnone <- haplo.em(Genone,locus.label=Snpnamesnone,control=haplo.em.control(min.posterior=1e-3))
HaploResnone

HaploResnone$hap.prob
nHaplonone = nrow(HaploResnone$haplotype)
haploListnone = HaploResnone$haplotype
indexMaxnone = which(HaploResnone$hap.prob == max(HaploResnone$hap.prob))
mostCommonHaplotypenone = HaploResnone$haplotype[indexMaxnone,]

