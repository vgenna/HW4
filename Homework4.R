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
length(HaploRes$hap1code)
length(HaploRes$hap2code)
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
which(listIds==20)

effectiveListID = c()
for(i in 1:length(listIds)) {
  if(i>which(listIds==20)) {
    effectiveListID[i]=listIds[i]+2
  }
  else {
    effectiveListID[i]=listIds[i]
  }
}


mostlikely = t(cbind(HaploRes$hap1code[effectiveListID], HaploRes$hap2code[effectiveListID]))
mostlikely=paste0(mostlikely[1,],"/",mostlikely[2,])
mostlikely = as.data.frame(t(mostlikely))
colnames(mostlikely)=listIds
xtable(mostlikely[,37:54])

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

locus = function(x){
  if (length(names(table(haploListnone[,x]))) > 1){
    geno = paste0(names(table(haploListnone[,x]))[1],names(table(haploListnone[,x]))[2])
  }
  else{
    geno = paste0(names(table(haploListnone[,x]))[1],names(table(haploListnone[,x]))[1])
  }
  return(geno)
} 

new_locus = sapply(1:ncol(haploListnone), locus)
unique(new_locus)

#7
mostlikelyAll = t(cbind(HaploRes$hap1code, HaploRes$hap2code))
mostlikelyAll=paste0(mostlikelyAll[1,],"/",mostlikelyAll[2,])
mostlikelyAll = mostlikelyAll[-22]
mostlikelyAll = mostlikelyAll[-22]
mostlikelyAll = (t(mostlikelyAll))

newLocusAlleles = c()

for(i in 1:nrow(Y)) {
  first = as.numeric(substr(mostlikelyAll[i],1,1))
  second = as.numeric(substr(mostlikelyAll[i],3,3))
  firstHap = HaploRes$haplotype[first, ]
  secondHap = HaploRes$haplotype[second, ]
  indivMostProb = paste0(firstHap,secondHap)
  newLocusAlleles[i]=names(sort(table(indivMostProb),decreasing=TRUE))[1]
}

mostLikelyGen = names(sort(table(newLocusAlleles),decreasing=TRUE))[1]
mostLikelyGen_prob = 100*sort(table(newLocusAlleles),decreasing=TRUE)[1]/sum(table(newLocusAlleles))

#8

maf = function(data){
  freqA = (data$AA * 2 + data$AB)/((data$AA + data$AB + data$BB)*2)
  freqB = (data$BB * 2 + data$AB)/((data$AA + data$AB + data$BB)*2)
  out = pmin(freqA, freqB)
  return(out)
}

mat = matrix(nrow = ncol(Y), ncol = 3)
alleles = matrix(nrow = ncol(Y), ncol = 3)
for(i in 1:ncol(Y)) {
  if(length(names(table(Y[,i])))==1) {
    line = c(table(Y[,i])[1], 0, 0)
    lineAL=c("x","x",names(table(Y[,i]))[1])
  } 
  else if (length(names(table(Y[,i])))==2) {
    if(substr(names(table(Y[,i]))[2],1,1) != substr(names(table(Y[,i]))[2],2,2)) {
      line = c(table(Y[,i])[1], table(Y[,i])[2], 0)
      lineAL=c("x",names(table(Y[,i]))[2],names(table(Y[,i]))[1])
    } else if(substr(names(table(Y[,i]))[1],1,1) != substr(names(table(Y[,i]))[1],2,2)) {
        line = c(0,table(Y[,i])[1] , table(Y[,i])[2])
        lineAL=c("x",names(table(Y[,i]))[1],names(table(Y[,i]))[2])
      } 
      else {
          line = c(table(Y[,i])[1], 0, table(Y[,i])[2])
          lineAL=c(names(table(Y[,i]))[1],"x",names(table(Y[,i]))[2])
      }
    }
  else if (length(names(table(Y[,i])))==3) {
    line = c(table(Y[,i])[1], table(Y[,i])[2], table(Y[,i])[3])
    lineAL=c(names(table(Y[,i]))[1],names(table(Y[,i]))[2],names(table(Y[,i]))[3])
  }
  mat[i,]=line 
  alleles[i,]=lineAL
}

mat = as.data.frame(mat)
alleles = as.data.frame(alleles)
colnames(mat) = c("AA", "AB", "BB")

MAF = maf(mat)

p = MAF[1]
q=1-p
simulatedCounts = rmultinom(1,139, c(p^2,2*p*q,q^2))

for(i in 2:length(MAF)){
  p = MAF[i]
  q=1-p
  simulatedCounts = cbind(simulatedCounts, rmultinom(1,139, c(p^2,2*p*q,q^2)) )
}

simulatedCounts = as.data.frame(t(simulatedCounts))

alleleMat=as.data.frame(matrix(nrow = 139, ncol = 28))
for(i in 1:nrow(simulatedCounts)) {
  line = c(as.character(rep(alleles[i,1], simulatedCounts[i,1])), as.character(rep(alleles[i,2], simulatedCounts[i,2])), as.character(rep(alleles[i,3], simulatedCounts[i,3])))
  alleleMat[,i]=line
}

alleleMat[1,13]="AA"
alleleMat[1:2,14]="CC"
alleleMat[1:3,15]="AA"
alleleMat[1:2,19]="CC"
alleleMat[1,23]="CC"

alleleMat = alleleMat[sample(1:nrow(alleleMat), replace = FALSE),]

Geno8 = cbind(substr(alleleMat[,1],1,1),substr(alleleMat[,1],2,2))

for (i in 2:snps) {
  Geno8 = cbind(Geno8, substr(alleleMat[,i],1,1) ,substr(alleleMat[,i],2,2))
}

Snpnames8 <- paste("SNP",1:snps,sep="")
Snpnames8
HaploRes8 <- haplo.em(Geno8,locus.label=Snpnames,control=haplo.em.control(min.posterior=1e-3))
HaploRes8








