#s='HP:3000050.txt'

HPO <- read.csv('/SAN/vyplab/UCLex/phenogenon/hpo.csv')
rownames(HPO) <- HPO$id


# TP=number of literature genes reported as significant by Phenogenon
# FP=number of non-literature genes reported as significant
# TN=number of non-literature genes reported as non-significant
# FN=number of literature genes not reported as significant


setwd('/SAN/vyplab/UCLex/phenogenon/recessive_genes/')
X <- data.frame()
for (s in list.files(pattern='*.txt')) {
print(s)
d <- read(s)
d <- d[1:200,]
X <- rbind(X,data.frame(hpo=hpo[gsub('.txt','',s),'name'],true=length(which(d$known=='TRUE')),false=length(which(d$known=='FALSE'))))
}
print(prop.table(table(X$true>0)))


setwd('/SAN/vyplab/UCLex/phenogenon/dominant_genes/')
X <- data.frame()
for (s in list.files(pattern='*.txt')) {
print(s)
d <- read(s)
d <- d[1:200,]
X <- rbind(X,data.frame(hpo=hpo[gsub('.txt','',s),'name'],true=length(which(d$known=='TRUE')),false=length(which(d$known=='FALSE'))))
}
print(prop.table(table(X$true>0)))


setwd('/SAN/vyplab/UCLex/phenogenon/literature_genes/')
X <- data.frame()
for (s in list.files(pattern='*.txt')) {
print(s)
d <- read(s)
d <- d[-which(d$phenogenon.dominant_pvalue=='None'|d$phenogenon.recessive_pvalue=='None'),]
X <- rbind(X,data.frame(hpo=hpo[gsub('.txt','',s),'name'], prop.signif=length(which(d$phenogenon.recessive_pvalue<0.05|d$phenogenon.dominant_pvalue<0.05))/nrow(d),n=nrow(d)))
}
X<-X[-which(X$n==0),]


print(prop.table(table(X$true>0)))

X <- data.frame()
for (s in list.files(pattern='*.txt')) {
print(s)
d <- read(s)
FP <- length(which(!d$known&d$p_val<0.05))/rnow(d)
#length(which(d$known&d$p_val<0.05))
TN <- length(which(!d$known&d$p_val>0.05))/nrow(d)
data.frame(hpo=,gene_name=gene_name,FP=FP,TN=TN)
}
print(prop.table(table(X$true>0)))
