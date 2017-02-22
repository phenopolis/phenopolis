#`<` -> silly.lt
#`<` <- function(e1,e2) { if (class(e1)!=class(e2)) { stop('love <3') } else { return(silly.lt(e1,e2)) }}


ROC <- data.frame()
results <- data.frame()

for ( pvalue.threshold in c(0.01,5/100,5/1000,5/10000,5/100000) ) {
HPO <- read.csv('/SAN/vyplab/UCLex/phenogenon/hpo.csv')
rownames(HPO) <- HPO$id
# TP and FN from literature genes
X <- data.frame()
for (s in list.files(pattern='*.txt',path='/SAN/vyplab/UCLex/phenogenon/literature_genes/')) {
print(s)
d <- read(file.path('/SAN/vyplab/UCLex/phenogenon/literature_genes/',s))
print(dim(d <- d[-which(d$phenogenon.dominant_pvalue=='None'|d$phenogenon.recessive_pvalue=='None'),]))
print(length(i <- which(as.numeric(d$phenogenon.recessive_pvalue)<pvalue.threshold|as.numeric(d$phenogenon.dominant_pvalue)<pvalue.threshold)))
hpo.id <- gsub('.txt','',s)
X <- rbind(X, data.frame( TP.N=length(i), TP=length(i)/nrow(d), FN.N=nrow(d[-i,]), FN=nrow(d[-i,])/nrow(d), hpo.id=hpo.id, hpo.name=HPO[hpo.id,'name'], TP.FN.N=nrow(d)))

}
X<-X[-which(X$TP.FN.N==0),]

# FP and TN from all genes
X2 <- data.frame()
for (s in list.files(pattern='*.txt',path='/SAN/vyplab/UCLex/phenogenon/recessive_genes/')) {
    print(s)
    hpo.id <- gsub('.txt','',s)
    if (!hpo.id %in% X$hpo.id) next
    d <- rbind(read(file.path('/SAN/vyplab/UCLex/phenogenon/dominant_genes/',s)),read(file.path('/SAN/vyplab/UCLex/phenogenon/recessive_genes/',s)))
    d$p_val <- as.numeric(d$p_val)
    d <- do.call('rbind', by(d,d$gene_id,function(x) data.frame(known=x$known,gene_name=x$gene_name,gene_id=x$gene_id,p_val=min(x$p_val))))
    d <- unique(d)
    X2 <- rbind(X2,data.frame(
    hpo.id=hpo.id,hpo.name=HPO[hpo.id,'name'],
    FP=length(which(!as.logical(d$known)&as.numeric(d$p_val)<pvalue.threshold))/nrow(d),
    FP.N=length(which(!as.logical(d$known)&as.numeric(d$p_val)<pvalue.threshold)),
    TN=length(which(!as.logical(d$known)&as.numeric(d$p_val)>pvalue.threshold))/nrow(d),
    TN.N=length(which(!as.logical(d$known)&as.numeric(d$p_val)>pvalue.threshold)),
    FP.TN.N=nrow(d)))
}

rownames(X) <- X$hpo.id
rownames(X2) <- X2$hpo.id

i <- intersect( rownames(X), rownames(X2) )
X3 <- cbind(X[i,],X2[i,])
# recall, TPR, sensitivity
X3$TPR <- with(X3,TP/(FN+TP))
#X3$TPR <- with(X3, TP/(TP+FN))
X3$precision <- with(X3,TP/(FP+TP))
X3$specificity <- with(X3,TN/(FP+TN))
X3$FPR <- with(X3, FP/(TN+FP))
X3 <- X3[,unique(colnames(X3))]

# number of individuals with HPO term
for (i in 1:nrow(X3)) {
    hpo_id <- X3[i,'hpo.id']
    print(f <- sprintf( '/SAN/vyplab/UCLex/phenogenon/hpo_individuals/%s.txt', hpo_id))
    individuals <- read.table(f)[,1]
    X3[i,'num_individuals'] <- length(individuals)
}

X3$pvalue.threshold <- pvalue.threshold

results <- rbind(results,X3)

#results[[pvalue.threshold]] <- list(results,list(X3))
#ROC <- rbind(ROC,data.frame(pvalue=pvalue.threshold, mean.TPR=mean(X3$TPR),mean.FPR=mean(X3$FPR),N=nrow(X3)))
ROC <- rbind(ROC,data.frame(pvalue=pvalue.threshold, mean.TPR=mean(X3$TPR,na.rm=T),mean.FPR=mean(X3$FPR,na.rmm=T),N=nrow(X3)))
}
