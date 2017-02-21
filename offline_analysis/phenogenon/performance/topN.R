# How many known genes are in the top N for both recessive and dominant?
N <- 10
top <- list()
# FP and TN from all genes
for (s in list.files(pattern='*.txt',path='/SAN/vyplab/UCLex/phenogenon/recessive_genes/')) {
    hpo.id <- gsub('.txt','',s)
    #if (!hpo.id %in% results$hpo.id) next
    print(s)
    d <- rbind(read(file.path('/SAN/vyplab/UCLex/phenogenon/dominant_genes/',s)),read(file.path('/SAN/vyplab/UCLex/phenogenon/recessive_genes/',s)))
    d$p_val <- as.numeric(d$p_val)
    d <- do.call('rbind', by(d,d$gene_id,function(x) data.frame(known=x$known,gene_name=x$gene_name,gene_id=x$gene_id,p_val=min(x$p_val))))
    d <- unique(d)
    d <- d[order(d$p_val),]
    d <- d[1:N,]
    top[[hpo.id]] <- length(which(d$known))>0
}
print(table(as.logical(top)))
