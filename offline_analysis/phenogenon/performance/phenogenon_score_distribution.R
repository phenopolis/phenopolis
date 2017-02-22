X.recessive <- data.frame()
for (s in list.files(pattern='*.txt',path='/SAN/vyplab/UCLex/phenogenon/recessive_genes/')) {
    print(s)
    hpo.id <- gsub('.txt','',s)
    #d <- rbind(read(file.path('/SAN/vyplab/UCLex/phenogenon/dominant_genes/',s)),read(file.path('/SAN/vyplab/UCLex/phenogenon/recessive_genes/',s)))
    d <- read(file.path('/SAN/vyplab/UCLex/phenogenon/recessive_genes/',s))
    #d <- do.call('rbind', by(d,d$gene_id,function(x) data.frame(known=x$known,gene_name=x$gene_name,gene_id=x$gene_id,p_val=min(x$p_val))))
    #d <- unique(d)
    d$p_val <- as.numeric(d$p_val)
    d <- d[order(d$p_val),]
    d$HPO <- hpo.id
    #print(f <- sprintf( '/SAN/vyplab/UCLex/phenogenon/hpo_individuals/%s.txt', hpo.id))
    #if (!file.exists(f)) next
    #individuals <- read.table(f)[,1]
    #d$num_individuals <- length(individuals)
    X.recessive <- rbind(X.recessive,d)
}

X.recessive$p_val <- as.numeric(X.recessive$p_val)


pdf('~/plot.pdf'); boxplot(-log10(p_val) ~ known, data=X.recessive); dev.off()
