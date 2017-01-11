#!/usr/bin/env Rscript

suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(rmongodb , quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))

mongo <- mongo.create(host='localhost',db='uclex')
mongo.drop(mongo,ns='uclex.simreg')
d <- readRDS('all-results.rds')

for (mode in c('dom','rec')) {
    for (gene in  names(d[[mode]])) {
        print(gene)
        x <- d[[mode]][[gene]]
        l <- list(gene=gene,phi=as.list(as.data.frame(t(x$phi))),p=x$p,mode=mode)
        names(l$phi) <- sapply(l$phi, '[[', 1)
        l$phi <- lapply(l$phi,function(x) list(hpo_id=x[[1]],desc=x[[2]],prob=as.numeric(x[[3]])))
        print(l)
        b=mongo.bson.from.list(l)
        mongo.insert(mongo,"uclex.simreg", b)
    }
}




