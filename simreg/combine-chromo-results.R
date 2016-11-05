sr_results <- lapply(setNames(nm=c("dom","rec")), function(moi) { 
    yl <- lapply(setNames(nm=dir(pattern=paste0(moi, "\\.rds$"))), readRDS)
    nm <- do.call(what=c, lapply(yl, names))
    if (any(duplicated(nm))) stop("same genes on two chromosomes?")
    setNames(nm=nm, do.call(what=c, yl))
})

saveRDS(sr_results, file=paste0("/SAN/vyplab/UCLex/SimReg/all-results.rds"))
