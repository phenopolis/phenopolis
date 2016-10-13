
library(bit)

base.dir <- '/slms/UGI/vm_exports/vyp/phenotips/uclex_files/mainset_February2016/'
files <- list.files(base.dir, pattern='*.bin')

variant.names <- as.character(read.csv(file.path(base.dir, 'variants.csv'))[,1])

variants <- lapply( files, function(f) readRDS( file.path(base.dir, f) ) )
names(variants) <- gsub('.bin','',files)


common.variants <- function(individuals) {
    n <- intersect(individuals,names(variants))
    return(variant.names[which(as.logical(Reduce(function(x,y) x&y, variants[n])))])
}



