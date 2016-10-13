
#source('patients-hpo.R')
library(Rserve)
library(magrittr)
library(ontologyIndex)
library(ontologySimilarity)
library(rmongodb)
#library(ontologyPlot)
#library(SimReg)

hpo <- ontologyIndex::get_ontology('/slms/UGI/vm_exports/vyp/phenotips/HPO/hp.obo')
# our hpo terms
d <- read.table('/slms/UGI/vm_exports/vyp/phenotips/HPO/hpo.txt',header=TRUE,sep='\t',stringsAsFactors=FALSE)
x <- strsplit(d$hpo,',')
names(x) <- d$eid

#map everyone's phenotype to non-redundant set
x <- lapply(x, function(x_i) minimal_set(hpo, c("HP:0000001", x_i)))

# search by hpo term
query.hpo <- function(hpo_term) {
    has_term <- sapply(x, function(phenotype) hpo_term %in% get_ancestors(hpo, phenotype))
    return(names(x)[has_term])
}


#source('variants.R')
library(bit)

base.dir <- '/slms/UGI/vm_exports/vyp/phenotips/uclex_files/current/bin/'
files <- list.files(base.dir, pattern='*.bin')

variant.names <- as.character(read.csv(file.path(base.dir, 'variants.csv'))[,1])

variants <- lapply( files, function(f) readRDS( file.path(base.dir, f) ) )
names(variants) <- gsub('.bin','',files)

lapply(variants, function(x) sum(as.logical(x)))->variant_counts

# for speed of access update mongo
for (n in names(variant_counts)) { mongo.update(mongo,'uclex-old.patients',list(external_id=n),list('$set'=list(variant_count=variant_counts[[n]]))) }

# write variant,individual pair to a file
#library(RPostgreSQL)
# loads the PostgreSQL driver
#drv <- dbDriver("PostgreSQL")
#con <- dbConnect(drv, dbname = "postgres", host = "localhost", port = 5432, user = "rmhanpo")

#for (n in names(variants)) {
for ( n in grep('IRDC',names(variants),value=TRUE) ) {
 for (v in variant.names[as.logical(variants[[n]])]) {
     v<-gsub('_','-',v)
     cat(v,n,'\n')
     res <- mongo.find.one(mongo,'uclex-old.variant_patient',list(variant_id=v,patient_id=n))
     if (is.null(res)) {
         print(mongo.insert(mongo,'uclex-old.variant_patient',list(variant_id=v,patient_id=n)))
     }
}
}

# seen in all our individuals
common_variants <- function(individuals) {
    n <- intersect(individuals,names(variants))
    return(variant.names[which(as.logical(Reduce(function(x,y) x&y, variants[n])))])
}

# seen only in all our individuals (not seen in anyone else)
# might still be seen in exac though so need to check
# this is unlikely to be useful for rare disease
private_variants <- function(individuals) {
    n <- intersect(individuals,names(variants))
    # everyone else
    others <- setdiff(names(variants),n)
    # seen in at least one other
    j <- (Reduce(function(x,y) x|y, variants[others]))
    if (length(n)==1) {
        i <- variants[n]
    } else {
        # seen in all of our individuals
        i <- (Reduce(function(x,y) x&y, variants[n]))
    }
    return(variant.names[which(as.logical(i & !j))])
}

lapply(names(variants), function(x) length(private_variants(x)))->private_variant_counts

for (n in names(private_variant_counts)) {
    mongo.update(mongo,'uclex-old.patients',list(external_id=n),list('$set'=list(private_variant_counts=private_variant_counts[[n]])))
}

run.Rserve()


