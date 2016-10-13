
library(magrittr)
library(ontologyIndex)
library(ontologySimilarity)
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


