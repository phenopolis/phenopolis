library(parallel)
library(data.table)
library(Matrix)
library(snpStats)
library(survival)
library(ontologyIndex)
library(ontologySimilarity)
library(jsonlite)
library(SimReg)
library(reshape2)

cat("Loading up data...\n")
data(hpo)
chr <- if (interactive()) "20" else commandArgs(trailingOnly=TRUE)[1]

independent <- read.table(file="unrelated.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)[[2]]

tab <- fread(sep=",", header=TRUE, input=paste0("zcat mainset_July2016_chr", chr, "_filtered3-genotypes.csv.gz"),stringsAsFactors=FALSE)
class(tab) <- "data.frame"
vxs <- as.matrix(tab[,-1])
rownames(vxs) <- tab[[1]]
rm(tab)

batches <- melt(lapply(dir(path="batches/", pattern="\\.samples$", full.name=TRUE), readLines))
names(batches) <- c("sample","batch")
batches$sample <- as.character(batches$sample)

anno <- read.table(check.names=FALSE, sep=",", header=TRUE, file=paste0("VEP_chr",chr,".csv.gz"), stringsAsFactors=FALSE)
pheno_df <- local({
        z <- pipe("python hpo_db.py")
        on.exit(close(z))
        lines <- strsplit(readLines(z),split=",")
        do.call(what=data.frame, c(list(check.names=FALSE, stringsAsFactors=FALSE), lapply(setNames(seq(length(lines[[1]])), nm=lines[[1]]), function(x) sapply(lines[-1], "[", x))))
})

cat("Processing data...\n")

phenotypes <- lapply(lapply(setNames(nm=pheno_df$eid, strsplit(pheno_df$observed_features, split=";")), intersect, hpo$id[!hpo$obsolete]), function(x) if (length(x) == 0) "HP:0000001" else x)


maf_mat <- as.matrix(anno[,c("AFR_MAF", "AMR_MAF", "EAS_MAF", "EUR_MAF", "SAS_MAF", "AA_MAF", "EA_MAF", "ExAC_MAF")])

anno_inds <- which(
        (is.na(anno$CADD) | anno$CADD > 15) &
        # MAF at least 0.01 in each population
        apply(ifelse(is.na(maf_mat), 0, maf_mat), 1, max) < 0.01 &
        (anno$IMPACT == "MODERATE" | anno$IMPACT == "HIGH") &
        anno$`#Uploaded_variation` %in% rownames(vxs)
)

vars_by_gene <- split(anno$`#Uploaded_variation`[anno_inds], anno$Gene[anno_inds])

mat <- local({
        raw <- vxs[unlist(use.names=FALSE, vars_by_gene),intersect(independent, intersect(colnames(vxs), names(phenotypes))),drop=FALSE]
        ifelse(is.na(raw), 0, raw)
})

using_phenos <- phenotypes[colnames(mat)]
at_least2 <- names(which(table(unlist(lapply(using_phenos, get_ancestors, ontology=hpo))) > 2))
using_phenos <- lapply(phenotypes, function(x) { t <- minimal_set(hpo, intersect(x, at_least2)); if (length(t) == 0) "HP:0000001" else t })

sample_batches <- ifelse(colnames(mat) %in% batches$sample, batches$batch[match(colnames(mat), batches$sample)], "missing")

dom_ys <- Filter(f=function(x) chisq.test(x=factor(x, levels=c(FALSE,TRUE)), y=factor(sample_batches))$p.value > 1e-3, x=Filter(x=lapply(vars_by_gene, function(x) apply(mat[x,,drop=FALSE], 2, sum) > 0), f=function(x) { sum(x) >= 2 & mean(x) < 0.1 }))
rec_ys <- Filter(f=function(x) chisq.test(x=factor(x, levels=c(FALSE,TRUE)), y=factor(sample_batches))$p.value > 1e-3, x=Filter(x=lapply(vars_by_gene, function(x) apply(mat[x,,drop=FALSE], 2, sum) > 1), f=function(x) { sum(x) >= 2 & mean(x) < 0.1 }))


saveRDS(dom_ys, file=paste0("/SAN/vyplab/UCLex/SimReg/chr",chr, "-dom-y.rds"))
saveRDS(rec_ys, file=paste0("/SAN/vyplab/UCLex/SimReg/chr",chr, "-rec-y.rds"))

rm(mat)
ic <- get_term_info_content(hpo, using_phenos)
tsm <- get_term_sim_mat(hpo, ic)

cat("Running SimReg...\n")
print(system.time(sr_results <- mclapply(
        mc.cores=16L,
        mc.preschedule=TRUE,
        X=dom_ys,
        FUN=function(y) {
                sr <- sim_reg(ontology=hpo, y=y, x=using_phenos, information_content=ic, sim_params=list(term_sim_mat=tsm))
                phi <- Filter(f=function(x) x > 0, x=round(digits=2, sort(decreasing=TRUE, get_term_marginals(sr))))
                bf <- sum_log_probs(sr$ml)

                list(phi=data.frame(row.names=NULL, check.names=FALSE, stringsAsFactors=FALSE, code=names(phi), name=hpo$name[names(phi)], prob=phi), p=0.01 * exp(bf) / (0.99 + 0.01 * exp(bf)))
        }
)))

saveRDS(sr_results, file=paste0("chr",chr, "-dom.rds"))

cat("Running SimReg recessive...\n")

print(system.time(sr_results_rec <- mclapply(
        mc.cores=16L,
        mc.preschedule=TRUE,
        X=rec_ys,
        FUN=function(y) {
                sr <- sim_reg(ontology=hpo, y=y, x=using_phenos, information_content=ic, sim_params=list(term_sim_mat=tsm))
                phi <- Filter(f=function(x) x > 0, x=round(digits=2, sort(decreasing=TRUE, get_term_marginals(sr))))
                bf <- sum_log_probs(sr$ml)

                list(phi=data.frame(row.names=NULL, check.names=FALSE, stringsAsFactors=FALSE, code=names(phi), name=hpo$name[names(phi)], prob=phi), p=0.01 * exp(bf) / (0.99 + 0.01 * exp(bf)))
        }
)))

saveRDS(sr_results_rec, file=paste0("chr",chr, "-rec.rds"))

