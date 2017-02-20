library(parallel)
library(data.table)
library(Matrix)
library(snpStats)
library(survival)
library(ontologyIndex)
library(ontologySimilarity)
library(jsonlite)
library(BeviMed)
library(reshape2)

data(hpo)

bevimed_chromosome <- function(chr) {

	independent <- read.table(file="/SAN/vyplab/UCLex/KING/independent_noOneKG.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)[[2]]

	tab <- fread(sep=",", header=TRUE, input=paste0("zcat /SAN/vyplab/UCLex/mainset_July2016/mainset_July2016_chr", chr, "_filtered3-genotypes.csv.gz"),stringsAsFactors=FALSE)
	class(tab) <- "data.frame"
	vxs <- as.matrix(tab[,-1])
	rownames(vxs) <- tab[[1]]
	rm(tab)

	batches <- melt(lapply(dir(path="/SAN/vyplab/UCLex/combinedGVCFs/batches/", pattern="\\.samples$", full.name=TRUE), readLines))
	names(batches) <- c("sample","batch")
	batches$sample <- as.character(batches$sample)

	anno <- read.table(check.names=FALSE, sep=",", header=TRUE, file=paste0("/SAN/vyplab/UCLex/mainset_July2016/mainset_July2016_VEP/VEP_chr",chr,".csv.gz"), stringsAsFactors=FALSE)
	pheno_df <- local({
			z <- pipe("python /SAN/vyplab/UCLex/scripts/hpo_db.py")
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

	unique_phenos <- unique(unlist(use.names=FALSE, using_phenos))
	has_pheno_matrix <- sapply(using_phenos, function(x) setNames(nm=unique_phenos, unique_phenos %in% get_ancestors(hpo, x)))
	cases_by_phenotype <- split(has_pheno_matrix, rownames(has_pheno_matrix))

	at_least2 <- cases_by_phenotype[sapply(cases_by_phenotype, function(x) length(x) > 1)]

	sample_batches <- ifelse(colnames(mat) %in% batches$sample, batches$batch[match(colnames(mat), batches$sample)], "missing")

	mats <- lapply(vars_by_gene, function(x) mat[x,,drop=FALSE])
	OK <- sapply(mats, function(x) if (sum(x) == 0) FALSE else chisq.test(x=factor(factor(apply(x > 0, 2, sum) > 0, levels=c(FALSE,TRUE))), y=factor(sample_batches))$p.value > 1e-4)

	rm(mat)

	cat("Running BeviMed...\n")
	print(system.time(bevimed_results <- mclapply(
		mc.cores=16L,
		mc.preschedule=TRUE,
		X=mats[OK],
		FUN=function(ac_mat) {
			Filter(f=function(x) length(x) > 0, x=lapply(
				at_least2,
				function(y) {
					Filter(f=function(x) !is.null(x), x=lapply(
						list(dominant=1, recessive=2),
						function(min_ac) {
							if (sum(apply(ac_mat[,y,drop=FALSE], 2, sum) >= min_ac) > 1) {
								bv <- summary(bevimed(G=t(ac_mat), y=y, min_ac=min_ac))
								list(
									p=bv$prob_gamma1,
									variants=setNames(nm=rownames(ac_mat), bv$conditional_prob_pathogenic)
								)
							} else {
								NULL
							}
						}
					))
				}
			))
		}
	)))
									
	saveRDS(bevimed_results, file=paste0("/SAN/vyplab/UCLex/BeviMed/chr",chr, ".rds"))
}

if (!interactive()) {
	for (chr in as.character(c(1:22, "X"))) {
		cat("Starting chromsome ", chr, "\n", sep="")
		bevimed_chromosome(chr)
	}
}


