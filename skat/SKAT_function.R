## how to use me. 
#  Rscript='/share/apps/R/bin/Rscript'
#  SKAT=/SAN/vyplab/UCLex/scripts/DNASeq_pipeline/ciangene/Gene_Based_Tests/SKAT/SKAT_function.R
#  CaseFile = path to one column file containing list of cases
#  ControlFile = path to one column file containing list of cases. If not specified, will use all non cases. 
#  $Rscript $SKAT --case.list $CaseFile --oDir outputDirectory --control.list $ControlFile

### Changes ####
# fixed bug with SKAT covariates
# added hwe filter - specify pval for cut off from exact test in controls.
# Fixed bug with zero count matrices
# Added HWE filter for variants - controls only
# added top 2 techPCs as covariates.
# Added CADD and ExAC maf filter. 
# Added ability to specify minReadDepth for filtering - default zero is fastest. 
# Added ability to specify genes you're interested in testing: --TargetGenes for a file with names or --SampleGene for a single gene name
# Added basic compoundHeterozygote Function and pvalue.  
# changed missingess to missing_cutoff=0.2. So variants with missingess >20% will be removed. 
# added fisher test and odds ratio. 
# changed nb.cases/ctrls to nb non NA calls (so nb.cases now is nb.patients*nb.variants)
# changed case.count/ctlr.count to allele counts. 1 het and 1 hom is a count of 3 for nb.variants.cases 
# now not removing non caucasians - adding top two PCs as covariates into SKATO instead. 
# removing related individuals. ReEdit: ## Turned off for now. 
# Fixed snp/gene matching
# Keeping only damaging variants
########

ldak='/cluster/project8/vyp/cian/support/ldak/ldak'

suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("SKAT"))
suppressPackageStartupMessages(library("snpStats"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("HardyWeinberg"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("data.table"))

option_list <- list(
	make_option(c("--chrom"), default=NULL,help="Chromosome"),
	make_option(c("-v", "--verbose"), action="store_true", default=TRUE,help="Print extra output [default]"),
 	make_option(c("--case.list"),  help="one column file containing list of cases",type='character'),
 	make_option(c("--control.list"), default=NULL, help="one column file containing list of controls",type='character'),
 	make_option(c("--SampleGene"), default=NULL, help="Gene Symbol eg ABCA4",type='character'),
 	make_option(c("--TargetGenes"), default=NULL, help="Gene File",type='character'),
 	make_option(c("--MinReadDepth"), default=10, help="Specify MinReadDepth",type='character'),
 	make_option(c("--SavePrep"), default=TRUE, help="Do you want to save an image of setup?",type='character'),
 	make_option(c("--minCadd"), default=30, help="minimum CADD score for retained variants",type='character'),
 	make_option(c("--maxExac"), default=0.001, help="Max EXAC maf for retained variants",type='character'),
 	make_option(c("--oDir"),default='SKATtest',type='character', help='Before per gene testing save of environment'),
 	make_option(c("--MinSNPs"),default=2,type='character',help='Minimum number of SNPs needed in gene.'),
 	make_option(c("--PlotPCA"),default=TRUE,type='character',help='If yes, the PCA graphs highlighting cases will be included'), 
 	make_option(c("--homozyg.mapping"),default=FALSE,type='character',help='Right now, homozygous mapping does not work...'),  
 	make_option(c("--MaxMissRate"),default=10,type='character',help='maximum per snp missingess rate'), 
 	make_option(c("--HWEp"),default=0 ,type='character',help='what hardy weinberg pvalue cut off for control disequilibrium to use to remove SNPs '), 
 	make_option(c("--compoundHets"),default=NULL,type='character',help='added function to look for compound Heterozygotes'),
 	make_option(c("--Release"),default='July2016',type='character',help='what release of UCLex do you want to use') ,
 	make_option(c("--MaxCtrlMAF"),default=0.1,type='character',help='Remove snps with maf above this in controls.') ,
 	make_option(c("--qcPREP"),default=FALSE,type='character',help='not used much, but handy for troubleshooting to save environment later on in function') 
 )


opt <- parse_args(OptionParser(option_list=option_list))
if ( opt$verbose ) {
 write("Starting Argument checks...\n", stderr())
}

MaxCtrlMAF<-as.numeric(opt$MaxCtrlMAF) 
qcPREP <- opt$qcPREP
chrom <- opt$chrom
case.list<-opt$case.list
control.list<-opt$control.list
outputDirectory<-opt$oDir
PlotPCA<-opt$PlotPCA
compoundHets<-opt$compoundHets
homozyg.mapping<-opt$homozyg.mapping
MaxMissRate<-as.numeric(opt$MaxMissRate)/100
min.depth<-as.numeric(opt$MinReadDepth) 
SavePrep<-opt$SavePrep
maxExac<-as.numeric(opt$maxExac)
minCadd<-as.numeric(opt$minCadd)
MinSNPs<-as.numeric(opt$MinSNPs)
HWEp<-as.numeric(opt$HWEp)
UCLex_Release<-opt$Release

if( (!is.null(opt$SampleGene)) && ( !is.null(opt$TargetGenes)) ) stop("Please supply either a single gene to SampleGene or a file of gene names to TargetGenes. Not both.")
if(!is.null(opt$SampleGene))
{
 	TargetGenes<-opt$SampleGene
 	message(paste("Setting up environment to test for gene:",TargetGenes))
} 
if(!is.null(opt$TargetGenes))
{
	TargetGenes<-read.table(opt$TargetGenes,header=FALSE)[,1]
	message(paste("Read",length(TargetGenes),'genes from file',paste0('--',opt$TargetGenes,'--')))
} 


## Check cases
if(!file.exists(case.list))stop("Case file doesn't exist. ")
message(paste('Reading cases from',paste0('--',case.list,'--')))
case.list<-read.table(case.list,header=FALSE)
message(paste('Found',nrow(case.list),'cases'))
case.list<-case.list[,1]


## Set up controls
if(is.null(control.list))message("Controls not specified, so will use all non case samples as controls.")
if(!is.null(control.list))
{
	if(!file.exists(control.list))message("Control file doesn't exist, so will use all non case samples as controls.")
	message(paste('Reading controls from',paste0('--',control.list,'--')))
	control.list<-read.table(control.list,header=FALSE,stringsAsFactors=FALSE)
	message(paste('Found',nrow(control.list),'controls'))
	control.list<-control.list[,1]
}

## make or break oDir
if(outputDirectory=='SKATtest')message('Output directory not specified so will use default folder of ./SKATtest')
if(outputDirectory!='SKATtest')print(paste("Outputting Results to:",outputDirectory)) 

message("Finished Argument checks.\n")
message("Starting test setup.\n")

doSKAT<-function(case.list=case.list,control.list=control.list,outputDirectory=outputDirectory,TargetGenes=TargetGenes,
	min.depth=min.depth,minCadd=minCadd,maxExac=maxExac,MinSNPs=MinSNPs,compoundHets=compoundHets,MaxMissRate=MaxMissRate,
	chrom=chrom,PlotPCA=PlotPCA,homozyg.mapping=homozyg.mapping,HWEp=HWEp,
	MaxCtrlMAF=MaxCtrlMAF,SavePrep=SavePrep)
{
	outputDirectory<-paste0(outputDirectory,'/')
	rootODir<-paste0('/SAN/vyplab/UCLex/mainset_',UCLex_Release,'/cian/') 
	if(file.exists(outputDirectory))file.remove(outputDirectory)
	dir.create(outputDirectory,recursive=TRUE)
	qc<-paste0(outputDirectory,'/qc/')
	if(!file.exists(qc))dir.create(qc)

	oFile<-paste0(outputDirectory,'SKAT_results_by_SNP.tab')
	if(file.exists(oFile)) file.remove(oFile)
	caseFile<-paste0(outputDirectory,'case_carriers')
	if(file.exists(caseFile)) file.remove(caseFile)
	ctrlFile<-paste0(outputDirectory,'ctrl_carriers')
	if(file.exists(ctrlFile)) file.remove(ctrlFile)
	compoundFileCases<-paste0(outputDirectory,'CompoundHets_cases')
	compoundFileCtrls<-paste0(outputDirectory,'CompoundHets_ctrls')
	if(file.exists(compoundFileCases)) file.remove(compoundFileCases)
	if(file.exists(compoundFileCtrls)) file.remove(compoundFileCtrls)

	message("Reading in snp/gene database")
	snp.gene.base<-as.data.frame(fread(paste0(rootODir,'snp_gene_id'),header=F))
	colnames(snp.gene.base)<-c("SNP",'ENSEMBL')
    print(dim(snp.gene.base))
    if (!is.null(chrom)) {
        snp.gene.base <- snp.gene.base[grep(sprintf('^%s_',chrom),snp.gene.base[,'SNP']),]

    }
    print(dim(snp.gene.base))
	gene.dict<-as.data.frame(fread(paste0(rootODir,'gene_dict_skat')))
	gene.dict$ENST<-NULL
	gene.dict<-unique(gene.dict)
	print(dim(snp.gene<-merge(gene.dict,snp.gene.base,by="ENSEMBL",all.y=T)))
	
	## This anno file contains Exac mafs, CADD scores etc. Filter SNP list based on user input. 
	snp.annotations<-as.data.frame(fread(paste0(rootODir,'Annotations/func.tab')))

    print(dim(snp.annotations))
    if (!is.null(chrom)) {
        snp.annotations <- snp.annotations[grep(sprintf('^%s_',chrom),snp.annotations[,'SNP']),]

    }
    print(dim(snp.annotations))

	filtered.snp.list<-subset(snp.annotations,snp.annotations$CADD>=minCadd & snp.annotations$ExAC_MAF<=maxExac)$SNP
	message(paste(length(filtered.snp.list),'SNPs kept after CADD and Exac filters of',minCadd,'and',maxExac,'respectively.')) 

	snplist.file<-paste0(outputDirectory,'qc/SNPlist')
	ldak.subset.file<-paste0(outputDirectory,'qc/SNPlist_ldak')
	if(!is.null(TargetGenes)) ## if we're testing a couple genes only, Im subsetting SNPmatrix to make it faster to read in. 
	{
		data<-paste0(rootODir,'allChr_snpStats_out') 
		#data<-'/SAN/vyplab/UCLex/mainset_July2016/cian/UCLex_phased/UCLex_phased_out'
		#data<-paste0(rootODir,'allChr_snpStats')

		target.snp.info<- snp.gene[ snp.gene$Symbol %in% TargetGenes,]
		target.snps<- target.snp.info$SNP[target.snp.info$SNP %in% filtered.snp.list]
		write.table(target.snps,snplist.file,col.names=F,row.names=F,quote=F,sep='\t')
		message(paste(length(target.snps),'SNPs found in target genes'))
		message("Extracting these from dataset...")
		run<-paste(ldak,'--make-sp',ldak.subset.file,'--bfile',data,'--extract',snplist.file)
		system(run)
		message("Reading in SNP data")
		snp.data<-paste0(outputDirectory,TargetGenes[1],'_out')
		snp.sp<-as.data.frame(fread(paste0(ldak.subset.file,'_out.sp'),header=F))
		snp.bim<-as.data.frame(fread(paste0(ldak.subset.file,'_out.bim'),header=F))
		snp.fam<-as.data.frame(fread(paste0(ldak.subset.file,'_out.fam'),header=F))
		colnames(snp.sp)<-snp.fam[,1]
		rownames(snp.sp)<-snp.bim[,2]
	}

	if(is.null(TargetGenes))
	{
		write.table(filtered.snp.list,snplist.file,col.names=F,row.names=F,quote=F,sep='\t')
		message("Extracting from dataset...")
		run<-paste(ldak,'--make-sp',ldak.subset.file,'--sp',data,'--extract',snplist.file)
		system(run)
		rgh<-'/SAN/vyplab/UCLex/mainset_July2016/cian/UCLex_phased'
		snp.sp<-read.table(paste0(rgh,'UCLex_phased_out.sp'),header=F)
		snp.bim<-read.table(paste0(rgh,'UCLex_phased_out.bim'),header=F)
		snp.fam<-read.table(paste0(rgh,'UCLex_phased_out.fam'),header=F)
		#snp.sp<-read.table(paste0(rootODir,'allChr_snpStats.sp'),header=F)
		#snp.bim<-read.table(paste0(rootODir,'allChr_snpStats.bim'),header=F)
		#snp.fam<-read.table(paste0(rootODir,'allChr_snpStats.fam'),header=F)
		#snp.sp<-read.table(paste0(rootODir,'gene_test_variants_out.sp'),header=F)
		#snp.bim<-read.table(paste0(rootODir,'gene_test_variants_out.bim'),header=F)
		#snp.fam<-read.table(paste0(rootODir,'gene_test_variants_out.fam'),header=F)
		colnames(snp.sp)<-snp.fam[,1]
		rownames(snp.sp)<-snp.bim[,2]
	}

		
	current.pheno<-rep(NA,ncol(	snp.sp))
	current.pheno[colnames(snp.sp)%in%case.list]<-1

	## ancestry matching
	ancestry<-read.table(paste0(rootODir,'UCLex_samples_ancestry'),header=T)
	techPCs<-read.table(paste0(rootODir,'TechPCs.vect'),header=F)
	depthPCs<-read.table(paste0(rootODir,'DepthPCs.vect'),header=F)

	#caucasians<-ancestry$V1[ancestry$Caucasian]
	#snps<-snp.sp[,colnames(snp.sp)%in%caucasians]
	## these are teh snps that have been filterd by func and maf
	unrelated<-read.table(paste0(dirname(rootODir),"/kinship/UCL-exome_unrelated.txt")) ## Keeping only unrelated individiuals
	#snp.sp<-snp.sp[,colnames(snp.sp)%in%unrelated[,1]]
	message(paste('Read Depth filter set at:',min.depth ) ) 

	if(min.depth>0)
	{
		message("Reading in Read Depth data")
		read.depth<-as.data.frame(fread(paste0(rootODir,'Average.read.depth.by.gene.by.sample.tab'),header=T))
		read.depth.vals<-read.depth[,5:ncol(read.depth)]
		sample.means<-colMeans(read.depth.vals)
		good.samples<-which(sample.means>=min.depth)
		##################################################################
		## should change this to only calculate mean based on genes in targetgenes/samplegenes. 
		##################################################################
		
		bad.samples<-names(which(sample.means<min.depth))
		bad.cases<-bad.samples[bad.samples %in% case.list]
		if(length(bad.cases)>0)print(bad.cases)
		write.table(bad.samples,paste0(qc,'samples_removed_because_of_low_read_depth.tab'),col.names=F,row.names=F,quote=F,sep='\t')
		write.table(bad.cases,paste0(qc,'cases_removed_because_of_low_read_depth.tab'),col.names=F,row.names=F,quote=F,sep='\t')

		gene.means<-rowMeans(read.depth.vals)
		good.genes<-which(gene.means>=min.depth)
		bad.genes<-read.depth$Gene[which(gene.means<min.depth)]
		write.table(bad.genes,paste0(qc,'genes_removed_because_of_low_read_depth.tab'),col.names=F,row.names=F,quote=F,sep='\t')

		good.genes.snps<-snp.gene$SNP[snp.gene$ENSEMBL %in% read.depth$Gene[good.genes]]
		clean.snp.data<- snp.sp[rownames(snp.sp) %in% good.genes.snps ,  colnames(snp.sp) %in% colnames(read.depth)[colnames(read.depth) %in% names(good.samples) ] ]

		snp.gene.clean<-snp.gene[snp.gene$SNP %in%rownames(clean.snp.data),]
		write.table(snp.gene.clean,paste0(qc,'read.depth.ancestry.func.clean.snps.tab'),col.names=F,row.names=F,quote=F,sep='\t')

		good.genes.data<-snp.gene[snp.gene$ENSEMBL %in% read.depth$Gene[good.genes],]
		uniq.genes<-unique(good.genes.data$ENSEMBL)
		nb.genes<-length(uniq.genes)
		rm(read.depth)
		message("Finished dealing with read depth")
	} else
	{
		message("No read depth filter specified. Skipping step.")
		clean.snp.data<- snp.sp
		good.genes.data<-snp.gene
		uniq.genes<-unique(good.genes.data$ENSEMBL)
		nb.genes<-length(uniq.genes)
	}

	###### clean
	rm(snp.sp,snp.gene.base)
	######

	if(homozyg.mapping)
	{
	    data<-paste0(rootODir,'read_depth/Depth_Matrix') 
	    #data<-'/SAN/vyplab/UCLex/mainset_July2016/cian/UCLex_phased/UCLex_phased_out'
	    #data<-paste0(rootODir,'allChr_snpStats')
	    ldak.subset.file.RD<-paste0(outputDirectory,'qc/SNPlist_ldak_depth')
	    message("Extracting these from dataset...")
	    run<-paste(ldak,'--make-sp', ldak.subset.file.RD,'--sp',data,'--extract',snplist.file)
	    system(run)
	    message("Reading in SNP data")
	    snp.sp<-as.data.frame(fread(paste0(ldak.subset.file.RD,'_out.sp'),header=F))
	    snp.bim<-as.data.frame(fread(paste0(ldak.subset.file.RD,'_out.bim'),header=F))
	    snp.fam<-as.data.frame(fread(paste0(ldak.subset.file.RD,'_out.fam'),header=F))
	    colnames(snp.sp)<-snp.fam[,1]
	    rownames(snp.sp)<-snp.bim[,2]
	    positions<-data.frame(t(data.frame(strsplit(snp.annotations$SNP,'_'))))
      	snp.annotations$Start<-positions[,2]
      	snp.annotations$Chr<-positions[,1]
#      	snp.sp<-snp.sp[rownames(snp.sp) %in% snp.annotations$SNP,colnames(snp.sp)%in%colnames(clean.snp.data)]
		snp.read.depth<-rowMeans(snp.sp,na.rm=T)
	}


	current.pheno<-rep(NA,ncol(clean.snp.data))
	current.pheno[colnames(clean.snp.data)%in%case.list]<-1

	my.cases<-colnames(clean.snp.data)[colnames(clean.snp.data)%in%case.list]
	print(my.cases)
	print(length(my.cases))
	if(!is.null(control.list))
	{
		current.pheno[colnames(clean.snp.data)%in%control.list]<-0
	} else current.pheno[!colnames(clean.snp.data)%in%case.list]<-0
	my.controls<-colnames(clean.snp.data)[current.pheno==0]
	bk.pheno<-current.pheno
	
	write.table(data.frame(colnames(clean.snp.data),current.pheno),paste0(outputDirectory,'PhenotypeFile'),col.names=FALSE,row.names=FALSE,quote=FALSE)
	write.table(my.cases,paste0(outputDirectory,'case.list'),col.names=FALSE,row.names=FALSE,quote=FALSE)
	write.table(my.controls,paste0(outputDirectory,'control.list'),col.names=FALSE,row.names=FALSE,quote=FALSE)

	if(PlotPCA)
	{
		prepPlot<-function(pcaData,caseList=case.list)
		{
			pcaData$is.case<-FALSE
			pcaData$is.case[pcaData[,1] %in% caseList]<-TRUE
			pcaData$is.case<-as.factor(pcaData$is.case)
			colnames(pcaData)[3:4]<-c('PC1','PC2')
			pcaData
		}
		ancestry.plot<-prepPlot(ancestry)
		tech.plot<-prepPlot(techPCs)
		depth.plot<-prepPlot(depthPCs)
		pca.plots<-paste0(outputDirectory,'qc/case_pca.plot.pdf')
		message(paste('PCA plot of cases is in',pca.plots))
		pdf(pca.plots) 

		pca.plot<- ggplot(ancestry.plot, aes(x=PC1, y=PC2,  
					color=is.case))+ geom_point(alpha=0.5) +ggtitle("Ancestry PCA")
		pca.plot<- pca.plot +annotate("point",ancestry.plot$PC1[ancestry.plot$is.case], ancestry.plot$PC2[ancestry.plot$is.case])
  		print(pca.plot)

		pca.plot<- ggplot(tech.plot, aes(x=PC1, y=PC2,  
					color=is.case))+ geom_point(alpha=0.5) +ggtitle("SNP Missingness PCA")
		pca.plot<- pca.plot +annotate("point",tech.plot$PC1[tech.plot$is.case], tech.plot$PC2[tech.plot$is.case])
  		print(pca.plot)

		pca.plot<- ggplot(depth.plot, aes(x=PC1, y=PC2,  
					color=is.case))+ geom_point(alpha=0.5) +ggtitle("Read Depth PCA")
		pca.plot<- pca.plot +annotate("point",depth.plot$PC1[depth.plot$is.case], depth.plot$PC2[depth.plot$is.case])
  		print(pca.plot)

		dev.off()
	}

	## Make the outptu dataframe
	cols<-c("Gene",'SKATO','nb.snps','nb.cases','nb.ctrls','nb.alleles.cases','nb.alleles.ctrls','case.maf','ctrl.maf','total.maf','nb.case.homs',
		'nb.case.hets','nb.ctrl.homs','nb.ctrl.hets','Chr','Start','End','FisherPvalue','OddsRatio','CompoundHetPvalue','minCadd','maxExac','min.depth',
		'MeanCallRateCases','MeanCallRateCtrls','MaxMissRate','HWEp','MinSNPs','SNPs','MaxCtrlMAF'
		)
	results<-data.frame(matrix(nrow=nb.genes,ncol=length(cols)))
	colnames(results)<-cols
	results$Gene<-uniq.genes
	results$minCadd<-minCadd
	results$maxExac<-maxExac
	results$min.depth<-min.depth
	results$HWEp<-HWEp
	results$MinSNPs<-MinSNPs
	results$MaxCtrlMAF<-MaxCtrlMAF
	results<-merge(gene.dict,results,by.y='Gene',by.x='ENSEMBL',all.y=T)
	srt<-data.frame(1:length(uniq.genes),uniq.genes)
	results<-merge(results,srt,by.y='uniq.genes',by.x='ENSEMBL')
	results<-results[order(results[,(ncol(results))]),]
	#results<-results[,1:(ncol(results)-1)]
	uniq.genes<-unique(results[,1])
	nb.genes<-length(uniq.genes)

	if(!is.null(TargetGenes))
	{
		results<-results[,1:(ncol(results)-1)]
		results<-results[results$Symbol %in% TargetGenes,]
		if(nrow(results)==0)stop('Specified genes not found ffs.')
		uniq.genes<-unique(results[,1])
		nb.genes<-length(uniq.genes)
	}

	if(SavePrep)
	{
		robj<-paste0(outputDirectory,'qc/test_setup.RData')
		message(paste('Saving workspace image to', robj))
		save(list=ls(environment()),file=robj)
	}	
	###################################
	message(paste("Starting tests on", nb.genes, 'genes')) 
	###################################

	for(gene in 1:nb.genes)
	{
		current.pheno<-bk.pheno
		print(paste('Gene xx',uniq.genes[gene], gene))
		gene.snps<-good.genes.data$SNP[ grep(uniq.genes[gene],good.genes.data$ENSEMBL) ]

		gene.data<- data.frame(t(data.frame(strsplit(gene.snps,'_')))) 
		gene.chr<-as.numeric(unique(gene.data[,1]))
		gene.start<-min(as.numeric(unique(gene.data[,2])))
		gene.end<-max(as.numeric(unique(gene.data[,2])))

		gene.snp.data<-clean.snp.data[ rownames(clean.snp.data) %in% gene.snps ,]
		nb.snps.in.gene<-nrow(gene.snp.data)
		#print(paste(nb.snps.in.gene,'snps in', uniq.genes[gene],'before full filtering complete')) 

		results$Chr[gene]<-gene.chr
		results$Start[gene]<-gene.start
		results$End[gene]<-gene.end
		results$MaxMissRate<-MaxMissRate
		print(paste('nb.snps.in.gene=',nb.snps.in.gene))
		if(nb.snps.in.gene>=MinSNPs)
		{
			test.gene<-unique(snp.gene$ENSEMBL[snp.gene$SNP %in% rownames(gene.snp.data)])  ## match to uniq. genes as a check, 
			if(length(unique(test.gene))>1) message ("SNPs span multiple genes")
			if(test.gene!=uniq.genes[gene]) stop ("Genes not sorted correctly")
			if(results$ENSEMBL[gene]!=uniq.genes[gene]) stop ("Genes not sorted correctly")

			case.snps<-gene.snp.data[,which(current.pheno==1)]
			ctrl.snps<-gene.snp.data[,which(current.pheno==0)]

			if(sum(as.matrix(case.snps),na.rm=T)>0)
			{
				case.snpStats<-invisible(new("SnpMatrix",t(case.snps+1)))
				case.summary<-col.summary(case.snpStats)
				if(!is.null(MaxMissRate)) #snps well covered in cases
				{
					good.qual.snps<-rownames(case.summary[case.summary$Call.rate>=as.numeric(MaxMissRate),]) 
				} else good.qual.snps<-rownames(ctrl.summary) 
				if(length(good.qual.snps)>0)
				{
					case.snps<-case.snps[rownames(case.snps)%in%good.qual.snps,]
					ctrl.snps<-ctrl.snps[rownames(ctrl.snps)%in%good.qual.snps,]
				}
			}

			if(sum(as.matrix(ctrl.snps),na.rm=T)>0)
			{
				ctrl.snpStats<-invisible(new("SnpMatrix",t(ctrl.snps+1)))
				ctrl.summary.by.sample<-row.summary(ctrl.snpStats)

				if(!is.null(MaxMissRate)) ## first select good qual controls
				{
					good.ctrls<-rownames(ctrl.summary.by.sample)[ctrl.summary.by.sample$Call.rate>=as.numeric(MaxMissRate)]
					ctrl.snps<-ctrl.snps[,colnames(ctrl.snps)%in%good.ctrls]
				}

				ctrl.snpStats<-invisible(new("SnpMatrix",t(ctrl.snps+1)))
				ctrl.summary<-col.summary(ctrl.snpStats)


				if(!is.null(MaxMissRate)) # then for these get good snps
				{
					good.qual.snps<-rownames(ctrl.summary[ctrl.summary$Call.rate>=as.numeric(MaxMissRate),]) 
				} else good.qual.snps<-rownames(ctrl.summary) 

				if(length(good.qual.snps)>0)
				{
					case.snps<-case.snps[rownames(case.snps)%in%good.qual.snps,]
					ctrl.snps<-ctrl.snps[rownames(ctrl.snps)%in%good.qual.snps,]

					if(HWEp>0)
					{
						## Hardy weinberg filtering in ctonrols.
						ctrl.snpStats<-invisible(new("SnpMatrix",t(ctrl.snps+1)))
						ctrl.summary<-col.summary(ctrl.snpStats)
						hwe.pvals<-pnorm(-abs(ctrl.summary$z.HWE))
						hwe.pvals[is.na(hwe.pvals)]<-1
						non.hwe.ctrl.snps<- rownames(ctrl.snps)[which(hwe.pvals <= HWEp) ]
						case.snps<-data.frame(case.snps[!rownames(case.snps)%in%non.hwe.ctrl.snps,]) 
						ctrl.snps<-data.frame(ctrl.snps[!rownames(ctrl.snps)%in%non.hwe.ctrl.snps,]) 
						message(paste(length(non.hwe.ctrl.snps),'SNPs removed because they are out of HWE in controls'))
					} else message('No HWE filter applied.')
					if(nrow(ctrl.snps) >=MinSNPs)
					{ 			
						ctrl.snpStats<-invisible(new("SnpMatrix",t(ctrl.snps+1)))
						ctrl.summary<-col.summary(ctrl.snpStats)
						results$MeanCallRateCtrls[gene]<-signif(mean(ctrl.summary$Call.rate,na.rm=T),2) 
					} else return
				}
			} else results$MeanCallRateCtrls[gene]<-0

			if(sum(as.matrix(case.snps),na.rm=T)>0)
			{
				case.snpStats<-invisible(new("SnpMatrix",t(case.snps+1)))
				case.summary<-col.summary(case.snpStats)
				results$MeanCallRateCases[gene]<-signif(mean(case.summary$Call.rate,na.rm=T),2)
			} else results$MeanCallRateCases[gene]<-0

			if(nrow(case.snps)>=MinSNPs)
			{ 
				maf.ctrl.set<- colnames(sample(ctrl.snps, ncol(ctrl.snps)/10)) # Separate 10% of ctrls for MAf filter
				maf.ctrl.set<-ctrl.snps
				maf.ctrls<- ctrl.snps[, colnames(ctrl.snps) %in% maf.ctrl.set ]
				maf.snp.cases<-apply(case.snps,1, function(x) signif(maf(as.numeric(unlist(table(unlist(x))))),2)  )
				maf.snp.ctrls<-apply(maf.ctrls,1, function(x) signif(maf(as.numeric(unlist(table(unlist(x))))),2)  )
				maf.snp.ctrls[is.na(maf.snp.ctrls)]<-0
				case.snps<-case.snps[rownames(case.snps) %in% names(which(maf.snp.ctrls<=MaxCtrlMAF)) ,]
				ctrl.snps<-ctrl.snps[rownames(ctrl.snps) %in% names(which(maf.snp.ctrls<=MaxCtrlMAF)),! colnames(ctrl.snps) %in% maf.ctrl.set]
			}
			
			nb.snps.in.gene2<-nrow(case.snps)
			print(paste('nb.snps.in.gene2=',nb.snps.in.gene2))
			if(nb.snps.in.gene2>=MinSNPs)
			{ 
				#case.snps<-final.snp.set[,which(current.pheno==1)]
				#ctrl.snps<-final.snp.set[,which(current.pheno==0)]

				ctrl.snps<-ctrl.snps[,colnames(ctrl.snps)%in%good.ctrls]
				final.snp.set<-as.matrix(cbind(case.snps,ctrl.snps))

				if(sum(as.matrix(ctrl.snps),na.rm=T)>0)
				{
					ctrl.snpStats<-invisible(new("SnpMatrix",t(ctrl.snps+1)))
					ctrl.summary<-col.summary(ctrl.snpStats)
					results$MeanCallRateCtrls[gene]<-signif(mean(ctrl.summary$Call.rate,na.rm=T),2) 
				} else results$MeanCallRateCtrls[gene]<-0
				if(sum(as.matrix(case.snps),na.rm=T)>0)
				{
					case.snpStats<-invisible(new("SnpMatrix",t(case.snps+1)))
					case.summary<-col.summary(case.snpStats)
					results$MeanCallRateCases[gene]<-signif(mean(case.summary$Call.rate,na.rm=T),2)
				} else results$MeanCallRateCases[gene]<-0

				results$nb.cases[gene]<-ncol(case.snps)
				results$nb.ctrls[gene]<-ncol(ctrl.snps)
		#		results$nb.cases[gene]<-length(which(!is.na(unlist(case.snps))) )
		#		results$nb.ctrls[gene]<-length(which(!is.na(unlist(ctrl.snps))) )

				results$nb.alleles.cases[gene]<-(length(grep(1,unlist(case.snps))))+ (length(grep(2,unlist(case.snps)))*2)
				results$nb.alleles.ctrls[gene]<-(length(grep(1,unlist(ctrl.snps))))+ (length(grep(2,unlist(ctrl.snps)))*2)

				##these counts are for gene total
				if(sum(as.matrix(case.snps),na.rm=T)>0)results$case.maf[gene]<-signif(maf(as.numeric(unlist(table(unlist(case.snps))))),2) 
				
				case.homs<-length(grep(2,case.snps))
				if(case.homs>0) results$nb.case.homs[gene]<-case.homs

				case.hets<-length(grep(1,case.snps))
				if(case.hets>0)results$nb.case.hets[gene]<-case.hets

				### now do ctrls	
				if(sum(as.matrix(ctrl.snps),na.rm=T)>0)results$ctrl.maf[gene]<-signif(maf(as.numeric(unlist(table(unlist(ctrl.snps))))),2) else results$ctrl.maf[gene]<-0

			#	if(results$ctrl.maf[gene]>0)
			#	{ 
					ctrl.homs<-length(grep(2,ctrl.snps))
					if(ctrl.homs>0)results$nb.ctrl.homs[gene]<-ctrl.homs
				
					ctrl.hets<-length(grep(1,ctrl.snps)) 
					if(ctrl.hets>0) results$nb.ctrl.hets[gene]<-ctrl.hets

					total.snps<-data.frame(cbind(case.snps,ctrl.snps))
					if(sum(as.matrix(total.snps),na.rm=T)>0)results$total.maf[gene]<-signif(maf(as.numeric(unlist(table(unlist(total.snps))))),2) 

					## these counts are for each snp in gene separately
					case.snp.hets<-apply(case.snps,1,function(x) length(grep(1,x)))
					case.snp.homs<-apply(case.snps,1,function(x) length(grep(2,x)))
					maf.snp.cases<-apply(case.snps,1, function(x) signif(maf(as.numeric(unlist(table(unlist(x))))),2)  )

					ctrl.snp.hets<-apply(ctrl.snps,1,function(x) length(grep(1,x)))
					ctrl.snp.homs<-apply(ctrl.snps,1,function(x) length(grep(2,x)))
					maf.snp.ctrls<-apply(ctrl.snps,1, function(x) signif(maf(as.numeric(unlist(table(unlist(x))))),2)  )

					current.pheno<-(c(rep(1,ncol(case.snps)),rep(0,ncol(ctrl.snps))))
					ancestry<-read.table(paste0(rootODir,'UCLex_samples_ancestry'),header=T)
					techPCs<-read.table(paste0(rootODir,'TechPCs.vect'),header=F)
					depthPCs<-read.table(paste0(rootODir,'DepthPCs.vect'),header=F)

					ancestry.pcs<-ancestry[match(colnames(final.snp.set),ancestry$V1),]
					techPCs<-techPCs[match(colnames(final.snp.set),techPCs$V1),]
					depthPCs<-depthPCs[match(colnames(final.snp.set),depthPCs$V1),]

					obj<-SKAT_Null_Model(current.pheno ~ 
						ancestry.pcs$V3+ancestry.pcs$V4#+ancestry.pcs$V5++ancestry.pcs$V6+ancestry.pcs$V7
						+techPCs$V3+techPCs$V4#+techPCs$V5+techPCs$V6+techPCs$V7
						#+depthPCs$V3+depthPCs$V4#+depthPCs$V5+depthPCs$V6+depthPCs$V7
						, out_type="D")
					if(sum(as.matrix(case.snps),na.rm=T)>0)
					{
						results$SKATO[gene] <- SKAT(t(as.matrix(final.snp.set)) , obj, missing_cutoff=.2, estimate_MAF=2,method="optimal.adj",impute.method="bestguess")$p.value
					}else results$SKATO[gene]<-1

					results$nb.snps[gene] <- nrow(case.snps)
					results$SNPs[gene]<-paste(rownames(case.snps),collapse=';')

					fisher.nb.cases<-length(which(!is.na(unlist(case.snps))) )
					fisher.nb.ctrls<-length(which(!is.na(unlist(ctrl.snps))) )
		       		mat<-matrix(c(fisher.nb.ctrls*2 - results$nb.alleles.ctrls[gene],
		       						results$nb.alleles.ctrls[gene],
		       						fisher.nb.cases*2 - results$nb.alleles.cases[gene],
		       						results$nb.alleles.cases[gene])
		       						, nrow = 2, ncol = 2)
					if (length(which(is.na(unlist(mat))))==0)
		       		{
		       			testy<-fisher.test(mat)
		       			results$FisherPvalue[gene]<-signif(testy$p.value,4) 
		       			results$OddsRatio[gene]<-signif(testy$estimate,4) 
		       		}


					snp.out<-data.frame( rownames(case.snps),case.snp.hets,case.snp.homs, maf.snp.cases,ctrl.snp.hets,ctrl.snp.homs, maf.snp.ctrls,results[gene,]) 
					write.table( snp.out, oFile, col.names=!file.exists(oFile),row.names=F,quote=F,sep='\t',append=T)

					GetCarriers<-function(snps)
					{
						car<- rownames( data.frame(unlist(apply(snps,1,function(x) which(x>0 )))))
						if(nrow(snps)==1)
						{
							carriers<-colnames(snps)[apply(snps,1,function(x) which(x>0 ))]
							variants<-str_extract(rownames(snps),"[0-9]{1,2}_[0-9]+_[A-Z]_[A-Z]")
							dat<-data.frame(cbind(variants,carriers))

						} else
						{
							carriers<-car
							carriers.clean<-gsub(carriers,pattern="[0-9]{1,2}_[0-9]+_[A-Z]_[A-Z]\\.",replacement="")
							variants<-str_extract(car,"[0-9]{1,2}_[0-9]+_[A-Z]_[A-Z]")
							dat<-data.frame(cbind(variants,carriers.clean))
						}
						
						if( ( identical(dat[,1],dat[,2])  | is.na(dat[,1])) && nrow(snps)>1)
						{
							index<-which(snps>0, arr.ind=TRUE)
							tt<-data.frame(matrix(nrow=nrow(index),ncol=2))
							tt[,1]<- rownames(snps)[index[1:nrow(index)]]
							tt[,2]<- colnames(snps)[index[nrow(index)+(1:nrow(index))]]

							carriers<-tt[,2]
							variants<-tt[,1]
							dat<-data.frame(cbind(variants,carriers))
						} 
						if( ( identical(dat[,1],dat[,2])  | is.na(dat[,1]) ) && nrow(snps)==1)
						{
							carriers<-car
							carriers.clean<-gsub(carriers,pattern="[0-9]{1,2}_[0-9]+_[A-Z]_[A-Z]\\.",replacement="")
							variants<-str_extract(car,"[0-9]{1,2}_[0-9]+_[A-Z]_[A-Z]")
							dat<-data.frame(cbind(variants,carriers.clean))
						} 
					return(dat)
					}
					case.calls<-FALSE
					ctrl.calls<-FALSE
					if(sum(as.matrix(case.snps),na.rm=T)>0)
					{
						case.dat<-GetCarriers(case.snps)
						case.dat<-data.frame(case.dat,uniq.genes[gene])
						write.table(case.dat, caseFile, col.names=FALSE,row.names=F,quote=F,sep='\t',append=T)
						case.calls<-TRUE
					}

					if(sum(as.matrix(ctrl.snps),na.rm=T)>0)
					{
						ctrl.dat<-GetCarriers(ctrl.snps)
						ctrl.dat<-data.frame(ctrl.dat,uniq.genes[gene])
						write.table(ctrl.dat,ctrlFile, col.names=FALSE,row.names=F,quote=F,sep='\t',append=T)
						ctrl.calls<-TRUE
					}
		
					if(!is.null(compoundHets))
					{

						GetCompoundHets<-function(snps,snp.dat,outFile)
						{
							compound.hets.names<-names(which(table(snp.dat[,grep("carriers",colnames(snp.dat))] )>1)) # who has more than one snp
							nb.hets<-0
							for(i in 1:length(compound.hets.names))
							{
								compound.snps<-t( snp.dat$variants[ grep(compound.hets.names[i],snp.dat[,grep("carriers",colnames(snp.dat))]) ] )  ## get names of snps seen in same individual
								compound.snp.calls<-paste(snps[rownames(snps)%in%compound.snps,colnames(snps)%in%compound.hets.names[i]],collapse=';')
								tt<-data.frame(uniq.genes[gene],compound.hets.names[i],paste(compound.snps,collapse=';'),t(compound.snp.calls) )
								if(length(tt)>0)
								{
							#		write.table(tt,temp,col.names=F,row.names=F,quote=F,sep='\t',append=T)
									write.table(tt,outFile,col.names=F,row.names=F,quote=F,sep='\t',append=T)
									nb.hets<-nb.hets+1
								}
							}
							return(nb.hets)
						} 

						caseTest<-FALSE
						ctrlTest<-FALSE

						if(case.calls && length(unique(case.dat[,grep("carriers",colnames(case.dat))]))<nrow(case.dat) )
						{
							case.compound.hets<-GetCompoundHets(snps=case.snps,snp.dat=case.dat,outFile=compoundFileCases)
							caseTest<-TRUE
						}
						if(ctrl.calls && length(unique(ctrl.dat[,grep("carriers",colnames(ctrl.dat))]))<nrow(ctrl.dat) )
						{
							ctrl.compound.hets<-GetCompoundHets(snps=ctrl.snps,snp.dat=ctrl.dat,outFile=compoundFileCtrls)
							ctrlTest<-TRUE
						}
						if(caseTest && ctrlTest)
						{
							nb.clean.cases<-ncol(case.snps)-case.compound.hets
							nb.clean.ctrls<-ncol(ctrl.snps)-ctrl.compound.hets

						fisher.nb.cases<-length(which(!is.na(unlist(case.snps))) )
						fisher.nb.ctrls<-length(which(!is.na(unlist(ctrl.snps))) )
			       		mat<-matrix(c(fisher.nb.ctrls*2 - results$nb.alleles.ctrls[gene],
			       						results$nb.alleles.ctrls[gene],
			       						fisher.nb.cases*2 - results$nb.alleles.cases[gene],
			       						results$nb.alleles.cases[gene])
			       						, nrow = 2, ncol = 2)
						results$CompoundHetPvalue[gene]<-fisher.test(mat,alternative='greater')$p.value
						}

						CheckHomozygosityRun<-function(compoundHetDat,snpDat,oFile)
						{
							homozyg.dir<- paste0(outputDirectory,'Homozyg_mapping') 
							if(!file.exists(homozyg.dir)) dir.create(homozyg.dir)
							pdf.out<-paste0(homozyg.dir,'Homozyg_mapping')

						}
					}
				print(results[gene,])
				if(qcPREP)
				{
					robj<-paste0(outputDirectory,'qc/test_setup.RData')
					message(paste('Saving workspace image to', robj))
					save(list=ls(environment()),file=robj)
				}	
			}
		}
	}
###################################
 		results <- results[order(results[,2]), ]
		results.out <- paste0(outputDirectory,'skat.csv')
		results<-results[order(results$SKATO),]
		qqplot.out <- paste0(outputDirectory,'skat_QQplot.png')

		write.table(results,results.out,col.names=T,row.names=F,quote=F,sep=',')
	
	#if(is.null(TargetGenes))
	#{
#		png(qqplot.out)
#		qq.chisq(-2*log(results$SKATO), df=2, x.max=30, main='SKAT',cex.main=0.7)	
#		dev.off() 
#	}

print("Finished testing.")
} #doSKAT



#### now run function.
doSKAT(case.list=case.list,control.list=control.list,outputDirectory=outputDirectory,TargetGenes=TargetGenes,
	min.depth=min.depth,minCadd=minCadd,maxExac=maxExac,MinSNPs=MinSNPs,compoundHets=compoundHets,MaxMissRate=MaxMissRate,
	chrom=chrom,PlotPCA=PlotPCA,homozyg.mapping=homozyg.mapping,HWEp=HWEp,
	MaxCtrlMAF=MaxCtrlMAF,SavePrep=SavePrep)