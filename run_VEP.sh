#!/bin/bash

####INPUT FILE
#vcfin=/home/zchads1/cluster/UCL-exomes_v2/variants/Levine_${chr}_single.vcf
####SET CHROMOSOME
#ch1r=1
#out dir

set -u
set -e
set -x

# prints to stderr in red
function error() { >&2 echo -e "\033[31m$*\033[0m"; }
function stop() { error "$*"; exit 1; }
try() { "$@" || stop "cannot $*"; }



function usage() {
    echo --input
    echo --chr
    echo --vcfout
    echo --reference
    echo --coding_only
    exit 1
}


coding_only=no
custom=

##
until [ -z "$1" ]
do
    # use a case statement to test vars. we always test $1 and shift at the end of the for block.
    case $1 in
    --input)
        shift
        input=$1;;
    --chr)
        shift
        chr=$1;;
    --vcfout)
        shift
        vcfout=$1;;
    --reference)
        shift
        reference=$1;;
    --coding_only)
        shift
        coding_only=$1;;
    --custom)
        shift
        custom=$1;;
    -* )
        echo "Unrecognized option: $1"
        usage
        exit 1;;
    esac
    shift
    if [ "$#" = "0" ]; then break; fi
done 

if [[ "$coding_only" == yes ]]
then
    coding_only="--coding_only"
else
    coding_only=
fi

####CONFIGURE SOFTWARE SHORTCUTS AND PATHS
ensembl=/cluster/project8/vyp/AdamLevine/software/ensembl/
#VEP=${ensembl}/src/ensembl-tools/scripts/variant_effect_predictor/variant_effect_predictor.pl
VEP=/cluster/project8/vyp/Software/ensembl-tools-release-82/scripts/variant_effect_predictor/variant_effect_predictor.pl
dir_cache=${ensembl}/cache/
perl=/share/apps/perl-5.14.2/bin/perl
PERL5LIB=${PERL5LIB}:${ensembl}/src/bioperl-1.6.1
PERL5LIB=${PERL5LIB}:${ensembl}/src/ensembl/modules
PERL5LIB=${PERL5LIB}:${ensembl}/src/ensembl-compara/modules
PERL5LIB=${PERL5LIB}:${ensembl}/src/ensembl-variation/modules
PERL5LIB=${PERL5LIB}:${ensembl}/src/ensembl-funcgen/modules
PERL5LIB=${PERL5LIB}:${ensembl}/Plugins
export PERL5LIB
export PATH=$PATH:/cluster/project8/vyp/vincent/Software/tabix-0.2.5/
condel_config=${ensembl}/Plugins/config/Condel/config

#fasta=/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta
#fasta=/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
#file:///scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
cat $input | grep '^##reference=' | cut -f2 -d'='

#annotations_dir=/cluster/project8/vyp/AdamLevine/annotations
annotations_dir=/cluster/project8/IBDAJE/VEP_custom_annotations/${reference}

# Custom annotations
####CADD http://cadd.gs.washington.edu/home
#This needs to be updated with the latest scores [ACTION!]
#They also now provide a script which is worth exploring
#custom_annotation="--custom ${annotations_dir}/CADD/chr${chr}.vcf.gz,CADD,vcf,exact"
custom_annotation="--custom ${annotations_dir}/CADD/chr${chr}.vcf.gz,CADD,vcf,exact"

####ExAC
for pop in AFR AMR Adj EAS FIN NFE OTH SAS
do
    shortname=EXAC_${pop}
    custom_annotation="${custom_annotation} --custom ${annotations_dir}/ExAC/0.3/chr${chr}_${pop}.vcf.gz,${shortname},vcf,exact"
done

####1kg
for pop in EUR AFR AMR ASN
do
    shortname=1KG_${pop}
    custom_annotation="${custom_annotation} --custom ${annotations_dir}/1kg/chr${chr}_${pop}.vcf.gz,${shortname},vcf,exact"
done

####ESP frequency annotations
for pop in EA AA
do
    shortname=ESP_${pop}
    custom_annotation="${custom_annotation} --custom ${annotations_dir}/esp/chr${chr}_${pop}.vcf.gz,${shortname},vcf,exact"
done

#### 
function UCLEX() {
    shortname=UCLEX
    custom_annotation="${custom_annotation} --custom ${annotations_dir}/UCLex/chr${chr}.vcf.gz,${shortname},vcf,exact"
}

###
function AJcontrols() {
    shortname=AJcontrols
    custom_annotation="${custom_annotation} --custom ${annotations_dir}/AJcontrols/chr${chr}.vcf.gz,${shortname},vcf,exact"
}

###
function AJcases() {
    shortname=AJcases
    custom_annotation="${custom_annotation} --custom ${annotations_dir}/AJcases/chr${chr}.vcf.gz,${shortname},vcf,exact"
}

###
function BroadAJcontrols() {
    shortname=BroadAJcontrols
    custom_annotation="${custom_annotation} --custom ${annotations_dir}/BroadAJcontrols/chr${chr}.vcf.gz,${shortname},vcf,exact"
}

###
function ImmunoBase() {
    for disease in CRO IBD UC
    do
        shortname=ImmunoBase_${disease}
        custom_annotation="${custom_annotation} --custom ${annotations_dir}/ImmunoBase/ImmunoBase_${disease}.bed.gz,${shortname},bed,overlap"
    done
}

###
function SZ_Curtis() {
    shortname=SZ_Curtis
    custom_annotation="${custom_annotation} --custom ${annotations_dir}/SZ_Curtis/chr${chr}.vcf.gz,${shortname},vcf,exact"
}



if [[ "$custom" != "" ]]
then
    for x in `echo $custom | tr ',' ' '`
    do
        $x
    done
fi

human_reference=/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/

if [[ "$reference" == "hg38_noAlt" ]]
then
    fasta=$human_reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
    chrPrefix='chr'
    assembly=GRCh38
    port=
elif [[ "$reference" == "1kg" ]]
then
    fasta=$human_reference/human_g1k_v37.fasta
    assembly=GRCh37
    chrPrefix=''
    port='--port 3337'
elif [[ "$reference" == "hg19" ]]
then
    fasta=$human_reference/hg19_UCSC.fa
    chrPrefix='chr'
    port=
    custom_annotation=
else
    stop Unsupported reference $reference
fi

# builtin VEP MAF
# from http://www.ensembl.org/info/docs/tools/vep/script/vep_options.html#opt_gmaf
# --gmaf
# Add the global minor allele frequency (MAF) from 1000 Genomes Phase 1 data for any existing variant to the output. Not used by default
# --maf_1kg
# Add allele frequency from continental populations (AFR,AMR,ASN,EUR) of 1000 Genomes Phase 1 to the output.
# Note the reported allele(s) and frequencies are for the non-reference allele from the original data, not necessarily the alternate allele from user input.
# Must be used with --cache Not used by default
# --maf_esp
# Include allele frequency from NHLBI-ESP populations. Note the reported allele(s) and frequencies are for the non-reference allele from the originial data, not necessarily the alternate allele from user input. Must be used with --cache Not used by default
# --old_maf
# For --maf_1kg and --maf_esp report only the frequency (no allele) and convert this frequency so it is always a minor frequency, i.e. < 0.5

maf="--maf_esp --gmaf --maf_1kg"
#fields="--fields SYMBOL,CLIN_SIG,AA_MAF,EA_MAF,SIFT,PolyPhen,CAROL,Condel"
fields=""

## Extra fields
# SYMBOL - the gene symbol
# SYMBOL_SOURCE - the source of the gene symbol
# STRAND - the DNA strand (1 or -1) on which the transcript/feature lies
# ENSP - the Ensembl protein identifier of the affected transcript
# SWISSPROT - UniProtKB/Swiss-Prot identifier of protein product
# TREMBL - UniProtKB/TrEMBL identifier of protein product
# UNIPARC - UniParc identifier of protein product
# HGVSc - the HGVS coding sequence name
# HGVSp - the HGVS protein sequence name
# SIFT - the SIFT prediction and/or score, with both given as prediction(score)
# PolyPhen - the PolyPhen prediction and/or score
# MOTIF_NAME - the source and identifier of a transcription factor binding profile aligned at this position
# MOTIF_POS - The relative position of the variation in the aligned TFBP
# HIGH_INF_POS - a flag indicating if the variant falls in a high information position of a transcription factor binding profile (TFBP)
# MOTIF_SCORE_CHANGE - The difference in motif score of the reference and variant sequences for the TFBP
# CELL_TYPE - List of cell types and classifications for regulatory feature
# CANONICAL - a flag indicating if the transcript is denoted as the canonical transcript for this gene
# CCDS - the CCDS identifer for this transcript, where applicable
# INTRON - the intron number (out of total number)
# EXON - the exon number (out of total number)
# DOMAINS - the source and identifer of any overlapping protein domains
# DISTANCE - Shortest distance from variant to transcript
# IND - individual name
# ZYG - zygosity of individual genotype at this locus
# SV - IDs of overlapping structural variants
# FREQS - Frequencies of overlapping variants used in filtering
# GMAF - Minor allele and frequency of existing variation in 1000 Genomes Phase 1
# AFR_MAF - Minor allele and frequency of existing variation in 1000 Genomes Phase 1 combined African population
# AMR_MAF - Minor allele and frequency of existing variation in 1000 Genomes Phase 1 combined American population
# ASN_MAF - Minor allele and frequency of existing variation in 1000 Genomes Phase 1 combined Asian population
# EUR_MAF - Minor allele and frequency of existing variation in 1000 Genomes Phase 1 combined European population
# AA_MAF - Minor allele and frequency of existing variant in NHLBI-ESP African American population
# EA_MAF - Minor allele and frequency of existing variant in NHLBI-ESP European American population
# CLIN_SIG - Clinical significance of variant from dbSNP
# BIOTYPE - Biotype of transcript or regulatory feature
# TSL - Transcript support level
# PUBMED - Pubmed ID(s) of publications that cite existing variant
# SOMATIC - Somatic status of existing variation(s)
# ALLELE_NUM - Allele number from input; 0 is reference, 1 is first alternate etc
# PICK - indicates if this block of consequence data was picked by --flag_pick or --flag_pick_allele
 
#output='--pick'
output='--vcf'
#plugins="--plugin Condel,${condel_config},b --plugin Carol --plugin LoF,human_ancestor_fa:/scratch2/vyp-scratch2/reference_datasets/loftee/human_ancestor.fa.rz,filter_position:0.05"
#plugins="--plugin Condel,${condel_config},b --plugin Carol --plugin CADD,${annotations_dir}/CADD/chr${chr}.vcf.gz --plugin LoF,human_ancestor_fa:/scratch2/vyp-scratch2/reference_datasets/loftee/human_ancestor.fa.rz,filter_position:0.05"
#plugins="--plugin Condel,${condel_config},b --plugin Carol --plugin CADD,${annotations_dir}/CADD/chr${chr}.vcf.gz --plugin GO --plugin ExAC,${annotations_dir}/ExAC/0.3/chr${chr}.vcf.gz"
plugins="--plugin Condel,${condel_config},b --plugin Carol --plugin CADD,${annotations_dir}/CADD/chr${chr}.vcf.gz --plugin LoF,human_ancestor_fa:/scratch2/vyp-scratch2/reference_datasets/loftee/human_ancestor.fa.rz,filter_position:0.05"

#$perl $VEP $port --verbose --ASSEMBLY $assembly --fasta $fasta --cache --dir_cache $dir_cache --input_file $vcfin --format vcf --sift b --polyphen b --symbol  --canonical --check_existing --check_alleles  --no_progress --output_file $vcfout  --force_overwrite $output --fork 2 $maf $fields $custom_annotation $plugins $coding_only --offline
$perl $VEP $port --verbose --ASSEMBLY $assembly --fasta $fasta --cache --dir_cache $dir_cache --input_file $input --sift b --polyphen b --symbol  --canonical --check_existing --check_alleles  --no_progress --output_file $vcfout  --force_overwrite $output --fork 2 $maf $fields $custom_annotation $plugins $coding_only --offline --allele_number --numbers --gene_phenotype
#$perl $VEP $port --verbose --ASSEMBLY $assembly --fasta $fasta --input_file $vcfin --cache --dir_cache $dir_cache --format vcf --sift b --polyphen b --symbol  --canonical --check_existing --check_alleles  --no_progress --output_file $vcfout  --force_overwrite $output --offline

