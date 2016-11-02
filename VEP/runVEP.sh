source env.sh

ensembl=
VEP_DIR=
DIR_CACHE=
DIR_PLUGINS=${DIR_CACHE}/Plugins

input=

VEP=${VEP_DIR}/variant_effect_predictor/variant_effect_predictor.pl 

reference=1kg

ensembl=/cluster/project8/vyp/AdamLevine/software/ensembl/
reference=/SAN/vyplab/UKIRDC/reference/human_g1k_v37.fasta 

export PERL5LIB=${PERL5LIB}:${ensembl}/src/bioperl-1.6.1::${ensembl}/src/ensembl/modules:${ensembl}/src/ensembl-compara/modules:${ensembl}/src/ensembl-variation/modules:${ensembl}/src/ensembl-funcgen/modules:${ensembl}/Plugins

    --custom /cluster/project9/IBDAJE/VEP_custom_annotations/1kg/1kg/chr${chr}_EUR.vcf.gz,1KG_EUR,vcf,exact \
    --custom /cluster/project9/IBDAJE/VEP_custom_annotations/1kg/1kg/chr${chr}_AFR.vcf.gz,1KG_AFR,vcf,exact \
    --custom /cluster/project9/IBDAJE/VEP_custom_annotations/1kg/1kg/chr${chr}_AMR.vcf.gz,1KG_AMR,vcf,exact \
    --custom /cluster/project9/IBDAJE/VEP_custom_annotations/1kg/1kg/chr${chr}_ASN.vcf.gz,1KG_ASN,vcf,exact \
    --custom /cluster/project9/IBDAJE/VEP_custom_annotations/1kg/esp/chr${chr}_EA.vcf.gz,ESP_EA,vcf,exact \
    --custom /cluster/project9/IBDAJE/VEP_custom_annotations/1kg/esp/chr${chr}_AA.vcf.gz,ESP_AA,vcf,exact \
    --custom /cluster/project9/IBDAJE/VEP_custom_annotations/1kg/CADD/chr${chr}.vcf.gz,CADD,vcf,exact \
    --custom /cluster/scratch3/vyp-scratch2/reference_datasets/Kaviar/Kaviar-160204-Public/hg19/VEP_annotation.vcf.gz,Kaviar,vcf,exact \


############### VEP_mongo chr${chr}
# split single lines
zcat $input | python multiallele_to_single_gvcf.py --headers CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO > chr${chr}_for_VEP.vcf
####CONFIGURE SOFTWARE SHORTCUTS AND PATHS
# RUN VEP
/share/apps/perl/bin/perl $VEP --offline \
    --ASSEMBLY GRCh37 --fasta $reference \
    --cache --dir_cache ${DIR_CACHE} \
    --dir_plugins ${DIR_PLUGINS} \
    --sift b --polyphen b --symbol --canonical --check_existing --check_alleles  \
    --fork 4 --maf_esp --gmaf --maf_1kg --maf_exac \
    --no_progress --quiet \
    --plugin Condel,${DIR_PLUGINS}/config/Condel/config,b \
    --plugin Carol \
    --plugin CADD,/cluster/project9/IBDAJE/VEP_custom_annotations/1kg/CADD/chr${chr}.vcf.gz \
    --no_stats \
    --hgvs \
    --pubmed \
    --plugin HGVSshift \
    --plugin SameCodon \
    --input_file chr${chr}_for_VEP.vcf \
    --json \
    --output_file STDOUT | python  postprocess_VEP_json.py | grep '^JSON:' | sed 's/^JSON://' > VEP_chr${chr}.json

bgzip -f -c chr${chr}_for_VEP.vcf > chr${chr}_for_VEP.vcf.gz
tabix -f -p vcf chr${chr}_for_VEP.vcf.gz
rm chr${chr}_for_VEP.vcf
