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


# split single lines
zcat $input | python multiallele_to_single_gvcf.py --headers CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO > chr${chr}_for_VEP.vcf

# RUN VEP
/share/apps/perl/bin/perl $VEP --offline \
    --input_file chr${chr}_for_VEP.vcf \
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
    --json \
    --output_file STDOUT | python  postprocess_VEP_json.py | grep '^JSON:' | sed 's/^JSON://' > VEP_chr${chr}.json

