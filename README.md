# Phenopolis



### Installation

Phenopolis requires a running mongo database and a running Phenotips server.

Clone the repository.

```
git clone git@github.com:pontikos/phenopolis.git
```

Download Phenotips.
```
https://phenotips.org/Download
```

Download the Exomiser stand alone.
```
https://1drv.ms/u/s!AnAWImk12qlQjo8FkToxDH5lxOS7Xw
```

Install latest version of mongo.
```
https://www.mongodb.com/download-center#community
```

### Creating database

First make sure mongoDB is running:
```
DBPATH=
mongod --dbpath $DBPATH --port 27017 --smallfiles
```

#### Importing data from JSON

The variants found in the VCF files are processed with VEP and the output is written to JSON.
This is further piped into another python script which adds further annotation and formatting and writes output to JSON.
The JSON is then imported with mongoimport.

```
git clone https://github.com/UCLGeneticsInstitute/DNASeq_pipeline
```
```
DNASeq_pipeline/annotation/postprocess_VEP_json.py
```
```
function VEP_mongo() {
    memo=30
    mkdir -p ${output}_VEP
    #VEP_output="--tab --output_file ${output}_VEP/VEP_${chr}.txt"
    ensembl=/cluster/project8/vyp/AdamLevine/software/ensembl/
    VEP_DIR=/cluster/project8/vyp/Software/ensembl-tools-release-82/scripts/
    DIR_CACHE=/SAN/vyplab/NCMD_raw/VEP/cache/
    DIR_PLUGINS=${DIR_CACHE}/Plugins
    for chr in `seq 1 22` X
    do
        #output_lines=`zcat ${output}_VEP/VEP_${chr}.json.gz | wc -l`
        #input_lines=`tail -n+2 ${output}_VEP/chr${chr}_for_VEP.vcf | wc -l`
        #echo ${output}_VEP/chr${chr}_for_VEP.vcf input lines: $input_lines
        #echo ${output}_VEP/VEP_${chr}.json.gz output lines: $output_lines
        #if [[ $lines -gt 0 ]]; then echo ${output}_VEP/VEP_${chr}.json.gz $lines gt than 0, skipping; continue; fi
        echo "
############### VEP_mongo chr${chr}
# split single lines
zcat ${output}_chr${chr}_for_annovar.vcf.gz | /share/apps/python/bin/python ${baseFolder}/annotation/multiallele_to_single_gvcf.py --headers CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO > ${output}_VEP/chr${chr}_for_VEP.vcf
####CONFIGURE SOFTWARE SHORTCUTS AND PATHS
reference=1kg
ensembl=/cluster/project8/vyp/AdamLevine/software/ensembl/
export PERL5LIB=${PERL5LIB}:${ensembl}/src/bioperl-1.6.1::${ensembl}/src/ensembl/modules:${ensembl}/src/ensembl-compara/modules:${ensembl}/src/ensembl-variation/modules:${ensembl}/src/ensembl-funcgen/modules:${ensembl}/Plugins
export PATH=$PATH:/cluster/project8/vyp/vincent/Software/tabix-0.2.5/
# RUN VEP
/share/apps/perl/bin/perl ${VEP_DIR}/variant_effect_predictor/variant_effect_predictor.pl --offline \
--ASSEMBLY GRCh37 --fasta /SAN/vyplab/UKIRDC/reference/human_g1k_v37.fasta \
--cache --dir_cache ${DIR_CACHE} \
--dir_plugins ${DIR_PLUGINS} \
--sift b --polyphen b --symbol --canonical --check_existing --check_alleles  \
--fork 4 --maf_esp --gmaf --maf_1kg --maf_exac \
--no_progress --quiet \
--custom /cluster/project9/IBDAJE/VEP_custom_annotations/1kg/1kg/chr${chr}_EUR.vcf.gz,1KG_EUR,vcf,exact \
--custom /cluster/project9/IBDAJE/VEP_custom_annotations/1kg/1kg/chr${chr}_AFR.vcf.gz,1KG_AFR,vcf,exact \
--custom /cluster/project9/IBDAJE/VEP_custom_annotations/1kg/1kg/chr${chr}_AMR.vcf.gz,1KG_AMR,vcf,exact \
--custom /cluster/project9/IBDAJE/VEP_custom_annotations/1kg/1kg/chr${chr}_ASN.vcf.gz,1KG_ASN,vcf,exact \
--custom /cluster/project9/IBDAJE/VEP_custom_annotations/1kg/esp/chr${chr}_EA.vcf.gz,ESP_EA,vcf,exact \
--custom /cluster/project9/IBDAJE/VEP_custom_annotations/1kg/esp/chr${chr}_AA.vcf.gz,ESP_AA,vcf,exact \
--custom /cluster/project9/IBDAJE/VEP_custom_annotations/1kg/CADD/chr${chr}.vcf.gz,CADD,vcf,exact \
--custom /cluster/scratch3/vyp-scratch2/reference_datasets/Kaviar/Kaviar-160204-Public/hg19/VEP_annotation.vcf.gz,Kaviar,vcf,exact \
--plugin Condel,${DIR_PLUGINS}/config/Condel/config,b \
--plugin Carol \
--plugin CADD,/cluster/project9/IBDAJE/VEP_custom_annotations/1kg/CADD/chr${chr}.vcf.gz \
--no_stats \
--hgvs \
--pubmed \
--plugin HGVSshift \
--plugin SameCodon \
--input_file ${output}_VEP/chr${chr}_for_VEP.vcf \
--json \
--output_file STDOUT | python ${baseFolder}/annotation/postprocess_VEP_json.py | grep '^JSON:' | sed 's/^JSON://' > ${output}_VEP/VEP_chr${chr}.json
/share/apps/genomics/htslib-1.1/bin/bgzip -f -c ${output}_VEP/chr${chr}_for_VEP.vcf > ${output}_VEP/chr${chr}_for_VEP.vcf.gz
/share/apps/genomics/htslib-1.1/bin/tabix -f -p vcf ${output}_VEP/chr${chr}_for_VEP.vcf.gz
rm ${output}_VEP/chr${chr}_for_VEP.vcf
" >> ${scripts_folder}/subscript_chr${chr}.sh
#--output_file STDOUT | python ${baseFolder}/annotation/postprocess_VEP_json.py | grep '^JSON:' | sed 's/^JSON://' | mongoimport --db uclex --collection variants --host phenotips
  done
}
```

Load individual for individual page:
```
 python views/load_individual.py --individual $ID --auth Admin:`cat ~/.pass/Admin`
```

#### Running pubmedbatch

The pubmedbatch, written by Jing, scores genes based on their pubmed relevance.

#### Running phenogenon

The phenogenon, written by Jing, does an enrichment test per gene and HPO term.

## Running server

Run Phenotips:
```
wget https://nexus.phenotips.org/nexus/content/repositories/releases/org/phenotips/phenotips-standalone/1.3-milestone-2/phenotips-standalone-1.3-milestone-2.zip
unzip phenotips-standalone-1.3-milestone-2.zip
cd phenotips-standalone-1.3-milestone-2
bash start.sh
```

Run Exomiser standalone:
```
EXOMISER_DATA=
cd $EXOMISER_DATA
wget ftp://ftp.sanger.ac.uk/pub/resources/software/exomiser/downloads/exomiser/exomiser-cli-7.2.1-data.zip
unzip exomiser-cli-7.2.1-data.zip
java -jar exomiser-rest-prioritiser-7.3.0-SNAPSHOT.jar --exomiser.data-directory=$EXOMISER_DATA
```

Run Phenopolis:
```
cd phenopolis
python run_server.py
```


### Acknowledgment

This code was originally forked from the ExAC browser but has since diverged considerably.


