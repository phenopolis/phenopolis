
# Mapping of variant properties from VEP JSON output to Pheno4J

| |  Value | Explanation | Pheno4J mapping |
| ------------- | ------------- | ----- | ----- |
|allele_string| A/G| | GeneticVariant.allele_string |
|start|38212762| | GeneticVariant.start |
|end|38212762| | GeneticVariant.end |
|seq_region_name| 22| chromosome | GeneticVariant.chrom |
|most_severe_consequence|intron_variant| worst consequence on a transcript | GeneticVariant.most_severe_consequence |
|strand|1| posi/nega tive strand | GeneticVariant.strand |


| custom_annotations |    Type |  Explanation | Pheno4J mapping |
| ------------- | ------------- | ----- | ----- |
| input[0].fields.name | "22-38212762-A-G" | | GeneticVariant.variantId |
| input[0].fields.AC | int | internal allele count | GeneticVariant.AC|
| input[0].fields.AF | 0<float<1 | internal allele frequency | GeneticVariant.allele_freq |
| input[0].fields.AN | int | internal allele number | GeneticVariant.AN | 
| input[0].fields.BaseQRankSum | 0<float<1 |
| input[0].fields.ClippingRankSum | 0<float<1 |
| input[0].fields.DP | int | | |
| input[0].fields.ExcessHet | int | Excess heterozygosity | GeneticVariant.ExcessHet | 
| input[0].fields.FS | float | Fisher strand | GeneticVariant.FS |
| input[0].fields.InbreedingCoeff  | float |  quality | GeneticVariant.InbreedingCoeff |
| input[0].fields.MLEAC | int | | GeneticVariant.MLEAC |
| input[0].fields.MLEAF | float | | GeneticVariant.MLEAF |
| input[0].fields.MQ | float | mapping quality | GeneticVariant.MQ |
| input[0].fields.MQ0 | float | | |
| input[0].fields.MQRankSum | float | | GeneticVariant.MQRankSum |
| input[0].fields.QD | float | | |
| input[0].fields.ReadPosRankSum | float | | GeneticVariant.ReadPosRankSum | 
| input[0].fields.SOR | float | | |
| input[0].fields.VQSLOD | float | variant quality score log odds | GeneticVariant.VQSLOD |
| input[0].fields.culprit | string | Culprit | GeneticVariant.Culprit |
| gnomad_(exomes/genomes)[0].fields.AC_AFR|int| allele count in africans | GeneticVariant.gnomad_(exomes/genomes)_AC_AFR |
| gnomad_(exomes/genomes)[0].fields.AC_AMR|int| allele count in americans | GeneticVariant.gnomad_(exomes/genomes)_AC_AMR |
| gnomad_(exomes/genomes)[0].fields.AC_ASJ|int| |  GeneticVariant.gnomad_(exomes/genomes)_AC_ASJ |
| gnomad_(exomes/genomes)[0].fields.AC_EAS|int| | |
| gnomad_(exomes/genomes)[0].fields.AC_FIN|int| | |
| gnomad_(exomes/genomes)[0].fields.AC_Female|int| | |
| gnomad_(exomes/genomes)[0].fields.AC_Male|int| | |
| gnomad_(exomes/genomes)[0].fields.AC_NFE|int| | |
| gnomad_(exomes/genomes)[0].fields.AC_OTH|int| | |
| gnomad_(exomes/genomes)[0].fields.AC_raw|int| raw allele count | GeneticVariant.gnomad_(exomes/genomes)_AC_raw |
| gnomad_(exomes/genomes)[0].fields.AF_AFR|0<float<1| | |
| gnomad_(exomes/genomes)[0].fields.AF_AMR|0<float<1| | |
| gnomad_(exomes/genomes)[0].fields.AF_ASJ|0<float<1| | |
| gnomad_(exomes/genomes)[0].fields.AF_EAS|0<float<1| | |
| gnomad_(exomes/genomes)[0].fields.AF_FIN|0<float<1| | |
| gnomad_(exomes/genomes)[0].fields.AF_Female|0<float<1|
| gnomad_(exomes/genomes)[0].fields.AF_Male|0<float<1|
| gnomad_(exomes/genomes)[0].fields.AF_NFE|0<float<1 | allele frequency in non-finish europeans | GeneticVariant.gnomad_(exomes/genomes)_AF_NFE |
| gnomad_(exomes/genomes)[0].fields.AF_OTH|0<float<1| |  GeneticVariant.gnomad_(exomes/genomes)_AF_OTH |
| gnomad_(exomes/genomes)[0].fields.AF_raw|0<float<1| | GeneticVariant.gnomad_(exomes/genomes)_AF_raw |
| gnomad_(exomes/genomes)[0].fields.AN_AFR|int| |  GeneticVariant.gnomad_(exomes/genomes)_AN_AFR |
| gnomad_(exomes/genomes)[0].fields.AN_AMR|int| | |
| gnomad_(exomes/genomes)[0].fields.AN_ASJ|int| | |
| gnomad_(exomes/genomes)[0].fields.AN_EAS|int| | |
| gnomad_(exomes/genomes)[0].fields.AN_FIN|int| | |
| gnomad_(exomes/genomes)[0].fields.AN_Female|int| | |
| gnomad_(exomes/genomes)[0].fields.AN_Male|int| | |
| gnomad_(exomes/genomes)[0].fields.AN_NFE|int| | |
| gnomad_(exomes/genomes)[0].fields.AN_OTH|int| | |
| gnomad_(exomes/genomes)[0].fields.AN_raw|int| | |
| gnomad_(exomes/genomes)[0].fields.Hom_AFR|int| | |
| gnomad_(exomes/genomes)[0].fields.Hom_AMR|int| | |
| gnomad_(exomes/genomes)[0].fields.Hom_ASJ|int| | |
| gnomad_(exomes/genomes)[0].fields.Hom_EAS|int| | |
| gnomad_(exomes/genomes)[0].fields.Hom_FIN|int| | |
| gnomad_(exomes/genomes)[0].fields.Hom_Female|int| | |
| gnomad_(exomes/genomes)[0].fields.Hom_Male|int| | |
| gnomad_(exomes/genomes)[0].fields.Hom_NFE|int| | |
| gnomad_(exomes/genomes)[0].fields.Hom_OTH|int| | |
| gnomad_(exomes/genomes)[0].fields.Hom_raw|int| | |
| gnomad_(exomes/genomes)[0].fields.Hom|int| | |
| kaviar[0].fields.AF | 0<float<1 | kaviar alelle frequency | GeneticVariant.kaviar_AF |
| kaviar[0].fields.AC | int | kaviar allele count | GeneticVariant.kaviar_AC |
| kaviar[0].fields.AN | int | kaviar allele number | GeneticVariant.kaviar_AN |
| dbsnp[0].name | string | dbsnp name | GeneticVariant.dbsnp_name | 


| transcript_consequences[*] | Value |  Explanation | Pheno4J mapping |
| ------------- | ------------- | ----- | ----- |
| exac_af|0<float<1| ExAC allele frequency| GeneticVariant.exac_AF |
| exac_af_adj|0<float<1| ExAC allele frequency in adjusted | GeneticVariant.exac_AF_ADJ |
| exac_af_afr|0<float<1| ExAC allele frequency in africans | GeneticVariant.exac_AF_AFR |
| exac_af_amr|0<float<1| ExAC allele frequency in americans | GeneticVariant.exac_AF_AMR |
| exac_af_consanguineous|0<float<1 | ExAC allele frequency in consanguineous | GeneticVariant.exac_AF_CONSANGUINEOUS |
| exac_af_eas|0<float<1| ExAC allele frequency in east asians | GeneticVariant.exac_AF_EAS |
| exac_af_female|0<float<1| ExAC allele frequency in females | GeneticVariant.exac_AF_FEMALE |
| exac_af_fin|0<float<1| ExAC allele frequency in finnish | GeneticVariant.exac_AF_FIN |
| exac_af_male|0<float<1| ExAC allele frequency in males | GeneticVariant.exac_AF_MALE |
| exac_af_nfe|0<float<1| ExAC allele frequency in non-finnish europeans | GeneticVariant.exac_AF_NFE |
| exac_af_oth|0<float<1| ExAC allele frequency in others | GeneticVariant.exac_AF_OTH |
| exac_af_popmax|0<float<1| Max ExAC allele frequency | GeneticVariant.exac_AF_POPMAX |
| exac_af_sas|0<float<1| ExAC allele frequency in south asians | GeneticVariant.exac_AF_SAS |
| cadd_phred | float | Phred CADD score | GeneticVariant.cadd_phred |
| cadd_raw | float | Raw CADD score | GeneticVariant.cadd_raw |
| hgvsc | "ENST00000323205.6:c.*7+30A>G" | cDNA mutation in HGVS format | TranscriptVariant.hgvsc |
| hgvsp |  | protein mutation in HGVS format | TranscriptVariant.hgvsp |
| impact | "MODIFIER" |  | TranscriptVariant.impact |
| consequence_terms | list of strings | consequence of mutation | TranscriptVariant.consequence_terms | 
| transcript_id | "ENST00000323205" | ensembl id of transcript | TranscriptVariant.trancript_id |
| intron/exon | "9/9" | which intron/exon is mutation in | TranscriptVariant.intron/exon |
| gene_id | "ENSG00000100116" | ensembl gene id | TranscriptVariant.gene_id, Gene.gene_id |
| gene_symbol | "GCAT" | | TranscriptVariant.gene_name, Gene.gene_name |
| hgnc_id | int | | Gene.hgnc_id |
| transcript_id | "ENST00000323205" | ensembl id of transcript | Transcript.trancript_id |
| canonical | bool | canonical transcript | Transcript.canonical |
| strand | -1 or 1 | positive/negative strand | Transcript.strand |
| protein_id | "ENSP00000371110" | ensembl protein id | Transcript.ensembl_protein_id |
| swissprot | ["O75600"] | swissprot protein id | Transcript.swissprot_protein_id |


## Variant quality

### FILTER
Values:
* PASS
* FAIL
* VQSR

### FS: Fisher Strand
Whether there is a bias in the fisher strand.

### BaseQRankSum

### MQ
Mapping quality

### QD
Quality depth

### ClippingRankSum

### MQRankSum

## Variant Frequency

### Kaviar

The [Kaviar db](db.systemsbiology.net/kaviar/)

Released: February 29, 2016 (version 160204-Public)
Kaviar contains 162 million SNV sites (including 25M not in dbSNP) and incorporates data from 35 projects encompassing 77,781 individuals (13.2K whole genome, 64.6K exome).
* Kaviar also contains 50 million short indels and substitutions from a subset of the data sources.
* Kaviar excludes cancer genomes but includes some data from cell lines and individuals affected by disease.
* An effort was made to exclude data from related individuals.

### Gnomad

The Genome Aggregation Database (gnomAD) is a resource developed by an international coalition of investigators, with the goal of aggregating and harmonizing both exome and genome sequencing data from a wide variety of large-scale sequencing projects, and making summary data available for the wider scientific community.

The data set provided on this website spans 123,136 exome sequences and 15,496 whole-genome sequences from unrelated individuals sequenced as part of various disease-specific and population genetic studies.


### ExAC

## Variant function impact

### CADD
The CADD is not transcript specific but depends on the DNA

### Transcript Consequences

The trancript consequences field contains the functional impact of the variant on all transcript on which it occurs.
Each transcript consequence contains:
* intron/exon affected : the number of the intron/exon affected
* hgvsc : change in cDNA (cDNA has introns spliced out)
* hgvsp: if the transcript is protein coding, what is the effect on the protein?

                  


