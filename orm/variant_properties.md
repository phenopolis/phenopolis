
# Mapping of variant properties from VEP JSON output to Pheno4J

| |  Value | Explanation | Pheno4J mapping |
| ------------- | ------------- | ----- | ----- |
|allele_string|A/G| | GeneticVariant.allele_string |
|start|38212762| | GeneticVariant.start |
|end|38212762| | GeneticVariant.end |
|seq_region_name| 22| chromosome | GeneticVariant.chrom |
|most_severe_consequence|intron_variant| worst consequence on a transcript | GeneticVariant.most_severe_consequence |
|strand|1| posi/nega tive strand | GeneticVariant.stand |


| custom_annotations |    Value |  Explanation | Pheno4J mapping |
| ------------- | ------------- | ----- | ----- |
| input[0].fields.name | "22-38212762-A-G" | | GeneticVariant.variantId |
| input[0].fields.AC | 6339 | internal allele count | GeneticVariant.AC|
| input[0].fields.AF | 0.651 | internal allele frequency | GeneticVariant.allele_freq |
| input[0].fields.AN | 9740 | internal allele number | GeneticVariant.AN | 
| input[0].fields.BaseQRankSum | 0.481 |
| input[0].fields.ClippingRankSum | 0.03 |
| input[0].fields.DP | 245471 | | |
| input[0].fields.ExcessHet | 0 | Excess heterozygosity | GeneticVariant.ExcessHet | 
| input[0].fields.FS | 1.189 | Fisher strand | GeneticVariant.FS |
| input[0].fields.InbreedingCoeff  | 0.0552 |  quality | GeneticVariant.InbreedingCoeff |
| input[0].fields.MLEAC | 6341 | | GeneticVariant.MLEAC |
| input[0].fields.MLEAF | 0.651 | | GeneticVariant.MLEAF |
| input[0].fields.MQ | 23.88 | mapping quality | GeneticVariant.MQ |
| input[0].fields.MQ0 | 0 | | |
| input[0].fields.MQRankSum | 0.032 | | GeneticVariant.MQRankSum |
| input[0].fields.QD | 20.93 | | |
| input[0].fields.ReadPosRankSum | 0.326 | | GeneticVariant.ReadPosRankSum | 
| input[0].fields.SOR | 0.829 | | |
| input[0].fields.VQSLOD | 3.37 | variant quality score log odds | GeneticVariant.VQSLOD |
| input[0].fields.culprit | FS | | Culprit | GeneticVariant.Culprit |
| gnomad_(exomes/genomes)[0].fields.AC_AFR|7224| allele count in africans | GeneticVariant.gnomad_(exomes/genomes)_AC_AFR |
| gnomad_(exomes/genomes)[0].fields.AC_AMR|540| americans | GeneticVariant.gnomad_(exomes/genomes)_AC_AMR |
| gnomad_(exomes/genomes)[0].fields.AC_ASJ|184| | |
| gnomad_(exomes/genomes)[0].fields.AC_EAS|1441| | |
| gnomad_(exomes/genomes)[0].fields.AC_FIN|2260| | |
| gnomad_(exomes/genomes)[0].fields.AC_Female|9802| | |
| gnomad_(exomes/genomes)[0].fields.AC_Male|12193| | |
| gnomad_(exomes/genomes)[0].fields.AC_NFE|9704| | |
| gnomad_(exomes/genomes)[0].fields.AC_OTH|642| | |
| gnomad_(exomes/genomes)[0].fields.AC_raw|22041| raw allele count | GeneticVariant.gnomad_(exomes/genomes)_AC_raw |
| gnomad_(exomes/genomes)[0].fields.AF_AFR|0.828821| | |
| gnomad_(exomes/genomes)[0].fields.AF_AMR|0.644391| | |
| gnomad_(exomes/genomes)[0].fields.AF_ASJ|0.609272| | |
| gnomad_(exomes/genomes)[0].fields.AF_EAS|0.891708| | |
| gnomad_(exomes/genomes)[0].fields.AF_FIN|0.647564| | |
| gnomad_(exomes/genomes)[0].fields.AF_Female|0.709159|
| gnomad_(exomes/genomes)[0].fields.AF_Male|0.713959|
| gnomad_(exomes/genomes)[0].fields.AF_NFE|0.64875| allele frequency in non-finish europeans | GeneticVariant.gnomad_(exomes/genomes)_AF_NFE |
| gnomad_(exomes/genomes)[0].fields.AF_OTH|0.655102| | |
| gnomad_(exomes/genomes)[0].fields.AF_raw|0.711184| | |
| gnomad_(exomes/genomes)[0].fields.AN_AFR|8716| | |
| gnomad_(exomes/genomes)[0].fields.AN_AMR|838| | |
| gnomad_(exomes/genomes)[0].fields.AN_ASJ|302| | |
| gnomad_(exomes/genomes)[0].fields.AN_EAS|1616| | |
| gnomad_(exomes/genomes)[0].fields.AN_FIN|3490| | |
| gnomad_(exomes/genomes)[0].fields.AN_Female|13822| | |
| gnomad_(exomes/genomes)[0].fields.AN_Male|17078| | |
| gnomad_(exomes/genomes)[0].fields.AN_NFE|14958| | |
| gnomad_(exomes/genomes)[0].fields.AN_OTH|980| | |
| gnomad_(exomes/genomes)[0].fields.AN_raw|30992| | |
| gnomad_(exomes/genomes)[0].fields.Hom_AFR|2989| | |
| gnomad_(exomes/genomes)[0].fields.Hom_AMR|171| | |
| gnomad_(exomes/genomes)[0].fields.Hom_ASJ|55| | |
| gnomad_(exomes/genomes)[0].fields.Hom_EAS|645| | |
| gnomad_(exomes/genomes)[0].fields.Hom_FIN|722| | |
| gnomad_(exomes/genomes)[0].fields.Hom_Female|3561| | |
| gnomad_(exomes/genomes)[0].fields.Hom_Male|4413| | |
| gnomad_(exomes/genomes)[0].fields.Hom_NFE|3177| | |
| gnomad_(exomes/genomes)[0].fields.Hom_OTH|215| | |
| gnomad_(exomes/genomes)[0].fields.Hom_raw|7990| | |
| gnomad_(exomes/genomes)[0].fields.Hom|7974| | |
| kaviar[0].fields.AF | 0.6665359 | kaviar alelle frequency | GeneticVariant.kaviar_AF |
| kaviar[0].fields.AC | 103649 | kaviar allele count | GeneticVariant.kaviar_AC |
| kaviar[0].fields.AN | 155504 | kaviar allele number | GeneticVariant.kaviar_AN |
| dbsnp[0].name | "rs13911" | dbsnp name | GeneticVariant.dbsnp_name | 


| transcript_consequences[*] | Value |  Explanation | Pheno4J mapping |
| ------------- | ------------- | ----- | ----- |
| exac_af|0.682| ExAC allele frequency| GeneticVariant.exac_AF |
| exac_af_adj|0.685| ExAC allele frequency in adjusted | GeneticVariant.exac_AF_ADJ |
| exac_af_afr|0.824| ExAC allele frequency in africans | GeneticVariant.exac_AF_AFR |
| exac_af_amr|0.624| ExAC allele frequency in americans | GeneticVariant.exac_AF_AMR |
| exac_af_consanguineous|0.73| ExAC allele frequency in consanguineous | GeneticVariant.exac_AF_CONSANGUINEOUS |
| exac_af_eas|0.878| ExAC allele frequency in east asians | GeneticVariant.exac_AF_EAS |
| exac_af_female|0.692| ExAC allele frequency in females | GeneticVariant.exac_AF_FEMALE |
| exac_af_fin|0.634| ExAC allele frequency in finnish | GeneticVariant.exac_AF_FIN |
| exac_af_male|0.679| ExAC allele frequency in males | GeneticVariant.exac_AF_MALE |
| exac_af_nfe|0.643| ExAC allele frequency in non-finnish europeans | GeneticVariant.exac_AF_NFE |
| exac_af_oth|0.713| ExAC allele frequency in others | GeneticVariant.exac_AF_OTH |
| exac_af_popmax|0.878| Max ExAC allele frequency | GeneticVariant.exac_AF_POPMAX |
| exac_af_sas|0.728| ExAC allele frequency in south asians | GeneticVariant.exac_AF_SAS |
| cadd_phred | 1.563 | Phred CADD score | GeneticVariant.cadd_phred |
| cadd_raw | -0.123251 | Raw CADD score | GeneticVariant.cadd_raw |
| hgvsc | "ENST00000323205.6:c.*7+30A>G" | cDNA mutation in HGVS format | TranscriptVariant.hgvsc |
| hgvsp |  | protein mutation in HGVS format | TranscriptVariant.hgvsp |
| impact | "MODIFIER" |  | TranscriptVariant.impact |
| consequence_terms | ["intron_variant"] | consequence of mutation | TranscriptVariant.consequence_terms | 
| transcript_id | "ENST00000323205" | ensembl id of transcript | TranscriptVariant.trancript_id |
| intron/exon | "9/9" | which intron/exon is mutation in | TranscriptVariant.intron/exon |
| gene_id | "ENSG00000100116" | ensembl gene id | TranscriptVariant.gene_id, Gene.gene_id |
| gene_symbol | "GCAT" | | TranscriptVariant.gene_name, Gene.gene_name |
| hgnc_id | 4188 | | Gene.hgnc_id |
| transcript_id | "ENST00000323205" | ensembl id of transcript | Transcript.trancript_id |
| canonical | 1 | canonical transcript | Transcript.canonical |
| strand | 1 | positive/negative strand | Transcript.strand |
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

                  
