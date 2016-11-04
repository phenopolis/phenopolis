# Mongodb database documentation

## hpo

### hpo
Stores the HPO ontology downloaded from the HPO website:
```
{
	"comment" : [
		"Root of all terms in the Human Phenotype Ontology."
	],
	"id" : [
		"HP:0000001"
	],
	"name" : [
		"All"
	]
}
```

## patients

### patients
Sync of the Phenotips database.  All patiens are accessed in JSON format from the API:
```
```

### hpo_cache
Stores all the HPO to person associations to allow fast lookup.
No need to traverse the HPO graph.
```
{
	
	"external_id" : ["individual_id"],
	"hpo_id" : "HP:0000003"
}
```

### hpo_freq



## pubmedbatch

## uclex

### variants
```
{
	"FS" : 0,
	"MISS_COUNT" : 4090,
	"WT_COUNT" : 711,
	"allele_count" : 2,
	"canonical_cadd" : [
		8.785
	],
	"variant_id" : "10-5674924-C-A",
	"hom_samples" : [
		""
	],
	"REF" : "C",
	"id" : "10-5674924-C-A",
	"DP" : 13544,
	"canonical_gene_name_upper" : [
		"RN7SL445P"
	],
	"allele_string" : "C/A",
	"allele_freq" : 0.0014064697609001407,
	"end" : 5674924,
	"culprit" : "FS",
	"synonym_variant_id" : "10-5674924-C-A",
	"index" : 1,
	"AF" : 0.001404,
	"POS" : 5674924,
	"ID" : "10-5674924-C-A",
	"start" : 5674924,
	"SOR" : 0.693,
	"hgvs" : "chr10:g.5674924C>A",
	"wt_samples" : [ ],
	"CHROM" : 10,
	"het_samples" : [ ],
	"AC" : 2,
	"EXAC" : null,
	"ExcessHet" : 0.0046,
	"VQSLOD" : 4.49,
	"HET_COUNT" : 0,
	"genes" : [
		"ENSG00000240577"
	],
	"AN" : 1424,
	"QUAL" : 35.62,
	"HOM_COUNT" : 1,
	"MLEAF" : 0.0007022,
	"transcript_consequences" : [
		{
			"impact" : "MODIFIER",
			"distance" : 1208,
			"hgnc_id" : 46461,
			"gene_id" : "ENSG00000240577",
			"gene_symbol_source" : "HGNC",
			"gene_symbol" : "RN7SL445P",
			"transcript_id" : "ENST00000495609",
			"cadd" : 8.785,
			"consequence_terms" : [
				"downstream_gene_variant"
			],
			"variant_allele" : "A",
			"strand" : 1,
			"canonical" : 1
		}
	],
	"most_severe_consequence" : "downstream_gene_variant",
	"strand" : 1,
	"allele_num" : 1424,
	"assembly_name" : "GRCh37",
	"InbreedingCoeff" : 0.1189,
	"MLEAC" : 1,
	"canonical_hgvsp" : [
		""
	],
	"canonical_transcript" : [
		"ENST00000495609"
	],
	"FILTER" : "PASS",
	"seq_region_name" : 10,
	"MQ" : 35,
	"QD" : 17.81,
	"ALT" : "A",
	"canonical_hgvsc" : [
		""
	]
}
```
### causal_variants
Variants which have been identified as causal

### ensembl_entrez

### exons

### gene_hpo

### genes
Genes table, loaded from Gencode file.

### hpo_freq

### hpo_gene

### lof

### patients

### read_depth

### retnet

### simreg

### solved_patients

### transcripts
