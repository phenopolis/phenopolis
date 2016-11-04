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

### genes
Genes table, loaded from Gencode file.

### hpo_freq

### gene_hpo
Phenogenon collection.
```
{
	"_id" : ObjectId("581257ef2ebb946cc56d9a77"),
	"het" : {

	},
	"description" : "add online-calculated CADD phred scores",
	"hom_comp" : {

	},
	"gene_id" : "ENSG00000015568",
	"version" : 4,
	"release" : "2016_Aug"
}
```

### hpo_gene
Phenogenon collection.
```
{
	"cadd_cutoff" : 15,
	"exac_dominant_cutoff" : 0.001,
	"exac_recessive_cutoff" : 0.01,
	"version" : 4,
	"hpo_id" : "HP:0000001",
	"release" : "2016_Aug",
	"data" : {
		"unrelated" : {
			"recessive" : [ ],
			"dominant" : [ ]
		},
		"related" : {
			"recessive" : [ ],
			"dominant" : [ ]
		}
	}
}
```

### lof

### patients

### read_depth

### retnet
```
{
	"omim" : [
		"268000",
		"608172",
		"613861"
	],
	"mode" : "r",
	"gene_name" : "DHDDS",
	"disease" : "recessive retinitis pigmentosa; protein: dehydrodolichyl diphosphate synthetase"
}
```

### simreg
```
{
	"_id" : ObjectId("580fa7c22ebb946cc56d5f19"),
	"gene" : "ENSG00000012779",
	"phi" : {
		"HP:0012045" : {
			"hpo_id" : "HP:0012045",
			"desc" : "Retinal flecks",
			"prob" : 0.45
		},
		"HP:0030506" : {
			"hpo_id" : "HP:0030506",
			"desc" : "Yellow/white lesions of the retina",
			"prob" : 0.45
		},
		"HP:0002597" : {
			"hpo_id" : "HP:0002597",
			"desc" : "Abnormality of the vasculature",
			"prob" : 0.12
		},
		"HP:0011276" : {
			"hpo_id" : "HP:0011276",
			"desc" : "Vascular skin abnormality",
			"prob" : 0.03
		},
		"HP:0001009" : {
			"hpo_id" : "HP:0001009",
			"desc" : "Telangiectasia",
			"prob" : 0.03
		},
		"HP:0000005" : {
			"hpo_id" : "HP:0000005",
			"desc" : "Mode of inheritance",
			"prob" : 0.02
		},
		"HP:0030466" : {
			"hpo_id" : "HP:0030466",
			"desc" : "Abnormal full-field electroretinogram",
			"prob" : 0.02
		},
		"HP:0000007" : {
			"hpo_id" : "HP:0000007",
			"desc" : "Autosomal recessive inheritance",
			"prob" : 0.02
		},
		"HP:0003745" : {
			"hpo_id" : "HP:0003745",
			"desc" : "Sporadic",
			"prob" : 0.02
		},
		"HP:0030453" : {
			"hpo_id" : "HP:0030453",
			"desc" : "Abnormal visual electrophysiology",
			"prob" : 0.01
		},
		"HP:0000479" : {
			"hpo_id" : "HP:0000479",
			"desc" : "Abnormality of the retina",
			"prob" : 0.01
		},
		"HP:0001626" : {
			"hpo_id" : "HP:0001626",
			"desc" : "Abnormality of the cardiovascular system",
			"prob" : 0.01
		}
	},
	"p" : 0.0002740577430736998,
	"mode" : "dom"
}
```

### solved_patients

### transcripts
```
{
	"_id" : ObjectId("57a4ba17cda06e21ec540df7"),
	"start" : 11870,
	"transcript_id" : "ENST00000456328",
	"strand" : "+",
	"stop" : 14410,
	"xstart" : 1000011870,
	"chrom" : "1",
	"gene_id" : "ENSG00000223972",
	"xstop" : 1000014410
}
```
