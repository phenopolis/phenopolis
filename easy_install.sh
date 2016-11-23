
git clone git@github.com:pontikos/phenopolis.git

# For local install without Phenotips, approve all logins by uncommenting line in views/__init__.py

sed -i '' 's/#NO_PHENOTIPS_INSTALLATION: //' phenopolis/views/__init__.py

# Basic build of db

wget --no-check-certificate https://uclex.cs.ucl.ac.uk/static/demo/uclex-genes.json
mongoimport --db uclex --collection genes --file uclex-genes.json --drop
mongo uclex --eval "db.variants.createIndex({'gene_id' : 1})"
mongo uclex --eval "db.variants.createIndex({'gene_name_upper' : 1})"
mongo uclex --eval "db.variants.createIndex({'gene_name' : 1})"
mongo uclex --eval "db.variants.createIndex({'other_names' : 1})"
mongo uclex --eval "db.variants.createIndex({'xstart' : 1})"
mongo uclex --eval "db.variants.createIndex({'xstop' : 1})"

wget --no-check-certificate https://uclex.cs.ucl.ac.uk/static/demo/uclex-simreg-TTLL5.json
mongoimport --db uclex --collection simreg --file uclex-simreg-TTLL5.json --drop
mongo uclex --eval "db.variants.createIndex({'gene' : 1})"
mongo uclex --eval "db.variants.createIndex({'mode' : 1})"
mongo uclex --eval "db.variants.createIndex({'p' : 1})"

wget --no-check-certificate https://uclex.cs.ucl.ac.uk/static/demo/uclex-variant-TTLL5.json
mongoimport --db uclex --collection variants --file uclex-variant-TTLL5.json --drop
mongo uclex --eval "db.variants.createIndex({'variant_id':1})"
mongo uclex --eval "db.variants.createIndex({'CHROM' : 1})"
mongo uclex --eval "db.variants.createIndex({'canonical_cadd' : 1})"
mongo uclex --eval "db.variants.createIndex({'EXAC' : 1})"
mongo uclex --eval "db.variants.createIndex({'FILTER' : 1})"
mongo uclex --eval "db.variants.createIndex({'canonical_transcript' : 1})"
mongo uclex --eval "db.variants.createIndex({'hom_samples' : 1})"
mongo uclex --eval "db.variants.createIndex({'het_samples' : 1})"
mongo uclex --eval "db.variants.createIndex({'canonical_gene_name_upper' : 1})"

wget --no-check-certificate https://uclex.cs.ucl.ac.uk/static/demo/uclex-gene_hpo-TTLL5.json
mongoimport --db uclex --collection gene_hpo --file uclex-gene_hpo-TTLL5.json --drop
mongo uclex --eval "db.variants.createIndex({'gene_id' : 1})"

wget --no-check-certificate https://uclex.cs.ucl.ac.uk/static/demo/hpo-gene_hpo-TTLL5.json
mongoimport --db hpo --collection gene_hpo --file hpo-gene_hpo-TTLL5.json --drop
mongo hpo --eval "db.variants.createIndex({'gene_id' : 1})"
mongo hpo --eval "db.variants.createIndex({'gene_name' : 1})"


wget --no-check-certificate https://uclex.cs.ucl.ac.uk/static/demo/hpo-genes_pheno-TTLL5.json
mongoimport --db hpo --collection genes_pheno --file hpo-genes_pheno-TTLL5.json --drop
mongo hpo --eval "db.variants.createIndex({'gene' : 1})"
mongo hpo --eval "db.variants.createIndex({'hpo' : 1})"


wget --no-check-certificate https://uclex.cs.ucl.ac.uk/static/demo/patients-patients-hidden.json
mongoimport --db patients --collection patients --file patients-patients-hidden.json --drop
mongo patients --eval "db.variants.createIndex({'external_id' : 1})"
mongo patients --eval "db.variants.createIndex({'report_id' : 1})"
mongo patients --eval "db.variants.createIndex({'features.id' : 1})"
mongo patients --eval "db.variants.createIndex({'sex' : 1})"
mongo patients --eval "db.variants.createIndex({'genes.gene' : 1})"
mongo patients --eval "db.variants.createIndex({'solved' : 1})"
mongo patients --eval "db.variants.createIndex({'clinicalStatus.clinicalStatus' : 1})"
mongo patients --eval "db.variants.createIndex({'specificity.score' : 1})"

# Run server

cd phenopolis 
python run_server.py


