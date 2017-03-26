# Get phenopolis and its submodule.
git clone https://github.com/pontikos/phenopolis.git
cd phenopolis
git submodule update --init --recursive
cd varnorm
python setup.py install --user
cd ../..

# python packages required
pip install Cython==0.25.2 --install-option="--no-cython-compile" --user 
pip install -r phenopolis/requirements.txt --user

# Make sure mongodb is running
DBPATH=db
mkdir -p $DBPATH
mongod --dbpath $DBPATH --port 27017 &

# Basic build of db
# download minimal files
# import them
# create indexes
cd phenopolis

#wget --no-check-certificate https://uclex.cs.ucl.ac.uk/static/demo/hpo-hpo.json -O hpo-hpo.json
mongoimport --db hpo --collection hpo --file tests/data/hpo-hpo.json --drop
mongo hpo --eval "db.hpo.createIndex({'id' : 1})"
mongo hpo --eval "db.hpo.createIndex({'name' : 1})"

#wget --no-check-certificate https://uclex.cs.ucl.ac.uk/static/demo/uclex-genes.json -O uclex-genes.json
mongoimport --db uclex --collection genes --file tests/data/uclex-genes.json --drop
mongo uclex --eval "db.genes.createIndex({'gene_id' : 1})"
mongo uclex --eval "db.genes.createIndex({'gene_name_upper' : 1})"
mongo uclex --eval "db.genes.createIndex({'gene_name' : 1})"
mongo uclex --eval "db.genes.createIndex({'other_names' : 1})"
mongo uclex --eval "db.genes.createIndex({'xstart' : 1})"
mongo uclex --eval "db.genes.createIndex({'xstop' : 1})"

#wget --no-check-certificate https://uclex.cs.ucl.ac.uk/static/demo/uclex-simreg-TTLL5.json -O uclex-simreg-TTLL5.json
mongoimport --db uclex --collection simreg --file tests/data/uclex-simreg-TTLL5.json --drop
mongo uclex --eval "db.simreg.createIndex({'gene' : 1})"
mongo uclex --eval "db.simreg.createIndex({'mode' : 1})"
mongo uclex --eval "db.simreg.createIndex({'p' : 1})"

#wget --no-check-certificate https://uclex.cs.ucl.ac.uk/static/demo/uclex-variant-TTLL5.json -O uclex-variant-TTLL5.json
mongoimport --db uclex --collection variants --file tests/data/uclex-variant-TTLL5.json --drop
mongo uclex --eval "db.variants.createIndex({'variant_id':1})"
mongo uclex --eval "db.variants.createIndex({'CHROM' : 1})"
mongo uclex --eval "db.variants.createIndex({'canonical_cadd' : 1})"
mongo uclex --eval "db.variants.createIndex({'EXAC' : 1})"
mongo uclex --eval "db.variants.createIndex({'FILTER' : 1})"
mongo uclex --eval "db.variants.createIndex({'canonical_transcript' : 1})"
mongo uclex --eval "db.variants.createIndex({'hom_samples' : 1})"
mongo uclex --eval "db.variants.createIndex({'het_samples' : 1})"
mongo uclex --eval "db.variants.createIndex({'canonical_gene_name_upper' : 1})"

#wget --no-check-certificate https://uclex.cs.ucl.ac.uk/static/demo/uclex-gene_hpo-TTLL5.json -O uclex-gene_hpo-TTLL5.json
mongoimport --db uclex --collection gene_hpo --file tests/data/uclex-gene_hpo-TTLL5.json --drop
mongo uclex --eval "db.gene_hpo.createIndex({'gene_id' : 1})"

#wget --no-check-certificate https://uclex.cs.ucl.ac.uk/static/demo/hpo-gene_hpo-TTLL5.json -O hpo-gene_hpo-TTLL5.json
mongoimport --db hpo --collection gene_hpo --file tests/data/hpo-gene_hpo-TTLL5.json --drop
mongo hpo --eval "db.gene_hpo.createIndex({'gene_id' : 1})"
mongo hpo --eval "db.gene_hpo.createIndex({'gene_name' : 1})"

#wget --no-check-certificate https://uclex.cs.ucl.ac.uk/static/demo/hpo-genes_pheno-TTLL5.json -O hpo-genes_pheno-TTLL5.json
mongoimport --db hpo --collection genes_pheno --file tests/data/hpo-genes_pheno-TTLL5.json --drop
mongo hpo --eval "db.genes_pheno.createIndex({'gene' : 1})"
mongo hpo --eval "db.genes_pheno.createIndex({'hpo' : 1})"

#wget --no-check-certificate https://uclex.cs.ucl.ac.uk/static/demo/patients-patients-hidden.json -O patients-patients-hidden.json
mongoimport --db patients --collection patients --file tests/data/patients-patients-hidden.json --drop
mongo patients --eval "db.patients.createIndex({'external_id' : 1})"
mongo patients --eval "db.patients.createIndex({'report_id' : 1})"
mongo patients --eval "db.patients.createIndex({'features.id' : 1})"
mongo patients --eval "db.patients.createIndex({'sex' : 1})"
mongo patients --eval "db.patients.createIndex({'genes.gene' : 1})"
mongo patients --eval "db.patients.createIndex({'solved' : 1})"
mongo patients --eval "db.patients.createIndex({'clinicalStatus.clinicalStatus' : 1})"
mongo patients --eval "db.patients.createIndex({'specificity.score' : 1})"

mongoimport --db users --collection users --file tests/data/users.json --drop


# Run server 
python runserver.py

exec $SHELL

