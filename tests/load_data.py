import json
import pymongo
from pymongo import MongoClient
from flask import Flask, current_app
import views
from config import config

# Load the test data set.
def load_data():
    app = Flask(__name__)
    with app.app_context():
        indexes = ['id', 'name']
        import_data('test_hpo', 'hpo', "./tests/data/hpo-hpo.json", indexes)
        indexes = ['gene_id', 'gene_name_upper', 'gene_name', 'other_names', 'xstart', 'xstop']
        import_data('test_uclex', 'genes', "./tests/data/uclex-genes.json", indexes)
        indexes = ['gene', 'mode', 'p']
        import_data('test_uclex', 'simreg', "./tests/data/uclex-simreg-TTLL5.json", indexes)
        # 'EXAC' should also be made an index but it throws an error - key too long.
        indexes = ['variant_id', 'CHROM', 'canonical_cadd', 'FILTER', 'canonical_transcript', 'hom_samples', 'het_samples', 'canonical_gene_name_upper']
        #indexes = ['variant_id', 'CHROM', 'canonical_cadd', 'EXAC', 'FILTER', 'canonical_transcript', 'hom_samples', 'het_samples', 'canonical_gene_name_upper']
        import_data('test_uclex', 'variants', "./tests/data/uclex-variant-TTLL5.json", indexes)
        indexes = ['gene_id']
        import_data('test_uclex', 'gene_hpo', "./tests/data/uclex-gene_hpo-TTLL5.json", indexes)
        indexes = ['gene_id', 'gene_name']
        import_data('test_hpo', 'gene_hpo', "./tests/data/hpo-gene_hpo-TTLL5.json", indexes)
        indexes = ['gene', 'hpo']
        import_data('test_hpo', 'genes_pheno', "./tests/data/hpo-genes_pheno-TTLL5.json", indexes)
        indexes = ['external_id', 'report_id', 'features.id', 'sex', 'genes.gene', 'solved', 'clinicalStatus.clinicalStatus', 'specificity.score']
        import_data('test_patients', 'patients', "./tests/data/patients-patients-hidden.json", indexes)
        import_data('test_users', 'users', "./tests/data/users.json")

# Load the test data set needed for test_login.
def load_user_data():
    app = Flask(__name__)
    with app.app_context():
        import_data('test_users', 'users', "./tests/data/users.json")

# Create a collection in the db, drop any existing data, add data from file.
# Create indexes in to the data.
def import_data(db_name, collection_name, file_location, indexes=None):
    db = views.get_db(db_name)
    collection = db.get_collection(collection_name)
    collection.drop()
    with open(file_location, 'r') as json_data:
        for line in json_data:
            dataset = json.loads(line.replace("$oid", "ObjectId").replace("$numberLong", "NumberLong"))
            collection.insert(dataset,check_keys=False)
    if indexes:
        for item in indexes:
            collection.create_index(item)
    
