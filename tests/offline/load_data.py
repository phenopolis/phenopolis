
import json
import pymongo
from pymongo import MongoClient
import sys
sys.path.append('offline_analysis/commons')
from phenopolis_utils import *

dbs = get_mongo_collections({
    'hpo_db':'test_hpo',
    'phenopolis_db':'test_uclex',
    'patient_db':'test_patients',
})
# Load the test data set.
def load_data():
    print '!!!!!!!'
    indexes = ['id', 'name']
    import_data('hpo_db', 'hpo', "./tests/data/hpo-hpo.json", indexes)

    indexes = ['gene_id', 'gene_name_upper', 'gene_name', 'other_names', 'xstart', 'xstop']
    import_data('phenopolis_db', 'genes', "./tests/data/uclex-genes.json", indexes)

    indexes = ['gene', 'mode', 'p']
    import_data('hpo_db', 'simreg', "./tests/data/uclex-simreg-TTLL5.json", indexes)

    # 'EXAC' should also be made an index but it throws an error - key too long.
    indexes = ['variant_id', 'CHROM', 'canonical_cadd', 'FILTER', 'canonical_transcript', 'hom_samples', 'het_samples', 'canonical_gene_name_upper']
    #indexes = ['variant_id', 'CHROM', 'canonical_cadd', 'EXAC', 'FILTER', 'canonical_transcript', 'hom_samples', 'het_samples', 'canonical_gene_name_upper']
    import_data('phenopolis_db', 'variants', "./tests/data/uclex-variant-TTLL5.json", indexes)

    indexes = ['gene_id']
    import_data('phenopolis_db', 'gene_hpo', "./tests/data/uclex-gene_hpo-TTLL5.json", indexes)

    indexes = ['gene_id', 'gene_name']
    import_data('hpo_db', 'gene_hpo', "./tests/data/hpo-gene_hpo-TTLL5.json", indexes)

    indexes = ['gene', 'hpo']
    import_data('hpo_db', 'genes_pheno', "./tests/data/hpo-genes_pheno-TTLL5.json", indexes)

    indexes = ['external_id', 'report_id', 'features.id', 'sex', 'genes.gene', 'solved', 'clinicalStatus.clinicalStatus', 'specificity.score']
    import_data('patient_db', 'patients', "./tests/data/patients-patients-hidden.json", indexes)


# Create a collection in the db, drop any existing data, add data from file.
# Create indexes in to the data.
def import_data(db_name, collection_name, file_location, indexes=None):
    db = dbs[db_name]
    collection = db[collection_name]
    collection.drop()
    with open(file_location, 'r') as json_data:
        for line in json_data:
            dataset = json.loads(line.replace("$oid", "ObjectId").replace("$numberLong", "NumberLong"))
            collection.insert(dataset)

    if indexes:
        for item in indexes:
            collection.create_index(item)
    
