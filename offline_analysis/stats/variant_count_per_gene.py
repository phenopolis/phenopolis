
import csv
from collections import Counter
import pymongo

conn = pymongo.MongoClient(host='localhost', port=27017)
db=conn['uclex-old']

pli_file='/slms/UGI/vm_exports/vyp/phenotips/ExAC/0.3.1/functional_gene_constraint/fordist_cleaned_exac_nonTCGA_z_pli_rec_null_data.txt'
pli=csv.DictReader(file(pli_file,'r'),delimiter='\t')
for l in pli:
    nvar_exp=sum([int(l[x]) for x in ('n_syn','n_lof','n_mis')])
    gene_name=l['gene']
    gene=db.genes.find_one({'gene_name_upper':gene_name})
    if gene is None:
        #print(gene,'missing',)
        continue
    gene_id=gene['gene_id']
    nvar_obs=db.variants.find({'genes':gene_id,'filter':'PASS'}).count()
    if nvar_obs>nvar_exp:
        print(gene_name, gene_id, nvar_exp,nvar_obs)
        
    


