

import pickle
import pandas
import numpy as np
import pymongo
from lookups import *

POPS=pandas.read_csv('pop.csv')
POPS=dict(zip(POPS.code,POPS.phenotype))

def connect_db():
    client = pymongo.MongoClient(host='localhost', port=27017)
    return client['exac']

db=connect_db()


# all genes
genes=[g for g in db.genes.find()]
gene_ids=[g['gene_id'] for g in genes]

pickle.dump(gene_ids,open('gene_ids.pkl','wb'))

hom_per_gene=dict()
for code in POPS: 
    print(code,pheno)
    pheno=POPS[code]
    hom_per_gene[code]=[ len([v['pop_homs'][pheno] for v in db.variants.find({'genes': gid}, fields={'_id': False})]) for gid in gene_ids ]

hom_per_gene=pickle.load(open('hom_per_gene.pkl','rb'))



hom_per_gene2=dict()
df=pandas.DataFrame()
for code in POPS: 
    print(code,pheno)
    pheno=POPS[code]
    hom_per_gene2[code]=zip(gene_ids, np.array(hom_per_gene[code]))
    df2=pandas.DataFrame(hom_per_gene2[code])
    df2.columns=['gene_id',code]
    df=pandas.concat([df,df2],axis=1)




code='KC'
df.columns=['gene_id',code]

code2='EYES'
df2=pandas.DataFrame(hom_per_gene2[code2])
df2.columns=['gene_id',code2]

df=pandas.concat([df, df2],axis=1)

#pandas.to_csv(




print(get_metrics(db, '7-86416220-G-A'))

print(get_variants_in_gene(db, 'ENSG00000167207'))

get_gene_by_name(db, 'NOD2')

get_gene_by_name(db, 'TITIN')

print(genes[0])

variants_per_gene=[ len([v for v in db.variants.find({'genes': gid}, fields={'_id': False})]) for gid in gene_ids ]


x = np.array(variants_per_gene)
y = np.bincount(x)

ii = np.nonzero(y)[0]
xx=zip(ii,y[ii]) 


#MUC4
print(len([v for v in db.variants.find({'genes': 'ENSG00000145113'})]))
#TTN
print(len([v for v in db.variants.find({'genes': 'ENSG00000155657'})]))


muc4_variants=[v for v in db.variants.find({'genes': 'ENSG00000145113'})]

muc4_variants[0]



'pop_homs'
'pop_ans'
'pop_acs'


