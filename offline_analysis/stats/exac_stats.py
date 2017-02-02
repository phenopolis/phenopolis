
from collections import Counter
import pymongo

conn = pymongo.MongoClient(host='localhost', port=27017)
db=conn['uclex-old']

#print(Counter([f for f in db.variants.find({'in_exac':False},{'filter':1,'_id':0})]))
#print(db.variants.find({'in_exac':False},{'filter':1,'_id':0}).count())
print(db.variants.find({'in_exac':False,'filter':'PASS'}).count())

#print(Counter([f for f in db.variants.find({'in_exac':True},{'filter':1,'_id':0})]))
#print(db.variants.find({'in_exac':True},{'filter':1,'_id':0}).count())
print(db.variants.find({'in_exac':True,'filter':'PASS'}).count())


# what is the distribution of passed variants no in ExAC
cumsum=[ db.variants.find({'in_exac':False,'filter':'PASS','allele_count':{'$gt':n}}).count() for n in range(0,5000) ]


HET=[]
for v in db.variants.find({'in_exac':False,'filter':'PASS','allele_count':2}):
    print(v['variant_id'], v.get('HET_INDIVIDUALS',[]), v.get('HOM_INDIVIDUALS',[]) )
    HET.append(v.get('HET_INDIVIDUALS',[]))

c=Counter([tuple(h) for h in HET])

sorted_x=sorted(c.items(),key=operator.itemgetter(1),reverse=True)

sorted_x[1:20]

import pickle
pickle.dump(sorted_x,file=file('shared_private_variants.pyd','wb'))


file('cumsum.txt','w').writelines(['%s,%s\n'%(i,x,) for i,x, in enumerate(cumsum)])

