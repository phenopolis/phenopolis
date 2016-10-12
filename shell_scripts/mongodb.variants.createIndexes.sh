
mongo
use uclex
db.variants.createIndex({'VARIANT_ID':1})
db.variants.createIndex({'Transript':1})
db.variants.createIndex({'Gene':1})
db.variants.createIndex({'SYMBOL':1})
db.variants.createIndex({'HET':1})
db.variants.createIndex({'HOM':1})

