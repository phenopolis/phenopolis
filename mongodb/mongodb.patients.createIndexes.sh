

mongo --port 27016
use patients
db.patients.createIndex({'external_id':1})
db.patients.createIndex({'features':1})

