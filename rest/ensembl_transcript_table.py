
import MySQLdb

# Open database connection
#mysql -u anonymous -h ensembldb.ensembl.org
db = MySQLdb.connect("ensembldb.ensembl.org","anonymous","","homo_sapiens_core_75_37" )

# prepare a cursor object using cursor() method
cursor = db.cursor()

# execute SQL query using execute() method.
cursor.execute("SELECT VERSION()")
# Fetch a single row using fetchone() method.
data = cursor.fetchone()
print "Database version : %s " % data

cursor.execute( """ SELECT transcript.stable_id, xref.display_label FROM transcript, object_xref, xref,external_db WHERE transcript.transcript_id = object_xref.ensembl_id AND object_xref.ensembl_object_type = 'Transcript' AND object_xref.xref_id = xref.xref_id AND xref.external_db_id = external_db.external_db_id AND external_db.db_name = 'RefSeq_mRNA'; """ )
cursor.fetchone()

# disconnect from server
db.close()


