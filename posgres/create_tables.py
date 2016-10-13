
import utils
import re
import traceback
#from utils import *
import copy
import sys
from collections import Counter
from collections import OrderedDict
import pysam
import itertools
import psycopg2

conn = psycopg2.connect(database="uclex")
cur = conn.cursor()

variants_headers=OrderedDict( [('variant_id','varchar PRIMARY KEY '), ('chrom','smallint'), ('pos','bigint'), ('xstart','bigint'), ('xstop','bigint'), ('ref','varchar'), ('alt','varchar'), ('filter','varchar'), ('site_quality','decimal'), ('miss_count','int'), ('wt_count','int'), ('het_count','int'), ('hom_count','int')] )
sql_create_table="CREATE TABLE variants ({});".format(','.join([' '.join(x) for x in variants_headers.items()]))
print(sql_create_table)
cur.execute( sql_create_table )

sql_create_table="CREATE TABLE het_variants_patients (variant_id varchar references variants(variant_id), patient_id varchar);"
print(sql_create_table)
cur.execute( sql_create_table )

sql_create_table="CREATE TABLE hom_variants_patients (variant_id varchar references variants(variant_id), patient_id varchar);"
print(sql_create_table)
cur.execute( sql_create_table )

conn.commit()


