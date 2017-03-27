import sys 

if len(sys.argv)>1 and sys.argv[1]=='SERVER':
    LOCAL=False
else:
    LOCAL=True

IMPORT_PYSAM_PRIMER3 = True # Set this to False to develop on Windows, where pysam and primer3 fail to install.

USE_ARGON2_AUTH = False # To test, populate USER_DB by entering  >>mongoimport --db users --collection users --file tests/data/users_argon2.json --drop