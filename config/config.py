import sys 

if len(sys.argv)>1 and sys.argv[1]=='SERVER':
    LOCAL=False
    LOCAL_WITH_PHENOTIPS = False 
else:
    LOCAL=True
    LOCAL_WITH_PHENOTIPS = True # Set to True to use Phenotips while in local mode.

IMPORT_PYSAM_PRIMER3 = False # Set this to False to develop on Windows, where pysam and primer3 fail to install.