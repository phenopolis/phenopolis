[![Build Status](https://travis-ci.org/phenopolis/phenopolis.svg?branch=master)](https://travis-ci.org/phenopolis/phenopolis)
[![Coverage Status](https://coveralls.io/repos/github/phenopolis/phenopolis/badge.svg?branch=master)](https://coveralls.io/github/phenopolis/phenopolis?branch=master)


# Phenopolis: an open platform for harmonization and analysis of sequencing and phenotype data

![alt tag](https://github.com/phenopolis/phenopolis/blob/master/static/phenopolis-pipeline.png)

Preprint on [biorxiv](http://biorxiv.org/content/early/2016/10/31/084582).

Phenopolis is used for research into the molecular diagnosis of rare genetic diseases by clinicians, geneticists and bioinformaticians at UCL, University of Leeds, University of Manchester and University of Oxford.

## The Phenopolis Website 
You can access a demo version of the server at:
https://phenopolis.org
username/password:
demo/demo123


# Contributors
We are especially interested in contributions to the UI (html, css, js) which could be greatly refactored and vastly improved.
Also any performance improvements to the db queries would be also greatly appreciated.
Let us know if you run into difficulties getting the code running!  Our goal is to make it easy for you to contribute so the project continues to grow!
​
# Installation
This section includes guides to a quick install and a full installation.
​
## Prerequisites
* Python 2 - you will need to use python2 as we are not python3 compatible since packages such as pygr which use the old ```print``` syntax are not compatible with python3. https://www.python.org/downloads/
* MongoDB - https://www.mongodb.com/download-center#community
* Neo4j - https://neo4j.com/download/

## Quick Install Demo for Coders

This quick install is for people who want to get a local version up and running quickly to contribute to the codebase of the project.

I have a written a shell script for quick installation, [easy_install.sh](https://github.com/phenopolis/phenopolis/blob/master/easy_install.sh), on some example data that is downloadable from our website.  This will only take ~256M of disk space. 

Start Neo4j and, if this is the first time you've run Neo4j, log in and change the password. Set your Neo4j uri and password in easy_install.sh.

Clone this repo and run easy_install.sh. This will install packages and get data. When complete, you should be able to browse to:
[http://localhost:8000/gene/TTLL5](http://localhost:8000/gene/TTLL5)

The example dataset covers only gene TTLL5. Web pages for other genes will show no information.


### Windows - additional steps
Phenopolis can be developed under Windows but requires some additional steps and some lesser-used functionality will not be available.
* pip - make sure it is up to date by running ```python -m pip install -U pip```
* VCForPython27 - install this from http://aka.ms/vcpython27
* scipy and biopython - get and install the .whl files below from http://www.lfd.uci.edu/~gohlke/pythonlibs/
 * pip install numpy-1.11.3+mkl-cp27-cp27m-win32.whl –user
 * pip install scipy-0.18.1-cp27-cp27m-win32.whl –user
 * pip install biopython-1.68-cp27-cp27m-win32.whl –user
* Execute [the shell script](https://github.com/phenopolis/phenopolis/blob/master/easy_install.sh) 
* pysam and primer3 - disable the install, these packages won't install on Windows.
* In [config.py](https://github.com/phenopolis/phenopolis/blob/master/config/config.py) set ```IMPORT_PYSAM_PRIMER3 = False```
* Rerun [the shell script](https://github.com/phenopolis/phenopolis/blob/master/easy_install.sh) (you may disable the commands ```git clone```, ```wget```, ```mongoimport``` and ```mongo```).

To debug in Visual Studio, first turn off the Flask debug by setting ```app.run(..,..,..,debug=False)``` in ```runserver.py```.

### Post-installation
When this is installed you should be able to browse to:
[http://localhost:8000/gene/TTLL5](http://localhost:8000/gene/TTLL5)

The example dataset covers only gene TTLL5. Web pages for other genes will show no information. 

### Full Installation

Phenopolis requires:
* a running mongo database
* a running neo4j database
* (optionally) a running Exomiser stand-alone server, which can be obtained on request as it being developed separately by [Julius Jacobsen](https://github.com/julesjacobsen).

You will then be able to run ```phenopolis.py```, the python Flask server.

The first step is to clone the repository.

```
git clone git@github.com:phenopolis/phenopolis.git
```

If you wish to download the Exomiser stand-alone server, please get in touch with [Julius Jacobsen](https://github.com/julesjacobsen).

### Creating database, importing data

First make sure mongoDB is running:
```
DBPATH=<path to db>
mongod --dbpath $DBPATH --port 27017 --smallfiles
```

#### Creating and importing data from JSON

The variants found in the VCF files are processed with [Variant Effect Predictor (VEP)](http://www.ensembl.org/info/docs/tools/vep/) and the output is written to JSON standard output.
The standard output is piped into another python script, ```VEP/postprocess_VEP_json.py```, which adds further annotation, formatting and writes output to JSON, which is then imported with mongoimport into the variants collection.

The bash command to run the VEP, assuming your variant files are in VCF format:
```
bash VEP/runVEP.sh --input <infile> --output <outfile>
```
Importing of the variants can then be done:
```
mongoimport --db $DBNAME --collection variants --host $HOST < <infile>
```
Load individual for individual page (this is currently tedious, we are going to streamline this):
```
 python views/load_individual.py --individual $ID --auth Admin:$PASSWORD
```

#### Running pubmedScore

The pubmedscore, written by [Jing Yu](https://github.com/logust79), scores genes based on their pubmed relevance.

The scripts can be found in [pubmedScore](pubmedScore):

Before running the script, it is preferable to write patients ids in `patients.txt`, which `pubmedScore.py` takes by default.

```
python pubmedScore.py
    -i patients.txt
    -g ABCA4 (if specified, will ignore -i and -p)
    -p patientID_1 (if specified, will ignore -i)
    -k retina,retinal,retinitis,blindness,macula,macular,stargardt,pigmentosa (Keywords to search on pubmed. Displayed is default)
```

#### Running phenogenon

The phenogenon, written by [Jing Yu](https://github.com/logust79), does an enrichment test (Fisher test) per gene and HPO term.

The scripts can be found in [phenogenon](phenogenon):

First, the user has to run `python snapshot_patient_hpo.py` to take a snapshot of patients' HPO at the time. Since the phenogenon analysis will take some time, this is to avoid any inconsistency that might be introduced by editing patients' HPO in the database when phenogenon is running.

Second, `python get_hpo_freq.py` will produce an HPO frequency file that phenogenon will use for its analysis.

Phenogenon can then be run as `python gene_hpo_analysis --chrom X` per chromosome. This feature can be utilised to parallelise the jobs on chromosomes. It uses `ExAC_freq` and `CADD_phred` scores to help filter the variants. The defaults are `ExAC_freq <= 0.01 and CADD_phred >= 15` for _recessive_ inheritance mode, and `ExAC_freq <= 0.001 and CADD_phred >= 15` for _dominant_ inheritance mode. It will produce a JSON file for each gene.

If one wishes to change the cutoffs to filter the variants after phenogenon is done, one can use `python recalculate_p.py --chrom X` to do the job quickly, without having to re-extracting info using the slow `gene_hpo_analysis.py`

After this, `python hpo_gene_anlaysis.py` will extract all genes with significant p values for each valid HPO term, and write to a JSON file for each HPO term.

## Running servers

Run Exomiser standalone:
```
EXOMISER_DATA=
cd $EXOMISER_DATA
wget ftp://ftp.sanger.ac.uk/pub/resources/software/exomiser/downloads/exomiser/exomiser-cli-7.2.1-data.zip
unzip exomiser-cli-7.2.1-data.zip
java -jar exomiser-rest-prioritiser-7.3.0-SNAPSHOT.jar --exomiser.data-directory=$EXOMISER_DATA
```
There is also an online version available at the [Monarch Initiative](https://monarchinitiative.org).
```
https://monarch-exomiser-prod.monarchinitiative.org/exomiser/api/prioritise/
```

Run Phenopolis:
```
cd phenopolis
python runserver.py
```


### Acknowledgment

This code was originally forked from the [ExAC browser](https://github.com/konradjk/exac_browser) but has since diverged considerably.


