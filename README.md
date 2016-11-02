# Phenopolis: an open platform for harmonization and analysis of sequencing and phenotype data

Preprint on [biorxiv](http://biorxiv.org/content/early/2016/10/31/084582).

You can access a demo version of the server at:
https://phenopolis.github.io
username/password:
demo/demo123

### Installation

Phenopolis requires:
* a running mongo database
* a running Phenotips server (needed to provide the HPO phenotypes per patient)
* (optionally) a running Exomiser stand-alone server, which can be obtained on request as it being developed separately by [Julius Jacobsen](https://github.com/julesjacobsen).

You will then be able to run the python Flask server.

The first step is to clone the repository.

```
git clone git@github.com:pontikos/phenopolis.git
```
Download Phenotips:
```
https://phenotips.org/Download
```
Install latest version of mongo:
```
https://www.mongodb.com/download-center#community
```
If you wish to download the Exomiser stand-alone server, please get in touch with [Julius Jacobsen](https://github.com/julesjacobsen).

### Creating database, importing data

First make sure mongoDB is running:
```
DBPATH=
mongod --dbpath $DBPATH --port 27017 --smallfiles
```

#### Creating and importing data from JSON

The variants found in the VCF files are processed with [Variant Effect Predictor (VEP)](http://www.ensembl.org/info/docs/tools/vep/) and the output is written to JSON standard output.
The standard output is piped into another python script, ```postprocess_VEP_json.py``` (available from another [repository](https://github.com/UCLGeneticsInstitute/DNASeq_pipeline)), which adds further annotation, formatting and writes output to JSON, which is then imported with mongoimport into the variants collection.

The bash command to run the VEP, assuming your variant files are:
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

The scripts can be found in `./pubmedScore`

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

The scripts can be found in `./phenogenon`

First, the user has to run `python snapshot_patient_hpo.py` to take a snapshot of patients' HPO at the time. Since the phenogenon analysis will take some time, this is to avoid any inconsistency that might be introduced by editing patients' HPO in the database when phenogenon is running.

Second, `python get_hpo_freq.py` will produce an HPO frequency file that phenogenon will use for its analysis.

Phenogenon can then be run as `python gene_hpo_analysis --chrom X` per chromosome. This feature can be utilised to parallelise the jobs on chromosomes. It uses `ExAC_freq` and `CADD_phred` scores to help filter the variants. The defaults are `ExAC_freq <= 0.01 and CADD_phred >= 15` for _recessive_ inheritance mode, and `ExAC_freq <= 0.001 and CADD_phred >= 15` for _dominant_ inheritance mode. It will produce a JSON file for each gene.

If one wishes to change the cutoffs to filter the variants after phenogenon is done, one can use `python recalculate_p.py --chrom X` to do the job quickly, without having to re-extracting info using the slow `gene_hpo_analysis.py`

After this, `python hpo_gene_anlaysis.py` will extract all genes with significant p values for each valid HPO term, and write to a JSON file for each HPO term.

## Running servers

Run Phenotips:
```
wget https://nexus.phenotips.org/nexus/content/repositories/releases/org/phenotips/phenotips-standalone/1.3-milestone-2/phenotips-standalone-1.3-milestone-2.zip
unzip phenotips-standalone-1.3-milestone-2.zip
cd phenotips-standalone-1.3-milestone-2
bash start.sh
```

Run Exomiser standalone:
```
EXOMISER_DATA=
cd $EXOMISER_DATA
wget ftp://ftp.sanger.ac.uk/pub/resources/software/exomiser/downloads/exomiser/exomiser-cli-7.2.1-data.zip
unzip exomiser-cli-7.2.1-data.zip
java -jar exomiser-rest-prioritiser-7.3.0-SNAPSHOT.jar --exomiser.data-directory=$EXOMISER_DATA
```

Run Phenopolis:
```
cd phenopolis
python run_server.py
```


### Acknowledgment

This code was originally forked from the [ExAC browser](https://github.com/konradjk/exac_browser) but has since diverged considerably.


