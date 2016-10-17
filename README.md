# Phenopolis



### Installation

Phenopolis requires a running mongo database and a running Phenotips server.

Clone the repository.

```
git@github.com:pontikos/phenopolis.git
````

Download Phenotips.
```
https://phenotips.org/Download
```

Download the Exomiser stand alone.
```
```

Install latest version of mongo.
```
https://www.mongodb.com/download-center#community
```

### Creating database

First make sure mongoDB is running:
```
DBPATH=
mongod --dbpath $DBPATH --port 27017 --smallfiles
```
The variants found in the VCF files are processed with VEP and the output is written to JSON.
This is further piped into another python script which adds further annotation and formatting and writes output to JSON.
The JSON is then imported with mongoimport.

#### Importing data from JSON

#### Running pubmedbatch

The pubmedbatch, written by Jing, scores genes based on their pubmed relevance.

#### Running phenogenon

The phenogenon, written by Jing, does an enrichment test per gene and HPO term.

## Running server

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

This code was originally forked from the ExAC browser but has since diverged considerably.


