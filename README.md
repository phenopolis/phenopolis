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
```

The variants found in the VCF files are processed with VEP and the output is written to JSON.
This is further piped into another python script which adds further annotation and formatting and writes output to JSON.
The JSON is then imported with mongoimport.


## Running server

Run Phenotips standalone:
```
```

Run webserver:
```
cd phenopolis
python run_server.py
```


### Acknowledgment

This code was originally forked from the ExAC browser.

