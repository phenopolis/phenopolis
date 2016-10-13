# Phenopolis

This code has been forked from the ExAC browser for the purpose of the exomes hosted at UCL (UCLex).


Installation
=======


Clone the repository.

Download Phenotips.

Download the Exomiser stand alone.

Install latest version of mongo.

### Creating database

The variants found in the VCF files are processed with VEP and the output is written to JSON.
This is further piped into another python script which adds further annotation and formatting and writes output to JSON.


## Running server

Run phenotips:
```
```
Run mongoDB:
```
```

Run webserver:
```
cd phenopolis
python run_server.py
```

