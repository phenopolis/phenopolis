This code has been forked from the ExAC browser for the purpose of the exomes hosted at UCL (UCLex).

Usage
=======

[http://phenotips.cs.ucl.ac.uk:8000](http://phenotips.cs.ucl.ac.uk:8000)


Access outside of UCL
======

The website is currently only accessible if you are on the UCL wired network but you can access it from outside the network by using an ssh tunnel:

```
ssh -f <username>@jb-gate1.biol.ucl.ac.uk -L1234:phenotips.cs.ucl.ac.uk:8000 -N
```

Then point your browser to [http://localhost:1234](http://localhost:1234).

Installation
=======

### Getting the code

Create a directory to put all this stuff in. This will serve as the parent directory of the actual uclex_browser repository 

    mkdir uclex
    cd uclex

First (as this can run in parallel), get the datasets that the browser uses and put them into an 'uclex_data' directory:


Now clone the repo: 

    git clone https://github.com/pontikos/uclex_browser.git

### Dependencies

Follow these instructions to get Python and Homebrew installed on your Mac:
http://docs.python-guide.org/en/latest/starting/install/osx/

Install MongoDB:

    brew install mongodb
    # or
    sudo port install mongodb

Create a directory to hold your mongo database files: 

    mkdir database

In a separate tab, start the mongo database server:

    mongod --dbpath database

This local server needs to be running at all times when you are working on the site.
You could do this in the background if you want or set up some startup service,
but I think it's easier just to open a tab you can monitor.

Finally, you may want to keep the system in a virtualenv:

    sudo port install py27-virtualenv # Or whatever version

If so, you can create a python virtual environment where the browser will live:

    mkdir uclex_env
    virtualenv uclex_env
    source uclex_env/bin/activate

### Installation

Install the python requirements:

    pip install -r requirements.txt

Note that this installs xBrowse too.

### Setup

At this point, it's probably worth quickly checking out the code structure if you haven't already :)

Now we must load the database from those flat files.
This is a single command, but it can take a while (can take advantage of parallel loads by modifying LOAD\_DB\_PARALLEL\_PROCESSES in uclex.py):

    python manage.py load_db

You won't have to run this often - most changes won't require rebuilding the database.
That said, this is (and will remain) idempotent,
so you can run it again at any time if you think something might be wrong - it will reload the database from scratch.
You can also reload parts of the database using any of the following commands:

    python manage.py load_variants_file
    python manage.py load_dbsnp_file
    python manage.py load_base_coverage
    python manage.py load_gene_models

Then, you need to create a cache for autocomplete and large gene purposes:

    python manage.py create_cache

### Running the site

Note that if you are revisiting the site after a break, make sure your virtualenv is `activate`'d.

You can run the development server with:

    python uclex.py

And visit on your browser:

    http://localhost:5000
    http://localhost:5000/gene/ENSG00000237683
    http://localhost:5000/variant/20-76735-A-T


For testing, you can open up an interactive shell with:

    python manage.py shell

