from views import *
from lookups import *
import requests
import re
from utils import *
import itertools
from config import config
if config.IMPORT_PYSAM_PRIMER3:
    import pysam
import csv
#hpo lookup
import random
from flask import Response, request
import os
from werkzeug.datastructures import Headers
import re

@app.route('/bam_viewer/')
def bam_viewer():
    return render_template('igv_viewer.html')

@app.route('/read_viz/bam/<sample>')
def read_viz(sample):
    BAM_FILES=app.config['BAM_FILES']
    print(request.method)
    headers=Headers()
    #headers.add('Content-Type','application/octet-stream')
    headers.add('Content-Transfer-Encoding','binary')
    #Date:Wed, 06 Jul 2016 17:19:52 GMT
    #ETag:"flask-1446310274.0-12661331-649139018"
    #Expires:Thu, 07 Jul 2016 05:19:52 GMT
    #Keep-Alive:timeout=5, max=93
    #Last-Modified:Sat, 31 Oct 2015 16:51:14 GMT
    headers.add('Accept-Ranges', 'bytes')
    #Server:Apache/2.4.12 (Red Hat) mod_wsgi/3.4 Python/2.7.8
    headers.add('X-Frame-Options','SAMEORIGIN')
    if sample=='gencode.v19.sorted.bed':
        bamfile=BAM_FILES+'/gencode.v19.sorted.bed'
    elif sample=='gencode.v19.sorted.bed.idx':
        bamfile=BAM_FILES+'/gencode.v19.sorted.bed.idx'
    elif sample.endswith('.bai'):
        bamfile=BAM_FILES+'/%s.bam.bai' % sample
    else:
        bamfile=BAM_FILES+'/%s.bam' % sample
    size = os.path.getsize(bamfile)
    print(size)
    status = 200
    begin  = 0
    end = size-1
    if request.headers.has_key("Range") and request.method=='GET':
        print(request.headers['Range'])
        headers.add('Accept-Ranges','bytes')
        ranges = re.findall(r"\d+", request.headers["Range"])
        begin  = int( ranges[0] )
        if len(ranges)>1: end = int( ranges[1] )
        headers.add('Content-Range','bytes %s-%s/%s' % (str(begin),str(end),size) )
        headers.add('Content-Length',str((end-begin)+1))
        with file(bamfile,'rb') as f:
            f.seek(begin)
            data=f.read(end-begin)
        print(len(data))
        response = Response( data, status=206, mimetype="application/octet-stream", headers=headers, direct_passthrough=True)
    else:
        if request.method=='HEAD':
            headers.add('Content-Length',size)
            response = Response( '', status=200, mimetype="application/octet-stream", headers=headers, direct_passthrough=True)
        elif request.method=='GET':
            response = Response( file(bamfile), status=200, mimetype="application/octet-stream", headers=headers, direct_passthrough=True)
    #Add mimetype   
    response.cache_control.public  = True
    response.make_conditional(request)
    return response




def read_viz2():
    print(sample)
    print(region)
    from subprocess import call
    tmpfile=subprocess.Popen('mktemp', shell=True, stdout=subprocess.PIPE).stdout.read().strip()+'.bam'
    print(tmpfile)
    print(subprocess.Popen("samtools view -b %s/%s_sorted_unique.bam %s > %s" % (BAM_FILES,sample,region, tmpfile), shell=True, stdout=subprocess.PIPE).stdout.read())
    subprocess.Popen('samtools index %s'%tmpfile).stdout.read()
    

