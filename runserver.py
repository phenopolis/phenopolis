import sys
from views import *

# Load default config and override config from an environment variable
global LOCAL

if len(sys.argv)>1 and sys.argv[1]=='SERVER':
    LOCAL=False
    app.config.from_pyfile('../phenopolis.cfg')
else:
    LOCAL=True
    app.config.from_pyfile('../local.cfg')



if __name__ == "__main__":
    # use ssl
    # add some common url. Would be good if can generate the url in real time
    home = ''
    #from OpenSSL import SSL
    # altnerative
    #context = SSL.Context(SSL.SSLv23_METHOD)
    #context = ssl.SSLContext(ssl.PROTOCOL_SSLv23)
    #context = ssl.SSLContext(ssl.PROTOCOL_TLSv1_2)
    #context.load_cert_chain( 'phenotips.crt', 'phenotips.key' )
    #context.use_privatekey_file('ssl/phenotips.cs.ucl.ac.uk.key')
    #context.use_privatekey_file('ssl/phenotips.cs.ucl.ac.uk.key')
    #context.use_privatekey_file('phenotips.key')
    #context.use_certificate_file('phenotips.crt')
    # this is now handled by Apache
    #app.run(host='0.0.0.0',port=8000,threaded=True,debug=True)
    app.run(host='0.0.0.0',port=8000,threaded=True)
    # threaded
    #app.run(threaded=True)
    #app.run(host='127.0.0.1',port=8000, debug = True, ssl_context=context)
    #app.run(host='0.0.0.0', port=8000, ssl_context=context)
    #app.run(host='0.0.0.0', port=8000, debug=True, ssl_context='adhoc')
    #app.run(host='0.0.0.0', port=8000, debug=True)
    #app.run(host='127.0.0.1', port=8000, debug=True)
    #app.run(host='0.0.0.0', port=8000, debug=True, ssl_context=('phenotips.key', 'phenotips.crt'))
    #toolbar=DebugToolbarExtension(app)
    #runner = Runner(app)  # adds Flask command line options for setting host, port, etc.
    #runner.run()





