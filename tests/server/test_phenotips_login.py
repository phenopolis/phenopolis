
import unittest

# TODO LMTW remove - for debugging
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))
import phenotips_python_client
# end TODO LMTW remove
from phenotips_python_client import PhenotipsClient
import httpie
from binascii import b2a_base64, a2b_base64
import requests


class PhenotipsLoginTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_login(self):
       

        print('doing session')
        url = 'http://localhost:8080/rest/patients/P0000001/permissions'
        #s = requests.Session()
        #s.auth('Admin','admin')
        ##r1 = s.get('http://localhost:8080/rest/patients/P0000001/permissions', 'Admin:admin')
        #r2 = s.get('http://localhost:8080/rest/patients/P0000001/permissions')

        s = requests.Session()
        s.auth = requests.auth.HTTPDigestAuth('Admin', 'admin')
        r2 = s.get('http://localhost:8080/rest/patients/P0000001/permissions')



        encoded_auth=b2a_base64('Admin:admin').strip()
        headers={'Authorization':'Basic %s'%encoded_auth, 'Accept':'application/json'}
        r=requests.get(url, headers=headers)
        assert(r.status_code == 200) 

        s = requests.Session()
        r = s.get(url, headers=headers)
        assert(r.status_code == 200) 

        small_headers={'Accept':'application/json'}
        r = s.get(url, headers=small_headers)
        status_code = r.status_code 

        new_url='http://localhost:8080/rest/patients?start=0&number=5'
        r = s.get(new_url, headers=small_headers)
        status_code = r.status_code 

        temp =''


    def dummy(self):
        conn=PhenotipsClient()

        # Test unauthorised credentials.
        #auth='%s:%s' % ('Adminx','admin') 
        #rv=conn.get_patient(auth=auth)
        #assert(rv == 401)

        #auth='%s:%s' % (session['user'],session['password2'],)
        auth='%s:%s' % ('Admin','admin') 
        all_patients=conn.get_patient(auth=auth).get('patientSummaries',[])
        all_eids=[p['eid'] for p in all_patients if p['eid']]
        total=len(all_eids)
        print('TOTAL NUMBER OF PATIENTS',total)
        assert(total == 1)

        ##
        #url = 'http://localhost:8080'
        #payload = 'Admin:admin'
        #session = httpie(url, payload)
        #http --session=Admin -a Admin:admin http://localhost:8080

#       $ http --session=logged-in -a username:password httpbin.org/get API-Key:123
#       $ http --session=logged-in httpbin.org/headers
        

        ##
        host='localhost'
        port='8080',
        site='localhost:8080'#'%s:%s'%(host,port,)
        auth=b2a_base64('Admin:admin').strip()
        headers={'Authorization':'Basic %s'%auth, 'Accept':'application/json'}
        start = 0
        number = 5
        url='http://%s/rest/patients?start=%d&number=%d' % (site,start,number)
        r=requests.get(url, headers=headers)
        print(r.status_code)

        #headers = {'--session=logged-in', '-a Admin:admin'}
        #r=requests.get(url, headers=headers)
        #print(r.status_code)

        print('doing session')
        s = requests.session()
        s.post(url, 'Admin:admin')
        #logged in! cookies saved for future requests.
        r2 = s.get(url)
        print(r2.status_code)
        try:
            r2=r2.json()
            all_patients=r2.get('patientSummaries',[])
            all_eids=[p['eid'] for p in all_patients if p['eid']]
            total=len(all_eids)
            print('TOTAL NUMBER OF PATIENTS',total)
        except:
            print('fail')

        headers={'Accept':'application/json; application/xml'}
        ID = 'P0000001'
        s = requests.session()
        s.post(url, 'Admin:admin')
        #r=s.get('http://%s/rest/patients/%s/permissions' % (site,ID), headers=headers) 
        r=s.get('http://%s/rest/patients/%s/permissions' % (site,ID))
        print(r.status_code)
        #assert(r.status_code == 200) #401

        # Invalid credentials
        s.post(url, 'Adminx:admin')
        r=s.get('http://%s/rest/patients/%s/permissions' % (site,ID), headers=headers) 
        #assert(r.status_code == 404)
        temp=''

        #http://localhost:8080/rest/patients/P0000001/permissions
        #http --session=JohnDoe -a JohnDoe:johndoespassword http://localhost:8080
        #http --session=AdminSession -a Admin:admin http://localhost:8080/rest/patients/P0000001/permissions

        #httpie('--session=AdminSession -a Admin:admin http://localhost:8080/rest/patients/P0000001/permissions')

        try:
            #version_number = subprocess.check_output(['git', 'describe', '--exact-match'])
            s2 = subprocess.check_output(['httpie', '-a Admin:admin http://localhost:8080/rest/patients/P0000001/permissions', '--session=AdminSession'])
        except:
            temp = None

        try:
            s_invalid = subprocess.check_output(['httpie', '-a Adminx:admin http://localhost:8080/rest/patients/P0000001/permissions', '--session=AdminSession'])
        except:
            temp = None

            
        s2.get('http://localhost:8080/rest/patients/P0000001/permissions')


if __name__ == '__main__':
    unittest.main()
