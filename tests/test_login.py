from __future__ import print_function
# Uncomment to run this module directly. TODO comment out.
#import sys, os
#sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
# End of uncomment.

import unittest
import subprocess
import runserver
from views import neo4j_driver
import helper

from passlib.hash import argon2

class Neo4jTestCase(unittest.TestCase):

    def setUp(self):
        helper.load_neo4j_test_data()
        runserver.app.config['TESTING'] = True
        self.app = runserver.app.test_client()


    def tearDown(self):
        helper.delete_neo4j_test_data()


    def test_users_data(self):
       with neo4j_driver.session() as neo4j_session:
            # Test unknown user
            results = neo4j_session.run("MATCH (u:User {user : 'xxx'}) RETURN u.user AS user, u.argon_password AS argon_password")
            result = results.single()
            assert(not result)

            # Test known user
            results = neo4j_session.run("MATCH (u:User {user : 'testSuite'}) RETURN u.user AS user, u.argon_password AS argon_password")
            result = results.single()
            assert(result)
            assert result['user'] == 'testSuite'
            assert(argon2.verify('demo123', result['argon_password']))


    def test_login_logout(self):
        rv = self.login('Testx', 'demo123')
        assert rv.status_code == 401
        assert 'Invalid Credentials. Please try again.' in rv.data
        rv = self.login('testSuite', 'demo123x')
        assert rv.status_code == 401
        assert 'Invalid Credentials. Please try again' in rv.data
        rv = self.login('testSuite', 'demo123')
        assert rv.status_code == 200
        assert 'Authenticated' in rv.data
        rv = self.logout()
        assert rv.status_code == 200
        assert 'Please login' and 'username' and 'password' in rv.data


    def test_change_password(self):
        rv = self.login('testSuite', 'demo123')
        assert rv.status_code == 200
        assert 'Authenticated' in rv.data

        rv = self.change_password('testSuite', 'demo123', 'demo456')
        assert rv.status_code == 200
        print(rv.data)
        assert 'Password for username \'testSuite\' changed' in rv.data

        rv = self.login('testSuite', 'demo456')
        assert rv.status_code == 200

        rv = self.login('testSuite', 'demo123')
        assert rv.status_code == 401

        rv = self.change_password('testSuite', 'demo456', 'demo123')
        assert rv.status_code == 200

        rv = self.change_password('x', 'demo123', 'demo456')
        assert rv.status_code == 401

        rv = self.change_password('testSuite', 'x', 'demo456')
        assert rv.status_code == 401


    def login(self, username, password):
        return self.app.post('/login', data=dict(
            name=username,
            password=password
        ), follow_redirects=True)


    def logout(self):
        return self.app.get('/logout', follow_redirects=True)


    def change_password(self, username, password, new_pass_1):
        return self.app.post('/change_password', data=dict(
            change_pwd_name=username,
            current_password=password,
            new_password_1=new_pass_1,
        ), follow_redirects=True)

if __name__ == '__main__':
    unittest.main()
