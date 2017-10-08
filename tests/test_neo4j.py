from __future__ import print_function
# Uncomment to run this module directly. TODO comment out.
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
# End of uncomment.

import unittest
import subprocess
import ConfigParser
import os
from io import StringIO
import runserver

from neo4j.v1 import GraphDatabase, basic_auth, CypherError
from passlib.hash import argon2

class Neo4jTestCase(unittest.TestCase):

    def setUp(self):
        self.setup_driver()
        self.load_neo4j_test_data()
        runserver.app.config['TESTING'] = True
        self.app = runserver.app.test_client()

    def tearDown(self):
        self.delete_neo4j_test_data()
        self.driver.close()

    def setup_driver(self):
        default_password = 'neo4j' # Travis will use a fresh Neo4j, with the default password.

        # Get Neo4j password, host and port from file.
        abs_file_path = os.path.dirname(__file__) + '/../local.cfg'
        vfile = StringIO(u'[Section1]\n%s'  % file(abs_file_path).read())
        parser = ConfigParser.ConfigParser(allow_no_value=True)
        parser.readfp(vfile)
        local_password = parser.get('Section1','NEO4J_PWD').strip("'")
        host = parser.get('Section1','NEO4J_HOST').strip("'")
        port = str(parser.get('Section1','NEO4J_PORT'))
        uri = "bolt://"+host+":"+port

        # Try local_password.
        try:
            self.driver = GraphDatabase.driver(uri, auth=basic_auth("neo4j", local_password))
            return
        except:
            pass

        # Try default_password.
        # Password handling from https://github.com/robinedwards/django-neomodel
        self.driver = GraphDatabase.driver(uri, auth=basic_auth("neo4j", default_password))
        with self.driver.session() as session:
            try:
                result = session.run("MATCH (a:Person) WHERE a.name = {name} RETURN a", {"name": "Crick"})
            except CypherError as ce:
                if 'The credentials you provided were valid, but must be changed before you can use this instance' in str(ce):
                    session.run("CALL dbms.changePassword({password})", {'password': local_password})
                    print("New database with no password set, setting password to '", local_password, "'.")
                    session.close()
                else:
                    raise ce

    def test_users_data(self):
       with self.driver.session() as session:
            # Test unknown user
            results = session.run("MATCH (u:User {user : 'xxx'}) RETURN u.user AS user, u.argon_password AS argon_password")
            result = results.single()
            assert(not result)

            # Test known user
            results = session.run("MATCH (u:User {user : 'testSuite'}) RETURN u.user AS user, u.argon_password AS argon_password")
            result = results.single()
            assert(result)
            assert result['user'] == 'testSuite'
            assert result['argon_password'] == '$argon2i$v=19$m=512,t=2,p=2$n1PK2XvvXcs5h/Aewzjn3A$PxZq8Hwae2EZ4RZX204qsQ'
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

    # We won't load from csv file because Neo4j is set up by default to load only from 
    # folder <neo4j-home>\import and we don't have access to change this on Travis-CI.
    def load_neo4j_test_data(self):  
        with self.driver.session() as session:
            session.run("CREATE (a:User {user: {username}, argon_password: {hash}})",
                        {"username": "testSuite", "hash": '$argon2i$v=19$m=512,t=2,p=2$n1PK2XvvXcs5h/Aewzjn3A$PxZq8Hwae2EZ4RZX204qsQ'})

    def delete_neo4j_test_data(self):
        with self.driver.session() as session:
            session.run("MATCH (a:User) WHERE a.user = {username} "
                        "DETACH DELETE a",
                        {"username": "testSuite"})

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
