from __future__ import print_function
# Uncomment to run this module directly. TODO comment out.
# import sys, os
# sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
# End of uncomment.

import unittest
import subprocess
import ConfigParser
import os
from io import StringIO

from neo4j.v1 import GraphDatabase, basic_auth, CypherError
from passlib.hash import argon2

class Neo4jTestCase(unittest.TestCase):

    def setUp(self):
        self.setup_driver()
        self.load_neo4j_test_data()

    def tearDown(self):
        self.delete_neo4j_test_data()

    def setup_driver(self):
        default_password = 'neo4j' # Travis will use a fresh Neo4j, with the default password.

        # Get Neo4j password from file.
        abs_file_path = os.path.dirname(__file__) + '/../local.cfg'
        vfile = StringIO(u'[Section1]\n%s'  % file(abs_file_path).read())
        parser = ConfigParser.ConfigParser(allow_no_value=True)
        parser.readfp(vfile)
        local_password = parser.get('Section1','NEO4J_PWD')

        # Try local_password.
        try:
            self.driver = GraphDatabase.driver("bolt://localhost:7687", auth=basic_auth("neo4j", local_password))
            return
        except:
            pass

        # Try default_password.
        # Password handling from https://github.com/robinedwards/django-neomodel
        self.driver = GraphDatabase.driver("bolt://localhost:7687", auth=basic_auth("neo4j", default_password))
        session = self.driver.session()
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
            result = session.run("MATCH (u:User) WHERE u.user = {username} "
                            "RETURN u.user AS user, u.argon_password AS argon_password",
                            {"username": "testSuite"})
        for record in result:
            assert record['user'] == 'testSuite'
            assert record['argon_password'] == '$argon2i$v=19$m=512,t=2,p=2$n1PK2XvvXcs5h/Aewzjn3A$PxZq8Hwae2EZ4RZX204qsQ'
            assert(argon2.verify('demo123', record['argon_password']))

    # We won't load from csv file because Neo4j is set up by default to load only from 
    # folder <neo4j-home>\import and we don't have access to change this on Travis-CI.
    def load_neo4j_test_data(self):  
        with self.driver.session() as session:
            session.run("CREATE (a:User {user: {username}, argon_password: {hash}})",
                        {"username": "Test Suite", "hash": '$argon2i$v=19$m=512,t=2,p=2$n1PK2XvvXcs5h/Aewzjn3A$PxZq8Hwae2EZ4RZX204qsQ'})

    def delete_neo4j_test_data(self):
        with self.driver.session() as session:
            session.run("MATCH (a:User) WHERE a.user = {username} "
                        "DETACH DELETE a",
                        {"username": "Test Suite"})

if __name__ == '__main__':
    unittest.main()
