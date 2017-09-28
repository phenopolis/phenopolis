from __future__ import print_function
# Uncomment to run this module directly. TODO comment out.
# import sys, os
# sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
# End of uncomment.

import unittest
import subprocess

from neo4j.v1 import GraphDatabase, basic_auth
from passlib.hash import argon2

class Neo4jTestCase(unittest.TestCase):

    def setUp(self):
        self.setup_driver()

    def tearDown(self):
        pass

    def setup_driver(self):
        local = True # TODO LMTW make this automatic
        password = "test" if local else "neo4j"
        self.driver = GraphDatabase.driver("bolt://localhost:7687", auth=basic_auth("neo4j", password))
        session = self.driver.session()

        # Password handling taken from https://github.com/robinedwards/django-neomodel
        try:
            result = session.run("MATCH (a:Person) WHERE a.name = {name} "
                                "RETURN a",
                                {"name": "Crick"})
        except CypherError as ce:
            if 'The credentials you provided were valid, but must be changed before you can use this instance' in str(ce):
                new_password = 'test'
                session.run("CALL dbms.changePassword({password})", {'password': new_password})
                print("New database with no password set, setting password to 'test'")
            else:
                raise ce

    def test_neo4j(self):

        session = self.driver.session()
        session.run("CREATE (a:Person {name: {name}, title: {title}})",
                    {"name": "Crick", "title": "Professor"})
        result = session.run("MATCH (a:Person) WHERE a.name = {name} "
                            "RETURN a.name AS name, a.title AS title",
                            {"name": "Crick"})
        for record in result:
            print("%s %s" % (record["title"], record["name"]))
        assert record["title"] == "Professor"
        assert record["name"] == "Crick"

        # clean up
        session.run("MATCH (a:Person) WHERE a.name = {name} "
                    "DETACH DELETE a",
                    {"name": "Crick"})
        session.close()

    def test_users_data(self):
        with self.driver.session() as session:
            result = session.run("MATCH (u:User) WHERE u.user = {username} "
                            "RETURN u.user AS user, u.argon_password AS argon_password",
                            {"username": "demo"})
        for record in result:
            print('%s %s' % (record['user'], record['argon_password']))
        assert record['user'] == 'demo'
        assert record['argon_password'] == '$argon2i$v=19$m=512,t=2,p=2$n1PK2XvvXcs5h/Aewzjn3A$PxZq8Hwae2EZ4RZX204qsQ'
        assert(argon2.verify('demo123', record['argon_password']))

    def test_commit_date(self): # TODO LMTW add to home.py as alternative to tag version no.
        try:
            commit_date = subprocess.check_output(['git', 'log', '-1', '--format=%cd', '--date=local'])
        except:
            commit_date = None
        print('Commit date is:- ', commit_date)

    def do_not_test_neo4j_data(self):  # TODO LMTW move to load data and load from .csv
        session = self.driver.session()
        #session.run("CREATE (a:User {user: {username}, argon_password: {hash}})",
        #            {"username": "neo4j", "hash": '$argon2i$v=19$m=512,t=2,p=2$bS1FyNn7fy9FyBnj/F/LGQ$55VPLCPvU4VgKa7rOneppA'})
        session.run("CREATE (a:User {user: {username}, argon_password: {hash}})",
                    {"username": "demo", "hash": '$argon2i$v=19$m=512,t=2,p=2$n1PK2XvvXcs5h/Aewzjn3A$PxZq8Hwae2EZ4RZX204qsQ'})




# TODO LMTW 28Sept17 where I am up to - 
# These tests pass - I can get user out of my neo4j db.
# Get check_auth working. Need to add a demo user to my neo4j first.

if __name__ == '__main__':
    unittest.main()
