import sys
import ConfigParser
import os
from io import StringIO
from neo4j.v1 import GraphDatabase, basic_auth, CypherError
from passlib.hash import argon2


def setup_neo4j_driver(host, port, password):
    default_password = 'neo4j' # Travis will use a fresh Neo4j, with the default password.
    local_password = password
    uri = "bolt://"+host+":"+str(port)

    # Try local_password.
    try:
        driver = GraphDatabase.driver(uri, auth=basic_auth("neo4j", local_password))
        return driver
    except:
        pass

    # Try default_password.
    # Password handling from https://github.com/robinedwards/django-neomodel
    driver = GraphDatabase.driver(uri, auth=basic_auth("neo4j", default_password))
    with driver.session() as neo4j_session:
        try:
            result = neo4j_session.run("MATCH (a:Person) WHERE a.name = {name} RETURN a", {"name": "Crick"})

        except CypherError as ce:
            if 'The credentials you provided were valid, but must be changed before you can use this instance' in str(ce):
                neo4j_session.run("CALL dbms.changePassword({password})", {'password': local_password})
                print("New database with no password set, setting password to '", local_password, "'.")
                neo4j_session.close()
            else:
                raise ce
    return driver

def create_demo_user(neo4j_session):
    results = neo4j_session.run("MATCH (u:User {user : 'demo'}) RETURN u")
    result = results.single()
    if not result:
        print("Adding user 'demo' to the neo4j database.")
        neo4j_session.run("CREATE (a:User {user: {username}, argon_password: {hash}})",
                {"username": "demo", "hash": argon2.hash("demo123")}) 


# For use in easy_install.sh to set up a demo user.
if __name__ == '__main__':
    uri = sys.argv[1] 
    password = sys.argv[2] 
    driver = GraphDatabase.driver(uri, auth=basic_auth("neo4j", password))
    with driver.session() as neo4j_session:  
        create_demo_user(neo4j_session)
