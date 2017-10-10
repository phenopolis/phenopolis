import ConfigParser
import os
from io import StringIO
from neo4j.v1 import GraphDatabase, basic_auth, CypherError


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
    with driver.session() as session:
        try:
            result = session.run("MATCH (a:Person) WHERE a.name = {name} RETURN a", {"name": "Crick"})

        except CypherError as ce:
            if 'The credentials you provided were valid, but must be changed before you can use this instance' in str(ce):
                session.run("CALL dbms.changePassword({password})", {'password': local_password})
                print("New database with no password set, setting password to '", local_password, "'.")
                session.close()
            else:
                raise ce
    return driver
