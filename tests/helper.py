from views import neo4j_driver
from passlib.hash import argon2


def login(app):
    return app.post('/login', data=dict(
        name='demo',
        password='demo123'
    ), follow_redirects=True)


# We won't load from csv file because Neo4j is set up by default to load only from 
# folder <neo4j-home>\import and we don't have access to change this on Travis-CI.
def load_neo4j_test_data():  
    with neo4j_driver.session() as session:
        session.run("CREATE (a:User {user: {username}, argon_password: {hash}})",
                    {"username": "testSuite", "hash": argon2.hash("demo123")})


def delete_neo4j_test_data():
    with neo4j_driver.session() as session:
        session.run("MATCH (a:User) WHERE a.user = {username} "
                    "DETACH DELETE a",
                    {"username": "testSuite"})


def create_neo4j_demo_user():  
    with neo4j_driver.session() as session:
        results = session.run("MATCH (u:User {user : 'demo'}) RETURN u")
        result = results.single()
        if not result:
            print("Adding user 'demo' to the neo4j database.")
            session.run("CREATE (a:User {user: {username}, argon_password: {hash}})",
                        {"username": "demo", "hash": argon2.hash("demo123")})

