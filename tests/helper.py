from views import neo4j_driver
from views.neo4j_setup import create_demo_user
from passlib.hash import argon2


def login(app):
    return app.post('/login', data=dict(
        name='demo',
        password='demo123'
    ), follow_redirects=True)


# We won't load from csv file because Neo4j is set up by default to load only from 
# folder <neo4j-home>\import and we don't have access to change this on Travis-CI.
def load_neo4j_test_data():  
    with neo4j_driver.session() as neo4j_session:
        neo4j_session.run("CREATE (a:User {user: {username}, argon_password: {hash}})",
                    {"username": "testSuite", "hash": argon2.hash("demo123")})


def delete_neo4j_test_data():
    with neo4j_driver.session() as neo4j_session:
        neo4j_session.run("MATCH (a:User) WHERE a.user = {username} "
                    "DETACH DELETE a",
                    {"username": "testSuite"})


def create_neo4j_demo_user():  
    with neo4j_driver.session() as neo4j_session:
        create_demo_user(neo4j_session)

def my_patients_neo4j_data():
    user='demo'
    with neo4j_driver.session() as neo4j_session:
        s="""
        MATCH (u:User {user:'%s'})
        MERGE (u)-[r:WRITES]->(p:Person {personId:"person1", gender:"M", score:0.69})
        MERGE (t:Term {termId:"HP:0000505", name:"Visual impairment"})
        MERGE (p)-[:PersonToObservedTerm]->(t)
        MERGE (p)-[:CandidateGene]->(g:Gene {gene_name:"TTLL5"})
        MERGE (gv:GeneticVariant {variantId:"22-38212762-A-G"})
        MERGE (gv14:GeneticVariant {variantId:"14-76201609-C-G"})
        MERGE (p)<-[:HomVariantToPerson]-(gv)
        MERGE (p)<-[:HetVariantToPerson]-(gv14);
        """ % (user)
        print(s)
        result = neo4j_session.run(s)

        # person2
        s="""
        MATCH (u:User {user:'%s'})
        MERGE (u)-[r:WRITES]->(p:Person {personId:"person2", gender:"F", score:0.69})
        MERGE (t:Term {termId:"HP:0000505", name:"Visual impairment"})
        MERGE (p)-[:PersonToObservedTerm]->(t)
        MERGE (p)-[:CandidateGene]->(g:Gene {gene_name:"DRAM2"})
        MERGE (gv:GeneticVariant {variantId:"22-38212762-A-G"})
        MERGE (gv14:GeneticVariant {variantId:"14-76201609-C-G"})
        MERGE (p)<-[:HomVariantToPerson]-(gv)
        MERGE (p)<-[:HetVariantToPerson]-(gv14);
        """ % (user)
        print(s)
        result = neo4j_session.run(s)


 

