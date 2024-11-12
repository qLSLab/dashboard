from neo4j import GraphDatabase
import streamlit as st

# @st.cache_data(suppress_st_warning=True)
def get_neo4j_data(uri='bolt://localhost:11007', user='neo4j', password='1', query='MATCH (n) RETURN n Limit 2'):
    driver = GraphDatabase.driver(uri, auth=(user, password))
    nodes = set()
    edges = set()
    with driver.session() as session:
        result = session.run(query)
        # Transforming the result into the required format
        data = []
        for record in result:
            path = record["relacion"]
            for node in path.nodes:
                node_properties = tuple(node.items())  # Convert dictionary to tuple
                nodes.add((node.id, node.labels, node_properties))
            for relationship in path.relationships:
                relationship_properties = tuple(relationship.items())  # Convert dictionary to tuple
                edges.add((relationship.id, relationship.type, relationship_properties, relationship.start_node.id, relationship.end_node.id))
    driver.close()
    return list(nodes), list(edges)
