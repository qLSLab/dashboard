import streamlit as st
from pyvis.network import Network
import pandas as pd
from py2neo import Graph, Node, Relationship, Path
import re

class GraphVisualizer:
    def __init__(self, neo4j_query):
        if neo4j_query:
            # print("Loading data from Neo4j... phase 1")
            self.nodes_data, self.edges_data = self.load_data_from_neo4j_query(query=neo4j_query)
        else:
            self.nodes_data, self.edges_data = self.create_example_graph()

    def load_data_from_neo4j_query(self, uri='bolt://localhost:11007', user='neo4j', password='1', query='MATCH (n) RETURN n Limit 2'):
        graph = Graph(uri, user=user, password=password)
        nodes = []
        edges = []
        node_ids = set()  # A set to keep track of node IDs and avoid duplicates

        results = graph.run(query)
        for record in results:
            # print("record ", record)
            # Iterate over each part of the record
            for part in record.values():
                # Check if the part is a Node
                if isinstance(part, Node):
                    node_data = self.extract_node_data(part)
                    if node_data['id'] not in node_ids:
                        nodes.append(node_data)
                        node_ids.add(node_data['id'])
                # Check if the part is a Relationship
                elif isinstance(part, Relationship):
                    self.extract_edge_data(part, edges)
                # Check if the part is a Path
                elif isinstance(part, Path):
                    for node in part.nodes:
                        node_data = self.extract_node_data(node)
                        if node_data['id'] not in node_ids:
                            nodes.append(node_data)
                            node_ids.add(node_data['id'])
                    for relationship in part.relationships:
                        self.extract_edge_data(relationship, edges)

        # Create Pandas DataFrames for nodes and edges
        nodes_data = pd.DataFrame(nodes).drop_duplicates(subset='id')
        edges_data = pd.DataFrame(edges)

        return nodes_data, edges_data

    def extract_edge_data(self, part, edges):
        # Add edges to the list
        if isinstance(part, Relationship):
            edges.append(self.format_edge_data(part))
        elif isinstance(part, Path):
            for relationship in part.relationships:
                edges.append(self.format_edge_data(relationship))

    def format_edge_data(self, relationship):
        return {
            'source': relationship.start_node['id'],
            'target': relationship.end_node['id']
           # 'title': relationship.type()
        }

    def extract_node_data(self, node):
        # Extract node data
        node_data = {
            'id': node['id'],
            'label': node['name'],  # Assuming 'name' is a common property
            'title': self.format_properties(node),
            'group': next(iter(node.labels))  # Take one of the labels as the group
        }
        return node_data

    def is_url(self, value):
        # Use a regular expression to detect URLs
        url_pattern = re.compile(
            r'^(?:http|ftp)s?://'  # http:// or https://
            r'(?:(?:[A-Z0-9](?:[A-Z0-9-]{0,61}[A-Z0-9])?\.)+(?:[A-Z]{2,6}\.?|[A-Z0-9-]{2,}\.?)|'  # domain...
            r'localhost|'  # localhost...
            r'\d{1,3}\.\d{1,3}\.\d{1,3}\.\d{1,3}|'  # ...or ipv4
            r'\[?[A-F0-9]*:[A-F0-9:]+\]?)'  # ...or ipv6
            r'(?::\d+)?'  # optional port
            r'(?:/?|[/?]\S+)$', re.IGNORECASE)
        return re.match(url_pattern, value) is not None

    def format_properties(self, node):
        # Format properties to include them in the title
        properties_html = '<br>'.join(
            f'<b>{key}</b>: <a href="{value}" target="_blank">{value}</a>' if self.is_url(value) else f'<b>{key}</b>: {value}'
            for key, value in node.items()
        )
        
        # Check if "href" is present in the text string
        if 'href' not in properties_html:
            properties_html += '<br><a href=""></a>'

        return properties_html  # If "href" is found in the text, allow the tooltip to be formatted in HTML

    @staticmethod
    def create_example_graph():
        nodes_data = pd.DataFrame({
            'id': [1, 2, 3],
            'label': ['Node1', 'Node2', 'Node3'],
            'title': ['<b>Node 1 data</b><br><b>http://example.com</b><br>bqbiol: http://biomodels.net/biology-qualifiers/<br><b>upperFluxBound: cobra_default_ub</b>', '<b>Node 2 data</b><br><a href="http://example.com">Link 2</a>', '<b>Node 3 data</b><br><a href="http://example.com">Link 3</a>'],
            'group': [1, 2, 2]
        })

        edges_data = pd.DataFrame({
            'source': [1, 1, 2],
            'target': [2, 3, 3],
            'title': ['Relation 1-2', 'Relation 1-3', 'Relation 2-3']
        })

        return nodes_data, edges_data

    def create_graph(self, height=400):
        net = Network(width="100%", bgcolor="#ffffff", font_color="black", notebook=True)

        for index, row in self.nodes_data.iterrows():
            net.add_node(row['id'], label=row['label'], title=row['title'], group=row['group'])

        for index, row in self.edges_data.iterrows():
            net.add_edge(row['source'], row['target']) #, title=row['title']

        return net

    def render(self):
        height = 300
        net = self.create_graph(height)
        net.show("mygraph.html") 
        HtmlFile = open("mygraph.html", 'r', encoding='utf-8')
        source_code = HtmlFile.read()
        source_code = '<div style="height: {}px;">{}</div>'.format(height, source_code)
        st.components.v1.html(source_code, height=height, scrolling=False)
