# This is a chatbot app that, using prompt engineering, will implement the RAG technique to answer questions about metabolic models

# Import libraries
import numpy as np
import streamlit as st
from langchain.chat_models import ChatOpenAI
from langchain.memory import ConversationBufferMemory
from langchain.chains import LLMChain
from langchain.prompts import PromptTemplate
from decouple import config
from langchain.graphs import Neo4jGraph
from langchain.chains import GraphCypherQAChain

from graph_message import GraphVisualizer
from neo4j_utils import get_neo4j_data

class Message:
    def __init__(self, role, chart_data=None, content="", graph_viz=False, query=None, key=None):
        self.role = role
        self.content = content
        self.chart_data = chart_data
        self.graph_viz = graph_viz
        self.query = query
        self.key = key

    def render(self, key=None):
        if self.graph_viz is True:
            graph_message = GraphVisualizer(self.query)
            graph_message.render()
        else:
            st.write(self.content)

graph = Neo4jGraph(
    url="bolt://localhost:7687",
    username="neo4j",
    password="11111111",
)

CYPHER_GENERATION_TEMPLATE = """
You are an expert Neo4j Developer translating user questions into Cypher to answer questions about a specific metabolic model that is represented by a Neo4j graph.
Convert the user's question based on the schema.

Instructions:
Use only the provided relationship types and properties in the schema.
Do not use any other relationship types or properties that are not provided.
Do not use any other nodes that are not provided.
The answer should be a Cypher query that returns the answer to the question.
The specie label is also named as metabolic, metabolites, molecule, or compound.
The GPR label represents the boolean logic of how the genes are related to the molecules.
Limit all the Cypher queries to 300 results (This is important).

Some cofactors are also represented as metabolites, but they are not part of the metabolic network.
Some examples of cofactors are:
'44215', 'NAD', 'C00003', '15846', '57540', 'C00001', 'WATER', '15377', '16908', 'C00004', 'NADH', '57945', '25805', 'OXYGEN-MOLECULE', '15379', 'C00007', 'NADP', '58349', 'B-HEP-1:5', 'C00667', '18009', 'C00006', '16474', '77312', 'NADPH', 'CPD-16005', 'C20745', '77177', '57783', 'C00005', 'CPD-16005', '15996', 'C00044', '37565', 'GTP', 'C00010', '57287', 'CO-A', '15346', 'CARBON-DIOXIDE', 'C00011', 'CARBON-DIOXIDE', '16526', 'AMMONIUM', 'C01342', 'AMMONIUM', 'AMMONIA', 'C00014', 'AMMONIA', '16134', '28938', 'C00002', '15422', '22258', 'ATP', '30616', '456216', '73342', 'ADP', '22251', 'C00008', 'G11113', '22252', '16761'
When you make a query of a path, you may exclude the cofactors.

Schema: {schema}

Cypher examples:
# Find the nodes that are connected to the node with id 'R_TMDK1' by a reaction relationship:
MATCH relacion = (n)-[r]-(p) 
WHERE n.id = 'R_TMDK1'
RETURN relacion

# Count the more frequent genes, return the 10 most frequent:
match (n:Gene)-[arcs]->() 
return n,count(arcs) 
order by count(arcs) desc 
limit 10

# Starting from a gene named G_HGNC__58__1374, obtain all the adjacent nodes to the reactions related:
MATCH relacionGen = (n)-[]-(r) 
WHERE n.id = "G_HGNC__58__1374"
MATCH relationReaction = (r)-[]-()
RETURN relationReaction

# Provide the shortest path between the reactions R_ATPtm and R_OMPDC:
Match (r1:Reaction)-[*0..1]->()
WHERE r1.id = 'R_ATPtm'
Match (r2:Reaction )<-[*0..1]-()
WHERE r2.id = 'R_OMPDC'
match p =  shortestPath((r1)-[*]->(r2))
with nodes(p) as nodes
unwind nodes as n
with n
where 'Reaction' IN LABELS(n)
return (n)-[]-()
Limit 300

# Show me the gene interaction relative to the reaction named FBPALDO:
MATCH p = (:GPR)-[:GPR_edge]->(reaction:Reaction)<-[:is_a_gene_that_works_in]-(gene:Gene)-[:GPR_edge]->(:GPR)
WHERE reaction.name = 'FBPALDO'
RETURN p

Note: Do not include any explanations or apologies in your responses.
Do not respond to any questions that might ask anything else than for you to construct a Cypher statement.
Do not include any text except the generated Cypher statement.

The question is:
Question: {question}
"""

# Creating the prompt template that will guide the AI model's responses
cypher_generation_prompt = PromptTemplate(
    template=CYPHER_GENERATION_TEMPLATE,
    input_variables=["schema", "question"],
)

# Initializing the necessary objects to interact with OpenAI and manage conversation memory
llm = ChatOpenAI(model_name="gpt-4o", openai_api_key=config("OPENAI_API_KEY"))  # ChatOpenAI instance with OpenAI API key

cypher_chain = GraphCypherQAChain.from_llm(
    llm,
    graph=graph,
    cypher_prompt=cypher_generation_prompt,
    verbose=True,
    return_intermediate_steps=True
)

# Configuring the Streamlit page
st.set_page_config(
    page_title="RISE",
    page_icon="🍂🍃",
    layout="wide"
)

# Initializing the session state if it does not exist, to store message history
if "messages" not in st.session_state:
    st.session_state.messages = [
        Message("assistant", content="Hello, how can I help you?"),
    ]

if "counter" not in st.session_state:
    st.session_state["counter"] = 0

# Render all messages in the stack each time the script runs
for message in st.session_state.messages:
    with st.chat_message(message.role):
        message.render()

# Initializing a list to store chat history
chat_history = []

# Collecting user input
user_prompt = st.chat_input()

# If there is user input, process and display the AI model's response
if user_prompt:
    user_message = Message("user", content=user_prompt)
    st.session_state.messages.append(user_message) # Adding the user's message to the session state
    with st.chat_message("user"):
        st.write(user_prompt)  # Displaying the user's message
    
    # Adding the user's input to the chat history
    chat_history.append(f'Human: {user_prompt}')
    # Incrementing the counter in the session state
    st.session_state["counter"] += 1
    with st.chat_message("assistant"):
        with st.spinner("Loading..."):  # Displaying a spinner while waiting for the AI model's response
            print(user_prompt)

            # Feeding the chat history and the user's question to llm_chain to get the AI model's response
            # ai_response = llm_chain.predict(chat_history=" ".join(chat_history), question=user_prompt)
            
            result = cypher_chain(user_prompt)
            ai_response = result["result"]

            # Generates the subgraph visualization
            query = result['intermediate_steps'][0]['query']
            chart_message = Message(role="assistant", graph_viz=True, query=query) # create the graph visualization
            st.session_state.messages.append(chart_message)
            chart_message.render()
            
            # Adding the assistant's response to the chat history
            predicted_message = Message("assistant", content=ai_response)
            st.session_state.messages.append(predicted_message)
            st.write(ai_response)  # Displaying the AI model's response
            
            # Adding the assistant's response to the chat history
            chat_history.append(f'AI: {ai_response}')
            print(st.session_state.messages)
