# Import libraries
import xml.dom.minidom as xmlFile
import xml.dom.minidom
from bs4 import BeautifulSoup
from neo4j import GraphDatabase
import shutil
import os
from decouple import config
import re
from lxml import etree
from neo4j import GraphDatabase, exceptions


database_url = config('DATABASE_IMPORT_FOLDER_URL')

class BatchProcessor:
    def __init__(self, session, batch_size=500):
        """
        Inicializa el procesador de lotes.

        Args:
            session (neo4j.Session): La sesión activa de Neo4j.
            batch_size (int): Número de operaciones por lote antes de ejecutar.
        """
        self.session = session
        self.batch_size = batch_size
        self.gpr_nodes = []
        self.gpr_edges = []
        self.gene_reaction_edges = []

    def add_gpr_node(self, node):
        self.gpr_nodes.append(node)
        if len(self.gpr_nodes) >= self.batch_size:
            self.execute_gpr_nodes()

    def add_gpr_edge(self, edge):
        self.gpr_edges.append(edge)
        if len(self.gpr_edges) >= self.batch_size:
            self.execute_gpr_edges()

    def add_gene_reaction_edge(self, edge):
        self.gene_reaction_edges.append(edge)
        if len(self.gene_reaction_edges) >= self.batch_size:
            self.execute_gene_reaction_edges()

    def execute_gpr_nodes(self):
        if not self.gpr_nodes:
            return
        query = '''
        UNWIND $nodes AS node
        CREATE (n:GPR {
            id: node.id,
            type: node.type,
            name: node.name,
            reactionId: node.reactionId
        })
        '''
        try:
            with self.session.begin_transaction() as tx:
                tx.run(query, nodes=self.gpr_nodes)
                tx.commit()
        except exceptions.Neo4jError as e:
            print(f"Error executing GPR nodes batch: {e}")
            # Manejo adicional si es necesario
        finally:
            self.gpr_nodes.clear()

    def execute_gpr_edges(self):
        if not self.gpr_edges:
            return
        query = '''
        UNWIND $edges AS rel
        MATCH (n {id: rel.from_id}), (p {id: rel.to_id})
        CREATE (n)-[:GPR_edge]->(p)
        '''
        try:
            with self.session.begin_transaction() as tx:
                tx.run(query, edges=self.gpr_edges)
                tx.commit()
        except exceptions.Neo4jError as e:
            print(f"Error executing GPR edges batch: {e}")
            # Manejo adicional si es necesario
        finally:
            self.gpr_edges.clear()

    def execute_gene_reaction_edges(self):
        if not self.gene_reaction_edges:
            return
        query = '''
        UNWIND $edges AS rel
        MATCH (n {id: rel.gene_id}), (p {id: rel.reaction_id})
        CREATE (n)-[:is_a_gene_that_works_in]->(p)
        '''
        try:
            with self.session.begin_transaction() as tx:
                tx.run(query, edges=self.gene_reaction_edges)
                tx.commit()
        except exceptions.Neo4jError as e:
            print(f"Error executing gene-reaction edges batch: {e}")
            # Manejo adicional si es necesario
        finally:
            self.gene_reaction_edges.clear()

    def execute_all(self):
        """
        Ejecuta cualquier operación restante que no haya sido procesada en lotes.
        """
        self.execute_gpr_nodes()
        self.execute_gpr_edges()
        self.execute_gene_reaction_edges()

# ---------------

def cargarLogicaGenes(path_to_xml_file, session):
    """
    Procesa un archivo SBML XML, extrayendo las reacciones y sus genes asociados,
    y actualiza la base de datos Neo4j en consecuencia.

    Args:
        session (neo4j.Session): La sesión activa de Neo4j.
        path_to_xml_file (str): La ruta al archivo SBML XML.
    """
    print(path_to_xml_file)
    # Parsear el archivo XML usando lxml para un mejor manejo de namespaces
    try:
        tree = etree.parse(path_to_xml_file)
    except Exception as e:
        print(f"Error parsing XML file: {e}")
        return
    root = tree.getroot()

    # Extraer namespaces y prefijos
    namespaces = root.nsmap  # Diccionario de prefix: URI
    # Mapeo inverso URI: prefix (excluyendo None)
    uri_to_prefix = {v: k for k, v in namespaces.items() if k is not None}

    # Función para identificar el prefijo para una etiqueta dada usando regex
    def get_prefix(tag_name):
        """
        Identifica el prefijo utilizado para una etiqueta dada en el archivo XML.

        Args:
            tag_name (str): El nombre local de la etiqueta (sin prefijo).

        Returns:
            str or None: El prefijo asociado con la etiqueta, o None si no se encuentra.
        """
        # Compilar un patrón regex para encontrar la etiqueta con cualquier prefijo
        pattern = re.compile(fr'<(\w+):{tag_name}')
        try:
            with open(path_to_xml_file, 'r', encoding='utf-8') as file:
                content = file.read()
                match = pattern.search(content)
                if match:
                    return match.group(1)
                else:
                    return None
        except Exception as e:
            print(f"Error reading XML file for prefix '{tag_name}': {e}")
            return None

    # Identificar prefijos para las etiquetas de interés
    reaction_prefix = get_prefix('reaction')
    gene_assoc_prefix = get_prefix('geneProductAssociation')
    annotation_prefix = get_prefix('annotation')
    gene_product_ref_prefix = get_prefix('geneProductRef')
    and_prefix = get_prefix('and')
    or_prefix = get_prefix('or')

    # Depuración: Imprimir prefijos identificados
    print(f"Identified prefixes:")
    print(f"reaction_prefix: {reaction_prefix}")
    print(f"gene_assoc_prefix: {gene_assoc_prefix}")
    print(f"annotation_prefix: {annotation_prefix}")
    print(f"gene_product_ref_prefix: {gene_product_ref_prefix}")
    print(f"and_prefix: {and_prefix}")
    print(f"or_prefix: {or_prefix}")

    # Función para construir el nombre completo de la etiqueta con URI de namespace
    def build_tag(prefix, tag):
        """
        Construye el nombre completo de la etiqueta con notación Clark.

        Args:
            prefix (str): El prefijo asociado con el namespace.
            tag (str): El nombre local de la etiqueta.

        Returns:
            str: La etiqueta en notación Clark.
        """
        if prefix and prefix in namespaces:
            return f'{{{namespaces[prefix]}}}{tag}'
        else:
            return tag  # Sin namespace

    # Definir nombres completos de etiquetas
    gene_product_association_tag = build_tag(gene_assoc_prefix, 'geneProductAssociation')
    annotation_tag = build_tag(annotation_prefix, 'annotation')

    # Definir etiquetas para geneProductRef, and, or
    gene_product_ref_tag = build_tag(gene_product_ref_prefix, 'geneProductRef')
    and_tag = build_tag(and_prefix, 'and')
    or_tag = build_tag(or_prefix, 'or')

    # Inicializar el procesador de lotes
    batch_size = 500  # Depends the memory size, adjust this value
    batch_processor = BatchProcessor(session, batch_size)

    # Función para verificar si una reacción tiene genes asociados
    def has_associate_genes(reaction):
        gene_assoc = reaction.find(f'.//{gene_product_association_tag}')
        return gene_assoc is not None

    # Función para verificar si una reacción tiene anotaciones
    def has_annotations(reaction):
        annotation = reaction.find(f'.//{annotation_tag}')
        return annotation is not None

    # Función para identificar el tipo de elemento de gen
    # Devuelve 'gene', 'and', 'or', o None
    def identify_element(element):
        tag = etree.QName(element).localname
        if tag == 'geneProductRef':
            return 'gene'
        elif tag == 'and':
            return 'and'
        elif tag == 'or':
            return 'or'
        else:
            return None

    # Función recursiva para procesar asociaciones de genes
    def add_elements(adjacent_gene_id, elements, reaction_id, batch_processor):
        cont = 0
        for elem in elements:
            element_type = identify_element(elem)
            if element_type == 'gene':
                # Obtener el atributo 'geneProduct'
                gene_product = elem.get(f'{{{namespaces.get(gene_product_ref_prefix)}}}geneProduct')
                if gene_product is None:
                    continue
                # Añadir a gene_reaction_edges
                batch_processor.add_gene_reaction_edge({
                    'gene_id': gene_product,
                    'reaction_id': reaction_id
                })
                # Añadir a gpr_edges
                batch_processor.add_gpr_edge({
                    'from_id': gene_product,
                    'to_id': adjacent_gene_id
                })
            elif element_type in ['and', 'or']:
                node_id = f"{adjacent_gene_id}_{element_type}{cont}"
                # Añadir a gpr_nodes
                batch_processor.add_gpr_node({
                    'id': node_id,
                    'type': element_type,
                    'name': element_type,
                    'reactionId': reaction_id
                })
                # Añadir a gpr_edges
                batch_processor.add_gpr_edge({
                    'from_id': node_id,
                    'to_id': adjacent_gene_id
                })
                # Procesar elementos hijos recursivamente
                child_elements = list(elem)
                add_elements(node_id, child_elements, reaction_id, batch_processor)
                cont += 1
            else:
                # Tipo de elemento desconocido, omitir
                continue

    # Función para procesar cada reacción
    def process_reaction(reaction, batch_processor):
        reaction_id = reaction.get('id')
        reaction_name = reaction.get('name')
        if has_associate_genes(reaction) and has_annotations(reaction):
            gene_assoc = reaction.find(f'.//{gene_product_association_tag}')
            if gene_assoc is not None:
                elements = list(gene_assoc)
                add_elements(reaction_id, elements, reaction_id, batch_processor)
        return

    # Encontrar todas las reacciones en el modelo
    if reaction_prefix:
        reactions = root.findall(f'.//{build_tag(reaction_prefix, "reaction")}')
    else:
        reactions = root.findall('.//reaction')  # Fallback a sin namespace

    print(f"{len(reactions)} Reactions found.")

    # Procesar cada reacción
    for reaction in reactions:
        process_reaction(reaction, batch_processor)
        reaction_id = reaction.get('id')
        reaction_name = reaction.get('name')
        print(f"Processed reaction {reaction_id}: {reaction_name}")

    # Ejecutar cualquier operación restante
    batch_processor.execute_all()


# ---------------- Uploading elements from the XML file ------------------------





def upload_elements(path_model, session):
    """
    This function uploads elements from an XML file into a Neo4j database. It handles species, reactions,
    and genes, and clears the database before uploading new data. The function processes the data in small
    batches and runs multiple separate queries.
    
    Args:
        path_model (str): The path to the XML model file.
        session (neo4j.Session): The active Neo4j session.
    """


    # Step 1: Delete all existing nodes and relationships
    delete_query = 'MATCH (n) DETACH DELETE n'
    session.run(delete_query)
 
    # Step 2: Process species in batches
    species_query = f"""
    CALL apoc.periodic.iterate(
        'CALL apoc.load.xml("file:///{path_model}") YIELD value 
         UNWIND value._children AS model
         UNWIND [item IN model._children WHERE item._type = "listOfSpecies"][0] AS ListOfSpecies
         UNWIND [item IN ListOfSpecies._children WHERE item._type = "species"] AS Specie
         RETURN DISTINCT Specie',
         
        'UNWIND [attr IN Specie._children WHERE attr._type = "annotation"] AS Annotation
         UNWIND [attr IN Annotation._children WHERE attr._type = "RDF"] AS RDF
         UNWIND [attr IN RDF._children WHERE attr._type = "Description"] AS Description
         UNWIND [attr IN Description._children WHERE attr._type = "is"] AS bqbiol
         UNWIND [attr IN bqbiol._children[0]._children WHERE attr._type = "li"] AS resource

         MERGE (s:Specie {{id: Specie.id}})
         SET s.constant = Specie.constant, 
             s.boundaryCondition = Specie.boundaryCondition, 
             s.hasOnlySubstanceUnits = Specie.hasOnlySubstanceUnits, 
             s.name = toLower(Specie.name),
             s.compartment = Specie.compartment,
             s.chemicalFormula = Specie["fbc:chemicalFormula"],
             s.type = Specie._type,
             s.Annotation = Annotation["xmlns:sbml"],
             s.RDF = RDF["xmlns:rdf"],
             s.Description = Description["rdf:about"],
             s.bqbiol = bqbiol["xmlns:bqbiol"],
             s.resource = resource["rdf:resource"]',
        {{batchSize: 1000, parallel: true}}
    );
    """
    session.run(species_query)

    # Step 3: Process reactions in batches
    reactions_query = f"""
    CALL apoc.periodic.iterate(
        'CALL apoc.load.xml("file:///{path_model}") YIELD value 
         UNWIND value._children AS model
         UNWIND [item IN model._children WHERE item._type = "listOfReactions"][0] AS ListOfReactions
         UNWIND [item IN ListOfReactions._children WHERE item._type = "reaction"] AS Reaction
         RETURN DISTINCT Reaction',

        'UNWIND [attr IN Reaction._children WHERE attr._type = "annotation"] AS Annotation_R
         UNWIND [attr IN Annotation_R._children WHERE attr._type = "RDF"] AS RDF_R
         UNWIND [attr IN RDF_R._children WHERE attr._type = "Description"] AS Description_R
         UNWIND [attr IN Description_R._children WHERE attr._type = "is"] AS bqbiol_R

         MERGE (r:Reaction {{id: Reaction.id}})
         SET r.name = toLower(Reaction.name),
             r.reversible = Reaction.reversible,
             r.lowerFluxBound = Reaction["fbc:lowerFluxBound"],
             r.fast = Reaction.fast,
             r.upperFluxBound = Reaction["fbc:upperFluxBound"],
             r.type = Reaction._type,
             r.sboTerm = Reaction.sboTerm,
             r.metaid = Reaction.metaid,
             r.Annotation = Annotation_R["xmlns:sbml"],
             r.RDF = RDF_R["xmlns:rdf"],
             r.Description = Description_R["rdf:about"],
             r.bqbiol = bqbiol_R["xmlns:bqbiol"]',
        {{batchSize: 1000, parallel: true}}
    );
    """
    session.run(reactions_query)

    # Step 4: Process reactants in batches
    reactants_query = f"""
    CALL apoc.periodic.iterate(
        'CALL apoc.load.xml("file:///{path_model}") YIELD value 
         UNWIND value._children AS model
         UNWIND [item IN model._children WHERE item._type = "listOfReactions"][0] AS ListOfReactions
         UNWIND [item IN ListOfReactions._children WHERE item._type = "reaction"] AS Reaction
         UNWIND [item IN Reaction._children WHERE item._type = "listOfReactants"] AS ListOfReactants
         UNWIND [item IN ListOfReactants._children WHERE item._type = "speciesReference"] AS Reactant
         RETURN DISTINCT Reactant, Reaction',

        'MATCH (n1:Specie {{id: Reactant.species}}), (n2:Reaction {{id: Reaction.id}})
         CREATE (n1)-[:Reactant_Of {{stoichiometry: Reactant.stoichiometry, constant: Reactant.constant}}]->(n2)',
        {{batchSize: 1000, parallel: true}}
    );
    """
    session.run(reactants_query)

    # Step 5: Process products in batches
    products_query = f"""
    CALL apoc.periodic.iterate(
        'CALL apoc.load.xml("file:///{path_model}") YIELD value 
         UNWIND value._children AS model
         UNWIND [item IN model._children WHERE item._type = "listOfReactions"][0] AS ListOfReactions
         UNWIND [item IN ListOfReactions._children WHERE item._type = "reaction"] AS Reaction
         UNWIND [item IN Reaction._children WHERE item._type = "listOfProducts"] AS ListOfProducts
         UNWIND [item IN ListOfProducts._children WHERE item._type = "speciesReference"] AS Product
         RETURN DISTINCT Product, Reaction',

        'MATCH (n3:Specie {{id: Product.species}}), (n4:Reaction {{id: Reaction.id}})
         CREATE (n3)<-[:Produce {{stoichiometry: Product.stoichiometry, constant: Product.constant}}]-(n4)',
        {{batchSize: 1000, parallel: true}}
    );
    """
    session.run(products_query)

    # Step 6: Process genes in batches
    gene_query = f"""
    CALL apoc.periodic.iterate(
        'CALL apoc.load.xml("file:///{path_model}") YIELD value 
         UNWIND value._children AS model
         UNWIND [item IN model._children WHERE item._type = "listOfGeneProducts"] AS ListOfGenes
         UNWIND [item IN ListOfGenes._children WHERE item._type = "geneProduct"] AS Gene
         RETURN DISTINCT Gene',

        'MERGE (g:Gene {{id: Gene["ns1:id"]}})
         SET g.label = Gene["ns1:label"],
             g.metaid = Gene.metaid,
             g.sboTerm = Gene.sboTerm,
             g.type = "gene"',
        {{batchSize: 1000, parallel: true}}
    );
    """
    session.run(gene_query)
    print("se procesaron los genes")







def copy_file(source, destination = database_url):
    """
    Copies a file from the source path to the destination path.

    :param source: Full path of the source file.
    :param destination: Full path of the destination file.
    :return: None
    """
    try:
        # Check if the source file exists
        if not os.path.isfile(source):
            print(f"The file {source} does not exist.")
            return

        # Check if the destination folder exists, if not, create it
        destination_folder = os.path.dirname(destination)
        if not os.path.exists(destination_folder):
            os.makedirs(destination_folder)

        # Copy the file
        shutil.copy2(source, destination)
        print(f"File copied from {source} to {destination}")
    except Exception as e:
        print(f"Error copying the file: {e}")



def get_import_folder(session):
    query_get_folder = """Call dbms.listConfig() YIELD name, value
WHERE name='server.directories.import'
RETURN value"""

    query_response = session.run(query_get_folder)

    # getting the import folder of the database
    for record in query_response:
        import_folder = record['value']
        return import_folder


# RUNING THE SCRIPT-------------------------------

def upload(pathModel= "models/yeast8.xml", dbImportFolderPath="", uri='bolt://localhost:7687', authUser='neo4j', authPass='11111111'):
    # Move db to db import folder
    # delete old db
    # drive query to upload reactions and metabolites
    # drive query to upload genes
    # drive program to upload logic genes


    # Import sbml model
    #pathModel = "models/yeast8.xml"

    # Connect with the graphdb
    graphdb = GraphDatabase.driver(uri=uri, auth=(authUser, authPass))
    session = graphdb.session()


    pathModel = "models/"+pathModel # ex. 'dir../yeast8.xml'
    ModelFileName = pathModel.split("/")[-1]
    dbImportFolderPath = pathModel
    copy_file(ModelFileName, get_import_folder(session))




    upload_elements(ModelFileName, session)
    cargarLogicaGenes(dbImportFolderPath, session)
    #uploadLogicGenes(pathModel, session)


    return






