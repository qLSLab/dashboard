# Import libraries
import xml.dom.minidom as xmlFile
import xml.dom.minidom
from bs4 import BeautifulSoup
from neo4j import GraphDatabase
import shutil
import os
from decouple import config


database_url = config('DATABASE_IMPORT_FOLDER_URL')



# This function asnwer if the reaction have genes to process, has a boolean output
def hasAssociateGenes(reaction):
    try:
        answ = reaction.find_all('ns1:geneProductAssociation')
#         print(answ)
        return bool(answ)
    except:
        print('nada encontrado')
        return False 

def hasAnnotations(reaction):
    try:
        answ = reaction.find_all('ns1:annotation')
#         print(answ)
#         print(bool(answ))
        return bool(answ)
    except:
        print('nada encontrado')
        return False 



def AddElements(adjacentGeneId, geneFile, reactionId, session):
    nodeList = []
    cont = 0
#     print('********')
    for i in geneFile:
#         print('imprimiendo i'+str(i))
#         print('identificativo'+ str(identifyElements(str(i))))
        # hacer una lista de elementos
        typeElement = identifyElements(str(i))
        if not typeElement: # es un gen
            nodeList.append(i)
#             print(i['fbc:geneProduct'])
            print('imprimiendo genes')
            _id = i['ns1:geneProduct']
            q2 =   '''match (n), (p)
            where n.id="{0}" and p.id="{1}"
            create (n)-[:GPR_edge]->(p)'''.format(_id, adjacentGeneId)
            session.run(q2)
            
            q3 =   '''match (n), (p)
            where n.id="{0}" and p.id="{1}"
            create (n)-[:is_a_gene_that_works_in]->(p)'''.format(_id, reactionId)
            session.run(q3)
            print(q3)
            cont = cont + 1
            
            
        elif typeElement == 1:
            nodeList.append('or')
            _id = adjacentGeneId + "_or" + str(cont)
            q1 = "create (N:GPR{id:'"+_id+"', type:'or', name:'or', reactionId:'"+reactionId+"' })"
#             print(q1)
            session.run(q1)
            print('imprmiendo or')
            
            q2 =   '''match (n), (p)
            where n.id="{0}" and p.id="{1}"
            create (n)-[:GPR_edge]->(p)'''.format(_id, adjacentGeneId)
            session.run(q2)

            cont = cont + 1
            AddElements(_id, cleanArray(i), reactionId, session)
            
        elif typeElement == 2:
            nodeList.append('and')
            _id = adjacentGeneId + "_and" + str(cont)
            q1 = "create (N:GPR{name:'and', id:'"+_id+"', type:'and', reactionId:'"+reactionId+"'  })"
#             print(q1)
            session.run(q1)
#             print('imprmiendo and')
            
            q2 =   '''match (n), (p)
            where n.id="{0}" and p.id="{1}"
            create (n)-[:GPR_edge]->(p)'''.format(_id, adjacentGeneId)
            session.run(q2)
            
            cont = cont + 1
            AddElements(_id, cleanArray(i), reactionId, session)
    


# This function indentify the kind of element (and/or/gene)
# input: the gene file
# output: and-2, or-1, gene-0
def identifyElements(geneFile):
    if(geneFile[0:9] == '<ns1:and>'):
        return 2
    elif(geneFile[0:8] == '<ns1:or>'):
        return 1
    elif(geneFile[0:19] == '<ns1:geneProductRef'):
        return 0
    else:
#         print(geneFile[0:19])
#         print(str(geneFile))
#         print('no se pudo identificar')
        return -3

def cleanArray(array):
    array = [x for x in array if x != '\n']
    return array


# In this function starts the process of each reaction to find and generate his genes in the Neo4j graph
def iteracionReaccion(reactionName, reactionId, bs_data, session):
    reactionFile = bs_data.find_all('reaction', {'id':reactionId})[0]
    if hasAssociateGenes(reactionFile) and hasAnnotations(reactionFile):
        # The first element
        geneFile = cleanArray(list(cleanArray(list(reactionFile))[3].children))
#        print(geneFile)
#         print(cleanArray(list(cleanArray(list(reactionFile))[3].children)))
        AddElements(reactionId, geneFile, reactionId, session)
        
        return
#     print('hola')
    return


def uploadLogicGenes(pathModel, session):

    doc = xmlFile.parse(pathModel);

    # Reading data from the xml file
    with open(pathModel, 'r') as f:
        data = f.read()

    # Passing the data of the xml
    # file to the xml parser of
    # beautifulsoup
    bs_data = BeautifulSoup(data, 'xml')


    # Obtain the list of reactions
    # ciclo donde se gira el codigo
    reactionList = doc.getElementsByTagName("ns0:reaction")  
    print("%d Reactions" % reactionList.length)
    # reactionName = 'Thymd Transport' #'PDH' # Prueba
    # reactionId = 'R_THYMDt1'# 'R_PDHm' # Prueba
    # iteracionReaccion(reactionName, reactionId) # Prueba

    for i in reactionList:
        reactionId = i.getAttribute("id")
        reactionName = i.getAttribute("name")
        iteracionReaccion(reactionName, reactionId, bs_data, session) 

        print(reactionId, reactionName)
    return



def uploadElements(pathModel, session):

    DelQuery = 'match (n) detach delete (n)'

    RMQuery = '''CALL apoc.load.xml("file:///{pathModel}") 
    YIELD value 
    UNWIND value._children AS model

    UNWIND [item  in model._children WHERE item._type = "listOfSpecies"][0] AS ListOfSpecies
    UNWIND [item  in ListOfSpecies._children WHERE item._type = "species"] AS Specie
    UNWIND  [attr IN Specie._children WHERE attr._type = 'annotation' ] as Annotation
    UNWIND  [attr IN Annotation._children WHERE attr._type = 'RDF' ] as RDF
    UNWIND  [attr IN RDF._children WHERE attr._type = 'Description' ] as Description
    UNWIND  [attr IN Description._children WHERE attr._type = 'is' ] as bqbiol
    UNWIND  [attr IN bqbiol._children[0]._children WHERE attr._type = 'li' ] as resource
    MERGE (s:Specie {{id: Specie.id}})
    SET s.constant = Specie.constant, 
        s.boundaryCondition = Specie.boundaryCondition, 
        s.hasOnlySubstanceUnits = Specie.hasOnlySubstanceUnits, 
        s.name = toLower(Specie.name),
        s.compartment = Specie.compartment,
        s.chemicalFormula = Specie['fbc:chemicalFormula'],
        s.type = Specie._type,
        s.Annotation = Annotation["xmlns:sbml"],
        s.RDF = RDF["xmlns:rdf"],
        s.Description = Description["rdf:about"],
        s.bqbiol = bqbiol["xmlns:bqbiol"],
        s.resource = resource["rdf:resource"]

    WITH DISTINCT model as model
    UNWIND [item  in model._children WHERE item._type = "listOfReactions"][0] AS ListOfReactions
    UNWIND [item  in ListOfReactions._children WHERE item._type = "reaction"] AS Reaction
    UNWIND  [attr IN Reaction._children WHERE attr._type = 'annotation' ] as Annotation_R
    UNWIND  [attr IN Annotation_R._children WHERE attr._type = 'RDF' ] as RDF_R
    UNWIND  [attr IN RDF_R._children WHERE attr._type = 'Description' ] as Description_R
    UNWIND  [attr IN Description_R._children WHERE attr._type = 'is' ] as bqbiol_R

    MERGE (r:Reaction {{id: Reaction.id}})
    SET r.name = toLower(Reaction.name),
        r.reversible = Reaction.reversible,
        r.lowerFluxBound = Reaction['fbc:lowerFluxBound'],
        r.fast = Reaction.fast,
        r.upperFluxBound = Reaction['fbc:upperFluxBound'],
        r.type = Reaction._type,
        r.sboTerm = Reaction.sboTerm,
        r.metaid = Reaction.metaid,
        r.Annotation = Annotation_R["xmlns:sbml"],
        r.RDF = RDF_R["xmlns:rdf"],
        r.Description = Description_R["rdf:about"],
        r.bqbiol = bqbiol_R["xmlns:bqbiol"]    

    WITH DISTINCT Reaction
    UNWIND [item  in Reaction._children WHERE item._type = "listOfReactants"] AS ListOfReactants
    UNWIND [item  in ListOfReactants._children WHERE item._type = "speciesReference"] AS Reactant

    MATCH (n1 {{id:Reactant.species}}), (n2 {{id:Reaction.id}})
    CREATE (n1) -[:Reactant_Of {{stoichiometry:Reactant.stoichiometry, constant:Reactant.constant }}]-> (n2)

    WITH DISTINCT Reaction
    UNWIND [item  in Reaction._children WHERE item._type = "listOfProducts"] AS ListOfProducts
    UNWIND [item  in ListOfProducts._children WHERE item._type = "speciesReference"] AS Product

    MATCH (n3 {{id:Product.species}}), (n4 {{id:Reaction.id}})
    CREATE (n3) <-[:Produce {{stoichiometry:Product.stoichiometry, constant:Product.constant }}]- (n4)'''.format(pathModel=pathModel)

    print(RMQuery)

    gQuery = '''CALL apoc.load.xml("file:///{pathModel}") 
    YIELD value 
    UNWIND value._children AS model

    WITH DISTINCT model as model
    UNWIND [item  in model._children WHERE item._type = "listOfGeneProducts"] AS ListOfGenes
    UNWIND [item  in ListOfGenes._children WHERE item._type = "geneProduct"] AS Gene

    MERGE (g:Gene {{id: Gene["ns1:id"]}})
    SET g.label = Gene["ns1:label"],
        g.metaid = Gene.metaid,
        g.sboTerm = Gene.sboTerm
    '''.format(pathModel=pathModel)

    session.run(DelQuery) # Delete old db
    session.run(RMQuery)  # add reactions and metabolites to db
    session.run(gQuery)   # add genes to db




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





# RUNING THE SCRIPT-------------------------------

def upload(pathModel= "models/yeast8.xml", dbImportFolderPath="", uri='bolt://localhost:11007', authUser='neo4j', authPass='1'):
    # Move db to db import folder
    # delete old db
    # drive query to upload reactions and metabolites
    # drive query to upload genes
    # drive program to upload logic genes

    pathModel = "models/"+pathModel # ex. 'dir../yeast8.xml'
    ModelFile = pathModel #pathModel.split("/")[-1]
    dbImportFolderPath = database_url + '/{ModelFile}'
    copy_file(ModelFile, dbImportFolderPath)


    # Import sbml model
    #pathModel = "models/yeast8.xml"

    # Connect with the graphdb
    graphdb = GraphDatabase.driver(uri=uri, auth=(authUser, authPass))
    session = graphdb.session()

    uploadElements(ModelFile, session)
    #uploadLogicGenes(pathModel, session)


    return






