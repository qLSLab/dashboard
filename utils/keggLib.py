import bioservices.kegg as kegg
import re
import bioservices.kegg as kegg
from Bio.KEGG.REST import *
from fuzzywuzzy import process

def getKeggInfo(query):
    dizInfo = {}
    try:
        k = kegg.KEGG()
        info = k.get(query)
        dizInfo = k.parse(info)
    except dizInfo == 400 or dizInfo == 404:
            dizInfo = {}
    return dizInfo

def getAllKEGGItemsInKEGGdb(dbName):
    k = kegg.KEGG()
    return k.list(dbName)

def getKEGGorganismCodes():
    k = kegg.KEGG()
    return k.organismIds


def searchMetInKegg(met2Search):
    threshold = 85
    met2Search = re.sub('[^0-9a-zA-Z]+', ' ',met2Search)
    print("met2Search:\t", met2Search)
    met2Search_joined = re.sub('[\s]', '+', met2Search)
    request = kegg_find(database = "compound", query = met2Search_joined)
    lAllResults = request.readlines()
    # print("lAllResults\t", lAllResults)
    found = False
    i = 0
    while i < len(lAllResults) and found == False:
        r = lAllResults[i]
        lsComps = r.strip('\n').split("\t")
        alternativeNames = [x.lower() for x in lsComps[1].split(';')]
        matches = process.extract(met2Search, alternativeNames, limit=10)
        # print("matches\t", matches)
        lMatches_aboveThresold = [match for match in matches if match[1] > threshold]

        idx = 0
        keep = "N"
        while idx < len(lMatches_aboveThresold) and keep == "N":
            possibleMatch = lMatches_aboveThresold[idx]
            print("possible Match: ", possibleMatch)
            keep = input("Is a valid match? Y/N ---> ")
            if keep == "Y":
                found = True
                correspondingKeggId = lsComps[0].split(":")[1]
            idx += 1
        i += 1
    return correspondingKeggId

def getKeggMetId(metName, dfConversionMetName2MetKegg):
    if metName not in dfConversionMetName2MetKegg["Name"].values:
        correspondingKeggId = searchMetInKegg(metName)
        print("correspondingKeggId: ", correspondingKeggId)
        dfConversionMetName2MetKegg.loc[len(dfConversionMetName2MetKegg)] = [correspondingKeggId, metName]
    else:
        correspondingKeggId = dfConversionMetName2MetKegg[dfConversionMetName2MetKegg["Name"] == metName]["keggId"].values[0]
    return correspondingKeggId, dfConversionMetName2MetKegg
