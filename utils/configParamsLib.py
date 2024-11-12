import sys

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib
import utils.genericLib as gL
import cobra as cb
import os
import bioservices.kegg as kegg
import utils.keggLib as kL

MODELDIR = gL.dDirs["mods"]

def validateParameter(paramTag, paramValue, lValues, dParams):
    if type(paramValue) is str and paramValue in lValues:
        dParams[paramTag] = paramValue
    elif type(paramValue) is list and all(elem in lValues for elem in paramValue):
        dParams[paramTag] = paramValue
    else:
        print(f"{paramTag} not recognized. It should be one of: {lValues}")
        sys.exit()

def loadConfigParams(confFileName="experimentSettings.toml"):
    gL.beingOrExit(confFileName)
    with open(confFileName, "rb") as confFile:
        dConfigParams = tomllib.load(confFile)
        # print("dConfigParams\n", dConfigParams)

    # dValidParams = {
    #     "modelName": ["yeast8", "eciYali_sbmlValidated", "yeast9"],
    #     "lRxns2Add": lRxnsKEGG.copy(),
    #     "organismCode": kL.getKEGGorganismCodes(),
    #     "checkMassBalance": ["Y", "N"],
    #     "connectSameIsolatedMetsWithSimpleTransport": ["Y", "N"],
    #     "closeIsolatedMetabolitesWithExchange": ["Y", "N"],
    #     "checkSomethingIsProducedFromNothing": ["Y", "N"],
    #     "checkIfUptakeSomethingBlockingDemand": ["Y", "N"],
    #     "checkLoops": ["Y", "N"],
    #     "refineBoundariesInReactionsFormingLoops": ["Y", "N"],
    #     "curateBiologExperimentAddedRxns": ["Y", "N"],
    #     "exp": ["endogenous", "aox"],
    #     "mediumData": ["Y", "N"],
    #     "extracelProdConsData": ["Y", "N"],
    #     "condition": ["glucose", "EG", "acetate", "ethylene glycol"],
    #     "medium": ["YNB", "YNByarr"],
    #     "mediumBounds": ["wet", "free"],
    #     'lActiveUptake': lRxns_originalModel.copy(),
    #     "keggRxnsToAddFileName": ["endogenous.tsv", "aox.tsv"],
    #     "addExchRxnsForTheseMetsFileName": ["endogenous_addExchForThoseMets.tsv", "aox_addExchForThoseMets.tsv"],
    # }
    # dMod = {
    #     "modelName": "yeast8",
    #     "lRxns2Add": lRxnsKEGG.copy(),
    #     "organismCode": "sce",
    #     "checkMassBalance": "Y",
    #     "connectSameIsolatedMetsWithSimpleTransport": "Y",
    #     "closeIsolatedMetabolitesWithExchange": "Y",
    #     "checkSomethingIsProducedFromNothing": "Y",
    #     "checkIfUptakeSomethingBlockingDemand": "Y",
    #     "checkLoops": "Y",
    #     "refineBoundariesInReactionsFormingLoops": "Y",
    #     "curateBiologExperimentAddedRxns": "Y",
    #     "exp": "endogenous",
    #     "mediumData": "Y",
    #     "extracelProdConsData": "Y",
    #     "condition": ["glucose", "EG"],
    #     "medium": "YNB",
    #     "mediumBounds": "free",
    #     "lActiveUptake": ["r_1992", "r_1832"],
    #     "keggRxnsToAddFileName": "endogenous.tsv",
    #     "addExchRxnsForTheseMetsFileName": "endogenous_addExchForThoseMets.tsv",
    # }
    # # for param in dValidParams:
    #
    # for param in dConfigParams:
    #     try:
    #         validateParameter(param, dConfigParams[param], dValidParams[param], dMod)
    #     except:
    #         print(f"{param} is not defined in the config file")
    #         sys.exit()
    return dConfigParams

def toLogConfigParams(confFileName="experimentSettings.toml", logStream=None):
    if logStream is not None:
        gL.toLog(logStream, "\n>>>>>>>> start toml file params >>>>>>>>\n\n")
        with open(confFileName, "r") as confFile:
            tomlConfFile = str(confFile.read())
            confFile.close()
            logStream.write(tomlConfFile)
        gL.toLog(logStream, "\n\n<<<<<<<< end toml file params <<<<<<<<\n\n")
    else:
        print("There is no log stream")
        sys.exit()
