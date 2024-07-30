#!/usr/bin/env python
import os
import sys
import time
import cobra as cb
import pandas as pd
import numpy as np
from cobra.util.solver import linear_reaction_coefficients
MODULEDIR = os.path.split(os.path.dirname(__file__))[0]
sys.path.append(MODULEDIR)
import utils.genericLib as gL
import utils.configParamsLib as cpL
import utils.experimentalSetupLib as expSL
import utils.manipulateModelLib as mmL
import utils.fluxAnalysisLib as faL
import utils.figuresLib as figL
import utils.modelReconstructionLib as recL
import utils.keggLib as kL
import re
from cobra import Model


timeStamp = gL.getTimeStamp()
####################################################################################################
# Setting the parameters
RAWDIR = gL.dDirs["raw"]
OUTDIR = gL.dDirs["out"]
MODELDIR = gL.dDirs["mods"]
FIGUREDIR = gL.dDirs["figs"]
UTILSDIR =  gL.dDirs["utils"]
MAPDIR =  gL.dDirs["map"]

# setting create model params
configFileName = "curateModel.toml"
dConfParams = cpL.loadConfigParams(confFileName=configFileName)
print(dConfParams)

logStrm = gL.logFileOpen(gL.dDirs["logs"], timeStamp, gL.getBaseName(sys.argv[0]))
gL.toLog(logStrm, f"\ntimeStamp: {timeStamp}")
cpL.toLogConfigParams(confFileName=configFileName, logStream=logStrm)

## Setting experiment params
configFileName = "experimentSettings.toml"
print(configFileName)
dConfParams_exp = cpL.loadConfigParams(confFileName=configFileName)
cpL.toLogConfigParams(confFileName=configFileName, logStream=logStrm)

## Load or create file including correspondence between metabolites name and their KEGG identifier
if os.path.exists(os.path.join(OUTDIR, dConfParams_exp["nameToKeggId"])) is False:
    dfName2KeggId = pd.DataFrame({"keggId": [], "Name": []})
else:
    dfName2KeggId = pd.read_csv(os.path.join(OUTDIR, dConfParams_exp["nameToKeggId"]), sep = "\t")

## Load starting model
if "modelName" in dConfParams:
    startingModelName = dConfParams["modelName"]
    if os.path.exists(os.path.join(MODELDIR, startingModelName)) is False:
        print(f"Attention: {startingModelName} does not exist.")
        sys.exit()
    else:
        modelGW = cb.io.read_sbml_model(os.path.join(MODELDIR, startingModelName))
        if "curation" in dConfParams:
            modelGW = recL.curateModel(modelGW, dfName2KeggId, logFile = logStrm, dConfigParams=dConfParams)
        dKeggId2ModelMet = mmL.getMet2KeggId(modelGW) # dict e' fatto cosi': d = {'C06197': {'c': 's_4317'}}

else:
    modelGW = Model()

print("LEN RXNS: ", len(modelGW.reactions), "\n")

dKeggId2ModelMet = mmL.getMet2KeggId(modelGW) # dict e' fatto cosi': d = {'C06197': {'c': 's_4317'}}

lActiveMet = []
for activeMet in ["H+", "H2O", "Oxygen", "Phosphate"]:
    correspondingKeggId, dfName2KeggId = kL.getKeggMetId(activeMet.lower(), dfName2KeggId)
    lActiveMet.append(correspondingKeggId)

## Create model
if os.path.exists(os.path.join(MODELDIR, dConfParams_exp["experimentModel"])) is False:
    newModel, dKeggId2ModelMet, dfName2KeggId = recL.createNewModel(dfName2KeggId, modelGW, dizKeggCompound2ModelCompound = dKeggId2ModelMet, dConfigParams=dConfParams_exp)
else:
    newModel = cb.io.read_sbml_model(os.path.join(MODELDIR, dConfParams_exp["experimentModel"]))

modelGW_plus_newModel = modelGW.merge(newModel)
dKeggId2ModelMet = mmL.getMet2KeggId(modelGW_plus_newModel)
cb.io.write_sbml_model(modelGW_plus_newModel, os.path.join(MODELDIR, startingModelName.split(".")[0] + "_" + dConfParams_exp["experimentModel"]))

# model constraining
if "mediumFileName" in dConfParams_exp:
    mediumFileName = dConfParams_exp["mediumFileName"]
    dfMedium = pd.read_csv(
              os.path.join(RAWDIR, mediumFileName),
              sep="\t",
              dtype={"Concentration": float})

    conditionMetsFileName = dConfParams_exp["condition"]
    dfConditionMets = pd.read_csv(
              os.path.join(RAWDIR, conditionMetsFileName),
              sep="\t",
              dtype={"Lb": float})

    dfName2KeggId, modelGW_plus_newModel = expSL.addMediumData2Model(modelGW_plus_newModel, dfMedium, dfConditionMets, lActiveMet, dfName2KeggId, dKeggId2ModelMet, dConfigParams = dConfParams_exp)

    ## save model
    lConditionMets = [re.sub('[\s]+', '', met) for met in dfConditionMets.Met.tolist()]
    cb.io.write_sbml_model(modelGW_plus_newModel, os.path.join(MODELDIR, dConfParams_exp["modelName"].split(".")[0] + "_" + dConfParams_exp["experimentModel"].split(".")[0] + "_" + "_".join(lConditionMets) + "_" + dConfParams_exp["mediumFileName"].split(".")[0] + ".xml"))


if "extracelProdConsFileName" in dConfParams_exp:
    dfConsProd = pd.read_csv(
        os.path.join(RAWDIR, dConfParams_exp["extracelProdConsFileName"]),
        sep="\t",
        dtype= {"Time": str})
    dfConsProd.set_index("Time", inplace = True)

    ## define type columns
    dMeasuredMetsCols = {}
    for col in dfConsProd.columns:
        dMeasuredMetsCols[col] = float
    dfConsProd = dfConsProd.astype(dMeasuredMetsCols)

    gL.toLog(
        logStrm,
        "* Constrain extracellular fluxes with variation between two consecutive time points of experimental consumption and production data")

    dfConditionMets = pd.read_csv(
              os.path.join(RAWDIR, conditionMetsFileName),
              sep="\t")
    lConditionMets = [re.sub('[\s]+', '', met) for met in dfConditionMets.Met.tolist()]

    baseModelName = dConfParams_exp["modelName"].split(".")[0] + "_" + dConfParams_exp["experimentModel"].split(".")[0] + "_" + "_".join(lConditionMets) + "_" + dConfParams_exp["mediumFileName"].split(".")[0]
    expSL.setExtracFluxes(dfConsProd, list(dMeasuredMetsCols.keys()), modelGW_plus_newModel, dKeggId2ModelMet, dfName2KeggId, baseModelName)

dfName2KeggId.to_csv(os.path.join(OUTDIR, dConfParams_exp["nameToKeggId"]), sep = "\t", index = False)


sys.exit()

# if dConfParams_exp["exp"] == "endogenous" or dConfParams_exp["exp"] == "aox":
#     print("modifico EG boundary")
#     for filePath in os.listdir(os.path.join(MODELDIR)):
#         if os.path.isfile(os.path.join(MODELDIR, filePath)) and filePath.startswith(dConfParams_exp["modelName"] + "_" + dConfParams_exp["exp"] + "_" + "_".join(lCondition)):
#             print("WT MODEL: ", filePath)
#             # set EG consumption to 0 for wt models
#             wt = cb.io.read_sbml_model(os.path.join(MODELDIR, filePath))
#
#             ## convert name to kegg id
#             metLower = "ethylene glycol"
#             print("metLower: ", metLower)
#             correspondingKeggId, dfName2KeggId = kL.getKeggMetId(metLower, dfName2KeggId)
#
#             # search for metabolite in model corresponding to the retrieved kegg id
#             rxnId, modelMetId = mmL.getExchInModelForKeggId(wt, correspondingKeggId, dKeggId2ModelMet)
#             print("met: ", correspondingKeggId)
#             print("modelMetId: ", modelMetId)
#             print("rxnId: ", rxnId, "\n")
#             if rxnId is None:
#                 mmL.addExchRxnForMet(wt, wt.metabolites.get_by_id(modelMetId), newLB=0, newUB=0)
#                 rxnId = 'EX_' + modelMetId
#                 print("None rxnId: ", rxnId, "\n")
#
#             rxn = wt.reactions.get_by_id(rxnId)
#             rxn.lower_bound = 0
#             rxn.upper_bound = 0
#
#             cb.io.write_sbml_model(wt, os.path.join(MODELDIR, filePath))
#
# dfName2KeggId.to_csv(os.path.join(OUTDIR, "name2KeggId_" + startingModelName + ".tsv"), sep = "\t", index = False)


## add AOX
#setting experiment params
configFileName = "experimentSettings.toml" # HA SENSO CREARE UN FILE PARAMETRI PER OGNI PATH DA CREARE OPPURE MODIFICARE QUESTO?
dConfParams_exp = cpL.loadConfigParams(confFileName=configFileName)
logStrm = gL.logFileOpen(gL.dDirs["logs"], timeStamp, gL.getBaseName(sys.argv[0]))
gL.toLog(logStrm, f"\ntimeStamp: {timeStamp}")
cpL.toLogConfigParams(confFileName=configFileName, logStream=logStrm)

startingModelName = dConfParams_exp["modelName"]
modelFromName = startingModelName + "_endogenous"
wt = cb.io.read_sbml_model(os.path.join(MODELDIR, modelFromName + ".xml"))

##
if os.path.exists(os.path.join(OUTDIR, "name2KeggId_" + startingModelName + ".tsv")) is False:
    dfName2KeggId = pd.DataFrame({"keggId": [], "Name": []})
else:
    dfName2KeggId = pd.read_csv(os.path.join(OUTDIR, "name2KeggId_" + startingModelName + ".tsv"), sep = "\t")

lActiveMet = []
for activeMet in ["H+", "H2O", "Oxygen", "Phosphate"]:
    correspondingKeggId, dfName2KeggId = kL.getKeggMetId(activeMet.lower(), dfName2KeggId)
    lActiveMet.append(correspondingKeggId)

dKeggId2ModelMet = mmL.getMet2KeggId(wt) # dict e' fatto cosi': d = {'C06197': {'c': 's_4317'}}


if os.path.exists(os.path.join(MODELDIR, dConfParams_exp["exp"] + ".xml")) is False:
    newModel, dKeggId2ModelMet, dfName2KeggId = recL.createNewModel(dKeggId2ModelMet,dfName2KeggId, dConfigParams=dConfParams_exp, model2GetInfoFrom = wt)
else:
    newModel = cb.io.read_sbml_model(os.path.join(MODELDIR, dConfParams_exp["exp"] + ".xml"))

wt_aox = wt.merge(newModel)
dKeggId2ModelMet = mmL.getMet2KeggId(wt_aox)

cb.io.write_sbml_model(wt_aox, os.path.join(MODELDIR, modelFromName + "_" + dConfParams_exp["exp"] + ".xml"))
dfName2KeggId.to_csv(os.path.join(OUTDIR, "name2KeggId_" + startingModelName + ".tsv"), sep = "\t", index = False)


# model constraining
medium = dConfParams_exp["medium"]
condition = dConfParams_exp["condition"]
lCondition = [re.sub('[\s]+', '_',cond) for cond in condition]

if "mediumData" in dConfParams_exp:
    dfMedium = pd.read_csv(
              os.path.join(RAWDIR, 'medium' + medium + ".tsv"),
              sep="\t",
              dtype={"Concentration": float})

    dfName2KeggId, wt_aox = expSL.addMediumData2Model(wt_aox, dfMedium, lActiveMet, dfName2KeggId, dKeggId2ModelMet, dConfigParams = dConfParams_exp)
    ## save model
    cb.io.write_sbml_model(wt_aox, os.path.join(MODELDIR, modelFromName + "_" + dConfParams_exp["exp"] + "_" + "_".join(lCondition) + "_" + dConfParams_exp["medium"] + ".xml"))

if "extracelProdConsData" in dConfParams_exp:
    conditionFile = "extracFluxes_" +  medium + '_' + '_'.join(lCondition) + ".tsv"
    dfConsProd = pd.read_csv(
        os.path.join(RAWDIR, conditionFile),
        sep="\t",
        dtype= {"Time": str})
    dfConsProd.set_index("Time", inplace = True)

    ## define type columns
    dMeasuredMetsCols = {}
    for col in dfConsProd.columns:
        dMeasuredMetsCols[col] = float
    dfConsProd = dfConsProd.astype(dMeasuredMetsCols)

    gL.toLog(
        logStrm,
        "* Constrain extracellular fluxes with variation between two consecutive time points of experimental consumption and production data")

    expSL.setExtracFluxes(dfConsProd, list(dMeasuredMetsCols.keys()), wt_aox, dKeggId2ModelMet, dfName2KeggId, dConfigParams=dConfParams_exp)

if dConfParams_exp["exp"] == "endogenous" or dConfParams_exp["exp"] == "aox":
    print("modifico EG boundary")
    for filePath in os.listdir(os.path.join(MODELDIR)):
        if os.path.isfile(os.path.join(MODELDIR, filePath)) and filePath.startswith(modelFromName + "_" + dConfParams_exp["exp"] + "_" + "_".join(lCondition)):
            print("WT MODEL: ", filePath)
            # set EG consumption to 0 for wt models
            wt_aox = cb.io.read_sbml_model(os.path.join(MODELDIR, filePath))

            ## convert name to kegg id
            metLower = "ethylene glycol"
            print("metLower: ", metLower)
            correspondingKeggId, dfName2KeggId = kL.getKeggMetId(metLower, dfName2KeggId)

            # search for metabolite in model corresponding to the retrieved kegg id
            rxnId, modelMetId = mmL.getExchInModelForKeggId(wt_aox, correspondingKeggId, dKeggId2ModelMet)
            print("met: ", correspondingKeggId)
            print("modelMetId: ", modelMetId)
            print("rxnId: ", rxnId, "\n")
            if rxnId is None:
                mmL.addExchRxnForMet(wt_aox, wt_aox.metabolites.get_by_id(modelMetId), newLB=0, newUB=0)
                rxnId = 'EX_' + modelMetId
                print("None rxnId: ", rxnId, "\n")

            rxn = wt_aox.reactions.get_by_id(rxnId)
            rxn.lower_bound = 0
            rxn.upper_bound = 0

            cb.io.write_sbml_model(wt_aox, os.path.join(MODELDIR, filePath))

dfName2KeggId.to_csv(os.path.join(OUTDIR, "name2KeggId_" + startingModelName + ".tsv"), sep = "\t", index = False)

sys.exit()



#mediumFileName = 'medium' + medium + '_' + '_'.join(condition) + ".tsv"
mediumFileName = 'medium' + medium + ".tsv"
dfMedium = pd.read_csv(
          os.path.join(RAWDIR, mediumFileName),
          sep="\t",
          dtype={"Concentration": float})

conditionFile = "extracFluxes_" +  medium + '_' + '_'.join(condition) + "_consumption.tsv"
dfConsProd_consumption = pd.read_csv(
    os.path.join(RAWDIR, conditionFile),
    sep="\t",
    dtype= {"Time": str})

dMeasuredMetsCols_consumption = {}
for col in dfConsProd_consumption.columns:
    if col.split("_")[0] not in ["Time", "OD"]:
        dMeasuredMetsCols_consumption[col] = float

dfConsProd_consumption = dfConsProd_consumption.astype(dMeasuredMetsCols_consumption)
# print(dfConsProd_consumption.dtypes)

conditionFile = "extracFluxes_" +  medium + '_' + '_'.join(condition) + "_production.tsv"
dfConsProd_production = pd.read_csv(
    os.path.join(RAWDIR, conditionFile),
    sep="\t",
    dtype= {"Time": str})

dMeasuredMetsCols_production = {}
for col in dfConsProd_production.columns:
    if col.split("_")[0] not in ["Time", "OD"]:
        dMeasuredMetsCols_production[col] = float

dfConsProd_production = dfConsProd_production.astype(dMeasuredMetsCols_production)

for experiment in lexperiments:
    print(experiment)
    gL.toLog(logStrm, f"Experiment: {experiment}")
    modelName_wExpDetails = modelName + '_' + '_'.join(experiment) + '_' + medium + '_' + '_'.join(condition)
    gL.toLog(logStrm, f"Model name: {modelName_wExpDetails}\n")
    ## creare file conversione mets in kegg id
    model = recL.getModel(experiment, logFile = logStrm, isolatedMetsClosing=closeMetsIsol, recreate = recreateModel, curateBaseModel = curateModel, dConfigParams=dConfParams)

    for tPair in lTimePointsCouples:
        t1 = tPair[0]
        t2 = tPair[1]
        print(t1, t2)
        gL.toLog(logStrm, f"\n>>> Experimental time points: {t2}h and {t1}h")
        constrainedModel = recL.constrainWwetData(model, modelName_wExpDetails, dfConsProd_consumption, dfConsProd_production, dfMedium, tStamp = timeStamp, dConfigParams=dConfParams, logFile = logStrm, tPoint1 = t1, tPoint2 = t2, expName = experiment, curateBaseModel = curateModel)


logStrm.close()
