#!/usr/bin/env python
import os
import sys
import time
import cobra as cb
import pandas as pd
import numpy as np
import itertools as itt

MODULEDIR = os.path.split(os.path.dirname(__file__))[0]
sys.path.append(MODULEDIR)
import utils.genericLib as gL
import utils.configParamsLib as cpL
import utils.experimentalSetupLib as expSL
import utils.manipulateModelLib as mmL
import utils.fluxAnalysisLib as faL
import utils.figuresLib as figL
import utils.modelReconstructionLib as recL

timeStamp = gL.getTimeStamp()
####################################################################################################
# Setting the parameters
RAWDIR = gL.dDirs["raw"]
OUTDIR = gL.dDirs["out"]
MODELDIR = gL.dDirs["mods"]
FIGUREDIR = gL.dDirs["figs"]
UTILSDIR =  gL.dDirs["utils"]
MAPDIR =  gL.dDirs["map"]

# setting experiment params
dConfParams = cpL.loadConfigParams()
modelName = dConfParams["modelName"]
medium = dConfParams["medium"]
condition = dConfParams["condition"]
mediumSetting = dConfParams["mediumBounds"]
o2Setting = dConfParams["oxygen"]
experiment = dConfParams["kind"]
scalingFactor = dConfParams["sf"]

logStrm = gL.logFileOpen(gL.dDirs["logs"], timeStamp, gL.getBaseName(sys.argv[0]))
gL.toLog(logStrm, f"\ntimeStamp: {timeStamp}")
cpL.toLogConfigParams(logStream=logStrm)
####################################################################################################

gL.toLog(logStrm, f"\n-----Step 1: Load model or create it (if does not exist yet)\n")
modelName_wExpDetails = modelName + '_' + '_'.join(experiment) + '_' + medium + '_' + '_'.join(condition)
gL.toLog(logStrm, f"Model name: {modelName_wExpDetails}\n")

model = recL.getModel(dConfParams, modelName, modelName_wExpDetails,medium, condition, experiment, logFile = logStrm)

gL.toLog(logStrm, f"\nModel summary:")
gL.toLog(logStrm, f"Number of reactions:\t{len(model.reactions)}")
gL.toLog(logStrm, f"Number of metabolites:\t{len(model.metabolites)}")
gL.toLog(logStrm, f"Number of genes:\t{len(model.genes)}")
gL.toLog(logStrm, f"Boundary rections:\t{len(model.boundary)}")
gL.toLog(logStrm, f"Exchange rections:\t{len(model.exchanges)}")
gL.toLog(logStrm, f"Sinks rections:\t{len(model.sinks)}\n")

####################################################################################################

gL.toLog(logStrm, f"\n-----Step 2: Constrain model with experimental data\n")
mediumFileName = 'medium' + medium + '_' + '_'.join(condition) + ".tsv"
gL.toLog(logStrm, f"Medium concentrations are included into the file called {mediumFileName}")

dfMedium = pd.read_csv(
          os.path.join(RAWDIR, mediumFileName),
          sep="\t",
          dtype={"Concentration": float})

dictTypes = {
    "Time": str,
    "Glucose_rep1": float,
    "Glucose_rep2": float,
    "Glucose_rep3": float,
    "Ethanol_rep1": float,
    "Ethanol_rep2": float,
    "Ethanol_rep3": float,
    "Glycerol_rep1": float,
    "Glycerol_rep2": float,
    "Glycerol_rep3": float,
    "Acetate_rep1": float,
    "Acetate_rep2": float,
    "Acetate_rep3": float,
    "OD_rep1": float,
    "OD_rep2": float,
    "OD_rep3": float,
}

conditionFile = "extracFluxes_" +  medium + '_' + '_'.join(condition) + ".tsv"
gL.toLog(logStrm, f"Extracellular fluxes are included into the file called {conditionFile}")

dfConsProd = pd.read_csv(
    os.path.join(RAWDIR, conditionFile),
    sep="\t",
    dtype=dictTypes,
)

gL.toLog(logStrm, f"Scaling factor of concentration of medium components: {scalingFactor}")
recL.constrainWwetData(model, modelName_wExpDetails, dfConsProd, dfMedium, tStamp = timeStamp, dConfigParams = dConfParams, logFile = logStrm, scale = scalingFactor)

####################################################################################################

gL.toLog(logStrm, f"\n-----Step 3: Simulate models\n")

lBiomass = []
dFeasible = {}
for timePointRow in dfConsProd.itertuples():
    tWet = timePointRow.Time
    gL.toLog(logStrm, f"tWet: {tWet}")

    constrainedModel = cb.io.read_sbml_model(os.path.join(MODELDIR, modelName_wExpDetails + '_' + tWet + 'h' + "_" + timeStamp + ".xml"))

    ## Biomass production
    gL.toLog(logStrm, "Setting the biomass reaction as objective function")
    biomassRxn = expSL.setBiomassReaction(constrainedModel, dConfParams)

    if biomassRxn is None:
        gL.toLog(logStrm, "It was not possible to set a biomass reaction. EXIT\n")
        sys.exit()
    else:
        gL.toLog(logStrm, f"Biomass Rxn: {biomassRxn}\n")
        gL.toLog(
            logStrm,
            f"Obj Rxn: {constrainedModel.objective.expression}\nObj max/min: {constrainedModel.objective.direction}",
        )

    constrainedModel.solver = "glpk"
    constrainedModel.reactions.get_by_id(biomassRxn).objective_coefficient = 1

    try:
        fba = constrainedModel.optimize()
        biomassFlux = fba.fluxes[biomassRxn]
        gL.toLog(logStrm, f">>>biomassFlux\t{biomassFlux}")
    except:
        biomassFlux = 0

    outputFileNamepFBA = "pFBA_" +  modelName_wExpDetails + '_' + tWet + 'h' + "_" + timeStamp
    outputFileNameFBA = "FBA_" +  modelName_wExpDetails + '_' + tWet + 'h' + "_" + timeStamp
    outputFileNameFVA = "FVA_" +  modelName_wExpDetails + '_' + tWet + 'h' + "_" + timeStamp

    gL.toLog(logStrm, f"\nCompute fluxes")
    biomass_currentTimepoint, statusFBA = faL.computeAndExportFluxes(constrainedModel,
                                biomassRxn,
                                outputFileNamepFBA,
                                outputFileNameFBA,
                                outputFileNameFVA,
                                logFile = logStrm)
    lBiomass.append(biomass_currentTimepoint)
    dFeasible[tWet] = statusFBA

gL.toLog(logStrm, f"\nModels feasibility:\t{dFeasible}\n")
gL.toLog(logStrm, f"List of biomass fluxes vs. time: {lBiomass}\n")

lODs = []
lTimes = []
for timePointRow in dfConsProd.itertuples():
    lTimes.append(int(timePointRow.Time))
    ## Compute average wet OD (necessary for generate plot OD and silico biomass vs. time)
    wetOD_currentTimePoint = expSL.averageWetOD(timePointRow)
    lODs.append(wetOD_currentTimePoint)

gL.toLog(logStrm, f"List of ODs vs. time: {lODs}")
gL.toLog(logStrm, f"List of biomass fluxes vs. time: {lBiomass}")
dfODvsBIOM= pd.DataFrame({'OD': lODs, 'BIOM': lBiomass})
dfODvsBIOM = dfODvsBIOM/dfODvsBIOM.iloc[0]
figL.plotODvsBiomasslux(lTimes, dfODvsBIOM['OD'].tolist(), dfODvsBIOM['BIOM'].tolist(), "ODBIOMvsTIME_" + timeStamp)

logStrm.close()
