import os
import sys
import time
import cobra as cb
import pandas as pd
import numpy as np

MODULEDIR = os.path.split(os.path.dirname(__file__))[0]
sys.path.append(MODULEDIR)
import utils.genericLib as gL
import utils.configParamsLib as cpL
import utils.experimentalSetupLib as expSL
import utils.manipulateModelLib as mmL
import utils.fluxAnalysisLib as faL
import utils.figuresLib as figL

timeStamp = gL.getTimeStamp()
####################################################################################################
# Setting the parameters
RAWDIR = gL.dDirs["raw"]
OUTDIR = gL.dDirs["out"]
MODELDIR = gL.dDirs["mods"]
FIGUREDIR = gL.dDirs["figs"]

# setting experiment params
dConfParams = cpL.loadConfigParams()
modelName = dConfParams["modelName"]
medium = dConfParams["medium"]
condition = dConfParams["condition"]
mediumSetting = dConfParams["mediumBounds"]
o2Setting = dConfParams["oxygen"]
experiment = dConfParams["kind"]

logStrm = gL.logFileOpen(gL.dDirs["logs"], timeStamp, gL.getBaseName(sys.argv[0]))
gL.toLog(logStrm, f"\ntimeStamp: {timeStamp}")
cpL.toLogConfigParams(logStream=logStrm)
####################################################################################################

gL.toLog(logStrm, f"\n-----Step 1: Load model or create it (if does not exist yet)\n")
modelName_wExpDetails = modelName + '_' + '_'.join(experiment) + '_' + medium + '_' + '_'.join(condition)
gL.toLog(logStrm, f"Model name: {modelName_wExpDetails}\n")

if os.path.exists(os.path.join(MODELDIR, modelName_wExpDetails + ".xml")) is True:
    gL.toLog(logStrm, "Model already exist. It will be now loaded.\n")
    model = cb.io.read_sbml_model(os.path.join(MODELDIR, modelName_wExpDetails + ".xml"))
else:
    gL.toLog(logStrm, "Model does not exist yet. It will be now created.")
    ## get list of heterologous pathways; for each pathway merge its model with the original model
    model = cb.io.read_sbml_model(os.path.join(MODELDIR, modelName + ".xml"))

    ## Add exchange reactions for isolated metabolites
    gL.toLog(logStrm, "* Add exchange reactions for isolated metabolites")
    expSL.exchgRxns4IsolatedMet(model)

    gL.toLog(logStrm, "* Add heterologous pathways")
    for exp in experiment:
        if os.path.exists(os.path.join(MODELDIR, exp + ".xml")) is True:
            pathModel = cb.io.read_sbml_model(os.path.join(MODELDIR, exp + ".xml"))
            model = model.merge(pathModel)
        else:
            pathModel = expSL.createHeterologousModel(exp)
            model = model.merge(pathModel)

    # Add to model uptake reactions of wet medium not already included in the model
    newRxnsFileName = 'newMediumRxns_' + medium + '_' + '_'.join(condition) + ".tsv"
    newMetsFileName = 'newMediumMets_' + medium + '_' + '_'.join(condition) + ".tsv"

    gL.toLog(logStrm, f"* Add medium reactions not already included in the original model. Reactions are included into the file called {newRxnsFileName}, metabolites are included into the file called {newMetsFileName}")

    dictTypes = {
        "Id": str,
        "Name": str,
        "Equation": str,
        "Lb": int,
        "Ub": int,
        "GPR": str
    }

    dfNewMediumRxns = pd.read_csv(
              os.path.join(RAWDIR, newRxnsFileName),
              sep="\t",
              dtype=dictTypes)

    dfNewMediumRxns = dfNewMediumRxns.fillna('')

    dictTypes = {
        "Id": str,
        "Name": str,
        "Compartment": str
    }

    dfNewMediumMets = pd.read_csv(
              os.path.join(RAWDIR, newMetsFileName), sep="\t",
              dtype=dictTypes)

    for newMet in dfNewMediumMets.itertuples():
        mmL.addMets(model, newMet)

    for newRxn in dfNewMediumRxns.itertuples():
        mmL.addRxns(model, newRxn)

    gL.toLog(logStrm,
        "Model construction is now complete and the SBML file will be saved into the models directory with the name {modelName_wExpDetails}.xml\n")
    cb.io.write_sbml_model(model, os.path.join(MODELDIR, modelName_wExpDetails + ".xml"))

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

#### ???? come chiamare ile se ho un altro set di dati in stesseccondizioni
dfConsProd = pd.read_csv(
    os.path.join(RAWDIR, conditionFile),
    sep="\t",
    dtype=dictTypes,
)

dModels_timePoints = {}
for timePointRow in dfConsProd.itertuples():
    tWet = timePointRow.Time
    gL.toLog(logStrm, f"\n>>> Experimental time point: {tWet}h")
    with model:
        # r_4315 need to be irreversible
        model.reactions.get_by_id("r_4315").lower_bound = 0
        model.reactions.get_by_id("r_4315").upper_bound = 1000

        ## Change boundaries of exchange reactions from infinite to the default value (1000 in this case)
        gL.toLog(logStrm, "* Scale the exchange reactions boundaries")
        model = expSL.scaleExchangeRxnsBoundaries(model)

        # Block all uptake reactions
        gL.toLog(logStrm,
            "* Block all uptake reactions (to block any undesired uptake)")
        model = expSL.switchOffAllSinkRxns(model, dConfParams)

        # add medium concentrations
        gL.toLog(logStrm, "* Constrain model with experimental concentrations of medium components")
        model = expSL.addMediumData(model, dfMedium, dConfParams)
        for rxnMedium in dfMedium.itertuples():
            rxn = model.reactions.get_by_id(rxnMedium.RxnID)
            gL.toLog(
                logStrm, f"\t{rxn.id} - {rxn.name} - ({rxn.lower_bound}, {rxn.upper_bound})"
            )

        gL.toLog(logStrm, "\n")

        # Set extracellular fluxes
        gL.toLog(
            logStrm,
            "* Constrain extracellular fluxes with experimental consumption or production data")
        model = expSL.setExtracellularFluxes(timePointRow, model, dConfParams)
        lExtracFluxs = [dConfParams['glc'], dConfParams['etoh'], dConfParams['glycerol'], dConfParams['acetate'], dConfParams['eg']]

        for rxnExtFlux in lExtracFluxs:
            rxn = model.reactions.get_by_id(rxnExtFlux)
            gL.toLog(
                logStrm, f"\t{rxn.id} - {rxn.name} - ({rxn.lower_bound}, {rxn.upper_bound})"
            )
        gL.toLog(logStrm, "\n")

        ## Set O2 uptake flux: variazione dell'O2 per ottenere le diverse OD della popolazione di cellule oppure ragiono sulla singola cellula e lasciamo svincolato e attivo il consumo di O2
        gL.toLog(
            logStrm,
            "* Constrain oxygen consumption flux")
        if dConfParams['oxygen'] == 'free':
            model, o2LB = expSL.setO2Flux(model,tWet, dConfigParams=dConfParams)
        elif dConfParams['oxygen'] == 'fixed':
            if tWet == '0':
                model, biomassFluxTime0, o2LB = expSL.setO2Flux(model,tWet, dConfigParams=dConfParams)
                gL.toLog(lgStream, f"Fixed O2 lower bound: {o2LB}")
            else:
                model, biomassFluxCurrentTime, deltaWetSilicoBiomOD, o2LB = expSL.setO2Flux(model,tWet, biomassFluxT0 = biomassFluxTime0, dConfigParams=dConfParams)
                gL.toLog(lgStream, f"Delta between Wet OD and in Silico biomass synthesis flux: {deltaWetSilicoBiomOD}")
                gL.toLog(lgStream, f"Fixed O2 lower bound: {o2LB}")
        dModels_timePoints[tWet] = model
        cb.io.write_sbml_model(model, os.path.join(MODELDIR, modelName_wExpDetails + '_' + tWet + 'h' + "_" + timeStamp + ".xml")) #### timestam?


sys.exit()

####################################################################################################

gL.toLog(logStrm, f"\n-----Step 3: Simulate models\n")

lBiomass = []
for tPoint in dModels_timePoints:
    ## Biomass production
    gL.toLog(logStrm, "Setting the biomass reaction as objective function")
    biomassRxn = expSL.setBiomassReaction(dModels_timePoints[tPoint], dConfParams)

    if biomassRxn is None:
        gL.toLog(logStrm, "It was not possible to set a biomass reaction. EXIT\n")
        sys.exit()
    else:
        gL.toLog(logStrm, f"Biomass Rxn: {biomassRxn}\n")
        gL.toLog(
            logStrm,
            f"Obj Rxn: {model.objective.expression}\nObj max/min: {model.objective.direction}",
        )

    dModels_timePoints[tPoint].solver = "glpk"

    outputFileNamepFBA = "pFBA_" +  modelName_wExpDetails + '_' + tWet + 'h' + "_" + timeStamp
    outputFileNameFBA = "FBA_" +  modelName_wExpDetails + '_' + tWet + 'h' + "_" + timeStamp
    outputFileNameFVA = "FVA_" +  modelName_wExpDetails + '_' + tWet + 'h' + "_" + timeStamp

    gL.toLog(logStrm, f"\nCompute fluxes")
    biomass_currentTimepoint = faL.computeAndExportFluxes(dModels_timePoints[tPoint],
                                OUTDIR,
                                outputFileNamepFBA,
                                outputFileNameFBA,
                                outputFileNameFVA,
                                dConfParams)
    lBiomass.append(biomass_currentTimepoint)

lODs = []
lTimes = []
for timePointRow in dfConsProd.itertuples():
    lTimes.append(int(timePointRow.Time))
    ## Compute average wet OD (necessary for generate plot OD and silico biomass vs. time)
    wetOD_currentTimePoint = expSL.averageWetOD(timePointRow)
    lODs.append(wetOD_currentTimePoint)

dfODvsBIOM= pd.DataFrame({'OD': lODs, 'BIOM': lBiomass})
dfODvsBIOM = dfODvsBIOM/dfODvsBIOM.iloc[0]
figL.plotODvsBiomasslux(lX, dfODvsBIOM['OD'].tolist(), dfODvsBIOM['BIOM'].tolist(), FIGUREDIR, "ODBIOMvsTIME_" +  modelName_wExpDetails + '_' + tWet + 'h' + "_" + timeStamp)

logStrm.close()
