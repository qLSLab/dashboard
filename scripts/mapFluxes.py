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
dConfParams = cpL.loadConfigParams(confFileName="experimentSettings_maps.toml")
modelName = dConfParams["modelName"]
medium = dConfParams["medium"]
condition = dConfParams["condition"]
sf = dConfParams["sf"]

logStrm = gL.logFileOpen(gL.dDirs["logs"], timeStamp, gL.getBaseName(sys.argv[0]))
gL.toLog(logStrm, f"\ntimeStamp: {timeStamp}")
cpL.toLogConfigParams(logStream=logStrm)
####################################################################################################

lScalingFactors = [1] + gL.getPrimeNumberInRange(2, 50)



############################à

gL.toLog(logStrm, "\n------ Map fluxes - all Conditions - scaling factor 11 - scaled only medium components ------\n")
gL.toLog(logStrm, "List of maps:\n")

dScaledOnlyMedium = {'endogeno': '20230719165435',
    'endogeno_batterico': '20230719173043',
    'endogeno_ncs': '20230719175312',
    'endogeno_taco': '20230719180916',
    'endogeno_marino': '20230719184105'}

dictTypes = {
    "Time": str}
conditionFile = "extracFluxes_YNB_EG.tsv"
dfConsProd = pd.read_csv(
    os.path.join(RAWDIR, conditionFile),
    sep="\t",
    dtype=dictTypes)

highestFlux_fba = -float("Inf")
lowestFlux_fba = float("Inf")
low_fba = float("Inf")

for timePoint in dfConsProd.itertuples():
    tWet = timePoint.Time
    for caso in dScaledOnlyMedium:
        baseFileName = modelName + "_" + caso + "_" + medium + "_" + condition + "_" + tWet + 'h' + "_sf" + sf + "_" + dScaledOnlyMedium[caso]
        outputFileNameFBA = "FBA_" + baseFileName + ".tsv"
        outputFileNameFVA = "FVA_" + baseFileName + ".tsv"

        dfFBA = pd.read_csv(
                os.path.join(OUTDIR, outputFileNameFBA),
                sep = '\t',
                dtype=dictTypes)

        dictTypes = {
            "Min": float,
            "Max": float,
            "delta": float
        }

        dfFVA = pd.read_csv(
            os.path.join(OUTDIR, outputFileNameFVA),
            sep = '\t',
            dtype=dictTypes)

        df_fba_fva = pd.merge(dfFBA, dfFVA, on = 'Rxn')
        outputFileNameFBAFVA = "FBA_FVA_" + baseFileName + ".tsv"
        df_fba_fva.to_csv(os.path.join(OUTDIR, outputFileNameFBAFVA), sep="\t", index=False)

        lFluxes_fba = df_fba_fva['Flux'].tolist() + df_fba_fva['delta'].tolist()

        [mapBoundaryColValue_fba, currentMin, currentMax, low] = faL.findMapBound(lFluxes_fba)

        if currentMax > highestFlux_fba:
            highestFlux_fba = currentMax*1

        if currentMin < lowestFlux_fba:
            lowestFlux_fba = currentMin*1

        if low < low_fba:
            low_fba = low*1

gL.toLog(logStrm, f"\nHighest Flux: {highestFlux_fba}")
gL.toLog(logStrm, f"Lowest Flux: {lowestFlux_fba}\n")

for timePoint in dfConsProd.itertuples():
    tWet = timePoint.Time
    gL.toLog(logStrm, f"Time Point: {tWet}")
    for caso in dScaledOnlyMedium:
        baseFileName = modelName + "_" + caso + "_" + medium + "_" + condition + "_" + tWet + 'h' + "_sf" + sf + "_" + dScaledOnlyMedium[caso]
        outputFileNameFBAFVA = "FBA_FVA_" +baseFileName + ".tsv"
        gL.toLog(logStrm, f"outputFileName Flux: {outputFileNameFBAFVA}")


        fba_fva_flux = os.path.join(OUTDIR, outputFileNameFBAFVA)
        fbaFlux = str(2)
        fvaFlux = str(5)

        ### RAINBOW
        scriptPath = os.path.join(UTILSDIR, "colorFluxes_py3_centralMet.py")
        mapColouredPath = os.path.join(MAPDIR,"centralMet_FBA_FVA" + baseFileName + "_rainbow.svg")
        gL.toLog(logStrm, f"{mapColouredPath}")
        mapPath =  os.path.join(MAPDIR, "y8_centralMetabolism_rainbow.svg")
        command = "python " + scriptPath + " " +  mapPath + " " + fba_fva_flux + " " + fbaFlux + " " + fvaFlux + " " + mapColouredPath + " " + str(lowestFlux_fba) + " " +  str(highestFlux_fba) + " " + str(low_fba) + " " + str('gist_rainbow_r') + " " + modelName + " " + caso+ " " + medium + " " + condition + " " + tWet + " " + sf
        os.system(str(command))

        scriptPath = os.path.join(UTILSDIR, "colorFluxes_py3.py")
        mapColouredPath = os.path.join(MAPDIR,"hetPaths_FBA_FVA_" + baseFileName + "_rainbow.svg")
        gL.toLog(logStrm, f"{mapColouredPath}")
        mapPath =  os.path.join(MAPDIR, "heterologousPaths_rainbow.svg")
        command = "python " + scriptPath + " " +  mapPath + " " + fba_fva_flux + " " + fbaFlux + " " + fvaFlux + " " + mapColouredPath + " " + str(lowestFlux_fba) + " " +  str(highestFlux_fba) + " " + str(low_fba) + " " + str('gist_rainbow_r') + " " + modelName + " " + caso+ " " + medium + " " + condition + " " + tWet + " " + sf
        os.system(str(command))

        mapColouredPath = os.path.join(MAPDIR,"glucoseSynthesis_FBA_FVA" + baseFileName + "_rainbow.svg")
        gL.toLog(logStrm, f"{mapColouredPath}")
        mapPath =  os.path.join(MAPDIR, "y8_glucose_rainbow.svg")
        command = "python " + scriptPath + " " +  mapPath + " " + fba_fva_flux + " " + fbaFlux + " " + fvaFlux + " " + mapColouredPath + " " + str(lowestFlux_fba) + " " +  str(highestFlux_fba) + " " + str(low_fba) + " " + str('gist_rainbow_r') + " " + modelName + " " + caso+ " " + medium + " " + condition + " " + tWet + " " + sf
        os.system(str(command))

        ### VIRIDIS

        scriptPath = os.path.join(UTILSDIR, "colorFluxes_py3_centralMet.py")
        mapColouredPath = os.path.join(MAPDIR,"centralMet_FBA_FVA_" + baseFileName + "_viridis.svg")
        gL.toLog(logStrm, f"{mapColouredPath}")
        mapPath =  os.path.join(MAPDIR, "y8_centralMetabolism_viridis.svg")
        command = "python " + scriptPath + " " +  mapPath + " " + fba_fva_flux + " " + fbaFlux + " " + fvaFlux + " " + mapColouredPath + " " + str(lowestFlux_fba) + " " +  str(highestFlux_fba) + " " + str(low_fba) + " " + str('viridis_r') + " " + modelName + " " + caso+ " " + medium + " " + condition + " " + tWet + " " + sf
        os.system(str(command))

        scriptPath = os.path.join(UTILSDIR, "colorFluxes_py3.py")
        mapColouredPath = os.path.join(MAPDIR,"hetPaths_FBA_FVA" + baseFileName + "_viridis.svg")
        gL.toLog(logStrm, f"{mapColouredPath}")
        mapPath =  os.path.join(MAPDIR, "heterologousPaths_viridis.svg")
        command = "python " + scriptPath + " " +  mapPath + " " + fba_fva_flux + " " + fbaFlux + " " + fvaFlux + " " + mapColouredPath + " " + str(lowestFlux_fba) + " " +  str(highestFlux_fba) + " " + str(low_fba) + " " + str('viridis_r') + " " + modelName + " " + caso+ " " + medium + " " + condition + " " + tWet + " " + sf
        os.system(str(command))

        mapColouredPath = os.path.join(MAPDIR,"glucoseSynthesis_FBA_FVA_" + baseFileName + "_viridis.svg")
        gL.toLog(logStrm, f"{mapColouredPath}\n")
        mapPath =  os.path.join(MAPDIR, "y8_glucose_viridis.svg")
        command = "python " + scriptPath + " " +  mapPath + " " + fba_fva_flux + " " + fbaFlux + " " + fvaFlux + " " + mapColouredPath + " " + str(lowestFlux_fba) + " " +  str(highestFlux_fba) + " " + str(low_fba) + " " + str('viridis_r')+ " " + modelName + " " + caso+ " " + medium + " " + condition + " " + tWet + " " + sf
        os.system(str(command))


logStrm.close()
