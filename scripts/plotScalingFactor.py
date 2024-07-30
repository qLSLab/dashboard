#!/usr/bin/env python
import os
import sys
import time
import cobra as cb
import pandas as pd
import numpy as np
import itertools as itt
import matplotlib.pyplot as plt

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
sf = dConfParams["sf"]

logStrm = gL.logFileOpen(gL.dDirs["logs"], timeStamp, gL.getBaseName(sys.argv[0]))
gL.toLog(logStrm, f"\ntimeStamp: {timeStamp}")
cpL.toLogConfigParams(logStream=logStrm)
####################################################################################################

lScalingFactors_p1 = [1] + gL.getPrimeNumberInRange(2, 20)
lScalingFactors_p2 = gL.getPrimeNumberInRange(20, 50)

dScaledOnlyMedium = {
    'endogeno': ['20230719165435', '20230720123753'],
    'endogeno_batterico': ['20230719173043', '20230720073035'],
    'endogeno_ncs': ['20230719175312', '20230720075512'],
    'endogeno_taco': ['20230719180916', '20230720083016'],
    'endogeno_marino': ['20230719184105', '20230720085935']}

# gL.toLog(logStrm, "List of maps:\n")

dictTypes = {
    "Time": str}
conditionFile = "extracFluxes_YNB_EG.tsv"
dfConsProd = pd.read_csv(
    os.path.join(RAWDIR, conditionFile),
    sep="\t",
    dtype=dictTypes)

for caso in dScaledOnlyMedium:
    dSF = {}
    print('caso: ', caso)
    tStamp1 = dScaledOnlyMedium[caso][0]
    tStamp2 = dScaledOnlyMedium[caso][0]
    lTime =[]
    for timePoint in dfConsProd.itertuples():
        tWet = timePoint.Time
        lTime.append(int(tWet))
        for sf in lScalingFactors_p1:
            baseFileName = modelName + "_" + caso + "_" + medium + "_" + condition + "_" + tWet + 'h' + "_sf" + str(sf) + "_" + tStamp1
            outputFileNameFBA = "FBA_" + baseFileName + ".tsv"

            dictTypes = {
                "Flux": float,
            }
            dfFBA = pd.read_csv(
                    os.path.join(OUTDIR, outputFileNameFBA),
                    sep = '\t',
                    dtype=dictTypes)

            biomassFlux = round(dfFBA[dfFBA['Rxn'] == 'r_2111']['Flux'].values[0], 5)
            if sf in dSF:
                dSF[sf] += [biomassFlux]
            else:
                dSF[sf] = [biomassFlux]
        for sf in lScalingFactors_p2:
            print('sf: ', sf)
            baseFileName = modelName + "_" + caso + "_" + medium + "_" + condition + "_" + tWet + 'h' + "_sf" + str(sf) + "_" + tStamp2
            outputFileNameFBA = "FBA_" + baseFileName + ".tsv"

            dictTypes = {
                "Flux": float,
            }
            dfFBA = pd.read_csv(
                    os.path.join(OUTDIR, outputFileNameFBA),
                    sep = '\t',
                    dtype=dictTypes)

            biomassFlux = round(dfFBA[dfFBA['Rxn'] == 'r_2111']['Flux'].values[0], 5)
            print('biomassFlux\n', biomassFlux)
            if sf in dSF:
                dSF[sf] += [biomassFlux]
            else:
                dSF[sf] = [biomassFlux]

    #print('dSF\n', dSF)

    fig = plt.figure(figsize=(13.3, 10))
    ax = fig.add_subplot(111)
    for sf in dSF:
        ax.scatter(
            lTime, [el / int(sf) for el in dSF[sf]], alpha=0.7, marker="o", s=50, label="sf " + str(sf)
        )


    plt.xlabel("Time (h)", size=18)
    plt.ylabel("Biomass flux scaled", size=18)
    plt.title(caso, fontsize=18)
    plt.xticks(lTime,fontsize=18)
    plt.yticks(fontsize=18)
    plt.legend(fontsize=10)
    plt.savefig(os.path.join(FIGUREDIR, caso + "_sf_1_50.png"))
    plt.close()




logStrm.close()
