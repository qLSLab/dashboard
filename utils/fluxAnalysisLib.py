import os
import sys
import cobra as cb
import pandas as pd
from cobra.flux_analysis import flux_variability_analysis
import utils.genericLib as gL

OUTDIR = gL.dDirs["out"]

def computeAndExportFluxes(
    modelIn,
    biomassRxnID,
    outputpFBAname,
    outputFBAname,
    outputFVAname, logFile = None):
    try:
        fba = modelIn.optimize()
        biomassFlux = fba.fluxes[biomassRxnID]
        gL.toLog(logFile, f">>>biomassFlux\t{biomassFlux}")
    except:
        biomassFlux = 0

    statusFBA = computeFBAandSaveFlux(
        modelIn, 'FBA', outputFBAname, logFile)
    if statusFBA == 'optimal':
        statuspFBA = computeFBAandSaveFlux(
            modelIn,'pFBA', outputpFBAname, logFile)
        computeFVAandSaveFlux(modelIn, outputFVAname, logFile)
    return biomassFlux, statusFBA

def computeFBAandSaveFlux(modelIn, kindSim, fileName, logFile = None):
    try:
        if kindSim == 'pFBA':
            fba = cb.flux_analysis.pfba(modelIn)
        elif kindSim == 'FBA':
            fba = modelIn.optimize()
        status = fba.status
        fbaFile = open(os.path.join(OUTDIR, fileName + ".tsv"), "w")
        gL.writeLineByLineToFile(fbaFile, ["Rxn", "Name", "Flux"], "\t")
        for row in fba.fluxes.iteritems():
            gL.writeLineByLineToFile(
                fbaFile, [row[0], modelIn.reactions.get_by_id(row[0]).name, row[1]], "\t"
            )
        fbaFile.close()

        saveFluxes2Excel(
            modelIn,
            fba,
            kindSim,
            fileName,
            logFile
        )
    except:
        status = "Infeasible"
    return status


def computeFVAandSaveFlux(modelIn, fileName, logFile = None):
    fva = flux_variability_analysis(modelIn)
    fva["minimum"] = fva["minimum"].round(7)
    fva["maximum"] = fva["maximum"].round(7)
    fva["delta"] = fva["maximum"]-fva["minimum"]
    fva.index.names = ["Rxn"]
    fva.reset_index(level=0, inplace=True)
    fva.rename(columns={"minimum": "Min", "maximum": "Max"}, inplace=True)
    fva.to_csv(os.path.join(OUTDIR, fileName + ".tsv"), sep="\t", index=False)

    saveFluxes2Excel(
        modelIn,
        fva,
        "FVA",
        fileName,
        logFile
    )


def saveFluxes2Excel(
    modelIn, dfFlux, simulation, fileName, logFile = None):
    lFlux_wEquation = []
    if simulation in ["pFBA", "FBA"]:
        for row in dfFlux.fluxes.iteritems():
            reactionEquation = getReactionEquation(modelIn, row[0])
            lFlux_wEquation.append(
                [
                    row[0],
                    modelIn.reactions.get_by_id(row[0]).name,
                    reactionEquation,
                    row[1],
                ]
            )

        dfFlux_wEquation = pd.DataFrame(
            lFlux_wEquation, columns=["Rxn", "Name", "Equation", "Flux"]
        )
    elif simulation == "FVA":
        for row in dfFlux.itertuples():
            reactionEquation = getReactionEquation(modelIn, row.Rxn)
            lFlux_wEquation.append(
                [
                    row.Rxn,
                    modelIn.reactions.get_by_id(row.Rxn).name,
                    reactionEquation,
                    row.Min,
                    row.Max,
                    row.delta
                ]
            )

        dfFlux_wEquation = pd.DataFrame(
            lFlux_wEquation,
            columns=["Rxn", "Name", "Equation", "Min", "Max", "Delta"],
        )


    # gL.toLog(
    #     logFile, f"\n\n>>> {simulation} -> Print flux values of interesting reactions:"
    # )


    dfFlux_wEquation.to_excel(os.path.join(OUTDIR, fileName + ".xlsx"))


def getReactionEquation(modelIn, rxnId):
    rxnObj = modelIn.reactions.get_by_id(rxnId)
    reactionString = rxnObj.reaction
    for react in rxnObj.reactants:
        reactionString = reactionString.replace(react.id, react.name)

    for prod in rxnObj.products:
        reactionString = reactionString.replace(prod.id, prod.name)
    return reactionString

def findMapBound(valuesList):
    ## finding the min and max values excluding INFs
    while max(valuesList) == float('Inf'):
        valuesList.remove(float('Inf'))
    maxVal = max(valuesList)

    while min(valuesList) == -float('Inf'):
        valuesList.remove(-float('Inf'))
    minVal = min(valuesList)
    low = min(x for x in valuesList if x != 0)

    ## fixing the boundaries for the color map to be simmetric around 0
    mapBoundary = max(abs(minVal), maxVal)
    return [mapBoundary, minVal, maxVal, low]
