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
        # gL.toLog(logFile, f">>>biomassFlux\t{biomassFlux}")
    except:
        biomassFlux = 0

    statusFBA = computeFBAandSaveFlux(
        modelIn, 'FBA', outputFBAname)
    if statusFBA == 'optimal':
        statuspFBA = computeFBAandSaveFlux(
            modelIn,'pFBA', outputpFBAname)
        computeFVAandSaveFlux(modelIn, outputFVAname)
    return biomassFlux, statusFBA

def computeFBAandSaveFlux(modelIn, kindSim, fileName):
# def computeFBAandSaveFlux(modelIn, kindSim, fileName, logFile = None):
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
        )
    except:
        status = "Infeasible"
    return status


# def computeFVAandSaveFlux(modelIn, fileName, logFile = None):
def computeFVAandSaveFlux(modelIn, fileName):
    fva = flux_variability_analysis(modelIn)
    fva["minimum"] = fva["minimum"].round(7)
    fva["maximum"] = fva["maximum"].round(7)
    fva["delta"] = fva["maximum"]-fva["minimum"]
    fva.index.names = ["Rxn"]
    fva.reset_index(level=0, inplace=True)
    fva.rename(columns={"minimum": "Min", "maximum": "Max"}, inplace=True)
    fva.to_csv(os.path.join(OUTDIR, fileName + ".tsv"), sep="\t", index=False)
    try:
        fva = flux_variability_analysis(modelIn)
        fva["minimum"] = fva["minimum"].round(7)
        fva["maximum"] = fva["maximum"].round(7)
        fva["delta"] = fva["maximum"]-fva["minimum"]
        fva.index.names = ["Rxn"]
        fva.reset_index(level=0, inplace=True)
        fva.rename(columns={"minimum": "Min", "maximum": "Max"}, inplace=True)
        fva.to_csv(os.path.join(OUTDIR, fileName + ".tsv"), sep="\t", index=False)
    except:
        print('Error during FVA execution\n')
    saveFluxes2Excel(
        modelIn,
        fva,
        "FVA",
        fileName
    )


# def saveFluxes2Excel(modelIn, dfFlux, simulation, fileName, logFile = None):
def saveFluxes2Excel(modelIn, dfFlux, simulation, fileName):
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


# def getInterestingRxnsFlux(dConfigParams=None):
#     lHeterologousPaths = dConfigParams["kind"]
def getInterestingRxnsFlux(lHeterologousPaths, dConfigParams=None):
    condition = dConfigParams["condition"]
    lInterestingRxns_endogeno = [
        "R01781",
        "R01333_NAD_c",
        "R01333_NADP_c",
        "R01333_NAD_m",
        "R01333_NADP_m",
        "Transport_glycolate_c_m",
        "r_4315",
        "EX_s_3997_e"
    ]
    lInterestingRxns_endogeno2 = [
        "AO1",
        "AO2",
        "AO3",
        "AO4",
    ]
    lInterestingRxns_batterico = [
        "B4",
        "B5a_NAD",
        "B5a_NADP",
        "B5b1",
        "B5b2_NAD",
        "B5b2_NADP",
        "B6",
    ]
    lInterestingRxns_NCL = ["NCL2", "NCL3"]
    lInterestingRxns_TaCo = ["T2", "T3", "T4"]
    lInterestingRxns_marino = ["M1", "M2", "M3", "M4A", "M4B"]
    if '_'.join(condition) in ["Glc_EG", "EG_Glc"]:
        lConditionRxns = ["r_1714", "EX_C01380", "r_1992"]
    elif '_'.join(condition) in ["EG"]:
        lConditionRxns = ["EX_C01380", "r_1992"]

    lInterestingRxns_energyRouteA = ["r_0662", "r_0658", "r_0832", "r_0831", "r_1022"]
    lInterestingRxns_energyRouteB = [
        "r_0716",
        "r_0717",
        "r_0714",
        "r_1930",
        "r_0715",
        "r_1930",
        "r_0713",
        "r_4185",
    ]
    dInterestingRxns =  {"medium": lConditionRxns,
                        "energy Route A": lInterestingRxns_energyRouteA,
                        "energy Route B": lInterestingRxns_energyRouteB}

    # dInterestingRxns =  {"medium": lConditionRxns}
    for het in lHeterologousPaths:
        if het == "endogeno":
            dInterestingRxns.update({
                "Path endogeno": lInterestingRxns_endogeno,
            })
        elif het == "endogeno2":
            dInterestingRxns.update({
                "Path endogeno 2": lInterestingRxns_endogeno2,
            })
        elif het == "batterico":
            dInterestingRxns.update({
                "Path batterico": lInterestingRxns_batterico,
            })
        elif het == "ncl":
            dInterestingRxns.update({
                "Path NCL": lInterestingRxns_NCL,
            })
        elif het == "taco":
            dInterestingRxns.update({
                "Path TaCo": lInterestingRxns_TaCo,
            })
        elif het == "marino":
            dInterestingRxns.update({
                "Path marino": lInterestingRxns_marino,
            })

    return dInterestingRxns


def rxnDeletion(rxnObj, model, biomassRxn):
    print(rxnObj)
    with model:
        model.reactions.get_by_id(rxnObj).lower_bound = 0
        model.reactions.get_by_id(rxnObj).upper_bound = 0
        try:
            fba = model.optimize()
            biomassFlux = round(fba.fluxes[biomassRxn], 4)
            return [rxnObj,biomassFlux,"optimal"]
            # delta = ((biomassFlux- biomassFlux_initial)/biomassFlux_initial)*100
            # reactionEquation = faL.getReactionEquation(constrainedModel, rxn.id)
            # gL.toLog(logStrm, f"=== RXN: {rxn.id} - EQUATION: {reactionEquation}\nSTATUS: {fba.status} - NEW BIOMASS FLUX: {biomassFlux} - DIFF: {delta}% ===\n")
        except:
            return [rxnObj,0,"infeasible"]
            # reactionEquation = faL.getReactionEquation(constrainedModel, rxn.id)
            # gL.toLog(logStrm, f"=== RXN: {rxn.id} - EQUATION: {reactionEquation}\nSTATUS: infeasible - NEW BIOMASS FLUX: 0 - DIFF: -100% ===\n")
