import utils.genericLib as gL
import utils.keggLib as kL
import cobra as cb
import sys
from collections import Counter
import utils.fluxAnalysisLib as faL
import pandas as pd
import re
import itertools as itt
import os
import utils.genericLib as gL

OUTDIR = gL.dDirs["out"]

"""
Questa libreria contiene funzioni per tutto cio' che riguarda gli elementi del modello:
prendere info esistenti, modificare info esistenti, aggiungere nuove info/mets/reazioni/geni,
modificare boundaries modello.
"""

def addMissingChemicalFormulae(model):
    for m in model.metabolites:
        if m.formula is None and "kegg.compound" in m.annotation:
            keggId = m.annotation["kegg.compound"]
            dizMet = kL.getKeggInfo("cpd:" + keggId)
            if "FORMULA" in dizMet:
                chemFormula = dizMet["FORMULA"]
                m.formula = chemFormula

def getMet2Compartment(model):
    dMet2Compartment = {}
    for met in model.metabolites:
        dMet2Compartment[met.id] = met.compartment
    return dMet2Compartment

def getMet2Name(model):
    dMet2Name = {}
    for met in model.metabolites:
        dMet2Name[met.id] = met.name.lower()
    return dMet2Name

def isIsolatedMet(model, met):
    m = model.metabolites.get_by_id(met)
    isolatedMet = False
    if len(m.reactions) == 1:
        isolatedMet = True
    elif len(m.reactions) > 1:
        reazioni = list(m.reactions)
        lonlySub = []
        lonlyProd = []
        lForward = []
        lBackward = []
        for r in reazioni:
            lonlySub.append(m in r.reactants)
            lonlyProd.append(m in r.products)
            if r.lower_bound == 0 and r.upper_bound > 0:
                lForward.append(True)
                lBackward.append(False)
            elif r.lower_bound < 0 and r.upper_bound == 0:
                lForward.append(False)
                lBackward.append(True)
            else:
                lForward.append(False)
                lBackward.append(False)
        if (all(lonlySub) is True and all(lForward) is True) or (all(lonlySub) is True and all(lBackward) is True) or (all(lonlyProd) is True and all(lForward) is True) or (all(lonlyProd) is True and all(lBackward) is True):
            isolatedMet = True
    return isolatedMet

def connectSameMetaboliteDifferentCompartment(model, dfbyName):
    dMet2Compartment = getMet2Compartment(model)
    for el in dfbyName.iteritems():
        metName = re.sub('[^0-9a-zA-Z]+', '_',el[0])
        if len(el[1]) == 2:
            lCompartments = []
            for m in el[1]:
                lCompartments.append(dMet2Compartment[m])
            addTransportRxn(model, metName, el[1], lCompartments)
        elif len(el[1]) > 2:
            lMet2BeTransported_initial = el[1]
            lMet2BeTransported_after = el[1]
            for isolMet in el[1]:
                if dMet2Compartment[isolMet] == "e":
                    lMet2BeTransported_after.remove(isolMet)
                    extracMet = isolMet
            if len(lMet2BeTransported_after) == len(lMet2BeTransported_initial):
                lCombinations = list(itt.combinations(lMet2BeTransported_after, 2))
                for pair in lCombinations:
                    lCompartments = []
                    for m in pair:
                        lCompartments.append(dMet2Compartment[m])
                    addTransportRxn(model, metName, pair, lCompartments)
            else:
                lExtracMet = gL.difference(lMet2BeTransported_initial, lMet2BeTransported_after)[0]
                for internalMet in lMet2BeTransported_after:
                    lCompartments = ["e", dMet2Compartment[internalMet]]
                    addTransportRxn(model, metName, [lExtracMet, internalMet], lCompartments)

def addTransportRxn(model, nameMet, lTransportedMets, lCompartments, lb = -1000, ub = 1000):
    idRxn = "Transport_" + nameMet + "_" + "_".join(lCompartments)
    if idRxn not in model.reactions:
        reaction = cb.Reaction(idRxn)
        reaction.name = nameMet + " transport between " + " and ".join(lCompartments)
        reaction.lower_bound = lb
        reaction.upper_bound = ub
        dSubs = gL.dizReaProd(lTransportedMets[0])
        dProds = gL.dizReaProd(lTransportedMets[1])
        dEquation = {}
        for s in dSubs:
            if s in model.metabolites:
                dEquation[model.metabolites.get_by_id(s)] = -1 * dSubs[s]
            else:
                addMets(model, p, p, "c", chemFormula = "", fromKEGG = False)

        for p in dProds:
            if p in model.metabolites:
                dEquation[model.metabolites.get_by_id(p)] = 1 * dProds[p]
            else:
                addMets(model, p, p, "c", chemFormula = "", fromKEGG = False)

        reaction.add_metabolites(dEquation)
        model.add_reactions([reaction])
        model.reactions.get_by_id(idRxn).gene_reaction_rule = ""
    else:
        print(f"Attention: {idRxn} is already in the model!")
        sys.exit()

def exchgRxns4IsolatedMet(model, dfbyName):
    """Add exchange reactions for isolated metabolites"""
    print("num rxns prima addition of exc rxns for isolated mets: " , len(model.reactions))
    for el in dfbyName.iteritems():
        metName = re.sub('[^0-9a-zA-Z]+', '_',el[0])
        if len(el[1]) == 1:
            isolatedMet = el[1][0]
            mObj = model.metabolites.get_by_id(isolatedMet)
            if len(mObj.reactions) == 1:
                reazione = list(mObj.reactions)[0]
                if reazione.lower_bound < 0 and reazione.upper_bound > 0:
                    addExchRxnForMet(model, mObj, -1000.0, 1000.0)
                elif reazione.lower_bound == 0 and reazione.upper_bound > 0:
                    if mObj in reazione.reactants:
                        addExchRxnForMet(model, mObj, -1000.0, 0.0)
                    elif mObj in reazione.products:
                        addExchRxnForMet(model, mObj, 0.0, 1000.0)
                elif reazione.lower_bound < 0 and reazione.upper_bound == 0:
                    if mObj in reazione.reactants:
                        addExchRxnForMet(model, mObj, 0.0, 1000.0)
                    elif mObj in reazione.products:
                        addExchRxnForMet(model, mObj, -1000.0, 0.0)
            # elif len(mObj.reactions) > 1:
            #     reazioni = list(mObj.reactions)
            #     lonlySub = []
            #     lonlyProd = []
            #     lForward = []
            #     lBackward = []
            #     for r in reazioni:
            #         lonlySub.append(mObj in r.reactants)
            #         lonlyProd.append(mObj in r.products)
            #         if r.lower_bound == 0 and r.upper_bound > 0:
            #             lForward.append(True)
            #             lBackward.append(False)
            #         elif r.lower_bound < 0 and r.upper_bound == 0:
            #             lForward.append(False)
            #             lBackward.append(True)
            #         else:
            #             lForward.append(False)
            #             lBackward.append(False)
            #     if all(lonlySub) is True and all(lForward) is True:
            #         addExchRxnForMet(model, mObj, -1000.0, 0.0)
            #     elif all(lonlySub) is True and all(lBackward) is True:
            #         addExchRxnForMet(model, mObj, 0, 1000.0)
            #     elif all(lonlyProd) is True and all(lForward) is True:
            #         addExchRxnForMet(model, mObj, 0, 1000.0)
            #     elif all(lonlyProd) is True and all(lBackward) is True:
            #         addExchRxnForMet(model, mObj, -1000.0, 0.0)
    print("num rxns dopo addition of exc rxns for isolated mets: " , len(model.reactions))

def addExchRxnForMet(modelIn, metIn, newLB=-1000.0, newUB=1000.0):
    """Add an exchange reaction (exchange or sink) for a given metabolite

    Keyword arguments:
         modelIn - is the cobra model to which the rxn should be added
         metIn - is the metabolite for which the rxn should be added
         newLB - is the lower bound of the exchange rxn
         newUB - is the upper bound of the exchange rxn

    Return:
         the new boundary reaction
    """
    # met = modelIn.metabolites.get_by_id(metIn.id)
    # # mId = met.id[:met.id.find('[')]
    # mId = met.id
    try:
        newBoundaryRxn = modelIn.add_boundary(metIn, type="my-bound", reaction_id = 'EX_' + metIn.id, lb=newLB, ub=newUB)
    except:
        print("An exchange reaction already exists for ", mId)
    return newBoundaryRxn

def blockUptakeAndDemand(model):
    for rxn in model.reactions:
        if len(rxn.reactants) * len(rxn.products) == 0.0:
            rxn.lower_bound = 0
            rxn.upper_bound = 0

def blockUptake(model):
    for rxn in model.reactions:
        if len(rxn.reactants) * len(rxn.products) == 0.0:
            if len(rxn.reactants) != 0 and rxn.lower_bound < 0:
                rxn.lower_bound = 0
            elif len(rxn.products) != 0 and rxn.upper_bound > 0:
                rxn.upper_bound = 0

def getUptake(model):
    lUptake_wPositiveFlux = []
    lUptake_wNegativeFlux = []
    for rxn in model.reactions:
        if len(rxn.reactants) * len(rxn.products) == 0.0:
            if len(rxn.reactants) != 0 and rxn.lower_bound < 0:
                lUptake_wNegativeFlux.append(rxn.id)
            elif len(rxn.products) != 0 and rxn.upper_bound > 0:
                lUptake_wPositiveFlux.append(rxn.id)
    return lUptake_wPositiveFlux, lUptake_wNegativeFlux

def blockDemand(model):
    lDemand_wPositiveFlux = []
    lDemand_wNegativeFlux = []
    for rxn in model.reactions:
        if len(rxn.reactants) * len(rxn.products) == 0.0:
            if len(rxn.reactants) != 0 and rxn.upper_bound > 0:
                rxn.upper_bound = 0
            elif len(rxn.products) != 0 and rxn.lower_bound < 0:
                rxn.lower_bound = 0

def getDemand(model):
    lDemand_wPositiveFlux = []
    lDemand_wNegativeFlux = []
    for rxn in model.reactions:
        if len(rxn.reactants) * len(rxn.products) == 0.0:
            if len(rxn.reactants) != 0 and rxn.upper_bound > 0:
                lDemand_wPositiveFlux.append(rxn.id)
            elif len(rxn.products) != 0 and rxn.lower_bound < 0:
                lDemand_wNegativeFlux.append(rxn.id)
    return lDemand_wPositiveFlux, lDemand_wNegativeFlux

def scaleExchangeRxnsBoundaries(modelIn, dConfigParams=None):
    """Force the boundaries of exchange reaction to be -1000 or 1000 (not inf)"""
    for rxn in modelIn.reactions:
        if rxn.lower_bound == -float("inf"):
            rxn.lower_bound = -1000.0
        if rxn.upper_bound == float("inf"):
            rxn.upper_bound = 1000.0
    # modName = dConfigParams["modelName"]
    # if modName == "eciYali_sbmlValidated":
    #     for rxn in modelIn.reactions:
    #         if len(rxn.reactants) * len(rxn.products) == 0.0 and len(rxn.reactants) == 0 and  len(rxn.products) != 0:
    #             if rxn.lower_bound == -float("inf"):
    #                 rxn.lower_bound = -1000.0
    #             if rxn.upper_bound == float("inf"):
    #                 rxn.upper_bound = 1000.0
    # else:
    for rxn in modelIn.reactions:
        if len(rxn.reactants) * len(rxn.products) == 0.0:
            if rxn.lower_bound == -float("inf"):
                rxn.lower_bound = -1000.0
            if rxn.upper_bound == float("inf"):
                rxn.upper_bound = 1000.0
    return modelIn

def checkSomethingIsProducedFromNothing(modelIn, logFile, dConfigParams=None):
    modelName = dConfigParams["modelName"].split(".")[0]
    lDemand_wPositiveFlux, lDemand_wNegativeFlux = getDemand(modelIn)
    blockUptake(modelIn)
    outputFileNameFVA = "FVA_" + modelName + "_checkIfMakeSomethingFromNothing"
    faL.computeFVAandSaveFlux(modelIn, outputFileNameFVA, logFile)
    fva = pd.read_csv(os.path.join(OUTDIR, outputFileNameFVA + ".tsv"), sep = "\t")

    checkDemand_wPositiveFlux = fva[fva["Rxn"].isin(lDemand_wPositiveFlux)]
    mask = (checkDemand_wPositiveFlux["Min"] > 0) | (checkDemand_wPositiveFlux["Max"] > 0)

    lDfs = []
    if checkDemand_wPositiveFlux.loc[mask, :].empty is False:
        lDfs.append(checkDemand_wPositiveFlux.loc[mask, :])

    # Rxn  Min     Max   delta
    # 1274  r_1832  0.0  999.65  999.65
    # One can see an error about that H+ can be made even if no reactions were
    # unbalanced. Protons are particularly problematic since it is rather
    # arbitary at which pH the formulas are written for. For the purpose of this
    # analysis, the protons can be ignored and fixed later.

    checkDemand_wNegativeFlux = fva[fva["Rxn"].isin(lDemand_wNegativeFlux)]
    mask = (checkDemand_wNegativeFlux["Min"] < 0) | (checkDemand_wNegativeFlux["Max"] < 0)
    if checkDemand_wNegativeFlux.loc[mask, :].empty is False:
        lDfs.append(checkDemand_wNegativeFlux.loc[mask, :])


    pd.concat(lDfs).to_csv(os.path.join(OUTDIR, modelName + "_checkIfSomethingFromNothing_positiveReactions.tsv"), sep = "\t", index = False)

def getMet2KeggId(model):
    dMet2KeggId = {}
    for met in model.metabolites:
        if "kegg.compound" in met.annotation:
            if met.annotation["kegg.compound"] not in dMet2KeggId:
                dMet2KeggId[met.annotation["kegg.compound"]] = {met.compartment: met.id}
            elif met.annotation["kegg.compound"] in dMet2KeggId:
                dMet2KeggId[met.annotation["kegg.compound"]].update({met.compartment: met.id})
    return dMet2KeggId

def checkIfUptakeSomethingBlockingDemand(modelIn, logFile, dConfigParams=None):
    modelName = dConfigParams["modelName"].split(".")[0]
    lUptake_wPositiveFlux, lUptake_wNegativeFlux = getUptake(modelIn)
    blockDemand(modelIn)
    outputFileNameFVA = "FVA_" + modelName + "_checkIfUptakeSomethingBlockingDemand"
    faL.computeFVAandSaveFlux(modelIn, outputFileNameFVA, logFile)
    fva = pd.read_csv(os.path.join(OUTDIR, outputFileNameFVA + ".tsv"), sep = "\t")

    checkUpt_wPositiveFlux = fva[fva["Rxn"].isin(lUptake_wPositiveFlux)]
    lDfs = []
    mask = (checkUpt_wPositiveFlux["Min"] > 0) | (checkUpt_wPositiveFlux["Max"] > 0)
    if checkUpt_wPositiveFlux.loc[mask, :].empty is False:
        lDfs.append(checkUpt_wPositiveFlux.loc[mask, :])

    checkUpt_wNegativeFlux = fva[fva["Rxn"].isin(lUptake_wNegativeFlux)]
    mask= (checkUpt_wNegativeFlux["Min"] < 0) | (checkUpt_wNegativeFlux["Max"] < 0)
    if checkUpt_wNegativeFlux.loc[mask, :].empty is False:
        lDfs.append(checkUpt_wNegativeFlux.loc[mask, :])

    pd.concat(lDfs).to_csv(os.path.join(OUTDIR, modelName + "_checkIfUptakeSomethingBlockingDemand_positiveReactions.tsv"), sep = "\t", index = False)

def getReactionsFormingLoops(modelIn, logFile, dConfigParams=None):
    modelName = dConfigParams["modelName"].split(".")[0]
    blockUptakeAndDemand(modelIn)
    outputFileNameFVA = "FVA_" + modelName + "_checkLoops"
    faL.computeFVAandSaveFlux(modelIn, outputFileNameFVA, logFile)
    dfLoops = pd.read_csv(os.path.join(OUTDIR, outputFileNameFVA + ".tsv"), sep = "\t")
    mask = (dfLoops["Min"] > 0) | (dfLoops["Max"] > 0)
    outputFileNameFVA_loops = "FVA_" + modelName + "_rectionsFormingloops"
    dfLoops.loc[mask, :].to_csv(os.path.join(OUTDIR, outputFileNameFVA_loops + ".tsv"), sep = "\t", index = False)

def convertInEquationModelMetsIds2Formula(modelIn, reactionString, lb, ub):
    equationSplt = re.split('-->|<=>|<--|=', reactionString)

    if equationSplt[0] != '':
        dSubs = gL.dizReaProd(equationSplt[0].strip())
    else:
        dSubs = {}

    lSubs_converted2Formula = []
    for s in dSubs:
        mObj = modelIn.metabolites.get_by_id(s)
        if mObj.formula is not None:
            lSubs_converted2Formula.append(str(dSubs[s]) + " " + mObj.formula)

    ## fare stessa cosa sui prodotti
    if equationSplt[1] != '':
        dProds = gL.dizReaProd(equationSplt[1].strip())
    else:
        dProds = {}

    lProds_converted2Formula = []
    for p in dProds:
        mObj = modelIn.metabolites.get_by_id(p)
        if mObj.formula is not None:
            lProds_converted2Formula.append(str(dProds[p]) + " " + mObj.formula)

    if lb < 0 and ub > 0:
        formulaEquation = " <=> ".join([" + ".join(lSubs_converted2Formula), " + ".join(lProds_converted2Formula)])
    elif lb >= 0 and ub > 0:
        formulaEquation = " --> ".join([" + ".join(lSubs_converted2Formula), " + ".join(lProds_converted2Formula)])
    elif lb < 0 and ub <= 0:
        formulaEquation = " <-- ".join([" + ".join(lSubs_converted2Formula), " + ".join(lProds_converted2Formula)])
    return formulaEquation

def getExchInModelForKeggId(model, keggId, dNamesConversion):
    rId = None
    activeModelMet = None
    if keggId in dNamesConversion and "e" in dNamesConversion[keggId]:
        activeModelMet = dNamesConversion[keggId]["e"]
        lRxns = model.metabolites.get_by_id(activeModelMet).reactions
        for r in lRxns:
            if len(r.reactants) * len(r.products) == 0.0:
                rId = r.id
    return (rId, activeModelMet)
#####################
def curateRxnsAddedAfterBiologUpdate(modelIn,dConfigParams = None):
    ## spengo rxns con tag "added after the biolog update" quando queste rxns non sono in cerevisiae
    # lrxns2On = dConfigParams["lrxns2LeftOn"]
    # def isRxnInOrganism():
    organismCode = dConfigParams["organismCode"]
    for rxn in modelIn.reactions:
        if "NOTES" in rxn.notes and "added after the Biolog update" in rxn.notes["NOTES"] and rxn.gene_reaction_rule == "":
            if (len(rxn.reactants) * len(rxn.products) != 0.0) and ("transport" not in rxn.name.lower()):
                missingKeggRxnId = True
                missingEc = True

                if "kegg.reaction" in rxn.annotation:
                    missingKeggRxnId = False
                    keggId = rxn.annotation["kegg.reaction"]
                    try:
                        dizRxn = kL.getKeggInfo("rn:" + keggId)
                        insideTargetOrganism_ko = False
                        if "ORTHOLOGY" in dizRxn:
                            lKos = dizRxn["ORTHOLOGY"]
                            for ko in lKos:
                                dizKo = kL.getKeggInfo("ko:" + ko)
                                if "GENES" in dizKo and organismCode in dizKo["GENES"]:
                                    insideTargetOrganism = True

                        insideTargetOrganism_ec = False
                        if "ENZYME" in dizRxn:
                            lEcs = dizRxn["ENZYME"]
                            for ec in lEcs:
                                dizEc = kL.getKeggInfo("ec:" + ec)
                                if "GENES" in dizEc and organismCode in dizEc["GENES"]:
                                    insideTargetOrganism_ec = True

                        if insideTargetOrganism_ko is False and insideTargetOrganism_ec is False:
                            rxn.lower_bound = 0
                            rxn.upper_bound = 0
                    except:
                        print(f"KEGG rxn ID {keggId} not valid\n")
                elif "ec-code" in rxn.annotation:
                    ecNumber = rxn.annotation["ec-code"]
                    missingEc = False
                    try:
                        lEcs = rxn.annotation["ec-code"]
                        insideTargetOrganism_ec = False
                        for ec in lEcs:
                            dizEc = kL.getKeggInfo("ec:" + ec)
                            if "GENES" in dizEc and organismCode in dizEc["GENES"]:
                                insideTargetOrganism_ec = True

                        if insideTargetOrganism_ec is False:
                            rxn.lower_bound = 0
                            rxn.upper_bound = 0
                    except:
                        print(f"EC number {ecNumber} not valid\n")
                if missingKeggRxnId is False and missingEc is False:
                    rxn.lower_bound = 0
                    rxn.upper_bound = 0

def getGPRfromKEGG(orgId, lOrths = None, lEc = None):
    if lOrths is not None:
        lGenes = []
        for ko in lOrths:
            dizKo = kL.getKeggInfo("ko:" + ko)
            if "GENES" in dizKo and orgId.upper() in dizKo["GENES"] and dizKo["GENES"][orgId.upper()] not in lGenes:
                lGenes.append(dizKo["GENES"]["SCE"].split("(")[0].strip())
    elif lEc is not None:
        lGenes = []
        for ec in lEc:
            dizEc = kL.getKeggInfo("ec:" + ec)
            if "GENES" in dizEc and orgId.upper() in dizEc["GENES"] and dizEc["GENES"][orgId.upper()] not in lGenes:
                lGenes.append(dizEc["GENES"]["SCE"].split("(")[0].strip())
    gpr = " OR " .join(lGenes)
    return gpr

def getKeggRxnCompartmentInModel(lKeggId, dNamesConversion):
    lComps = []
    for met in lKeggId:
        met = met.strip()
        if met in dNamesConversion:
            lComps += list(dNamesConversion[met].keys())

    compartment = next(iter(Counter(lComps))) # take first element of the sorted list
    return compartment

def fromKeggRxn2ModelRxn(model, lKeggId, rxnCompartment, equation, dNamesConversion):
    # rxnCompartment = getKeggRxnCompartmentInModel(lKeggId, dNamesConversion)
    for met in lKeggId:
        if met in dNamesConversion and rxnCompartment in dNamesConversion[met]:
            equation = equation.replace(met, dNamesConversion[met][rxnCompartment])
        else:
            # crea nuovo metabolita
            dizMet = kL.getKeggInfo("cpd:" + met)
            nameMet = ""
            if "NAME" in dizMet:
                nameMet = dizMet["NAME"][0]
            formula = ""
            if "FORMULA" in dizMet:
                formula = dizMet["FORMULA"]
            addMets(model, met + "_" + rxnCompartment, nameMet, rxnCompartment, chemFormula = formula, fromKEGG = True)
            equation = equation.replace(met, met + "_" + rxnCompartment)
    return equation

def convert2ModelMet(keggMetId, compartment, modelFromWhichInfo, dNamesConversion = {}):
    print("keggMetId: ", keggMetId, " - compartment: ", compartment, "\n")
    if modelFromWhichInfo is not None and keggMetId in dNamesConversion and compartment in dNamesConversion[keggMetId]:
        print("OLD MET:\t", dNamesConversion[keggMetId][compartment])
        # print("dNamesConversion[keggMetId][compartment] ESISTE: ", dNamesConversion[keggMetId][compartment] in model2AddTo.reactions)
        # sys.exit()
        metObj = modelFromWhichInfo.metabolites.get_by_id(dNamesConversion[keggMetId][compartment])
        fromKEGG = False
        return metObj.id, metObj.name, metObj.formula, dNamesConversion, fromKEGG
    else:
        print("NEW MET")
        # crea nuovo metabolita
        dizMet = kL.getKeggInfo("cpd:" + keggMetId)
        nameMet = ""
        if "NAME" in dizMet:
            nameMet = dizMet["NAME"][0].strip(";").lower()
        formula = ""
        if "FORMULA" in dizMet:
            formula = dizMet["FORMULA"]
        fromKEGG = True
        # addMets(model2AddTo, keggMetId + "_" + compartment, nameMet, compartment, chemFormula = formula)
        # if keggMetId not in dNamesConversion:
        #     print(keggMetId, " not in the dict")
        #     dNamesConversion[keggMetId] = {compartment: keggMetId + "_" + compartment}
        # else:
        #     print(keggMetId, " already in the dict!")
        #     dNamesConversion[keggMetId].update({compartment: keggMetId + "_" + compartment})

        return keggMetId + "_" + compartment, nameMet, formula, dNamesConversion, fromKEGG

def fromKEGG2ModelIdAndAdd2Model(met, compartment, dNamesConversion, model2AddTo, modelTakenInfoFrom):
    metObjId, metObjName, metObjFormula, dNamesConversion, isKEGG = convert2ModelMet(met, compartment, modelTakenInfoFrom, dNamesConversion = dNamesConversion)
    print("RESULTS: ", metObjId,"\n", metObjName,"\n", metObjFormula,"\n", isKEGG)
    addMets(model2AddTo, metObjId, metObjName, compartment, chemFormula = metObjFormula, fromKEGG = isKEGG)
    print("ADDED MET\t", model2AddTo.metabolites.get_by_id(metObjId), "\n")
    return metObjId,dNamesConversion

def convertKeggRxns2ModelRxn(dizKeggRxn, model2AddTo, organismId,rxnIdentifier, compartment, modelTakenInfoFrom, dNamesConversion = None):
    reaction = cb.Reaction(rxnIdentifier + "_" + compartment)

    if "NAME" in dizKeggRxn:
        reaction.name = dizKeggRxn["NAME"][0]
    else:
        reaction.name = rxnIdentifier + "_" + compartment

    lEquation = dizKeggRxn["EQUATION"].split(' <=> ')
    dSubs = gL.dizReaProd(lEquation[0])
    dProds = gL.dizReaProd(lEquation[1])

    dKeggMets2ModelMets_subs = {}
    for met in dSubs:
        metObjId, dNamesConversion = fromKEGG2ModelIdAndAdd2Model(met, compartment, dNamesConversion, model2AddTo,  modelTakenInfoFrom)
        dKeggMets2ModelMets_subs[metObjId] = 1* dSubs[met]

    dKeggMets2ModelMets_prods = {}
    for met in dProds:
        metObjId, dNamesConversion = fromKEGG2ModelIdAndAdd2Model(met, compartment, dNamesConversion, model2AddTo,  modelTakenInfoFrom)
        dKeggMets2ModelMets_prods[metObjId] = 1* dProds[met]

    dKeggMets2ModelMets = {}

    for s in dKeggMets2ModelMets_subs:
        dKeggMets2ModelMets[model2AddTo.metabolites.get_by_id(s)] = -1 * dKeggMets2ModelMets_subs[s]
    for p in dKeggMets2ModelMets_prods:
        dKeggMets2ModelMets[model2AddTo.metabolites.get_by_id(p)] = 1 * dKeggMets2ModelMets_prods[p]

    reaction.add_metabolites(dKeggMets2ModelMets)
    if "ORTHOLOGY" in dizKeggRxn:
        gpr = getGPRfromKEGG(organismId, lOrths = list(dizKeggRxn['ORTHOLOGY'].keys()))
    elif "ENZYME" in dizKeggRxn:
        gpr = getGPRfromKEGG(organismId, lEc = dizKeggRxn['ENZYME'])

    return reaction, gpr, dNamesConversion

def convertEquation2ModelEquation(model2AddTo, rxnIdentifier,compartment, equation,dfConversionMetName2MetKegg, model2GetInfoFrom, dNamesConversion = {}):
    reaction = cb.Reaction(rxnIdentifier + "_" + compartment)
    reaction.name = rxnIdentifier

    lEquation = equation.split(' = ')
    dSubs = gL.dizReaProd(lEquation[0])
    dProds = gL.dizReaProd(lEquation[1])

    dMets2ModelMets_subs = {}
    for met in dSubs:
        metLower = met.lower()
        print("metLower: ", metLower)
        correspondingKeggId, dfConversionMetName2MetKegg= kL.getKeggMetId(metLower, dfConversionMetName2MetKegg)
        print("correspondingKeggId: \t", correspondingKeggId)
        metObjId, dNamesConversion = fromKEGG2ModelIdAndAdd2Model(correspondingKeggId,compartment, dNamesConversion, model2AddTo, model2GetInfoFrom)
        dMets2ModelMets_subs[metObjId] = 1* dSubs[met]

    print("\ndMets2ModelMets_subs\n", dMets2ModelMets_subs)

    dMets2ModelMets_prods = {}

    for met in dProds:
        metLower = met.lower()
        print("metLower: ", metLower)
        correspondingKeggId, dfConversionMetName2MetKegg= kL.getKeggMetId(metLower, dfConversionMetName2MetKegg)
        print("correspondingKeggId: \t", correspondingKeggId)
        metObjId, dNamesConversion = fromKEGG2ModelIdAndAdd2Model(correspondingKeggId, compartment, dNamesConversion, model2AddTo, model2GetInfoFrom)
        dMets2ModelMets_prods[metObjId] = 1* dProds[met]

    print("\dMets2ModelMets_prods\n", dMets2ModelMets_prods)

    dKeggMets2ModelMets = {}
    for s in dMets2ModelMets_subs:
        dKeggMets2ModelMets[model2AddTo.metabolites.get_by_id(s)] = -1 * dMets2ModelMets_subs[s]
    for p in dMets2ModelMets_prods:
        dKeggMets2ModelMets[model2AddTo.metabolites.get_by_id(p)] = 1 * dMets2ModelMets_prods[p]
    print("dKeggMets2ModelMets\t", dKeggMets2ModelMets)
    reaction.add_metabolites(dKeggMets2ModelMets)

    return reaction, dNamesConversion, dfConversionMetName2MetKegg


def addRxns(model,dNamesConversion, organismId, lb = -1000, ub = 1000, rxnIdentifier = None):
    if rxnIdentifier not in model.reactions:
        dizRxn = kL.getKeggInfo("rn:" + rxnIdentifier)
        reaction = cb.Reaction(rxnIdentifier)
        reaction.name = dizRxn["NAME"][0]
        reaction.lower_bound = lb
        reaction.upper_bound = ub
        equation = dizRxn["EQUATION"].split(' <=> ')
        if equation[0] != '':
            dSubs = gL.dizReaProd(equation[0])
        else:
            dSubs = {}
        if equation[1] != '':
            dProds = gL.dizReaProd(equation[1])
        else:
            dProds = {}

        equationConverted = fromKeggRxn2ModelRxn(model, list(dSubs.keys()) + list(dProds.keys()), dizRxn["EQUATION"], dNamesConversion)

        if equationConverted.split(' <=> ')[0] != '':
            dSubs = gL.dizReaProd(equationConverted.split(' <=> ')[0])
        else:
            dSubs = {}
        if equationConverted.split(' <=> ')[1] != '':
            dProds = gL.dizReaProd(equationConverted.split(' <=> ')[1])
        else:
            dProds = {}

        dEquation = {}
        for s in dSubs:
            dEquation[model.metabolites.get_by_id(s)] = -1 * dSubs[s]
        for p in dProds:
            dEquation[model.metabolites.get_by_id(p)] = 1 * dProds[p]
        reaction.add_metabolites(dEquation)
        model.add_reactions([reaction])
        if "ORTHOLOGY" in dizRxn:
            gpr = getGPRfromKEGG(organismId, lOrths = list(dizRxn['ORTHOLOGY'].keys()))
        elif "ENZYME" in dizRxn:
            gpr = getGPRfromKEGG(organismId, lEc = dizRxn['ENZYME'])
        # gpr = getGPRfromKEGG(model, list(dizRxn['ORTHOLOGY'].keys()))
        model.reactions.get_by_id(rxnIdentifier).gene_reaction_rule = gpr
    else:
        print(f"Attention: {rxnIdentifier} is already in the model!")
        sys.exit()

def addRxn2Model(model, rxnIdentifier, compartment, dNamesConversion, lb = -1000, ub = 1000):
    idRxn = rxnIdentifier + "_" + compartment
    if idRxn not in model.reactions:
        reaction = cb.Reaction(idRxn)
        reaction.name = nameMet + " transport between " + " and ".join(lCompartments)
        reaction.lower_bound = lb
        reaction.upper_bound = ub
        dSubs = gL.dizReaProd(lTransportedMets[0])
        dProds = gL.dizReaProd(lTransportedMets[1])
        dEquation = {}
        for s in dSubs:
            if s in model.metabolites:
                dEquation[model.metabolites.get_by_id(s)] = -1 * dSubs[s]
            else:
                addMets(model, p, p, "c", chemFormula = "", fromKEGG = False)

        for p in dProds:
            if p in model.metabolites:
                dEquation[model.metabolites.get_by_id(p)] = 1 * dProds[p]
            else:
                addMets(model, p, p, "c", chemFormula = "", fromKEGG = False)

        reaction.add_metabolites(dEquation)
        model.add_reactions([reaction])
        model.reactions.get_by_id(idRxn).gene_reaction_rule = ""
    else:
        print(f"Attention: {idRxn} is already in the model!")
        sys.exit()

    if rxnIdentifier not in model.reactions:
        dizRxn = kL.getKeggInfo("rn:" + rxnIdentifier)
        reaction = cb.Reaction(rxnIdentifier)
        reaction.name = dizRxn["NAME"][0]
        reaction.lower_bound = lb
        reaction.upper_bound = ub
        equation = dizRxn["EQUATION"].split(' <=> ')
        if equation[0] != '':
            dSubs = gL.dizReaProd(equation[0])
        else:
            dSubs = {}
        if equation[1] != '':
            dProds = gL.dizReaProd(equation[1])
        else:
            dProds = {}

        equationConverted = fromKeggRxn2ModelRxn(model, list(dSubs.keys()) + list(dProds.keys()), dizRxn["EQUATION"], dNamesConversion)

        if equationConverted.split(' <=> ')[0] != '':
            dSubs = gL.dizReaProd(equationConverted.split(' <=> ')[0])
        else:
            dSubs = {}
        if equationConverted.split(' <=> ')[1] != '':
            dProds = gL.dizReaProd(equationConverted.split(' <=> ')[1])
        else:
            dProds = {}

        dEquation = {}
        for s in dSubs:
            dEquation[model.metabolites.get_by_id(s)] = -1 * dSubs[s]
        for p in dProds:
            dEquation[model.metabolites.get_by_id(p)] = 1 * dProds[p]
        reaction.add_metabolites(dEquation)
        model.add_reactions([reaction])

        if "ORTHOLOGY" in dizRxn:
            gpr = getGPRfromKEGG(organismId, lOrths = list(dizRxn['ORTHOLOGY'].keys()))
        elif "ENZYME" in dizRxn:
            gpr = getGPRfromKEGG(orgId, lEc = dizRxn['ENZYME'])

        # gpr = getGPRfromKEGG(model, list(dizRxn['ORTHOLOGY'].keys()))
        model.reactions.get_by_id(rxnIdentifier).gene_reaction_rule = gpr
    else:
        print(f"Attention: {rxnIdentifier} is already in the model!")
        sys.exit()


# def addMets(model, row):
def addMets(model, metId, metName, metCompartment, chemFormula = "", fromKEGG = False):
    lMets_objects = []
    if metId not in model.metabolites:
        metObject = cb.Metabolite(metId, formula=chemFormula, name=metName, compartment=metCompartment)
        lMets_objects.append(metObject)
    model.add_metabolites(lMets_objects)
    if fromKEGG == True:
        print("metId KEGG: ", metId.split("_")[0])
        model.metabolites.get_by_id(metId).annotation["kegg.compound"] = metId.split("_")[0]
    print("ANNOTATION\n", model.metabolites.get_by_id(metId), " - ", model.metabolites.get_by_id(metId).annotation, "\n")
