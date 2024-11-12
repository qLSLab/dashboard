import os
import sys
import pandas as pd
import utils.manipulateModelLib as mmL
import utils.experimentalSetupLib as expSL
import utils.genericLib as gL
import cobra as cb
import itertools as itt
import re
import utils.fluxAnalysisLib as faL
import utils.keggLib as kL
from cobra import Model
import numpy as np

MODELDIR = gL.dDirs["mods"]
RAWDIR = gL.dDirs["raw"]
OUTDIR = gL.dDirs["out"]

def curateModel(modelIn, organismKeggCode, dfName2KeggId, modelPath, curateFile_newRxns, curateFile_bounds2Change,
            # curateFile_biologExp,
            curateFile_grp2Refine, dSpecies2Annotation):
    """
    Load base model or create it if it does not exist yet and, if necessary, curate it.
    """
    # modelFileName = dConfigParams["modelName"]
    # modelName = modelFileName.split(".")[0]

    # print("\nmodelPath: ", modelPath, "\n")
    modelFileName = os.path.split(modelPath)[1]
    # print("\nmodelFileName: ", modelFileName, "\n")
    modelName = os.path.splitext(modelFileName)[0]
    # print("\nmodelName: ", modelName, "\n")


    ## check if rxns are mass balanced
    mmL.addMissingChemicalFormulae(modelIn, dSpecies2Annotation)
    dUnbalancedRxns = cb.manipulation.validate.check_mass_balance(modelIn)
    # This function will return elements which violate mass balance
    lItems = []
    for unbRxn in dUnbalancedRxns:
        rxnObj = modelIn.reactions.get_by_id(unbRxn.id)
        lbOriginal = rxnObj.lower_bound
        ubOriginal = rxnObj.upper_bound
        reactionString = rxnObj.reaction
        formulaEquation = mmL.convertInEquationModelMetsIds2Formula(modelIn, reactionString, lbOriginal, ubOriginal)
        lItems.append([unbRxn.id, reactionString, formulaEquation, dUnbalancedRxns[unbRxn]])

    dfUnbRxns = pd.DataFrame(lItems, columns = ["RxnId", "Equation", "Equation_wFormula", "UnbalancedItems"])
    # gL.toLog(logFile, f"Reactions in {modelName}_unbalancedReactions.tsv need to be manually fixed\n")
    dfUnbRxns.to_csv(os.path.join(OUTDIR, modelName + "_unbalancedReactions.tsv"), sep = "\t", index = True)

    # Add exchange reactions for isolated metabolites
    dMet2Name = mmL.getMet2Name(modelIn)
    df = pd.DataFrame(dMet2Name.items(), columns=['Id', 'Name'])

    lIsolatedMets = []
    for met in modelIn.metabolites:
        isIsolatedMet = mmL.isIsolatedMet(modelIn, met.id)
        if isIsolatedMet is True:
            lIsolatedMets.append(met.id)

    # gL.toLog(logFile, f"Model includes {len(lIsolatedMets)} isolated metabolites")

    dfIsolatedMets = df[df["Id"].isin(lIsolatedMets)]
    dfIsolatedMetbyName = dfIsolatedMets.groupby("Name")["Id"].apply(list)

    ## Add simple transport reactions for isolated metabolites
    mmL.connectSameMetaboliteDifferentCompartment(modelIn, dfIsolatedMetbyName)

    dKeggId2ModelMet = mmL.getMet2KeggId(modelIn, dSpecies2Annotation)

    # print("curateFile_newRxns: ", curateFile_newRxns, "\n-----------\n")

    if curateFile_newRxns != "ATTENTION: You have not selected any file. Please fill and then export the table below.":
        curateFile_newRxns_fileName = curateFile_newRxns.split(":\t")[1].strip()
        # print("curateFile_newRxns: ", curateFile_newRxns_fileName, "\n-----------\n")
        if os.path.exists(curateFile_newRxns_fileName) is True:
            dfRxnsToAdd = pd.read_csv(curateFile_newRxns_fileName,
                sep="\t")
    elif os.path.exists(os.path.join(RAWDIR,"modelRefinement_newRxns2Add.tsv")) is True:
        dfRxnsToAdd = pd.read_csv(os.path.join(RAWDIR,"modelRefinement_newRxns2Add.tsv"), sep="\t")
    else:
        dfRxnsToAdd = pd.DataFrame()

    if dfRxnsToAdd.empty == False:
        for rowkeggRxnsToAdd in dfRxnsToAdd.itertuples():
            # print("\n------------\nrowkeggRxnsToAdd\n", rowkeggRxnsToAdd, "\n")
            # print("\n------------\organismKeggCode\n", organismKeggCode, "\n")
            rxnId = rowkeggRxnsToAdd.Rxn + "_" + rowkeggRxnsToAdd.Compartment
            if rxnId not in modelIn.reactions:
                dizRxn = kL.getKeggInfo("rn:" + rowkeggRxnsToAdd.Rxn)
                if len(dizRxn) != 0:
                    rxnObj, gprRule, dKeggId2ModelMet = mmL.convertKeggRxns2ModelRxn(dizRxn,dSpecies2Annotation, modelIn, organismKeggCode, rowkeggRxnsToAdd.Rxn, rowkeggRxnsToAdd.Compartment, modelIn, dNamesConversion = dKeggId2ModelMet)
                else:
                    rxnObj, dKeggId2ModelMet, dfName2KeggId = mmL.convertEquation2ModelEquation(modelIn, rowkeggRxnsToAdd.Rxn, rowkeggRxnsToAdd.Compartment, rowkeggRxnsToAdd.Equation, dfConversionMetName2MetKegg, modelIn, dNamesConversion= dKeggId2ModelMet)
                    if pd.isna(rowkeggRxnsToAdd.Gpr) is False:
                        gprRule = rowkeggRxnsToAdd.Gpr
                    else:
                        gprRule = ""

                modelIn.add_reactions([rxnObj])
                modelIn.reactions.get_by_id(rxnId).gene_reaction_rule = gprRule

    ## Close all remaining isolated mets with an exchange reaction
    # gL.toLog(logFile, "* Add exchange reactions for isolated metabolites")
    mmL.exchgRxns4IsolatedMet(modelIn, dfIsolatedMetbyName)

    # A first control is that the modelIn should not be able to produce any metabolites without uptake of metabolites.
    with modelIn:
        mmL.checkSomethingIsProducedFromNothing(modelIn, modelName)
        # gL.toLog(logFile, "Reactions in checkIfSomethingFromNothing_positiveReactions.tsv need to be manually fixed\n")

    # A second control is that the modelIn should not be able to consume any metabolites without demand of metabolites.
    with modelIn:
        mmL.checkIfUptakeSomethingBlockingDemand(modelIn,modelName)
        # gL.toLog(logFile, "Reactions in checkIfUptakeSomethingBlockingDemand_positiveReactions.tsv need to be manually fixed\n")

    ## check reactions forming loops
    with modelIn:
        mmL.getReactionsFormingLoops(modelIn, modelName)

    if curateFile_bounds2Change != "ATTENTION: You have not selected any file. Please fill and then export the table below.":
        curateFile_bounds2Change_fileName = curateFile_bounds2Change.split(":\t")[1].strip()
        if os.path.exists(curateFile_bounds2Change_fileName) is True:
            dfBounds = pd.read_csv(curateFile_bounds2Change_fileName,
                sep="\t",
                dtype= {"Lb": float, "Ub": float})

    elif os.path.exists(os.path.join(RAWDIR,"yeast8_refinement_bounds2Change.tsv")) is True:
        dfBounds = pd.read_csv(os.path.join(RAWDIR,"yeast8_refinement_bounds2Change.tsv"), sep="\t", dtype= {"Lb": float, "Ub": float})
    else:
        dfBounds = pd.DataFrame()

    # print("\n", dfBounds, "\n")
    if dfBounds.empty == False:
        for row in dfBounds.itertuples():
            modelIn.reactions.get_by_id(row.Id).lower_bound = row.Lb
            modelIn.reactions.get_by_id(row.Id).upper_bound = row.Ub

    ## aggiungo le gpr identificate
    if curateFile_grp2Refine != "ATTENTION: You have not selected any file. Please fill and then export the table below.":
        curateFile_grp2Refine_fileName = curateFile_grp2Refine.split(":\t")[1].strip()
        if os.path.exists(curateFile_grp2Refine_fileName) is True:
            dfnewGpr = pd.read_csv(
                curateFile_grp2Refine_fileName,
                sep="\t",
                dtype= {"Lb": int, "Ub": int})

    elif os.path.exists(os.path.join(RAWDIR,"yeast8_refinement_newGpr.tsv")) is True:
        dfnewGpr = pd.read_csv(os.path.join(RAWDIR,"yeast8_refinement_newGpr.tsv"), sep="\t", dtype= {"Lb": int, "Ub": int})
    else:
        dfnewGpr = pd.DataFrame()

    if dfnewGpr.empty == False:
        for row in dfnewGpr.itertuples():
            modelIn.reactions.get_by_id(row.Id).gene_reaction_rule = row.Gpr

    ## Change boundaries of exchange reactions from infinite to the default value (1000 in this case)
    # gL.toLog(logFile, "* Scale the exchange reactions boundaries")
    modelIn = mmL.scaleExchangeRxnsBoundaries(modelIn)

    print("\n\n-------------------FINE CURA--------------\n\n")

    # cb.io.write_sbml_model(modelIn, os.path.join(MODELDIR, modelName + "_postCuration.xml"))

    return modelIn


# def createNewModel(dfConversionMetName2MetKegg, model2GetInfoFrom, rxnsToAddFileName,exchangeRxnsFileName, experimentModel, keggOrgCode, dizKeggCompound2ModelCompound = {}):
# # def createNewModel(dfConversionMetName2MetKegg, model2GetInfoFrom, dizKeggCompound2ModelCompound = {}, dConfigParams=None):
#     # experimentModel = experimentModel.split(".")[0]
#     newModel = Model(experimentModel + " pathway")
#     # newModel = Model(dConfigParams["experimentModel"].split(".")[0] + " pathway")
#
#     # rxnsToAddFileName = dConfigParams["experimentRxsToAdd"]
#     # dfRxnsToAdd = pd.read_csv(os.path.join(RAWDIR, rxnsToAddFileName), sep = "\t")
#
#     if rxnsToAddFileName != "ATTENTION: You have not selected any file. Please fill and then export the table below.":
#         rxnsToAddFile = rxnsToAddFileName.split(":\t")[1].strip()
#         if os.path.exists(os.path.join(RAWDIR, rxnsToAddFile)) is True:
#             dfRxnsToAdd = pd.read_csv(os.path.join(RAWDIR, rxnsToAddFile),
#                 sep="\t")
#
#     elif os.path.exists(os.path.join(RAWDIR, experimentModel + "_internalRxns.tsv")) is True:
#         dfRxnsToAdd = pd.read_csv(os.path.join(RAWDIR, experimentModel + "_internalRxns.tsv"), sep="\t")
#     else:
#         dfRxnsToAdd = pd.DataFrame()
#
#     if dfRxnsToAdd.empty == False:
#         for rowkeggRxnsToAdd in dfRxnsToAdd.itertuples():
#             rxnId = rowkeggRxnsToAdd.Rxn + "_" + rowkeggRxnsToAdd.Compartment
#             if rxnId not in newModel.reactions:
#                 dizRxn = kL.getKeggInfo("rn:" + rowkeggRxnsToAdd.Rxn)
#                 if len(dizRxn) != 0:
#                     rxnObj, gprRule, dizKeggCompound2ModelCompound = mmL.convertKeggRxns2ModelRxn(dizRxn, newModel, keggOrgCode, rowkeggRxnsToAdd.Rxn, rowkeggRxnsToAdd.Compartment, model2GetInfoFrom, dNamesConversion = dizKeggCompound2ModelCompound)
#                     # rxnObj, gprRule, dizKeggCompound2ModelCompound = mmL.convertKeggRxns2ModelRxn(dizRxn, newModel, dConfigParams["organismCode"], rowkeggRxnsToAdd.Rxn, rowkeggRxnsToAdd.Compartment, model2GetInfoFrom, dNamesConversion = dizKeggCompound2ModelCompound)
#                 else:
#                     rxnObj, dizKeggCompound2ModelCompound, dfConversionMetName2MetKegg = mmL.convertEquation2ModelEquation(newModel, rowkeggRxnsToAdd.Rxn, rowkeggRxnsToAdd.Compartment, rowkeggRxnsToAdd.Equation, dfConversionMetName2MetKegg, model2GetInfoFrom, dNamesConversion= dizKeggCompound2ModelCompound)
#                     if pd.isna(rowkeggRxnsToAdd.Gpr) is False:
#                         gprRule = rowkeggRxnsToAdd.Gpr
#                     else:
#                         gprRule = ""
#
#                 newModel.add_reactions([rxnObj])
#                 newModel.reactions.get_by_id(rxnId).gene_reaction_rule = gprRule
#
#     # addExchRxnsForTheseMetsFileName = dConfigParams["exchangeRxns"]
#     # dfaddExchRxnsForTheseMets = pd.read_csv(os.path.join(RAWDIR, addExchRxnsForTheseMetsFileName), sep = "\t")
#     #####################################################
#     #####################################################
#     if exchangeRxnsFileName != "ATTENTION: You have not selected any file. Please fill and then export the table below.":
#         exchangeRxnsFile = exchangeRxnsFileName.split("\t")[1].strip()
#         if os.path.exists(os.path.join(RAWDIR, exchangeRxnsFile)) is True:
#             dfaddExchRxnsForTheseMets = pd.read_csv(os.path.join(RAWDIR, exchangeRxnsFile),
#                 sep="\t")
#
#     elif os.path.exists(os.path.join(RAWDIR, experimentModel + "_exchRxns.tsv")) is True:
#         dfaddExchRxnsForTheseMets = pd.read_csv(os.path.join(RAWDIR, experimentModel + "_exchRxns.tsv"), sep="\t")
#     else:
#         dfaddExchRxnsForTheseMets = pd.DataFrame()
#
#     if dfaddExchRxnsForTheseMets.empty == False:
#         for rowaddExch in dfaddExchRxnsForTheseMets.itertuples():
#             # metObjId, metObjName, metObjFormula, dizKeggCompound2ModelCompound, isKEGG = mmL.convert2ModelMet(rowaddExch.Met, rowaddExch.Compartment, dNamesConversion = dizKeggCompound2ModelCompound, modelFromWhichInfo = model2GetInfoFrom)
#             # mmL.addMets(newModel, metObjId, metObjName, rowaddExch.Compartment, chemFormula = metObjFormula, fromKEGG = isKEGG)
#             metObjId, dizKeggCompound2ModelCompound = mmL.fromKEGG2ModelIdAndAdd2Model(rowaddExch.Met, rowaddExch.Compartment, dizKeggCompound2ModelCompound,newModel, model2GetInfoFrom)
#             mmL.addExchRxnForMet(newModel, newModel.metabolites.get_by_id(metObjId), newLB=rowaddExch.Lb, newUB=rowaddExch.Ub)
#
#     # Add exchange reactions for isolated metabolites
#     dMet2Name = mmL.getMet2Name(newModel)
#     df = pd.DataFrame(dMet2Name.items(), columns=['Id', 'Name'])
#
#     lIsolatedMets = []
#     for met in newModel.metabolites:
#         isIsolatedMet = mmL.isIsolatedMet(newModel, met.id)
#         if isIsolatedMet is True:
#             lIsolatedMets.append(met.id)
#
#     # print(f"newModel includes {len(lIsolatedMets)} isolated metabolites\n")
#
#     dfIsolatedMets = df[df["Id"].isin(lIsolatedMets)]
#     dfIsolatedMetbyName = dfIsolatedMets.groupby("Name")["Id"].apply(list)
#
#     ## Add simple transport reactions for isolated metabolites
#     mmL.connectSameMetaboliteDifferentCompartment(newModel, dfIsolatedMetbyName)
#
#     ## export newModel
#     # cb.io.write_sbml_model(newModel, os.path.join(MODELDIR, dConfigParams["experimentModel"]))
#
#     # print("\n\nFINE CREATION!!!\n")
#
#     return newModel, dizKeggCompound2ModelCompound, dfConversionMetName2MetKegg

def createNewModel_part1_addInternalRxns(dfConversionMetName2MetKegg, model2GetInfoFrom, rxnsToAddFileName, experimentModel, keggOrgCode, dSpecies2Annotation, dizKeggCompound2ModelCompound = {}):
    newModel = Model(experimentModel + " pathway")

    if rxnsToAddFileName != "ATTENTION: You have not selected any file. Please fill and then export the table below.":
        rxnsToAddFile = rxnsToAddFileName.split(":\t")[1].strip()
        if os.path.exists(os.path.join(RAWDIR, rxnsToAddFile)) is True:
            dfRxnsToAdd = pd.read_csv(os.path.join(RAWDIR, rxnsToAddFile),
                sep="\t")

    elif os.path.exists(os.path.join(RAWDIR, experimentModel + "_internalRxns.tsv")) is True:
        dfRxnsToAdd = pd.read_csv(os.path.join(RAWDIR, experimentModel + "_internalRxns.tsv"), sep="\t")
    else:
        dfRxnsToAdd = pd.DataFrame()

    if dfRxnsToAdd.empty == False:
        for rowkeggRxnsToAdd in dfRxnsToAdd.itertuples():
            dAnnotation = {}
            rxnId = rowkeggRxnsToAdd.Rxn + "_" + rowkeggRxnsToAdd.Compartment
            if rxnId not in newModel.reactions:
                dizRxn = kL.getKeggInfo("rn:" + rowkeggRxnsToAdd.Rxn)
                if len(dizRxn) != 0:
                    rxnObj, gprRule, dizKeggCompound2ModelCompound = mmL.convertKeggRxns2ModelRxn(dizRxn, dSpecies2Annotation, newModel, keggOrgCode, rowkeggRxnsToAdd.Rxn, rowkeggRxnsToAdd.Compartment, model2GetInfoFrom, dNamesConversion = dizKeggCompound2ModelCompound)
                    dAnnotation["kegg.reaction"] = rowkeggRxnsToAdd.Rxn
                    # rxnObj, gprRule, dizKeggCompound2ModelCompound = mmL.convertKeggRxns2ModelRxn(dizRxn, newModel, dConfigParams["organismCode"], rowkeggRxnsToAdd.Rxn, rowkeggRxnsToAdd.Compartment, model2GetInfoFrom, dNamesConversion = dizKeggCompound2ModelCompound)
                else:
                    rxnObj, dizKeggCompound2ModelCompound, dfConversionMetName2MetKegg = mmL.convertEquation2ModelEquation(newModel, rowkeggRxnsToAdd.Rxn, rowkeggRxnsToAdd.Compartment, rowkeggRxnsToAdd.Equation, dfConversionMetName2MetKegg, model2GetInfoFrom, dNamesConversion= dizKeggCompound2ModelCompound)
                    if pd.isna(rowkeggRxnsToAdd.Gpr) is False:
                        gprRule = rowkeggRxnsToAdd.Gpr
                    else:
                        gprRule = ""

                newModel.add_reactions([rxnObj])
                newModel.reactions.get_by_id(rxnId).gene_reaction_rule = gprRule
                newModel.reactions.get_by_id(rxnId).annotation = dAnnotation

    return newModel, dizKeggCompound2ModelCompound, dfConversionMetName2MetKegg


def createNewModel_part2_addExchangeRxns(newModel):
    lExchangeRxns = []
    lColumns = ["MetId", "MetName", "Lb", "Ub"]
    lIsolatedMets = []
    for met in newModel.metabolites:
        isolatedMet, addDemand, addSink = mmL.isIsolatedMet_and_defineIsolatedMetType(newModel, met)
        if isolatedMet is True:
            lIsolatedMets.append(met.id)
            if addDemand is True:
                lb = 0
                ub = 1000.0
                # mmL.addExchRxnForMet(newModel, met, lb, ub)
                lExchangeRxns.append([met.id, met.name,lb, ub])
            elif addSink is True:
                lb = -1000.0
                ub = 0
                # mmL.addExchRxnForMet(newModel, met, lb, ub)
                lExchangeRxns.append([met.id, met.name, lb, ub])

    dfProposedExchangeRxns = pd.DataFrame(lExchangeRxns, columns = lColumns)
    dfProposedExchangeRxns.to_csv(os.path.join(RAWDIR, "exchangeReactions_newModel.tsv"), sep ="\t", index=False)

    ## Add simple transport reactions for isolated metabolites
    dMet2Name = mmL.getMet2Name(newModel)
    df = pd.DataFrame(dMet2Name.items(), columns=['Id', 'Name'])
    dfIsolatedMets = df[df["Id"].isin(lIsolatedMets)]
    dfIsolatedMetbyName = dfIsolatedMets.groupby("Name")["Id"].apply(list)
    mmL.connectSameMetaboliteDifferentCompartment(newModel, dfIsolatedMetbyName)

    return newModel, dfProposedExchangeRxns


def createNewModel_part3_addRevisedExchangeRxns(model, dfProposedExchangeRxns):
    for row in dfProposedExchangeRxns.itertuples():
        met = model.metabolites.get_by_id(row.MetId)
        mmL.addExchRxnForMet(model, met, float(row.Lb), float(row.Ub))

    return model

##############################
### OLD PART #################
##############################


def getModel(lHeterologousPaths, logFile = None, isolatedMetsClosing="no", recreate = "no", curateBaseModel = "no", dConfigParams=None):
    """
    Load model or create it if it does not exist yet.
    """
    originalModelName = dConfigParams["modelName"]
    mediumName = dConfigParams["medium"]
    lconditions = dConfigParams["condition"]
    outModelFileName = originalModelName + '_' + '_'.join(lHeterologousPaths) + '_' + mediumName + '_' + '_'.join(lconditions)

    if os.path.exists(os.path.join(MODELDIR, outModelFileName + ".xml")) is False:
        gL.toLog(logFile, "Model does not exist yet. It will be now created.")
        ## get list of heterologous pathways; for each pathway merge its model with the original model
        model = cb.io.read_sbml_model(os.path.join(MODELDIR, originalModelName + ".xml"))

        # if curateBaseModel == "yes":
        #     if os.path.exists(os.path.join(RAWDIR, originalModelName + "_refinement_newRxns2Add.tsv")) is True:
        #         dfRxns = pd.read_csv(
        #             os.path.join(RAWDIR, originalModelName + "_refinement_newRxns2Add.tsv"),
        #             sep="\t",
        #             dtype= {"Lb": int, "Ub": int, "GPR": str, "Id":str,	"Name": str, "Equation":str})
        #         # dfRxns = dfRxns.astype({"GPR":str})
        #         dfRxns["GPR"] = dfRxns["GPR"].astype(str)
        #
        #         for row in dfRxns.itertuples():
        #             mmL.addRxns(model, row)


        if curateBaseModel == "yes":
            if os.path.exists(os.path.join(RAWDIR, originalModelName + "_refinement_newRxns2Add_1.tsv")) is True:
                dfRxns = pd.read_csv(
                    os.path.join(RAWDIR, originalModelName + "_refinement_newRxns2Add.tsv"),
                    sep="\t",
                    dtype= {"Lb": int, "Ub": int, "GPR": str, "Id":str,	"Name": str, "Equation":str})
                # dfRxns = dfRxns.astype({"GPR":str})
                dfRxns["GPR"] = dfRxns["GPR"].astype(str)

                for row in dfRxns.itertuples():
                    mmL.addRxns(model, row)
        if originalModelName in ["yeast8", "yeast9"]:
            lRxns = [("R02235", "dihydrofolate:NAD+ oxidoreductase", "s_0719 + s_1203 + s_0794 = s_0625 + s_1198", -1000, 1000, "YOR236W"),
                    ("R02236","dihydrofolate:NADP+ oxidoreductase", "s_0719 + s_1212 + s_0794 = s_0625 + s_1207", -1000, 1000, "YOR236W"),
                    ("R00936","5,6,7,8-tetrahydrofolate:NAD+ oxidoreductase", "s_0625 + s_1203 + s_0794 = s_1487 + s_1198",-1000,1000,"YOR236W"),
                    ("hexadecanalTransport","hexadecanal transport", "s_0825 = s_0826",-1000,1000,"")]

            dfRxns = pd.DataFrame(lRxns, columns = ['Id', 'Name', 'Equation', 'Lb', 'Ub', 'GPR'])
            for r in dfRxns.itertuples():
                mmL.addRxns(model, r)

            # ### add lb e ub per queste rxns
            # model.reactions.get_by_id("r_4167").lower_bound = -1000
            # model.reactions.get_by_id("r_4167").upper_bound = 0

        # Add exchange reactions for isolated metabolites
        if isolatedMetsClosing == 'yes':
            gL.toLog(logFile, "* Add exchange reactions for isolated metabolites")
            expSL.exchgRxns4IsolatedMet(model)

        gL.toLog(logFile, "* Add heterologous pathways")
        for exp in lHeterologousPaths:
            if os.path.exists(os.path.join(MODELDIR, exp + ".xml")) is False or recreate == "yes":
                pathModel = expSL.createHeterologousModel(exp)
                model = model.merge(pathModel)
            elif os.path.exists(os.path.join(MODELDIR, exp + ".xml")) is True or recreate == "no":
                pathModel = cb.io.read_sbml_model(os.path.join(MODELDIR, exp + ".xml"))
                model = model.merge(pathModel)

        # if originalModelName in ["yeast8", "yeast-GEM_9"]:
        # Add to model uptake reactions of wet medium not already included in the model
        newRxnsFileName = 'newMediumRxns_' + mediumName + '_' + '_'.join(lconditions) + ".tsv"
        newMetsFileName = 'newMediumMets_' + mediumName + '_' + '_'.join(lconditions) + ".tsv"

        gL.toLog(logFile, f"* Add medium reactions not already included in the original model. Reactions are included into the file called {newRxnsFileName}, metabolites are included into the file called {newMetsFileName}")

        dictTypes = {
            "Id": str,
            "Name": str,
            "Compartment": str
        }

        if os.path.exists( os.path.join(RAWDIR, newMetsFileName)) is True:
            dfNewMediumMets = pd.read_csv(
                      os.path.join(RAWDIR, newMetsFileName), sep="\t",
                      dtype=dictTypes)

            for newMet in dfNewMediumMets.itertuples():
                mmL.addMets(model, newMet)

        dictTypes = {
            "Id": str,
            "Name": str,
            "Equation": str,
            "Lb": int,
            "Ub": int,
            "GPR": str
        }

        if os.path.exists( os.path.join(RAWDIR, newRxnsFileName)) is True:
            dfNewMediumRxns = pd.read_csv(
                      os.path.join(RAWDIR, newRxnsFileName),
                      sep="\t",
                      dtype=dictTypes)

            dfNewMediumRxns = dfNewMediumRxns.fillna('')
            for newRxn in dfNewMediumRxns.itertuples():
                mmL.addRxns(model, newRxn)


        gL.toLog(logFile,
            "Model construction is now complete and the SBML file will be saved into the models directory with the name {outModelFileName}.xml\n")
        cb.io.write_sbml_model(model, os.path.join(MODELDIR, outModelFileName + ".xml"))
    elif os.path.exists(os.path.join(MODELDIR, outModelFileName + ".xml")) is True:
        gL.toLog(logFile, "Model already exist. It will be now loaded.\n")
        model = cb.io.read_sbml_model(os.path.join(MODELDIR, outModelFileName + ".xml"))
    return model

def constrainWwetData(modelIn, outModelFileName, dfExtracFluxes_consumption, dfExtracFluxes_production, dfMediumComposition, tStamp = None, dConfigParams=None, logFile = None,  tPoint1 = 0, tPoint2 = 0, expName = None, curateBaseModel = "no"):
    originalModelName = dConfigParams["modelName"]

    rowT1_c = dfExtracFluxes_consumption[dfExtracFluxes_consumption.Time == tPoint1]
    rowT2_c = dfExtracFluxes_consumption[dfExtracFluxes_consumption.Time == tPoint2]

    rowT1_p = dfExtracFluxes_production[dfExtracFluxes_production.Time == tPoint1]
    rowT2_p = dfExtracFluxes_production[dfExtracFluxes_production.Time == tPoint2]

    lMeasuredMets_consumption = []
    for col in dfExtracFluxes_consumption.columns:
        measuredMet = col.split("_")
        if measuredMet[0] not in ["Time", "OD"] and measuredMet[0] not in lMeasuredMets_consumption:
            lMeasuredMets_consumption.append(measuredMet[0])

    lMeasuredMets_production = []
    for col in dfExtracFluxes_production.columns:
        measuredMet = col.split("_")
        if measuredMet[0] not in ["Time", "OD"] and measuredMet[0] not in lMeasuredMets_production:
            lMeasuredMets_production.append(measuredMet[0])

    gL.toLog(logFile, f"\n>>> Experimental time points: {tPoint2}h and {tPoint1}h")
    with modelIn:

        # if dConfigParams["modelName"] in ["yeast8", "yeast-GEM_9"]:
        if curateBaseModel == "yes":
            # if os.path.exists(os.path.join(RAWDIR, originalModelName + "_refinement_bounds2Change_1.tsv")) is True:
            #     dfBounds = pd.read_csv(
            #         os.path.join(RAWDIR, originalModelName + "_refinement_bounds2Change.tsv"),
            #         sep="\t",
            #         dtype= {"Lb": float, "Ub": float})
            #     for row in dfBounds.itertuples():
            #         print(type(row.Lb), type(row.Ub))
            #         modelIn.reactions.get_by_id(row.Id).lower_bound = row.Lb
            #         modelIn.reactions.get_by_id(row.Id).upper_bound = row.Ub
            #         print(row.Id, modelIn.reactions.get_by_id(row.Id).lower_bound, modelIn.reactions.get_by_id(row.Id).upper_bound)

            ### r_4315 solo backward da dati termodinamica
            # modelIn.reactions.get_by_id("r_4315").lower_bound = -1000
            # modelIn.reactions.get_by_id("r_4315").upper_bound = 0

            # ## spengo queste rxns che coinvolgono dead-end metabolites
            # lRxns2Off = ["r_1543", "r_1554", "r_0070", "r_1640", "r_1008", "r_4363", "r_4364", "r_4365",
            # "r_4382", "r_4383", "r_4499", "r_4384", "r_4385", "r_4462", "r_4477", "r_4568", "r_4569"]
            # for rxn in lRxns2Off:
            #     modelIn.reactions.get_by_id(rxn).lower_bound = 0
            #     modelIn.reactions.get_by_id(rxn).upper_bound = 0

            if os.path.exists(os.path.join(RAWDIR, originalModelName + "_refinement_bounds2Change.tsv")) is True:
                dfBounds = pd.read_csv(
                    os.path.join(RAWDIR, originalModelName + "_refinement_bounds2Change.tsv"),
                    sep="\t",
                    dtype= {"Lb": float, "Ub": float})
                for row in dfBounds.itertuples():
                    modelIn.reactions.get_by_id(row.Id).lower_bound = row.Lb
                    modelIn.reactions.get_by_id(row.Id).upper_bound = row.Ub

            ## spengo rxns con tag "added after the biolog update" quando queste rxns non sono in cerevisiae tranne quelle isolate manualmente nella lista lrxns2On
            lrxns2On = dConfigParams["lrxns2LeftOn"]
            if len(lrxns2On) != 0:
                for rxn in modelIn.reactions:
                    if "NOTES" in rxn.notes and "added after the Biolog update" in rxn.notes["NOTES"] and rxn.gene_reaction_rule == "" and rxn.id not in lrxns2On:
                        rxn.lower_bound = 0
                        rxn.upper_bound = 0

            ## aggiungo le gpr identificate
            if os.path.exists(os.path.join(RAWDIR, originalModelName + "_refinement_newGpr.tsv")) is True:
                dfnewGpr = pd.read_csv(
                    os.path.join(RAWDIR, originalModelName + "_refinement_newGpr.tsv"),
                    sep="\t",
                    dtype= {"Lb": int, "Ub": int})
                for row in dfnewGpr.itertuples():
                    modelIn.reactions.get_by_id(row.Id).gene_reaction_rule = row.Gpr

            # dRxn2GPR = {"r_4355": "YHR208W and YJR148W", "r_4372": "YGR247W and YKR034W", "r_4388": "YGR247W",
            # "r_4407": "YGR208W", "r_4412": "YGR247W"}
            # for reac in dRxn2GPR:
            #     modelIn.reactions.get_by_id(reac).gene_reaction_rule = dRxn2GPR[reac]

            # ## LOOPS
            # ## rxns curate manualmente nella direzione
            # lForwardDirections = ["r_2194", "r_2195", "r_2196", "r_2197", "r_2198", "r_2199", "r_2200", "r_2201", "r_2202", "r_2203", "r_2204", "r_2205", "r_2211", "r_2212", "r_2214", "r_2215", "r_2217", "r_2218", "r_2620", "r_2621", "r_2622", "r_2623", "r_2624", "r_2625", "r_2626", "r_2627", "r_2628", "r_2629", "r_2630", "r_2631", "r_2632", "r_2633", "r_2634", "r_2635", "r_2636", "r_2637", "r_2638", "r_2639", "r_2640", "r_2641", "r_2642", "r_2643", "r_2644", "r_2645", "r_2646", "r_2647", "r_2648", "r_2649", "r_2650", "r_2651", "r_2652", "r_2653", "r_2654", "r_2655", "r_2656", "r_2657", "r_2658", "r_2659", "r_2660", "r_2661", "r_2662", "r_2663", "r_2664", "r_2665", "r_2666", "r_2667", "r_2668", "r_2669", "r_2670", "r_2671", "r_2672", "r_2673", "r_2674", "r_2675", "r_2676", "r_2677", "r_2678", "r_2679", "r_2680", "r_2681", "r_2682", "r_2683", "r_2684", "r_2685", "r_2686", "r_2687", "r_2688", "r_2689", "r_2690", "r_2691", "r_2692", "r_2693", "r_2694", "r_2695", "r_2696", "r_2697", "r_2698", "r_2699", "r_2700", "r_2701", "r_2702", "r_2703", "r_2704", "r_2705", "r_2706", "r_2707", "r_2708", "r_2709", "r_2710", "r_2711", "r_2712", "r_2713", "r_2714", "r_2715", "r_2716", "r_2717", "r_2718", "r_2719", "r_2720", "r_2721", "r_2722", "r_2723", "r_2724", "r_2725", "r_2726", "r_2727", "r_2728", "r_2729", "r_2730", "r_2731", "r_2732", "r_2733", "r_2734", "r_2735", "r_2736", "r_2737", "r_2738", "r_2739", "r_2740", "r_2741", "r_2742", "r_2743", "r_2744", "r_2745", "r_2746", "r_2747", "r_2748", "r_2749", "r_2750", "r_2751", "r_2752", "r_2753", "r_2754", "r_2755", "r_2756", "r_2757", "r_2758", "r_2759", "r_2760", "r_2761", "r_2762", "r_2763", "r_2764", "r_2765", "r_2766", "r_2767", "r_2768", "r_2769", "r_2770", "r_2771", "r_2772", "r_2773", "r_2774", "r_2775", "r_2776", "r_2777", "r_2778", "r_2779", "r_2780", "r_2781", "r_2782", "r_2783", "r_2784", "r_2785", "r_2786", "r_2787", "r_2788", "r_2789", "r_2790", "r_2791", "r_2792", "r_2793", "r_2794", "r_2795", "r_2796", "r_2797", "r_2798", "r_2799", "r_2800", "r_2801", "r_2802", "r_2803", "r_2804", "r_2805", "r_2806", "r_2807", "r_2808", "r_2809", "r_2810", "r_2811", "r_2884", "r_2885", "r_2886", "r_2887", "r_2888", "r_2889", "r_2890", "r_2891", "r_2892", "r_2893", "r_2894", "r_2895", "r_2896", "r_2897", "r_2898", "r_2899", "r_2900", "r_2901", "r_2902", "r_2903", "r_2904", "r_2905", "r_2906", "r_2907", "r_2908", "r_2909", "r_2910", "r_2911", "r_2912", "r_2913", "r_2914", "r_2915", "r_2916", "r_2917", "r_2918", "r_2919", "r_2920", "r_2921", "r_2922", "r_2923", "r_2924", "r_2925", "r_2926", "r_2927", "r_2928", "r_2929", "r_2930", "r_2931", "r_2932", "r_2933", "r_2934", "r_2935", "r_2936", "r_2937", "r_2938", "r_2939", "r_2940", "r_2941", "r_2942", "r_2943", "r_2944", "r_2945", "r_2946", "r_2947", "r_2948", "r_2949", "r_2950", "r_2951", "r_2952", "r_2953", "r_2954", "r_2955", "r_2956", "r_2957", "r_2958", "r_2959", "r_2960", "r_2961", "r_2962", "r_2963", "r_2964", "r_2965", "r_2966", "r_2967", "r_2968", "r_2969", "r_2970", "r_2971", "r_2972", "r_2973", "r_2974", "r_2975", "r_2976", "r_2977", "r_2978", "r_2979", "r_2980", "r_2981", "r_2982", "r_2983", "r_2984", "r_2985", "r_2986", "r_2987", "r_2988", "r_2989", "r_2990", "r_2991", "r_2992", "r_2993", "r_2994", "r_2995", "r_2996", "r_2997", "r_2998", "r_2999", "r_3000", "r_3001", "r_3002", "r_3003", "r_3004", "r_3005", "r_3006", "r_3007", "r_3008", "r_3009", "r_3010", "r_3011", "r_4161", "r_4162", "r_4163",
            # "r_1148", "r_1667", "r_1202", "r_1208", "r_1220", "r_4280", "r_4281", "r_4323", "r_4324", "r_0302", "r_0280","r_0303", "r_2305", "r_0672", "r_4274",
            # "r_1084", "r_1665", "r_1000", "r_0454", "r_0455", "r_1112", "r_0714"]
            #
            # for rxn in lForwardDirections:
            #     modelIn.reactions.get_by_id(rxn).lower_bound = 0
            #     modelIn.reactions.get_by_id(rxn).upper_bound = 1000

            # lBackDirection = ["r_2119",
            # #"r_1063",
            # "r_0815", "r_4236", "r_4264", "r_1831"] # -1000,0
            #
            # for rxn in lBackDirection:
            #     modelIn.reactions.get_by_id(rxn).lower_bound = -1000
            #     modelIn.reactions.get_by_id(rxn).upper_bound = 0
            #
            # lOff = ["r_1758", "r_4226", "r_4262", "r_1652"] # 0,0 r_1760, r_4566 era in lista ma nel nuovo modello l'hanno tolta
            #
            # for rxn in lOff:
            #     modelIn.reactions.get_by_id(rxn).lower_bound = 0
            #     modelIn.reactions.get_by_id(rxn).upper_bound = 0
            #
            # # set constraints 0,1000 to r_4046 instead of 0.7 set in the original model
            # gL.toLog(logFile, "r_4046 with boundaries 0,1000")
            # modelIn.reactions.get_by_id("r_4046").lower_bound = 0
            # modelIn.reactions.get_by_id("r_4046").upper_bound = 1000

        ## Change boundaries of exchange reactions from infinite to the default value (1000 in this case)
        gL.toLog(logFile, "* Scale the exchange reactions boundaries")
        modelIn = expSL.scaleExchangeRxnsBoundaries(modelIn, dConfigParams)

        # Block all uptake reactions
        if dConfigParams["blockUptakes"] == "yes":
            gL.toLog(logFile,
                "* Block all uptake reactions (to block any undesired uptake)")
            modelIn = expSL.switchOffAllSinkRxns(modelIn, dConfigParams)

        # if dConfigParams["modelName"] in ["yeast8", "yeast-GEM_9"]:
        #     # r_4291 according to reactome needs to be irreversible
        #     gL.toLog(logFile, "r_4291 with boundaries 0,1000")
        #     modelIn.reactions.get_by_id("r_4291").lower_bound = 0
        #     modelIn.reactions.get_by_id("r_4291").upper_bound = 1000

        # add medium concentrations
        gL.toLog(logFile, "* Constrain modelIn with experimental concentrations of medium components")
        modelIn = expSL.addMediumData(modelIn, dfMediumComposition,dConfigParams = dConfigParams)
        for rxnMedium in dfMediumComposition.itertuples():
            rxn = modelIn.reactions.get_by_id(rxnMedium.RxnID)
            gL.toLog(
                logFile, f"\t{rxn.id} - {rxn.name} - ({rxn.lower_bound}, {rxn.upper_bound})"
            )


        gL.toLog(logFile, "\n")

        # Set extracellular fluxes
        gL.toLog(
            logFile,
            "* Constrain extracellular fluxes with variation between two consecutive time points of experimental consumption and production data")

        modelIn = expSL.setExtracFluxes(lMeasuredMets_consumption, lMeasuredMets_production, tPoint1, tPoint2, rowT1_c, rowT2_c, rowT1_p, rowT2_p, modelIn, dConfigParams=dConfigParams, experiment = expName)

        # if "EG" in dConfigParams['condition'] and dConfigParams["modelName"] in ["yeast8", "yeast-GEM_9"]:
        #     lExtracFluxs = [dConfigParams['glucose'], dConfigParams['ethanol'], dConfigParams['glycerol'], dConfigParams['acetate'], dConfigParams['eg']]
        # elif "EG" not in dConfigParams['condition'] and dConfigParams["modelName"] in ["yeast8", "yeast-GEM_9"]:
        #     lExtracFluxs = [dConfigParams['glucose'], dConfigParams['ethanol'], dConfigParams['glycerol'], dConfigParams['acetate']]
        # elif "EG" in dConfigParams['condition'] and dConfigParams["modelName"] == "eciYali_sbmlValidated":
        #     lExtracFluxs = [dConfigParams['glucose'],  dConfigParams['eg']]
        # elif "EG" not in dConfigParams['condition'] and dConfigParams["modelName"] == "eciYali_sbmlValidated":
        #     lExtracFluxs = [dConfigParams['glucose']]
        #
        # for rxnExtFlux in lExtracFluxs:
        for rxnExtFlux in lMeasuredMets_production + lMeasuredMets_consumption:
            rxn = modelIn.reactions.get_by_id(dConfigParams[rxnExtFlux])
            gL.toLog(
                logFile, f"\t{rxn.id} - {rxn.name} - ({rxn.lower_bound}, {rxn.upper_bound})"
            )

        gL.toLog(logFile, "\n")

        ## Set O2 uptake flux: variazione dell'O2 per ottenere le diverse OD della popolazione di cellule oppure ragiono sulla singola cellula e lasciamo svincolato e attivo il consumo di O2
        gL.toLog(
            logFile,
            "* Constrain oxygen consumption flux")

        if dConfigParams['oxygen'] == 'free':
            modelIn, o2LB = expSL.setO2Flux(modelIn, dConfigParams=dConfigParams)
        #
        # if dConfigParams["modelName"] in ["yeast8", "yeast-GEM_9"]:
        #     ## open constraints for certain medium components
        #     # Ammonium	r_1654	5
        #     # Inositol	r_1947	2e-3
        #     # Riboflavin	r_2038	200e-6
        #     # Folic acid	r_1792	2e-6
        #     modelIn.reactions.get_by_id("r_1654").lower_bound = -1000
        #     modelIn.reactions.get_by_id("r_1947").lower_bound = -1000
        #     modelIn.reactions.get_by_id("r_2038").lower_bound = -1000
        #     modelIn.reactions.get_by_id("r_1792").lower_bound = -1000
        #     # # Copper	r_4594	40e-6
        #     modelIn.reactions.get_by_id("r_4594").lower_bound = -1000
        #     # # Manganese	r_4595	400e-6
        #     modelIn.reactions.get_by_id("r_4595").lower_bound = -1000
        #     # Zinc	r_4596	400e-6
        #     modelIn.reactions.get_by_id("r_4596").lower_bound = -1000
        #     # # Sodium	r_2049	0.1002
        #     modelIn.reactions.get_by_id("r_2049").lower_bound = -1000
        #     # Pi
        #     modelIn.reactions.get_by_id("r_2005").lower_bound = -1000


        outModelName = outModelFileName + "_" + tPoint1 + "h_" + tPoint2 + "h_" + tStamp
        cb.io.write_sbml_model(modelIn, os.path.join(MODELDIR, outModelName + ".xml"))

    return modelIn
