import numpy as np
import pandas as pd
import utils.manipulateModelLib as mmL
from cobra import Model
from cobra.io import write_sbml_model
import os
import utils.genericLib as gL
import sys
import re
import utils.keggLib as kL
import itertools as itt
import cobra as cb

MODELDIR = gL.dDirs["mods"]

def setBiomassReaction(modelIn, dConfigParams=None):
    """Sets the biomass reaction name accoroding to the model or the given rxnName and sets its boundaries.

    Keyword arguments:
         dConfigParams - is the configuration parameters dictionary where it should be stored a validated model name
         rxnName - [None | string] is the reaction name to be used as biomass reaction (if not None, replaces the default one )

    Return:
         [biomassRxn | None] The biomass reaction name. If None somethong went wrong.
    """
    biomassRxn = dConfigParams["biomassID"]

    modelIn.reactions.get_by_id(biomassRxn).lower_bound = 0
    modelIn.reactions.get_by_id(biomassRxn).upper_bound = 1000
    return biomassRxn


def switchOffAllSinkRxns(modelIn, lActiveUptake):
    """Block all the sink reactions in order to include in the model only uptake reactions in line with the experimental medium.

    Keyword arguments:
        modelIn - is the cobra model where the sink reactions are blocked
    """
    # lActiveUptake = dConfigParams["lActiveUptake"]

    for rxn in modelIn.reactions:
        # if len(rxn.reactants) * len(rxn.products) == 0.0 and rxn.id not in lActiveUptake and rxn.id.startswith("prot_") == False:
        if len(rxn.reactants) * len(rxn.products) == 0.0 and rxn.id not in lActiveUptake:
            if len(rxn.reactants) != 0 and rxn.lower_bound < 0:
                rxn.lower_bound = 0
            elif len(rxn.products) != 0 and rxn.upper_bound > 0:
                rxn.upper_bound = 0
    return modelIn

def addMediumData(modelIn, dfMedium, dNamesConversion, dConfigParams=None):
    """Add to the model all the medium data in terms of their experimental concentration.

    Keyword arguments:
         modelIn - is the cobra model where the medium data are added
         dfMedium - is the dataframe including all experimental concentrations of medium actors
         dConfigParams - is the configuration parameters dictionary where it should be stored a validated model name
    """
    #mediumSetting = dConfigParams["mediumBounds"]

    #for rxnMedium in dfMedium.itertuples():
    for metMedium in dfMedium.itertuples():
        print("ROW\n", metMedium, "\n")
        rxnId, modelMetId = mmL.getExchInModelForKeggId(modelIn, metMedium.KeggMetId, dNamesConversion)
        # if rxnId == "":
        print("met: ", metMedium.Met)
        print("kegg met: ", metMedium.KeggMetId)
        print("modelMetId: ", modelMetId)
        print("rxnId: ", rxnId, "\n")
        if rxnId is None:
            mmL.addExchRxnForMet(modelIn, modelIn.metabolites.get_by_id(modelMetId), newLB=-1000, newUB=0)
            rxnId = 'EX_' + modelMetId
            print("None rxnId: ", rxnId, "\n")

        rxn = modelIn.reactions.get_by_id(rxnId)
        # rxn = modelIn.reactions.get_by_id(rxnMedium.RxnID)
        #if mediumSetting == "wet":
        rxn.lower_bound = -1 * (metMedium.Concentration)
        #elif mediumSetting == "free":
        #    rxn.lower_bound = -1000

    return modelIn



def setExtracFluxes(dfConsProd, lMeasuredMets, modelIn, dNamesConversion,dfConversionMetName2MetKegg, baseModelName):
    """Set lower and upper boundaries for extracellular consumptions or production variations between two time points according to experimental values.

    Keyword arguments:
        rowT1 - is the row for the first time point of the dataframe including experimental extracellular consumption rates
        rowT2 - is the row for the second time point of the dataframe including experimental extracellular consumption rates
        modelIn - is the cobra model where the sink reactions are blocked
        dConfigParams - is the configuration parameters dictionary
    """

    lTimePoints = dfConsProd.index.tolist()
    print("time points: ", lTimePoints)

    # sort list (essential for grouping)
    lMeasuredMets.sort()
    # using lambda + itertools.groupby() + split() to group similar substrings
    lMeasuredMets_grouped = [list(i) for j, i in itt.groupby(lMeasuredMets,
                      lambda a: a.split('_')[0])]
    print ("The grouped list is : " + str(lMeasuredMets_grouped))

    if len(lTimePoints) == 1:
        dValues = {}
        for metMeasured in lMeasuredMets_grouped:
            lvalues = dfConsProd.iloc[0][metMeasured].tolist()
            dValues[metMeasured[0].split("_")[0]] = [min(lvalues), max(lvalues)]
    else:
        # use itertools.tee to create two iterators from the list
        a, b = itt.tee(lTimePoints)
        # advance the iterator by one element
        next(b, None)
        # use zip to pair the elements from the two iterators
        lTimePointsCouples = list(zip(a, b))
        print("lTimePointsCouples: ", lTimePointsCouples)

        for pair in lTimePointsCouples:
            with modelIn:
                # create dictionary where for each measured metabolite there is the corresponding consumption or production rate in a given time interval
                # pair[1] + "_" + pair[0] ## key
                print(pair)
                deltaT = int(pair[1]) - int(pair[0])
                print("deltaT: ", deltaT)

                dfConsProd_currentTimePoints = dfConsProd[dfConsProd.index.isin(pair)]
                print("\nDF\n", dfConsProd_currentTimePoints)
                # dfConsProd_currentTimePoints.set_index('Time', inplace = True)
                print("\nDF_T\n", dfConsProd_currentTimePoints.T)

                dfConsProd_currentTimePoints_transpose = dfConsProd_currentTimePoints.T
                dfConsProd_currentTimePoints_transpose["variation"] =(dfConsProd_currentTimePoints_transpose[pair[1]] - dfConsProd_currentTimePoints_transpose[pair[0]])/deltaT

                dfConsProd_currentT_wVariation = dfConsProd_currentTimePoints_transpose[["variation"]].T
                print("dfConsProd_currentT_wVariation: \n", dfConsProd_currentT_wVariation)

                for metMeasured in lMeasuredMets_grouped:
                    print("metMeasured: ", metMeasured)
                    lvalues = dfConsProd_currentT_wVariation.loc["variation"][metMeasured].tolist()
                    print(lvalues)
                    print(min(lvalues), max(lvalues))

                    ## convert name to kegg id
                    metSplitted = metMeasured[0].split("_")
                    del metSplitted[-1]
                    metLower = " ".join(metSplitted).lower()
                    print("metLower: ", metLower)
                    correspondingKeggId, dfConversionMetName2MetKegg= kL.getKeggMetId(metLower, dfConversionMetName2MetKegg)

                    # search for metabolite in model corresponding to the retrieved kegg id
                    rxnId, modelMetId = mmL.getExchInModelForKeggId(modelIn, correspondingKeggId, dNamesConversion)
                    print("met: ", correspondingKeggId)
                    print("modelMetId: ", modelMetId)
                    print("rxnId: ", rxnId, "\n")
                    if rxnId is None:
                        mmL.addExchRxnForMet(modelIn, modelIn.metabolites.get_by_id(modelMetId), newLB=-1000, newUB=0)
                        rxnId = 'EX_' + modelMetId
                        print("None rxnId: ", rxnId, "\n")

                    rxn = modelIn.reactions.get_by_id(rxnId)

                    # use the dValues dictionary to set constraints for the exchange reactions corresponding to each measured metabolite

                    if max(lvalues) < 0:
                        if len(rxn.reactants) != 0 and len(rxn.products) == 0: # A = 0
                            try:
                                rxn.upper_bound = 0
                                rxn.lower_bound = 1 * min(lvalues)
                            except:
                                rxn.lower_bound = 1 * min(lvalues)
                                rxn.upper_bound = 0
                        elif len(rxn.reactants) == 0 and len(rxn.products) != 0: # 0 = A
                            try:
                                rxn.upper_bound = -1 * min(lvalues)
                                rxn.lower_bound = 0
                            except:
                                rxn.lower_bound = 0
                                rxn.upper_bound = -1 * min(lvalues)
                    else:
                        if len(rxn.reactants) != 0 and len(rxn.products) == 0: # A = 0
                            try:
                                rxn.upper_bound = 1 * max(lvalues)
                                rxn.lower_bound = 1 * min(lvalues)
                            except:
                                rxn.lower_bound = 1 * min(lvalues)
                                rxn.upper_bound = 1 * max(lvalues)
                        elif len(rxn.reactants) == 0 and len(rxn.products) != 0: # 0 = A
                            try:
                                rxn.upper_bound = -1 * min(lvalues)
                                rxn.lower_bound = -1 * max(lvalues)
                            except:
                                rxn.lower_bound = -1 * max(lvalues)
                                rxn.upper_bound = -1 * min(lvalues)

                    print("LB: ", rxn.lower_bound)
                    print("UB: ", rxn.upper_bound, "\n")

                # save model
                cb.io.write_sbml_model(modelIn, os.path.join(MODELDIR, baseModelName + "_" + "_".join(pair) + ".xml"))

def addMediumData2Model(modelwMediumData, dfMedium, dfConditionMets, lMetRemainActive, dfName2KeggId, dKeggId2ModelMet, dConfigParams=None):
    ldfMedium_wInfo = []
    for row in dfMedium.itertuples():
        correspondingKeggId, dfName2KeggId = kL.getKeggMetId(row.Met.lower(), dfName2KeggId)
        ldfMedium_wInfo.append([row.Met, correspondingKeggId, row.Concentration])
    dfMedium_wKeggId = pd.DataFrame(ldfMedium_wInfo, columns = ["Met", "KeggMetId", "Concentration"])

    ## aggiungo al file medium delle righe per definire cosa entra dalla variabile "medium"
    for metCondition in dfConditionMets.itertuples():
        correspondingKeggId, dfName2KeggId = kL.getKeggMetId(metCondition.Met.lower(), dfName2KeggId)
        dfMedium_wKeggId.loc[len(dfMedium_wKeggId)] = [metCondition.Met, correspondingKeggId, metCondition.Lb]

    lActiveUptake = []
    for activeMet in lMetRemainActive:
        rxnId, modelMetId = mmL.getExchInModelForKeggId(modelwMediumData, activeMet, dKeggId2ModelMet)
        lActiveUptake.append(rxnId)

    #lActiveUptake += dConfigParams["lActiveUptake"] # aggiungerei per yarrowia reazioni che nel modello iniziano per "prot_"
    modelwMediumData = switchOffAllSinkRxns(modelwMediumData, lActiveUptake)

    # add medium concentrations data
    modelwMediumData = addMediumData(modelwMediumData, dfMedium_wKeggId,dKeggId2ModelMet, dConfigParams = dConfigParams)
    return dfName2KeggId, modelwMediumData

########################################
# OLD PART
########################################

def createHeterologousModel(expKind):
    if expKind == "endogenous":
        model = addPathwayEndogeno()
    elif expKind == "aox":
        model = addPathwayEndogeno_alternative()
    elif expKind == "bacterial":
        model = addPathwayBatterico()
    elif expKind == "ncl":
        model = addPathwayNCL()
    elif expKind == "taco":
        model = addPathwayTaCo()
    elif expKind == "marine":
        model = addPathwayMarine()
    elif expKind == "wt":
        model = EGpathwayYarrowia()
    elif expKind == "wt_sceLike":
        model = EGpathwayYarrowia_sceLike()

    write_sbml_model(model, os.path.join(MODELDIR, expKind + ".xml"))
    return model

def addPathwayEndogeno(expName = 'Endogenous pathway'):
    model = Model(expName)
    lMets = [("s_3997_m", "glycolate [mitochondrion]", "m"),
    ("C01380_c", "ethylene glycol [cytoplasm]", "c"),
    ("s_1198", "NAD [cytoplasm]", "c"),
    ("s_0775", "glycolaldehyde [cytoplasm]", "c"),
    ("s_1203", "NADH [cytoplasm]", "c"),
    ("s_0794", "H+ [cytoplasm]", "c"),
    ("s_0803", "H2O [cytoplasm]", "c"),
    ("s_3997", "glycolate [cytoplasm]", "c"),
    ("s_3997_e", "glycolate [extracellular]", "e"),
    ("s_1207", "NADP(+) [cytoplasm]", "c"),
    ("s_1212", "NADPH [cytoplasm]", "c"),
    ("s_0777", "glycolaldehyde [mitochondrion]", "m"),
    ("s_1200", "NAD [mitochondrion]", "m"),
    ("s_0807", "H2O [mitochondrion]", "m"),
    ("s_1205", "NADH [mitochondrion]", "m"),
    ("s_0799", "H+ [mitochondrion]", "m"),
    ("s_1210", "NADP(+) [mitochondrion]", "m"),
    ("s_1214", "NADPH [mitochondrion]", "m")]

    dfMets = pd.DataFrame(lMets, columns = ['Id', 'Name', 'Compartment'])
    for m in dfMets.itertuples():
        mmL.addMets(model, m)

    lRxns = [
        (
            "R01781",
            "1,2-Ethanediol:NAD+ oxidoreductase",
            "C01380_c + s_1198 = s_0775 + s_1203 + s_0794",
            0,
            1000,
            "YLL056C or P0A9S1 or Q5FU50",
        ),
        (
            "R01333_NAD_c",
            "glycolaldehyde:NAD+ oxidoreductase",
            "s_0775 + s_1198 + s_0803 = s_3997 + s_1203 + 2 s_0794",
            0,
            1000,
            "YMR170C or YMR169C or B1N7J6",
        ),
        (
            "R01333_NADP_c",
            "glycolaldehyde:NADP+ oxidoreductase",
            "s_0775 + s_1207 + s_0803 = s_3997 + s_1212 + 2 s_0794",
            0,
            1000,
            "YMR170C or YMR169C",
        ),
        (
            "R01333_NAD_m",
            "glycolaldehyde:NAD+ oxidoreductase",
            "s_0777 + s_1200 + s_0807 = s_3997_m + s_1205 + 2 s_0799",
            0,
            1000,
            "YOR374W or YER073W",
        ),
        (
            "R01333_NADP_m",
            "glycolaldehyde:NADP+ oxidoreductase",
            "s_0777 + s_1210 + s_0807 = s_3997_m + s_1214 + 2 s_0799",
            0,
            1000,
            "YOR374W or YER073W",
        ),
        (
            "Transport_glycolate_c_m",
            "Glycolate transport cytoplasm_mitochondrion",
            "s_3997 = s_3997_m",
            -1000,
            1000,
            "",
        ),
        (
            "Transport_glycolate_c_e",
            "Glycolate transport cytoplasm_extracellular",
            "s_3997 = s_3997_e",
            0,
            1000,
            "",
        ),
    ]
    dfRxns = pd.DataFrame(lRxns, columns = ['Id', 'Name', 'Equation', 'Lb', 'Ub', 'GPR'])
    for r in dfRxns.itertuples():
        mmL.addRxns(model, r)

    metIsolatedExchgRxn(model, model.metabolites.get_by_id("s_3997_e"), newLB=0.0, newUB=1000.0)

    return model



def addPathwayEndogeno_alternative(expName = 'Endogenous alternative pathway'):
    model = Model(expName)
    lMets = [
    ("C01380_c", "ethylene glycol [cytoplasm]", "c"),
    ("s_1275", "oxygen", "c"),
    ("s_0775", "glycolaldehyde [cytoplasm]", "c"),
    ("s_0837", "hydrogen peroxide", "c"),
    ("C14448_c", "glyoxal", "c"),
    ("s_3997", "glycolate [cytoplasm]", "c"),
    ("s_3997_e", "glycolate [extracellular]", "e"),
    ("s_1207", "NADP(+) [cytoplasm]", "c"),
    ("s_1212", "NADPH [cytoplasm]", "c"),
    ("s_0779","glyoxylate [cytoplasm]", "c")]

    dfMets = pd.DataFrame(lMets, columns = ['Id', 'Name', 'Compartment'])
    for m in dfMets.itertuples():
        mmL.addMets(model, m)

    lRxns = [
        (
            "AO1",
            "AO1",
            "C01380_c + s_1275 = s_0775 + s_0837",
            0,
            1000,
            "",
        ),
        (
            "AO2",
            "AO2",
            "s_0775 + s_1275 = C14448_c + s_0837",
            0,
            1000,
            "",
        ),
        (
            "AO3",
            "AO3",
            "C14448_c + s_1207 = s_3997 + s_1212",
            0,
            1000,
            "C",
        ),
        (
            "Transport_glycolate_c_e",
            "Glycolate transport cytoplasm_extracellular",
            "s_3997 = s_3997_e",
            0,
            1000,
            "",
        ),
    ]
    dfRxns = pd.DataFrame(lRxns, columns = ['Id', 'Name', 'Equation', 'Lb', 'Ub', 'GPR'])
    for r in dfRxns.itertuples():
        mmL.addRxns(model, r)

    metIsolatedExchgRxn(model, model.metabolites.get_by_id("s_3997_e"), newLB=0.0, newUB=1000.0)
    return model



def addPathwayBatterico(expName = 'Bacterial pathway'):
    model = Model(expName)
    lMets = [
        ("C01146", "Tartronate semialdehyde [cytoplasm]", "c"),
        ("C00168", "3-Hydroxypyruvate [cytoplasm]", "c"),
        ("s_0794", "H+ [cytoplasm]", "c"),
        ("s_1203", "NADH [cytoplasm]", "c"),
        ("s_1198", "NAD [cytoplasm]", "c"),
        ("s_0779","glyoxylate [cytoplasm]", "c"),
        ("s_1207", "NADP(+) [cytoplasm]", "c"),
        ("s_1212", "NADPH [cytoplasm]", "c"),
        ("s_0394", "ADP [cytoplasm]", "c"),
        ("s_0188", "2-phospho-D-glyceric acid [cytoplasm]", "c"),
        ("s_0456", "carbon dioxide [cytoplasm]", "c"),
        ("s_3958", "D-Glycerate [cytoplasm]", "c"),
        ("s_0434", "ATP [cytoplasm]", "c")
    ]

    dfMets = pd.DataFrame(lMets, columns = ['Id', 'Name', 'Compartment'])
    for m in dfMets.itertuples():
        mmL.addMets(model, m)

    lRxns = [
        (
            "B4",
            "glc-catalyzed reaction",
            "2 s_0779 + s_0794 = C01146 + s_0456",
            0,
            1000,
            "P0AEP7",
        ),
        (
            "B5a_NAD",
            "glxR-catalyzed reaction NADH dependent",
            "C01146 + s_1203 + s_0794 = s_3958 + s_1198",
            0,
            1000,
            "P77161",
        ),
        (
            "B5a_NADP",
            "glxR-catalyzed reaction NADPH dependent",
            "C01146 + s_1212 + s_0794 = s_3958 + s_1207",
            0,
            1000,
            "P77161",
        ),
        (
            "B5b1",
            "hyi-catalyzed reaction",
            "C01146 = C00168",
            0,
            1000,
            "P30147",
        ),
        (
            "B5b2_NAD",
            "glxR-catalyzed reaction NADH dependent",
            "C00168 + s_0794 + s_1203 = s_3958 + s_1198",
            0,
            1000,
            "P77161",
        ),
        (
            "B5b2_NADP",
            "glxR-catalyzed reaction NADPH dependent",
            "C00168 + s_0794 + s_1212 = s_3958 + s_1207",
            0,
            1000,
            "P77161",
        ),
        (
            "B6",
            "garK or glxK-catalyzed reaction",
            "s_3958 + s_0434 = s_0794 + s_0394 + s_0188",
            0,
            1000,
            "P23524 or P77364",
        ),
    ]

    dfRxns = pd.DataFrame(lRxns, columns = ['Id', 'Name', 'Equation', 'Lb', 'Ub', 'GPR'])
    for r in dfRxns.itertuples():
        mmL.addRxns(model, r)
    return model


def addPathwayNCL(expName = 'NCL pathway'):
    model = Model(expName)
    lMets = [("C00227", "Acetyl phosphate [cytoplasm]", "c"),
            ("s_0794", "H+ [cytoplasm]", "c"),
            ("s_1203", "NADH [cytoplasm]", "c"),
            ("s_1198", "NAD [cytoplasm]", "c"),
            ("s_0775", "glycolaldehyde [cytoplasm]", "c"),
            ("s_1475", "TDP [cytoplasm]", "c"),
            ("s_1489", "thiamine [cytoplasm]", "c"),
            ("s_0529", "coenzyme A [cytoplasm]", "c"),
            ("s_0373", "acetyl-CoA [cytoplasm]", "c"),
            ("C01380_c", "ethylene glycol [cytoplasm]", "c")]

    dfMets = pd.DataFrame(lMets, columns = ['Id', 'Name', 'Compartment'])
    for m in dfMets.itertuples():
        mmL.addMets(model, m)

    lRxns = [
        (
            "NCL2",
            "NCL2",
            "s_0775 + s_1475 = C00227 + s_1489",
            0,
            1000,
            "P24224",
        ),
        (
            "NCL3",
            "NCL3",
            "C00227 + s_1489 + s_0529 = s_0373 + s_1475",
            0,
            1000,
            "P77218",
        ),
    ]

    dfRxns = pd.DataFrame(lRxns, columns = ['Id', 'Name', 'Equation', 'Lb', 'Ub', 'GPR'])
    for r in dfRxns.itertuples():
        mmL.addRxns(model, r)

    return model


def addPathwayTaCo(expName = 'TaCo pathway'):
    model = Model(expName)
    lMets = [
        ("CPD22847", "Glycolyl-CoA [cytoplasm]", "c"),
        ("tartronylCoA", "tartronyl-CoA [cytoplasm]", "c"),
        ("s_0775", "glycolaldehyde [cytoplasm]", "c"),
        ("s_0794", "H+ [cytoplasm]", "c"),
        ("s_1203", "NADH [cytoplasm]", "c"),
        ("s_1198", "NAD [cytoplasm]", "c"),
        ("s_1207", "NADP(+) [cytoplasm]", "c"),
        ("s_1212", "NADPH [cytoplasm]", "c"),
        ("s_0451", "biotin [cytoplasm]", "c"),
        ("s_0394", "ADP [cytoplasm]", "c"),
        ("s_0434", "ATP [cytoplasm]", "c"),
        ("s_3958", "D-Glycerate [cytoplasm]", "c"),
        ("s_0529", "coenzyme A [cytoplasm]", "c"),
        ("s_0445", "bicarbonate [cytoplasm]", "c"),
        ("s_1322", "phosphate [cytoplasm]", "c"),
        ("int1", "intermediate T3 reaction [cytoplasm]", "c")]

    dfMets = pd.DataFrame(lMets, columns = ['Id', 'Name', 'Compartment'])
    for m in dfMets.itertuples():
        mmL.addMets(model, m)

    lRxns = [
        (
            "T2",  # in kegg per il gene annotato per questa rxn c'e' solo la rxn generica R09097
            "T2",
            "s_0775 + s_0529 + s_1198 = CPD22847 + s_1203 + s_0794",
            0,
            1000,
            "Q21A49",
        ),
        (
            "T3_1",
            "T3 part1",
            "CPD22847 + s_0445 + s_0434 + s_0451 = int1",
            0,
            1000,
            "GCCA and GCCB",
        ),

        (
            "T3_2",
            "T3 part2",
            "int1 = tartronylCoA + s_0394 + s_1322 + s_0451",
            0,
            1000,
            "GCCA and GCCB",
        ),
        (
            "T4",
            "T4",
            "tartronylCoA + s_1212 + s_0794 = s_3958 + s_1207 + s_0529",
            0,
            1000,
            "TCR",
        ),
    ]

    dfRxns = pd.DataFrame(lRxns, columns = ['Id', 'Name', 'Equation', 'Lb', 'Ub', 'GPR'])
    for r in dfRxns.itertuples():
        mmL.addRxns(model, r)

    lMets2 = [("C01380_c", "ethylene glycol [cytoplasm]", "c"),
    ("s_0188", "2-phospho-D-glyceric acid [cytoplasm]", "c"),]
    dfMets2 = pd.DataFrame(lMets2, columns = ['Id', 'Name', 'Compartment'])
    for m2 in dfMets2.itertuples():
        mmL.addMets(model, m2)

    lRxns2 = [
        (
            "B6",
            "garK or glxK-catalyzed reaction",
            "s_3958 + s_0434 = s_0794 + s_0394 + s_0188",
            0,
            1000,
            "P23524 or P77364",
        ),
    ]

    dfRxns2 = pd.DataFrame(lRxns2, columns = ['Id', 'Name', 'Equation', 'Lb', 'Ub', 'GPR'])
    for r2 in dfRxns2.itertuples():
        mmL.addRxns(model, r2)

    return model


def addPathwayMarine(expName = 'Marine pathway'):
    model = Model(expName)
    lMets = [
        ("C03961", "beta-hydroxy-Aspartate [cytoplasm]", "c"),
        ("C05840", "Iminosuccinate [cytoplasm]", "c"),
        ("s_3997_m", "glycolate [mitochondrion]", "m"),
        ("s_0794", "H+ [cytoplasm]", "c"),
        ("s_1203", "NADH [cytoplasm]", "c"),
        ("s_1198", "NAD [cytoplasm]", "c"),
        ("s_1207", "NADP(+) [cytoplasm]", "c"),
        ("s_1212", "NADPH [cytoplasm]", "c"),
        ("s_0803", "H2O [cytoplasm]", "c"),
        ("s_0779", "glyoxylate [cytoplasm]", "c"),
        ("s_0973", "L-aspartate [cytoplasm]", "c"),
        ("s_1393", "pyridoxal 5&apos;-phosphate [cytoplasm]", "c"),
        ("s_1003", "L-glycine [cytoplasm]", "c"),
        ("s_1271", "oxaloacetate [cytoplasm]", "c"),
        ("s_0419", "ammonium [cytoplasm]", "c"),
        ("C01380_c", "ethylene glycol [cytoplasm]", "c"),
        ("s_0775", "glycolaldehyde [cytoplasm]", "c"),
        ("s_3997", "glycolate [cytoplasm]", "c"),
        ("s_0777", "glycolaldehyde [mitochondrion]", "m"),
        ("s_1200", "NAD [mitochondrion]", "m"),
        ("s_0807", "H2O [mitochondrion]", "m"),
        ("s_1205", "NADH [mitochondrion]", "m"),
        ("s_0799", "H+ [mitochondrion]", "m"),
        ("s_1210", "NADP(+) [mitochondrion]", "m"),
        ("s_1214", "NADPH [mitochondrion]", "m"),
        ("int2", "intermediate reaction M1 [cytoplasm]", "c")]

    dfMets = pd.DataFrame(lMets, columns = ['Id', 'Name', 'Compartment'])
    for m in dfMets.itertuples():
        mmL.addMets(model, m)

    lRxns = [
        (
            "M1_1",
            "M1 part1",
            "s_0779 + s_0973 + s_1393 = int2",
            0,
            1000,
            "A1B8Z3",
        ),

        (
            "M1_2",
            "M1 part2",
            "int2 = s_1003 + s_1271 + s_1393",
            0,
            1000,
            "A1B8Z3",
        ),
        (
            "M2",
            "M2",
            "s_1003 + s_0779 = C03961",
            0,
            1000,
            "A1B8Z1",
        ),
        (
            "M3",
            "M3",
            "C03961 = C05840 + s_0803",
            0,
            1000,
            "A1B8Z2",
        ),
        (
            "M4A",
            "M4A",
            "C05840 + s_1203 + s_0794 = s_0973 + s_1198",
            0,
            1000,
            "A1B8Z0",
        ),
        (
            "M4B",
            "M4B",
            "C05840 + s_0803 = s_1271 + s_0419",
            0,
            1000,
            "",
        ),
    ]

    dfRxns = pd.DataFrame(lRxns, columns = ['Id', 'Name', 'Equation', 'Lb', 'Ub', 'GPR'])
    for r in dfRxns.itertuples():
        mmL.addRxns(model, r)

    return model


def EGpathwayYarrowia_sceLike(expName = 'EG pathway like Sce'):
    model = Model(expName)
    lMets = [
        ("C01380[c]", "ethylene glycol [cytoplasm]", "c"),
        ("s_0359[c]", "acetaldehyde [cytoplasm]", "c"),
        ("s_0803[c]", "H2O [cytoplasm]", "c"),
        ("s_1198[c]", "NAD [cytoplasm]", "c"),
        ("s_0775[c]", "glycolaldehyde [cytoplasm]", "c"),
        ("s_1203[c]", "NADH [cytoplasm]", "c"),
        ("s_0794[c]", "H+ [cytoplasm]", "c"),
        ("C00160[c]", "Glycolate [cytoplasm]", "c"),
        ("ETR-Quinones[c]", "electron-transfer quinone [cytoplasm]", "c"),
        ("ETR-Quinols[c]", "electron-transfer quinol [cytoplasm]", "c"),
        ("s_0779[c]", "glyoxylate [cytoplasm]", "c"),
        ("s_0456[c]", "carbon dioxide [cytoplasm]", "c"),
        ("C01146[c]", "2-Hydroxy-3-oxopropanoate [cytoplasm]", "c"),
        ("C00258[c]", "D-Glycerate [cytoplasm]", "c"),
        ("C03961[c]", "erythro-3-Hydroxy-Ls-aspartate [cytoplasm]", "c"),
        ("s_1464[c]", "succinyl-CoA [cytoplasm]", "c"),
        ("s_1212[c]", "NADPH [cytoplasm]", "c"),
        ("s_1207[c]", "NADP(+) [cytoplasm]", "c"),
        ("s_0434[c]", "ATP [cytoplasm]", "c"),
        ("s_0188[c]", "2-phospho-D-glyceric acid [cytoplasm]", "c"),
        ("s_0394[c]", "ADP [cytoplasm]", "c"),
        ("s_0373[c]", "acetyl-CoA [cytoplasm]", "c"),
        ("s_0066[c]", "(S)-malate [cytoplasm]", "c"),
        ("s_0529[c]", "coenzyme A [cytoplasm]", "c"),
        ("s_1003[c]", "L-glycine [cytoplasm]", "c"),
        ("s_1271[c]", "oxaloacetate [cytoplasm]", "c"),
        ("s_0419[c]", "ammonium [cytoplasm]", "c"),
    ]

    dfMets = pd.DataFrame(lMets, columns = ['Id', 'Name', 'Compartment'])
    for m in dfMets.itertuples():
        mmL.addMets(model, m)

    lRxns = [
        (
            "R01781",
            "R01781",
            "C01380[c] + s_1198[c] = s_0775[c] + s_1203[c] + s_0794[c]",
            0,
            1000,
            "",
        ),
        (
            "R01333",
            "R01333",
            "s_0775[c] + s_1198[c] + s_0803[c] = C00160[c] + s_1203[c] + 2 s_0794[c]",
            0,
            1000,
            "",
        ),
            (
            "RXN0-7229",
            "RXN0-7229",
            "C00160[c] + s_1198[c] = s_0779[c] + s_0794[c] + s_1203[c]",
            -1000,
            0,
            "",
        ),
    ]

    dfRxns = pd.DataFrame(lRxns, columns = ['Id', 'Name', 'Equation', 'Lb', 'Ub', 'GPR'])
    for r in dfRxns.itertuples():
        mmL.addRxns(model, r)

    metIsolatedExchgRxn(model, model.metabolites.get_by_id("C01380[c]"), newLB=-1000.0, newUB=0.0)
    return model

def EGpathwayYarrowia(expName = 'EG pathway from literature'):
    model = Model(expName)
    lMets = [
        ("C01380[c]", "ethylene glycol [cytoplasm]", "c"),
        ("s_0359[c]", "acetaldehyde [cytoplasm]", "c"),
        ("s_0803[c]", "H2O [cytoplasm]", "c"),
        ("s_1198[c]", "NAD [cytoplasm]", "c"),
        ("s_0775[c]", "glycolaldehyde [cytoplasm]", "c"),
        ("s_1203[c]", "NADH [cytoplasm]", "c"),
        ("s_0794[c]", "H+ [cytoplasm]", "c"),
        ("C00160[c]", "Glycolate [cytoplasm]", "c"),
        ("ETR-Quinones[c]", "electron-transfer quinone [cytoplasm]", "c"),
        ("ETR-Quinols[c]", "electron-transfer quinol [cytoplasm]", "c"),
        ("s_0779[c]", "glyoxylate [cytoplasm]", "c"),
        ("s_0456[c]", "carbon dioxide [cytoplasm]", "c"),
        ("C01146[c]", "2-Hydroxy-3-oxopropanoate [cytoplasm]", "c"),
        ("C00258[c]", "D-Glycerate [cytoplasm]", "c"),
        ("C03961[c]", "erythro-3-Hydroxy-Ls-aspartate [cytoplasm]", "c"),
        ("s_1464[c]", "succinyl-CoA [cytoplasm]", "c"),
        ("s_1212[c]", "NADPH [cytoplasm]", "c"),
        ("s_1207[c]", "NADP(+) [cytoplasm]", "c"),
        ("s_0434[c]", "ATP [cytoplasm]", "c"),
        ("s_0188[c]", "2-phospho-D-glyceric acid [cytoplasm]", "c"),
        ("s_0394[c]", "ADP [cytoplasm]", "c"),
        ("s_0373[c]", "acetyl-CoA [cytoplasm]", "c"),
        ("s_0066[c]", "(S)-malate [cytoplasm]", "c"),
        ("s_0529[c]", "coenzyme A [cytoplasm]", "c"),
        ("s_1003[c]", "L-glycine [cytoplasm]", "c"),
        ("s_1271[c]", "oxaloacetate [cytoplasm]", "c"),
        ("s_0419[c]", "ammonium [cytoplasm]", "c"),
        # ("C06337[c]", "Terephthalate [cytoplasm]", "c"),
        # (
        #     "C06318[c]",
        #     "(3S,4R)-3,4-Dihydroxycyclohexa-1,5-diene-1,4-dicarboxylate [cytoplasm]",
        #     "c",
        # ),
        # ("C00230[c]", "3,4-Dihydroxybenzoate [cytoplasm]", "c"),
        # ("C01163[c]", "3-Carboxy-cis,cis-muconate [cytoplasm]", "c"),
        # ("C01278[c]", "2-Carboxy-2,5-dihydro-5-oxofuran-2-acetate [cytoplasm]", "c"),
        # ("C03586[c]", "2-Oxo-2,3-dihydrofuran-5-acetate [cytoplasm]", "c"),
        # ("C00846[c]", "3-Oxoadipate [cytoplasm]", "c"),
        # ("C02232[c]", "3-Oxoadipyl-CoA [cytoplasm]", "c"),
        # ("C04484[c]", "4-Carboxy-2-hydroxymuconate semialdehyde [cytoplasm]", "c"),
        # ("C05375[c]", "2-Hydroxy-2-hydropyrone-4,6-dicarboxylate [cytoplasm]", "c"),
        # ("C03671[c]", "2-Pyrone-4,6-dicarboxylate [cytoplasm]", "c"),
        # ("C04434[c]", "(1E)-4-Oxobut-1-ene-1,2,4-tricarboxylate [cytoplasm]", "c"),
        # ("C04115[c]", "4-Carboxy-4-hydroxy-2-oxoadipate [cytoplasm]", "c"),
        # (
        #     "MNXM9831[c]",
        #     "(2Z,4Z)-2-hydroxy-5-carboxymuconate-6-semialdehyde [cytoplasm]",
        #     "c",
        # ),
        # ("C00682[c]", "(2Z,4E)-2-hydroxy-6-oxohexa-2,4-dienoate [cytoplasm]", "c"),
        # ("C02501[c]", "2-Hydroxymuconate [cytoplasm]", "c"),
        # ("C03453[c]", "gamma-Oxalocrotonate [cytoplasm]", "c"),
        # ("C00596[c]", "2-Hydroxy-2,4-pentadienoate [cytoplasm]", "c"),
        # ("C03589[c]", "4-Hydroxy-2-oxopentanoate [cytoplasm]", "c"),
    ]

    dfMets = pd.DataFrame(lMets, columns = ['Id', 'Name', 'Compartment'])
    for m in dfMets.itertuples():
        mmL.addMets(model, m)

    lRxns = [
            # (
            #     "R05148",
            #     "R05148",
            #     "C06337[c] + s_1275[c] + s_1203[c] + s_0794[c] = C06318[c] + s_1198[c]",
            #     0,
            #     1000,
            #     "",
            # ),
            # (
            #     "R01633",
            #     "R01633",
            #     "C06318[c] + s_1198[c] = C00230[c] + s_0456[c] + s_1203[c] + s_0794[c]",
            #     0,
            #     1000,
            #     "",
            # ),
            # (
            #     "R01631",
            #     "R01631",
            #     "C00230[c]  + s_1275[c] = C01163[c]",
            #     0,
            #     1000,
            #     "",
            # ),
            # (
            #     "R03307",
            #     "R03307",
            #     "C01163[c] = C01278[c]",
            #     0,
            #     1000,
            #     "",
            # ),
            # (
            #     "R03470",
            #     "R03470",
            #     "C01278[c] = C03586[c] + s_0456[c]",
            #     -1000,
            #     1000,
            #     "",
            # ),
            # (
            #     "R02991",
            #     "R02991",
            #     "C03586[c] + s_0803[c] = C00846[c]",
            #     0,
            #     1000,
            #     "",
            # ),
            # (
            #     "R02990",
            #     "R02990",
            #     "s_1464[c] + C00846[c] = s_1458[c] + C02232[c]"),
            #     0,
            #     1000,
            #     "",
            # ),
            # (
            #     "R00829",
            #     "R00829",
            #     "s_0529[c] + C02232[c] = s_1464[c] + s_0373[c]",
            #     0,
            #     1000,
            #     "",
            # ),
            # (
            #     "R01632",
            #     "R01632",
            #     "C00230[c] + s_1275[c] = C04484[c]",
            #     0,
            #     1000,
            #     "",
            # ),
            # (
            #     "R04489",
            #     "R04489",
            #     "C04484[c] = C05375[c]",
            #     0,
            #     1000,
            #     "",
            # ),
            # (
            #     "R04279",
            #     "R04279",
            #     "C05375[c] + s_1207[c] = C03671[c] + s_1212[c] + s_0794[c]",
            #     0,
            #     1000,
            #     "",
            # ),
            # (
            #     "R04277",
            #     "R04277",
            #     "C03671[c] + s_0803[c] = C04434[c]",
            #     0,
            #     1000,
            #     "",
            # ),
            # (
            #     "R04478",
            #     "R04478",
            #     "C04434[c] + s_0803[c] = C04115[c]",
            #     0,
            #     1000,
            #     "",
            # ),
            # (
            #     "R00350",
            #     "R00350",
            #     "C04115[c] = s_1271[c] + s_1399[c]",
            #     0,
            #     1000,
            #     "",
            # ),
            # (
            #     "R00217",
            #     "R00217",
            #     "s_1271[c] = s_1399[c] + s_0456[c]",
            #     0,
            #     1000,
            #     "",
            # ),
            # (
            #     "RXN-10882",
            #     "RXN-10882",
            #     "C00230[c] + s_1275[c] = MNXM9831[c]",
            #     0,
            #     1000,
            #     "",
            # ),
            # (
            #     "RXN-10883",
            #     "RXN-10883",
            #     "MNXM9831[c] + s_0794[c] = C00682[c] + s_0456[c]",
            #     0,
            #     1000,
            #     "",
            # ),
            # (
            #     "R02762",
            #     "R02762",
            #     "C00682[c] + s_1198[c] + s_0803[c] = C02501[c] + s_1203[c] + s_0794[c]",
            #     0,
            #     1000,
            #     "",
            # ),
            # (
            #     "R03966",
            #     "R03966",
            #     "C02501[c] = C03453[c]",
            #     -1000,
            #     1000,
            #     "",
            # ),
            # (
            #     "R02602",
            #     "R02602",
            #     "C03453[c] = C00596[c] + s_0456[c]",
            #     0,
            #     1000,
            #     "",
            # ),
            # (
            #     "R02601",
            #     "R02601",
            #     "C00596[c] + s_0803[c] = C03589[c]",
            #     0,
            #     1000,
            #     "",
            # ),
            # (
            #     "R00750",
            #     "R00750",
            #     "C03589[c] = s_0359[c] + s_1399[c]",
            #     0,
            #     1000,
            #     "",
            # ),
            # (
            #     "R00228",
            #     "R00228",
            #     "s_0359[c] + s_0529[c] + s_1198[c] = s_0373[c] + s_1203[c] + s_0794[c]",
            #     0,
            #     1000,
            #     "",
            # ),
            (
            "RXN-21079",
            "RXN-21079",
            "C01380[c] = s_0359[c] + s_0803[c]",
            0,
            1000,
            "",
        ),
        (
            "R01781",
            "R01781",
            "C01380[c] + s_1198[c] = s_0775[c] + s_1203[c] + s_0794[c]",
            0,
            1000,
            "",
        ),
        (
            "R01333",
            "R01333",
            "s_0775[c] + s_1198[c] + s_0803[c] = C00160[c] + s_1203[c] + 2 s_0794[c]",
            0,
            1000,
            "",
        ),
        (
            "RXN0-7229",
            "RXN0-7229",
            "2 C00160[c] + 2 ETR-Quinones[c] = 2 s_0779[c] + 2 s_0794[c] + 2 ETR-Quinols[c]",
            0,
            1000,
            "",
        ),
        (
            "R00013",
            "R00013",
            "2 s_0779[c] = C01146[c] + s_0456[c]",
            0,
            1000,
            "",
        ),
        (
            "R01747_NADPH",
            "R01747_NADPH",
            "C01146[c] + s_1212[c] + s_0794[c] = C00258[c] + s_1207[c]",
            0,
            1000,
            "",
        ),
        (
            "R01747_NADH",
            "R01747_NADH",
            "C01146[c] + s_1203[c] + s_0794[c] = C00258[c] + s_1198[c]",
            0,
            1000,
            "",
        ),
        (
            "R08572",
            "R08572",
            "C00258[c] + s_0434[c] = s_0188[c] + s_0394[c]",
            0,
            1000,
            "",
        ),
        (
            "R00472",
            "R00472",
            "s_0373[c] + s_0803[c] + s_0779[c] = s_0066[c] + s_0529[c]",
            0,
            1000,
            "",
        ),
        (
            "R00478",
            "R00478",
            "s_1003[c] + s_0779[c] = C03961[c]",
            0,
            1000,
            "",
        ),
        (
            "R00347",
            "R00347",
            "C03961[c] = s_1271[c] + s_0419[c]",
            0,
            1000,
            "",
        )
    ]

    dfRxns = pd.DataFrame(lRxns, columns = ['Id', 'Name', 'Equation', 'Lb', 'Ub', 'GPR'])
    for r in dfRxns.itertuples():
        mmL.addRxns(model, r)

    #metIsolatedExchgRxn(model, model.metabolites.get_by_id("C06337[c]"), newLB=-1000.0, newUB=0.0)
    metIsolatedExchgRxn(model, model.metabolites.get_by_id("C01380[c]"), newLB=-1000.0, newUB=0.0)
    metIsolatedExchgRxn(model, model.metabolites.get_by_id("ETR-Quinones[c]"), newLB=-1000.0, newUB=0.0)
    metIsolatedExchgRxn(model, model.metabolites.get_by_id("ETR-Quinols[c]"), newLB=0, newUB=1000.0)
    return model


def addTPAEGpathway_yarrowia(model):
    lMets = [
        ("C06337[c]", "Terephthalate [cytoplasm]", "c"),
        (
            "C06318[c]",
            "(3S,4R)-3,4-Dihydroxycyclohexa-1,5-diene-1,4-dicarboxylate [cytoplasm]",
            "c",
        ),
        ("C00230[c]", "3,4-Dihydroxybenzoate [cytoplasm]", "c"),
        ("C01163[c]", "3-Carboxy-cis,cis-muconate [cytoplasm]", "c"),
        ("C01278[c]", "2-Carboxy-2,5-dihydro-5-oxofuran-2-acetate [cytoplasm]", "c"),
        ("C03586[c]", "2-Oxo-2,3-dihydrofuran-5-acetate [cytoplasm]", "c"),
        ("C00846[c]", "3-Oxoadipate [cytoplasm]", "c"),
        ("C02232[c]", "3-Oxoadipyl-CoA [cytoplasm]", "c"),
        ("C04484[c]", "4-Carboxy-2-hydroxymuconate semialdehyde [cytoplasm]", "c"),
        ("C05375[c]", "2-Hydroxy-2-hydropyrone-4,6-dicarboxylate [cytoplasm]", "c"),
        ("C03671[c]", "2-Pyrone-4,6-dicarboxylate [cytoplasm]", "c"),
        ("C04434[c]", "(1E)-4-Oxobut-1-ene-1,2,4-tricarboxylate [cytoplasm]", "c"),
        ("C04115[c]", "4-Carboxy-4-hydroxy-2-oxoadipate [cytoplasm]", "c"),
        (
            "MNXM9831[c]",
            "(2Z,4Z)-2-hydroxy-5-carboxymuconate-6-semialdehyde [cytoplasm]",
            "c",
        ),
        ("C00682[c]", "(2Z,4E)-2-hydroxy-6-oxohexa-2,4-dienoate [cytoplasm]", "c"),
        ("C02501[c]", "2-Hydroxymuconate [cytoplasm]", "c"),
        ("C03453[c]", "gamma-Oxalocrotonate [cytoplasm]", "c"),
        ("C00596[c]", "2-Hydroxy-2,4-pentadienoate [cytoplasm]", "c"),
        ("C03589[c]", "4-Hydroxy-2-oxopentanoate [cytoplasm]", "c"),
        ("C01380[c]", "ethylene glycol [cytoplasm]", "c"),
        ("C00160[c]", "Glycolate [cytoplasm]", "c"),
        ("ETR-Quinones[c]", "electron-transfer quinone [cytoplasm]", "c"),
        ("ETR-Quinols[c]", "electron-transfer quinol [cytoplasm]", "c"),
        ("C01146[c]", "2-Hydroxy-3-oxopropanoate [cytoplasm]", "c"),
        ("C00258[c]", "D-Glycerate [cytoplasm]", "c"),
        ("C03961[c]", "erythro-3-Hydroxy-Ls-aspartate [cytoplasm]", "c"),
        ("s_1464[c]", "succinyl-CoA [cytoplasm]", "c"),
    ]

    model = mmL.addMets(model, lMets)

    lRxns = [
        (
            "R05148",
            "R05148",
            {
                model.metabolites.get_by_id("C06337[c]"): -1,
                model.metabolites.get_by_id("s_1275[c]"): -1,
                model.metabolites.get_by_id("s_1203[c]"): -1,
                model.metabolites.get_by_id("s_0794[c]"): -1,
                model.metabolites.get_by_id("C06318[c]"): 1,
                model.metabolites.get_by_id("s_1198[c]"): 1,
            },
            0,
            1000,
            "",
        ),
        (
            "R01633",
            "R01633",
            {
                model.metabolites.get_by_id("C06318[c]"): -1,
                model.metabolites.get_by_id("s_1198[c]"): -1,
                model.metabolites.get_by_id("C00230[c]"): 1,
                model.metabolites.get_by_id("s_0456[c]"): 1,
                model.metabolites.get_by_id("s_1203[c]"): 1,
                model.metabolites.get_by_id("s_0794[c]"): 1,
            },
            0,
            1000,
            "",
        ),
        (
            "R01631",
            "R01631",
            {
                model.metabolites.get_by_id("C00230[c]"): -1,
                model.metabolites.get_by_id("s_1275[c]"): -1,
                model.metabolites.get_by_id("C01163[c]"): 1,
            },
            0,
            1000,
            "",
        ),
        (
            "R03307",
            "R03307",
            {
                model.metabolites.get_by_id("C01163[c]"): -1,
                model.metabolites.get_by_id("C01278[c]"): 1,
            },
            0,
            1000,
            "",
        ),
        (
            "R03470",
            "R03470",
            {
                model.metabolites.get_by_id("C01278[c]"): -1,
                model.metabolites.get_by_id("C03586[c]"): 1,
                model.metabolites.get_by_id("s_0456[c]"): 1,
            },
            -1000,
            1000,
            "",
        ),
        (
            "R02991",
            "R02991",
            {
                model.metabolites.get_by_id("C03586[c]"): -1,
                model.metabolites.get_by_id("s_0803[c]"): -1,
                model.metabolites.get_by_id("C00846[c]"): 1,
            },
            0,
            1000,
            "",
        ),
        (
            "R02990",
            "R02990",
            {
                model.metabolites.get_by_id("s_1464[c]"): -1,
                model.metabolites.get_by_id("C00846[c]"): -1,
                model.metabolites.get_by_id("s_1458[c]"): 1,
                model.metabolites.get_by_id("C02232[c]"): 1,
            },
            0,
            1000,
            "",
        ),
        (
            "R00829",
            "R00829",
            {
                model.metabolites.get_by_id("s_0529[c]"): -1,
                model.metabolites.get_by_id("C02232[c]"): -1,
                model.metabolites.get_by_id("s_1464[c]"): 1,
                model.metabolites.get_by_id("s_0373[c]"): 1,
            },
            0,
            1000,
            "",
        ),
        (
            "R01632",
            "R01632",
            {
                model.metabolites.get_by_id("C00230[c]"): -1,
                model.metabolites.get_by_id("s_1275[c]"): -1,
                model.metabolites.get_by_id("C04484[c]"): 1,
            },
            0,
            1000,
            "",
        ),
        (
            "R04489",
            "R04489",
            {
                model.metabolites.get_by_id("C04484[c]"): -1,
                model.metabolites.get_by_id("C05375[c]"): 1,
            },
            0,
            1000,
            "",
        ),
        (
            "R04279",
            "R04279",
            {
                model.metabolites.get_by_id("C05375[c]"): -1,
                model.metabolites.get_by_id("s_1207[c]"): -1,
                model.metabolites.get_by_id("C03671[c]"): 1,
                model.metabolites.get_by_id("s_1212[c]"): 1,
                model.metabolites.get_by_id("s_0794[c]"): 1,
            },
            0,
            1000,
            "",
        ),
        (
            "R04277",
            "R04277",
            {
                model.metabolites.get_by_id("C03671[c]"): -1,
                model.metabolites.get_by_id("s_0803[c]"): -1,
                model.metabolites.get_by_id("C04434[c]"): 1,
            },
            0,
            1000,
            "",
        ),
        (
            "R04478",
            "R04478",
            {
                model.metabolites.get_by_id("C04434[c]"): -1,
                model.metabolites.get_by_id("s_0803[c]"): -1,
                model.metabolites.get_by_id("C04115[c]"): 1,
            },
            0,
            1000,
            "",
        ),
        (
            "R00350",
            "R00350",
            {
                model.metabolites.get_by_id("C04115[c]"): -1,
                model.metabolites.get_by_id("s_1271[c]"): 1,
                model.metabolites.get_by_id("s_1399[c]"): 1,
            },
            0,
            1000,
            "",
        ),
        (
            "R00217",
            "R00217",
            {
                model.metabolites.get_by_id("s_1271[c]"): -1,
                model.metabolites.get_by_id("s_1399[c]"): 1,
                model.metabolites.get_by_id("s_0456[c]"): 1,
            },
            0,
            1000,
            "",
        ),
        (
            "RXN-10882",
            "RXN-10882",
            {
                model.metabolites.get_by_id("C00230[c]"): -1,
                model.metabolites.get_by_id("s_1275[c]"): -1,
                model.metabolites.get_by_id("MNXM9831[c]"): 1,
            },
            0,
            1000,
            "",
        ),
        (
            "RXN-10883",
            "RXN-10883",
            {
                model.metabolites.get_by_id("MNXM9831[c]"): -1,
                model.metabolites.get_by_id("s_0794[c]"): -1,
                model.metabolites.get_by_id("C00682[c]"): 1,
                model.metabolites.get_by_id("s_0456[c]"): 1,
            },
            0,
            1000,
            "",
        ),
        (
            "R02762",
            "R02762",
            {
                model.metabolites.get_by_id("C00682[c]"): -1,
                model.metabolites.get_by_id("s_1198[c]"): -1,
                model.metabolites.get_by_id("s_0803[c]"): -1,
                model.metabolites.get_by_id("C02501[c]"): 1,
                model.metabolites.get_by_id("s_1203[c]"): 1,
                model.metabolites.get_by_id("s_0794[c]"): 1,
            },
            0,
            1000,
            "",
        ),
        (
            "R03966",
            "R03966",
            {
                model.metabolites.get_by_id("C02501[c]"): -1,
                model.metabolites.get_by_id("C03453[c]"): 1,
            },
            -1000,
            1000,
            "",
        ),
        (
            "R02602",
            "R02602",
            {
                model.metabolites.get_by_id("C03453[c]"): -1,
                model.metabolites.get_by_id("C00596[c]"): 1,
                model.metabolites.get_by_id("s_0456[c]"): 1,
            },
            0,
            1000,
            "",
        ),
        (
            "R02601",
            "R02601",
            {
                model.metabolites.get_by_id("C00596[c]"): -1,
                model.metabolites.get_by_id("s_0803[c]"): -1,
                model.metabolites.get_by_id("C03589[c]"): 1,
            },
            0,
            1000,
            "",
        ),
        (
            "R00750",
            "R00750",
            {
                model.metabolites.get_by_id("C03589[c]"): -1,
                model.metabolites.get_by_id("s_0359[c]"): 1,
                model.metabolites.get_by_id("s_1399[c]"): 1,
            },
            0,
            1000,
            "",
        ),
        (
            "R00228",
            "R00228",
            {
                model.metabolites.get_by_id("s_0359[c]"): -1,
                model.metabolites.get_by_id("s_0529[c]"): -1,
                model.metabolites.get_by_id("s_1198[c]"): -1,
                model.metabolites.get_by_id("s_0373[c]"): 1,
                model.metabolites.get_by_id("s_1203[c]"): 1,
                model.metabolites.get_by_id("s_0794[c]"): 1,
            },
            0,
            1000,
            "",
        ),
        (
            "RXN-21079",
            "RXN-21079",
            {
                model.metabolites.get_by_id("C01380[c]"): -1,
                model.metabolites.get_by_id("s_0359[c]"): 1,
                model.metabolites.get_by_id("s_0803[c]"): 1,
            },
            0,
            1000,
            "",
        ),
        (
            "R01781",
            "R01781",
            {
                model.metabolites.get_by_id("C01380[c]"): -1,
                model.metabolites.get_by_id("s_1198[c]"): -1,
                model.metabolites.get_by_id("s_0775[c]"): 1,
                model.metabolites.get_by_id("s_1203[c]"): 1,
                model.metabolites.get_by_id("s_0794[c]"): 1,
            },
            0,
            1000,
            "",
        ),
        (
            "R01333",
            "R01333",
            {
                model.metabolites.get_by_id("s_0775[c]"): -1,
                model.metabolites.get_by_id("s_1198[c]"): -1,
                model.metabolites.get_by_id("s_0803[c]"): -1,
                model.metabolites.get_by_id("C00160[c]"): 1,
                model.metabolites.get_by_id("s_1203[c]"): 1,
                model.metabolites.get_by_id("s_0794[c]"): 2,
            },
            0,
            1000,
            "",
        ),
        (
            "RXN0-7229",
            "RXN0-7229",
            {
                model.metabolites.get_by_id("C00160[c]"): -2,
                model.metabolites.get_by_id("ETR-Quinones[c]"): -2,
                model.metabolites.get_by_id("s_0779[c]"): 2,
                model.metabolites.get_by_id("s_0794[c]"): 2,
                model.metabolites.get_by_id("ETR-Quinols[c]"): 2,
            },
            0,
            1000,
            "",
        ),
        (
            "R00013",
            "R00013",
            {
                model.metabolites.get_by_id("s_0779[c]"): -2,
                model.metabolites.get_by_id("C01146[c]"): 1,
                model.metabolites.get_by_id("s_0456[c]"): 1,
            },
            0,
            1000,
            "",
        ),
        (
            "R01747_NADPH",
            "R01747_NADPH",
            {
                model.metabolites.get_by_id("C01146[c]"): -1,
                model.metabolites.get_by_id("s_1212[c]"): -1,
                model.metabolites.get_by_id("s_0794[c]"): -1,
                model.metabolites.get_by_id("C00258[c]"): 1,
                model.metabolites.get_by_id("s_1207[c]"): 1,
            },
            0,
            1000,
            "",
        ),
        (
            "R01747_NADH",
            "R01747_NADH",
            {
                model.metabolites.get_by_id("C01146[c]"): -1,
                model.metabolites.get_by_id("s_1203[c]"): -1,
                model.metabolites.get_by_id("s_0794[c]"): -1,
                model.metabolites.get_by_id("C00258[c]"): 1,
                model.metabolites.get_by_id("s_1198[c]"): 1,
            },
            0,
            1000,
            "",
        ),
        (
            "R08572",
            "R08572",
            {
                model.metabolites.get_by_id("C00258[c]"): -1,
                model.metabolites.get_by_id("s_0434[c]"): -1,
                model.metabolites.get_by_id("s_0188[c]"): 1,
                model.metabolites.get_by_id("s_0394[c]"): 1,
            },
            0,
            1000,
            "",
        ),
        (
            "R00472",
            "R00472",
            {
                model.metabolites.get_by_id("s_0373[c]"): -1,
                model.metabolites.get_by_id("s_0803[c]"): -1,
                model.metabolites.get_by_id("s_0779[c]"): -1,
                model.metabolites.get_by_id("s_0066[c]"): 1,
                model.metabolites.get_by_id("s_0529[c]"): 1,
            },
            0,
            1000,
            "",
        ),
        (
            "R00478",
            "R00478",
            {
                model.metabolites.get_by_id("s_1003[c]"): -1,
                model.metabolites.get_by_id("s_0779[c]"): -1,
                model.metabolites.get_by_id("C03961[c]"): 1,
            },
            0,
            1000,
            "",
        ),
        (
            "R00347",
            "R00347",
            {
                model.metabolites.get_by_id("C03961[c]"): -1,
                model.metabolites.get_by_id("s_1271[c]"): 1,
                model.metabolites.get_by_id("s_0419[c]"): 1,
            },
            0,
            1000,
            "",
        ),
        (
            "EX_C06337",
            "EX_C06337",
            {model.metabolites.get_by_id("C06337[c]"): 1},
            0,
            1000,
            "",
        ),
        (
            "EX_C01380",
            "EX_C01380",
            {model.metabolites.get_by_id("C01380[c]"): 1},
            0,
            1000,
            "",
        ),
    ]

    model = mmL.addRxns(model, lRxns)
    return model


# def exchgRxns4IsolatedMet(modelIn):
#     """Add exchange reactions for isolated metabolites"""
#     for m in modelIn.metabolites:
#         if len(m.reactions) == 1:
#             reazione = list(m.reactions)[0]
#             if reazione.lower_bound < 0 and reazione.upper_bound > 0:
#                 metIsolatedExchgRxn(modelIn, m, -1000.0, 1000.0)
#             elif reazione.lower_bound == 0 and reazione.upper_bound > 0:
#                 if m in reazione.reactants:
#                     metIsolatedExchgRxn(modelIn, m, -1000.0, 0.0)
#                 elif m in reazione.products:
#                     metIsolatedExchgRxn(modelIn, m, 0.0, 1000.0)
#             elif reazione.lower_bound < 0 and reazione.upper_bound == 0:
#                 if m in reazione.reactants:
#                     metIsolatedExchgRxn(modelIn, m, 0.0, 1000.0)
#                 elif m in reazione.products:
#                     metIsolatedExchgRxn(modelIn, m, -1000.0, 0.0)
#         elif len(m.reactions) > 1:
#             reazioni = list(m.reactions)
#             lonlySub = []
#             lonlyProd = []
#             lForward = []
#             lBackward = []
#             for r in reazioni:
#                 lonlySub.append(m in r.reactants)
#                 lonlyProd.append(m in r.products)
#                 if r.lower_bound == 0 and r.upper_bound > 0:
#                     lForward.append(True)
#                     lBackward.append(False)
#                 elif r.lower_bound < 0 and r.upper_bound == 0:
#                     lForward.append(False)
#                     lBackward.append(True)
#                 else:
#                     lForward.append(False)
#                     lBackward.append(False)
#             if all(lonlySub) is True and all(lForward) is True:
#                 metIsolatedExchgRxn(modelIn, m, -1000.0, 0.0)
#             elif all(lonlySub) is True and all(lBackward) is True:
#                 metIsolatedExchgRxn(modelIn, m, 0, 1000.0)
#             elif all(lonlyProd) is True and all(lForward) is True:
#                 metIsolatedExchgRxn(modelIn, m, 0, 1000.0)
#             elif all(lonlyProd) is True and all(lBackward) is True:
#                 metIsolatedExchgRxn(modelIn, m, -1000.0, 0.0)
#     return modelIn




# def setExtracFluxes(lMeasuredMets_consumption, lMeasuredMets_production, t1, t2, rowT1_c, rowT2_c, rowT1_p, rowT2_p, modelIn, dConfigParams=None, experiment = None):
#     """Set lower and upper boundaries for extracellular consumptions or production variations between two time points according to experimental values.
#
#     Keyword arguments:
#         rowT1 - is the row for the first time point of the dataframe including experimental extracellular consumption rates
#         rowT2 - is the row for the second time point of the dataframe including experimental extracellular consumption rates
#         modelIn - is the cobra model where the sink reactions are blocked
#         dConfigParams - is the configuration parameters dictionary
#     """
#
#     condition = dConfigParams['condition']
#
#     dMeasuredMet2Rxn_consumption = {}
#     for measuredMet in lMeasuredMets_consumption:
#         dMeasuredMet2Rxn_consumption[measuredMet] = modelIn.reactions.get_by_id(dConfigParams[measuredMet])
#
#     dMeasuredMet2Rxn_production = {}
#     for measuredMet in lMeasuredMets_production:
#         dMeasuredMet2Rxn_production[measuredMet] = modelIn.reactions.get_by_id(dConfigParams[measuredMet])
#
#     # glc = dConfigParams['glc']
#     # glcRxn = modelIn.reactions.get_by_id(glc)
#     # eg = dConfigParams['eg']
#     # egRxn = modelIn.reactions.get_by_id(eg)
#     # acetate = dConfigParams['acetate']
#     # acetateRxn = modelIn.reactions.get_by_id(acetate)
#     # if dConfigParams["modelName"] in ["yeast8", "yeast-GEM_9"]:
#     #     etoh = dConfigParams['etoh']
#     #     glycerol = dConfigParams['glycerol']
#     #     etohRxn = modelIn.reactions.get_by_id(etoh)
#     #     glycerolRxn = modelIn.reactions.get_by_id(glycerol)
#
#     # create dictionary where for each measured metabolite there is the corresponding consumption or production rate in a given time interval
#     if t1 == '0' and t2 == '0':
#         # if dConfigParams["modelName"] in ["yeast8", "yeast-GEM_9"]:
#         #     lStrings  = ['Glucose', 'Ethanol', 'Glycerol', 'Acetate', 'EG']
#         # elif dConfigParams["modelName"] == "eciYali_sbmlValidated":
#         #     lStrings  = ['Glucose']
#         # for s in lStrings:
#         dValues_consumption = {}
#         for s in lMeasuredMets_consumption:
#             lvalues = [rowT2_c[s+'_rep1'].values[0],
#                         rowT2_c[s+'_rep2'].values[0],
#                         rowT2_c[s+'_rep3'].values[0]]
#             dfTmp = pd.DataFrame({'flux': lvalues})
#             dValues_consumption[s] = [min(dfTmp['flux'].tolist()), max(dfTmp['flux'].tolist())]
#
#         dValues_production = {}
#         for s in lMeasuredMets_production:
#             lvalues = [rowT2_p[s+'_rep1'].values[0],
#                         rowT2_p[s+'_rep2'].values[0],
#                         rowT2_p[s+'_rep3'].values[0]]
#             dfTmp = pd.DataFrame({'flux': lvalues})
#             dValues_production[s] = [min(dfTmp['flux'].tolist()), max(dfTmp['flux'].tolist())]
#
#     #elif t1 in ['0','2','4','6','8'] and t2 in ['0','2','4','6','8']:
#     else:
#         #lStrings  = ['Glucose', 'Ethanol', 'Glycerol', 'Acetate']
#         dValues_consumption = {}
#         for s in lMeasuredMets_consumption:
#         #for s in lStrings:
#             lvalues = [(rowT2_c[s+'_rep1'].values[0] - rowT1_c[s+'_rep1'].values[0])/(int(t2)-int(t1)),
#                     (rowT2_c[s+'_rep2'].values[0] - rowT1_c[s+'_rep2'].values[0])/(int(t2)-int(t1)),
#                     (rowT2_c[s+'_rep3'].values[0] - rowT1_c[s+'_rep3'].values[0])/(int(t2)-int(t1))]
#             # lvalues = [(rowT2[s+'_rep1'].values[0] - rowT1[s+'_rep1'].values[0]),
#             #         (rowT2[s+'_rep2'].values[0] - rowT1[s+'_rep2'].values[0]),
#             #         (rowT2[s+'_rep3'].values[0] - rowT1[s+'_rep3'].values[0])]
#             dfTmp = pd.DataFrame({'flux': lvalues})
#             dValues_consumption[s] = [min(dfTmp['flux'].tolist()), max(dfTmp['flux'].tolist())]
#
#         dValues_production = {}
#         for s in lMeasuredMets_production:
#         #for s in lStrings:
#             lvalues = [(rowT2_p[s+'_rep1'].values[0] - rowT1_p[s+'_rep1'].values[0])/(int(t2)-int(t1)),
#                     (rowT2_p[s+'_rep2'].values[0] - rowT1_p[s+'_rep2'].values[0])/(int(t2)-int(t1)),
#                     (rowT2_p[s+'_rep3'].values[0] - rowT1_p[s+'_rep3'].values[0])/(int(t2)-int(t1))]
#             # lvalues = [(rowT2[s+'_rep1'].values[0] - rowT1[s+'_rep1'].values[0]),
#             #         (rowT2[s+'_rep2'].values[0] - rowT1[s+'_rep2'].values[0]),
#             #         (rowT2[s+'_rep3'].values[0] - rowT1[s+'_rep3'].values[0])]
#             dfTmp = pd.DataFrame({'flux': lvalues})
#             dValues_production[s] = [min(dfTmp['flux'].tolist()), max(dfTmp['flux'].tolist())]
#
#         if "EG" in condition:
#             if '_'.join(experiment) in ["endogenous","endogenous_aox", "aox_endogenous", "aox"]:
#                 dValues_consumption["eg"] = [0, 0] # deciso di settarlo come non consumato
#             else:
#                 dValues_consumption["eg"] = [0, -1000]
#
#         if "Acetate" in condition:
#             dValues_consumption["eg"] = [0, -1000]
#
#     # use the dValues dictionary to set constraints for the exchange reactions corresponding to each measured metabolite
#     for extrMet in dValues_production:
#         if len(dMeasuredMet2Rxn_production[extrMet].reactants) != 0 and len(dMeasuredMet2Rxn_production[extrMet].products) == 0: # A = 0
#             try:
#                 dMeasuredMet2Rxn_production[extrMet].upper_bound = 1 * (dValues_production[extrMet][1])
#                 dMeasuredMet2Rxn_production[extrMet].lower_bound = 1 * (dValues_production[extrMet][0])
#             except:
#                 dMeasuredMet2Rxn_production[extrMet].lower_bound = 1 * (dValues_production[extrMet][0])
#                 dMeasuredMet2Rxn_production[extrMet].upper_bound = 1 * (dValues_production[extrMet][1])
#         elif len(dMeasuredMet2Rxn_production[extrMet].reactants) == 0 and len(dMeasuredMet2Rxn_production[extrMet].products) != 0: # 0 = A
#             try:
#                 dMeasuredMet2Rxn_production[extrMet].upper_bound = -1 * (dValues_production[extrMet][1])
#                 dMeasuredMet2Rxn_production[extrMet].lower_bound = -1 * (dValues_production[extrMet][0])
#             except:
#                 dMeasuredMet2Rxn_production[extrMet].lower_bound = -1 * (dValues_production[extrMet][0])
#                 dMeasuredMet2Rxn_production[extrMet].upper_bound = -1 * (dValues_production[extrMet][1])
#
#     for extrMet in dValues_consumption:
#         if dValues_consumption[extrMet][1] < 0:
#             if len(dMeasuredMet2Rxn_consumption[extrMet].reactants) != 0 and len(dMeasuredMet2Rxn_consumption[extrMet].products) == 0: # A = 0
#                 try:
#                     dMeasuredMet2Rxn_consumption[extrMet].upper_bound = 0
#                     dMeasuredMet2Rxn_consumption[extrMet].lower_bound = 1 * (dValues_consumption[extrMet][1])
#                 except:
#                     dMeasuredMet2Rxn_consumption[extrMet].lower_bound = 1 * (dValues_consumption[extrMet][1])
#                     dMeasuredMet2Rxn_consumption[extrMet].upper_bound = 0
#             elif len(dMeasuredMet2Rxn_consumption[extrMet].reactants) == 0 and len(dMeasuredMet2Rxn_consumption[extrMet].products) != 0: # 0 = A
#                 try:
#                     dMeasuredMet2Rxn_consumption[extrMet].upper_bound = -1 * (dValues_consumption[extrMet][1])
#                     dMeasuredMet2Rxn_consumption[extrMet].lower_bound = 0
#                 except:
#                     dMeasuredMet2Rxn_consumption[extrMet].lower_bound = 0
#                     dMeasuredMet2Rxn_consumption[extrMet].upper_bound = -1 * (dValues_consumption[extrMet][1])
#         else:
#             dMeasuredMet2Rxn_consumption[extrMet].upper_bound = 0
#             dMeasuredMet2Rxn_consumption[extrMet].lower_bound = 0
#
#     # if "EG" in condition:
#     #     try:
#     #         dMeasuredMet2Rxn["eg"].lower_bound = -1 * dValues["eg"][1]
#     #         dMeasuredMet2Rxn["eg"].upper_bound = 0
#     #     except:
#     #         dMeasuredMet2Rxn["eg"].upper_bound = 0
#     #         dMeasuredMet2Rxn["eg"].lower_bound = -1 * dValues["eg"][1]
#
#     # if dConfigParams["modelName"] in ["yeast8", "yeast-GEM_9"]:
#         # try:
#         #     dMeasuredMet2Rxn["ethanol"].upper_bound = 1 * (dValues["ethanol"][1])
#         #     dMeasuredMet2Rxn["ethanol"].lower_bound = 1 * (dValues["ethanol"][0])
#         # except:
#         #     dMeasuredMet2Rxn["ethanol"].lower_bound = 1 * (dValues["ethanol"][0])
#         #     dMeasuredMet2Rxn["ethanol"].upper_bound = 1 * (dValues["ethanol"][1])
#         #
#         # try:
#         #     dMeasuredMet2Rxn["glycerol"].lower_bound = 1 * (dValues["glycerol"][0])
#         #     dMeasuredMet2Rxn["glycerol"].upper_bound = 1 * (dValues["glycerol"][1])
#         # except:
#         #     dMeasuredMet2Rxn["glycerol"].upper_bound = 1 * (dValues["glycerol"][1])
#         #     dMeasuredMet2Rxn["glycerol"].lower_bound = 1 * (dValues["glycerol"][0])
#         #
#         # try:
#         #     dMeasuredMet2Rxn["acetate"].lower_bound = 1 * (dValues["acetate"][0])
#         #     dMeasuredMet2Rxn["acetate"].upper_bound = 1 * (dValues["acetate"][1])
#         # except:
#         #     dMeasuredMet2Rxn["acetate"].upper_bound = 1 * (dValues["acetate"][1])
#         #     dMeasuredMet2Rxn["acetate"].lower_bound = 1 * (dValues["acetate"][0])
#
#     # if '_'.join(condition) in ["EG"]:
#     #     try:
#     #         dMeasuredMet2Rxn["glucose"].lower_bound = 0
#     #         dMeasuredMet2Rxn["glucose"].upper_bound = 0
#     #     except:
#     #         dMeasuredMet2Rxn["glucose"].upper_bound = 0
#     #         dMeasuredMet2Rxn["glucose"].lower_bound = 0
#     #
#     #     try:
#     #         dMeasuredMet2Rxn["eg"].lower_bound = -1 * dValues["eg"][1]
#     #         dMeasuredMet2Rxn["eg"].upper_bound = 0
#     #     except:
#     #         dMeasuredMet2Rxn["eg"].upper_bound = 0
#     #         dMeasuredMet2Rxn["eg"].lower_bound = -1 * dValues["eg"][1]
#
#         # if dConfigParams["modelName"] in ["yeast8", "yeast-GEM_9"]:
#         #     try:
#         #         #egRxn.lower_bound = -1 * (dValues["EG"][1])
#         #         egRxn.lower_bound = -1000
#         #         egRxn.upper_bound = 0
#         #     except:
#         #         egRxn.upper_bound = 0
#         #         #egRxn.lower_bound = -1 * (dValues["EG"][1])
#         #         egRxn.lower_bound = -1000
#         # elif dConfigParams["modelName"] == "eciYali_sbmlValidated":
#         #     try:
#         #         egRxn.lower_bound = -1000
#         #         egRxn.upper_bound = 0
#         #     except:
#         #         egRxn.upper_bound = 0
#         #         egRxn.lower_bound = -1000
#
#     # elif '_'.join(condition) in ["Glc_EG", "EG_Glc"]:
#         # if dValues["glucose"][1] > 0:
#         #     try:
#         #         dMeasuredMet2Rxn["glucose"].lower_bound = -1 * dValues["glucose"][1]
#         #         dMeasuredMet2Rxn["glucose"].upper_bound = 0
#         #     except:
#         #         dMeasuredMet2Rxn["glucose"].upper_bound = 0
#         #         dMeasuredMet2Rxn["glucose"].lower_bound = -1 * dValues["glucose"][1]
#         # else:
#         #     try:
#         #         dMeasuredMet2Rxn["glucose"].lower_bound = 1 * dValues["glucose"][1]
#         #         dMeasuredMet2Rxn["glucose"].upper_bound = 0
#         #     except:
#         #         dMeasuredMet2Rxn["glucose"].upper_bound = 0
#         #         dMeasuredMet2Rxn["glucose"].lower_bound = 1 * dValues["glucose"][1]
#         #
#         # try:
#         #     dMeasuredMet2Rxn["eg"].lower_bound = -1 * dValues["eg"][1]
#         #     dMeasuredMet2Rxn["eg"].upper_bound = 0
#         # except:
#         #     dMeasuredMet2Rxn["eg"].upper_bound = 0
#         #     dMeasuredMet2Rxn["eg"].lower_bound = -1 * dValues["eg"][1]
#
#         # if dConfigParams["modelName"] in ["yeast8", "yeast-GEM_9"]:
#         #     if dValues["glucose"][1] > 0:
#         #         try:
#         #             dMeasuredMet2Rxn["glucose"].lower_bound = -1 * dValues["glucose"][1]
#         #             dMeasuredMet2Rxn["glucose"].upper_bound = 0
#         #         except:
#         #             dMeasuredMet2Rxn["glucose"].upper_bound = 0
#         #             dMeasuredMet2Rxn["glucose"].lower_bound = -1 * dValues["glucose"][1]
#         #     else:
#         #         try:
#         #             dMeasuredMet2Rxn["glucose"].lower_bound = 1 * dValues["glucose"][0]
#         #             dMeasuredMet2Rxn["glucose"].upper_bound = 0
#         #         except:
#         #             dMeasuredMet2Rxn["glucose"].upper_bound = 0
#         #             dMeasuredMet2Rxn["glucose"].lower_bound = 1 * dValues["glucose"][0]
#         #
#         #     if '_'.join(experiment) == "endogeno":
#         #         try:
#         #             dMeasuredMet2Rxn["eg"].lower_bound = 0
#         #             dMeasuredMet2Rxn["eg"].upper_bound = 0
#         #         except:
#         #             dMeasuredMet2Rxn["eg"].upper_bound = 0
#         #             dMeasuredMet2Rxn["eg"].lower_bound = 0
#         #     else:
#         #         try:
#         #             dMeasuredMet2Rxn["eg"].lower_bound = -1 * (dValues["eg"][1])
#         #             # dMeasuredMet2Rxn["eg"].lower_bound = 1 * (dValues["Glucose"][0] * 10)
#         #             dMeasuredMet2Rxn["eg"].upper_bound = 0
#         #         except:
#         #             dMeasuredMet2Rxn["eg"].upper_bound = 0
#         #             dMeasuredMet2Rxn["eg"].lower_bound = -1 * (dValues["eg"][1])
#         #             # egRxn.lower_bound = 1 * (dValues["Glucose"][0] * 10)
#         # elif dConfigParams["modelName"] == "eciYali_sbmlValidated":
#         #     try:
#         #         dMeasuredMet2Rxn["glucose"].lower_bound = 0
#         #         dMeasuredMet2Rxn["glucose"].upper_bound = 1 * dValues["glucose"][1]
#         #     except:
#         #         dMeasuredMet2Rxn["glucose"].upper_bound = 1 * dValues["glucose"][1]
#         #         dMeasuredMet2Rxn["glucose"].lower_bound = 0
#         #
#         #     try:
#         #         dMeasuredMet2Rxn["eg"].lower_bound = -1000
#         #         dMeasuredMet2Rxn["eg"].upper_bound = 0
#         #     except:
#         #         dMeasuredMet2Rxn["eg"].upper_bound = 0
#         #         dMeasuredMet2Rxn["eg"].lower_bound = -1000
#
#     # elif '_'.join(condition) in ["Glc"]:
#     #     if dConfigParams["modelName"] == "eciYali_sbmlValidated":
#     #         try:
#     #             dMeasuredMet2Rxn["glucose"].lower_bound = 0
#     #             dMeasuredMet2Rxn["glucose"].upper_bound = 1 * dValues["glucose"][1]
#     #         except:
#     #             dMeasuredMet2Rxn["glucose"].upper_bound = 1 * dValues["glucose"][1]
#     #             dMeasuredMet2Rxn["glucose"].lower_bound = 0
#     #
#     #     try:
#     #         dMeasuredMet2Rxn["eg"].lower_bound = 0
#     #         dMeasuredMet2Rxn["eg"].upper_bound = 0
#     #     except:
#     #         dMeasuredMet2Rxn["eg"].upper_bound = 0
#     #         dMeasuredMet2Rxn["eg"].lower_bound = 0
#     # elif '_'.join(condition) in ["EG_Acetate", "Acetate_EG"]:
#     #     try:
#     #         dMeasuredMet2Rxn["glucose"].lower_bound = 0
#     #         dMeasuredMet2Rxn["glucose"].upper_bound = 0
#     #     except:
#     #         dMeasuredMet2Rxn["glucose"].upper_bound = 0
#     #         dMeasuredMet2Rxn["glucose"].lower_bound = 0
#     #
#     #     try:
#     #         dMeasuredMet2Rxn["eg"].lower_bound = -1000
#     #         dMeasuredMet2Rxn["eg"].upper_bound = 0
#     #     except:
#     #         dMeasuredMet2Rxn["eg"].upper_bound = 0
#     #         dMeasuredMet2Rxn["eg"].lower_bound = -1000
#
#         # try:
#         #     dMeasuredMet2Rxn["acetate"].lower_bound = 0
#         #     dMeasuredMet2Rxn["acetate"].upper_bound = 1000
#         # except:
#         #     dMeasuredMet2Rxn["acetate"].upper_bound = 1000
#         #     dMeasuredMet2Rxn["acetate"].lower_bound = 0
#
#         # dMeasuredMet2Rxn["acetate"].lower_bound = -1000
#
#     return modelIn


def setO2Flux(modelIn, wetTime = None, biomassFluxT0 = None, dConfigParams=None):
    """Setting of O2 uptake flux: variation of O2 uptake rate to fit different ODs of the cell population (o2Setting is set to fixed) or working on an unique cell left O2 consumption rate as free (o2Setting is set to free).

    Keyword arguments:
        modelIn - is the cobra model where the sink reactions are blocked
        wetTime - is the current experimental time point
        biomassFluxT0 - is the biomass flux value computed at the initial time point
        dConfigParams - is the configuration parameters dictionary
    """
    o2Exch = dConfigParams['oxygenRxn']
    o2Setting = dConfigParams['oxygen']
    biomassRxn = setBiomassReaction(modelIn, dConfigParams = dConfigParams)
    o2 = modelIn.reactions.get_by_id(o2Exch)
    if o2Setting == "free":
        o2.lower_bound = -1000
        return modelIn, o2.lower_bound

    elif o2Setting == "fixed":
        if wetTime == '0':
            try:
                o2.upper_bound = 0
                o2.lower_bound = -1
            except:
                o2.lower_bound = -1
                o2.upper_bound = 0
            fba = modelIn.optimize()
            biomassFlux = fba.fluxes[biomassRxn]
            return modelIn, biomassFlux,o2.lower_bound
        else:
            # Parallelizing using Pool.apply().
            # Init multiprocessing.Pool()
            # pool = mp.Pool(mp.cpu_count())
            # 'pool.apply' the varyOxygen2fitOD function
            # deltaWetSilico, biomassFlux_wO2MinFlux, O2flux = pool.apply(varyOxygen2fitOD, args=(modelIn, biomassFluxT0, dConfigParams))
            # pool.close()
            dOx2Biomass = varyOxygen2fitOD(modelIn, biomassFluxT0, dConfigParams)
            try:
                o2.upper_bound = 0
                o2.lower_bound = -1 * O2flux
            except:
                o2.lower_bound = -1 * O2flux
                o2.upper_bound = 0
            return modelIn, biomassFlux_wO2MinFlux, deltaWetSilico, o2.lower_bound

def varyOxygen2fitOD(modelIn, biomassFluxT0 = None, dConfigParams=None):
    """Setting of O2 uptake flux: variation of O2 uptake rate to fit different ODs of the cell population (o2Setting is set to fixed) or working on an unique cell left O2 consumption rate as free (o2Setting is set to free).

    Keyword arguments:
        modelIn - is the cobra model where the sink reactions are blocked
        biomassFluxT0 - is the biomass flux value computed at the initial time point
        dConfigParams - is the configuration parameters dictionary
    """
    o2Exch = dConfigParams['oxygenRxn']
    biomassRxn = dConfigParams['biomassRxn']

    dOx2Biomass = {}
    o2flux = 1.01
    for o2flux in np.arange(1.01, 1000.01, 0.01):
        ratioSilico = changeyLbAndTestBiomass(modelId, o2flux, o2Exch, biomassRxn, biomassFluxT0)
        dOx2Biomass[o2flux]=ratioSilico

    return dOx2Biomass

def changeyLbAndTestBiomass(modelId, v, rxnId, biomassRxnId, biomassFluxT0 =None):
    rxn = modelId.reactions.get_by_id(rxnId)
    try:
        rxn.upper_bound = 0
        rxn.lower_bound = -1 * v
    except:
        rxn.lower_bound = -1 * v
        rxn.upper_bound = 0

    try:
        fba = modelId.optimize()
        biomassFlux = fba.fluxes[biomassRxnId]
        ratioSilico = biomassFlux / biomassFluxT0
    except:
        ratioSilico = None
    return ratioSilico

def averageWetOD(row):
    """Compute average wet OD (necessary for generate plot OD and silico biomass vs. time).

    Keyword arguments:
        row - is the row for a specific time point of the dataframe including experimental extracellular consumption rates and optical densities
    """
    return (row.OD_rep1 + row.OD_rep2 + row.OD_rep3) / 3

def setBiomassReaction(modelIn, dConfigParams=None, rxnName=None):
    """Sets the biomass reaction name accoroding to the model or the given rxnName and sets its boundaries.

    Keyword arguments:
         dConfigParams - is the configuration parameters dictionary where it should be stored a validated model name
         rxnName - [None | string] is the reaction name to be used as biomass reaction (if not None, replaces the default one )

    Return:
         [biomassRxn | None] The biomass reaction name. If None somethong went wrong.
    """

    modName = dConfigParams["modelName"]
    biomassRxn = None
    if rxnName is not None:
        biomassRxn = rxnName
    if modName in ["yeast8", "yeast-GEM_9"]:
        biomassRxn = "r_2111"
    elif modName == "eciYali_sbmlValidated":
        biomassRxn = "xBIOMASS"

    modelIn.reactions.get_by_id(biomassRxn).lower_bound = 0
    modelIn.reactions.get_by_id(biomassRxn).upper_bound = 1000
    return biomassRxn
