import os
import sys
import time
import cobra as cb
import pandas as pd
import re
import numpy as np
import matplotlib.pyplot as plt
import genericLib as gL
import configParamsLib as cpL
import fluxAnalysisLib as faL

timeStamp = gL.getTimeStamp()

print("----------timeStamp:\t", timeStamp, "\n")

####################################################################################################
# Setting the parameters

# setting working dirs
print(gL.dDirs, "\n")
RAWDIR = gL.dDirs["raw"]
OUTDIR = gL.dDirs["out"]
MODELDIR = gL.dDirs["mods"]
FIGUREDIR = gL.dDirs["figs"]

# setting experiment params
dConfParams = cpL.loadConfigParams()
modelName = dConfParams["modelName"]
condition = dConfParams["condition"]
mediumSetting = dConfParams["medium"]
o2Setting = dConfParams["oxygen"]
experiment = dConfParams["kind"]

logStrm = gL.logFileOpen(gL.dDirs["logs"], timeStamp, gL.getBaseName(sys.argv[0]))
gL.toLog(logStrm, f"timeStamp: {timeStamp}")
cpL.toLogConfigParams(logStream=logStrm)
print("\n-----Loaded settings:\n")
print("Model name:\t", modelName)
print("condition:\t", condition)
print("medium:\t", mediumSetting)
print("oxygen:\t", o2Setting)

####################################################################################################

####################################################################################################
print(f"\n-----Step 0: Load {modelName} model\n")
gL.toLog(logStrm, f"\n-----Step 0: Load {modelName} model\n")
####################################################################################################

model = cb.io.read_sbml_model(os.path.join(MODELDIR, modelName + ".xml"))
model.solver = "glpk"

# print("Number of reactions:\t", len(model.reactions))
# print("Number of metabolites:\t", len(model.metabolites))
# print("Number of genes:\t", len(model.genes))
# print("Obj:\t", model.objective, "\n")

gL.toLog(logStrm, f"Number of reactions:\t{len(model.reactions)}")
gL.toLog(logStrm, f"Number of metabolites:\t{len(model.metabolites)}")
gL.toLog(logStrm, f"Number of genes:\t{len(model.genes)}")
gL.toLog(logStrm, f"Obj: {model.objective}\n")

## Biomass production
if modelName == "yeast8":
    biomassRxn = "r_2111"
elif modelName == "eciYali_sbmlValidated":
    biomassRxn = "xBIOMASS"  ## Yarrowia

model.reactions.get_by_id(biomassRxn).lower_bound = 0
model.reactions.get_by_id(biomassRxn).upper_bound = 1000

## Change boundaries of exchange reactions from infinite to the default value (1000 in this case)
for rxn in model.reactions:
    if (len(rxn.reactants) == 0 and len(rxn.products) != 0) or (
        len(rxn.reactants) != 0 and len(rxn.products) == 0
    ):
        if rxn.lower_bound == float("inf"):
            rxn.lower_bound = 1000
        elif rxn.lower_bound == -float("inf"):
            rxn.lower_bound = -1000
        elif rxn.upper_bound == float("inf"):
            rxn.upper_bound = 1000
        elif rxn.upper_bound == -float("inf"):
            rxn.upper_bound = -1000

## Add exchange reaction for isolated metabolites
for m in model.metabolites:
    if len(m.reactions) == 1:
        reazione = list(m.reactions)[0]

        if reazione.lower_bound == -1000 and reazione.upper_bound == 1000:
            metNew = m.id.replace("[", "")
            metNew = metNew.replace("]", "")
            metNew = metNew.replace("_", "")
            lRxns = [
                (
                    "EX_" + metNew,
                    metNew + " exchange",
                    {model.metabolites.get_by_id(m.id): -1},
                    -1000,
                    1000,
                    "",
                )
            ]
            model = faL.addRxns(model, lRxns)

        elif reazione.lower_bound == 0 and reazione.upper_bound == 1000:
            if m in reazione.reactants:
                metNew = m.id.replace("[", "")
                metNew = metNew.replace("]", "")
                metNew = metNew.replace("_", "")
                lRxns = [
                    (
                        "EX_" + metNew,
                        metNew + " exchange",
                        {model.metabolites.get_by_id(m.id): -1},
                        -1000,
                        0,
                        "",
                    )
                ]
                model = faL.addRxns(model, lRxns)
            elif m in reazione.products:
                metNew = m.id.replace("[", "")
                metNew = metNew.replace("]", "")
                metNew = metNew.replace("_", "")
                lRxns = [
                    (
                        "EX_" + metNew,
                        metNew + " exchange",
                        {model.metabolites.get_by_id(m.id): -1},
                        0,
                        1000,
                        "",
                    )
                ]
                model = faL.addRxns(model, lRxns)

        elif reazione.lower_bound == -1000 and reazione.upper_bound == 0:
            if m in reazione.reactants:
                metNew = m.id.replace("[", "")
                metNew = metNew.replace("]", "")
                metNew = metNew.replace("_", "")
                lRxns = [
                    (
                        "EX_" + metNew,
                        metNew + " exchange",
                        {model.metabolites.get_by_id(m.id): -1},
                        0,
                        1000,
                        "",
                    )
                ]
                model = faL.addRxns(model, lRxns)
            elif m in reazione.products:
                metNew = m.id.replace("[", "")
                metNew = metNew.replace("]", "")
                metNew = metNew.replace("_", "")
                lRxns = [
                    (
                        "EX_" + metNew,
                        metNew + " exchange",
                        {model.metabolites.get_by_id(m.id): -1},
                        -1000,
                        0,
                        "",
                    )
                ]
                model = faL.addRxns(model, lRxns)


####################################################################################################
print(
    "\n-----Step 1: Block all uptake reactions (cosi che entra solo quello che troviamo nel medium wet)\n"
)
gL.toLog(
    logStrm,
    "\n-----Step 1: Block all uptake reactions (cosi che entra solo quello che troviamo nel medium wet)\n",
)
####################################################################################################
# print('Uptake reactions che sto bloccando:\n')
if modelName == "yeast8":
    for rxn in model.reactions:
        if len(rxn.reactants) != 0 and len(rxn.products) == 0 and rxn.lower_bound < 0:
            rxn.lower_bound = 0
            # print('rxn:\t', rxn.id, '\t', rxn.name)
elif modelName == "eciYali_sbmlValidated":
    for rxn in model.reactions:
        if (
            len(rxn.reactants) == 0
            and len(rxn.products) != 0
            and rxn.id.startswith("prot_") == False
        ):
            rxn.lower_bound = 0
            rxn.upper_bound = 0
    ## nel modello di yarrowia ogni exchange è splittata in due (se l'ID della prima e' 'E1', l'inversa si chimera' 'E1_REV'; le exchange di entrata sono quelle col suffisso _REV)


####################################################################################################
print("\n-----Step 2: Add medium composition data\n")
gL.toLog(logStrm, "\n-----Step 2: Add medium composition data\n")
####################################################################################################

if condition == "YNB_Glc":
    dfMedium = pd.read_csv(
        os.path.join(RAWDIR, "mediumComposition_YNB_Glc.tsv"),
        sep="\t",
        dtype={"Concentration": float},
    )

elif condition == "YNB_Glc_EG" or condition == "YNB_EG":
    if condition == "YNB_Glc_EG":
        dfMedium = pd.read_csv(
            os.path.join(RAWDIR, "mediumComposition_YNB_Glc_EG.tsv"),
            sep="\t",
            dtype={"Concentration": float},
        )
    elif condition == "YNB_EG":
        dfMedium = pd.read_csv(
            os.path.join(RAWDIR, "mediumComposition_YNB_EG.tsv"),
            sep="\t",
            dtype={"Concentration": float},
        )

    lMets = [
        ("C01380[e]", "Ethylene glycol [extracellular]", "e"),
        ("C01380[c]", "Ethylene glycol [cytoplasm]", "c"),
    ]

    model = faL.addMets(model, lMets)

    lRxns = [
        (
            "EX_C01380",
            "EG uptake",
            {model.metabolites.get_by_id("C01380[e]"): -1},
            -1000,
            0,
            "",
        ),
        (
            "Transport_EG_e_c",
            "EG transport extracellular_cytoplasm",
            {
                model.metabolites.get_by_id("C01380[e]"): -1,
                model.metabolites.get_by_id("C01380[c]"): 1,
            },
            -1000,
            1000,
            "",
        ),
    ]

    model = faL.addRxns(model, lRxns)


# add Ferric ion exchange reaction (e' incluso tra i metaboliti del medium ma non e' presente la reazione di exchange nel modello)
reaction = cb.Reaction("EX_s3936")
reaction.name = "Ferric ion exchange"
reaction.lower_bound = -1000
reaction.upper_bound = 0
reaction.add_metabolites({model.metabolites.get_by_id("s_3936[e]"): -1})
model.add_reactions([reaction])

for rxnMedium in dfMedium.itertuples():
    rxn = model.reactions.get_by_id(rxnMedium.RxnID)
    if mediumSetting == "wet":
        rxn.lower_bound = -1 * rxnMedium.Concentration
    elif mediumSetting == "free":
        rxn.lower_bound = -1000
        if condition == "YNB_EG":
            model.reactions.get_by_id("r_1714").lower_bound = 0

    gL.toLog(
        logStrm, f"Rxn:\t{rxn.id}\t{rxn.name}\t{rxn.lower_bound}\t{rxn.upper_bound}"
    )

ngam = model.reactions.get_by_id("r_4046")
ngam.lower_bound = 0
ngam.upper_bound = 1000
gL.toLog(logStrm, f"Rxn ngam:\t{ngam.id}\t{ngam.lower_bound}\t{ngam.upper_bound}")

# riattivo exchange di O2
o2 = model.reactions.get_by_id("r_1992")
o2.lower_bound = -1000
o2.upper_bound = 0

# riattivo exchange di H+
h = model.reactions.get_by_id("r_1832")
h.lower_bound = -1000

# check Glucose and EG constraints
gL.toLog(logStrm, "\nCheck Glucose and EG constraints:")
gL.toLog(
    logStrm,
    f"""Glc:\t{model.reactions.get_by_id("r_1714").lower_bound}\t{model.reactions.get_by_id("r_1714").upper_bound}""",
)
if condition == "YNB_Glc_EG" or condition == "YNB_EG":
    gL.toLog(
        logStrm,
        f"""EG:\t{model.reactions.get_by_id("EX_C01380").lower_bound}\t{model.reactions.get_by_id("EX_C01380").upper_bound}\n""",
    )
else:
    print("EG:\t 0, 0,\n")

####################################################################################################
print(
    "\n-----Step 3: Set extracellular fluxes and O2 uptake flux, and add heterologous pathways\n"
)
gL.toLog(
    logStrm,
    "\n-----Step 3: Set extracellular fluxes and O2 uptake flux, and add heterologous pathways\n",
)
####################################################################################################

glc = "r_1714"
etoh = "r_1761"
glycerol = "r_1808"
acetate = "r_1634"

dictTypes = {
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

if condition == "YNB_Glc":
    conditionFile = "consumptionProduction_YNB_Glc.tsv"
elif condition == "YNB_Glc_EG":
    conditionFile = "consumptionProduction_YNB_Glc_EG.tsv"
elif condition == "YNB_EG":
    conditionFile = "consumptionProduction_YNB_EG.tsv"

dfConsProd = pd.read_csv(
    os.path.join(RAWDIR, conditionFile),
    sep="\t",
    dtype=dictTypes,
)

## variazione dell'O2 per ottenere le diverse OD della popolazione di cellule oppure ragiono sulla singola cellula e lasciamo svincolato e attivo il consumo di O2

## Experiment 1: add all pathways
lGrowth_wet = []
lGrowth_silico = []

fileSuffixName = (
    condition
    + "_medium_"
    + mediumSetting
    + "_O2_"
    + o2Setting
    + "_experiment_"
    + experiment
)

gL.toLog(logStrm, "\n>>> Wet time: 0h")
row0 = dfConsProd.iloc[0]
biomassFlux0, wetOD0 = faL.doExperiment(
    o2Setting,
    "y8_pFBA_" + fileSuffixName + "_t0_" + timeStamp,
    "y8_FBA_" + fileSuffixName + "_t0_" + timeStamp,
    "y8_FVA_" + fileSuffixName + "_t0_" + timeStamp,
    condition,
    experiment,
    0,
    row0,
    row0,
    model,
    lgStream=logStrm
)
lGrowth_wet.append(wetOD0 / wetOD0)
lGrowth_silico.append(biomassFlux0 / biomassFlux0)

gL.toLog(logStrm, "\n>>> Wet time: 2h")
row = dfConsProd.iloc[1]
biomassFlux, wetOD = faL.doExperiment(
    o2Setting,
    "y8_pFBA_" + fileSuffixName + "_t2_" + timeStamp,
    "y8_FBA_" + fileSuffixName + "_t2_" + timeStamp,
    "y8_FVA_" + fileSuffixName + "_t2_" + timeStamp,
    condition,
    experiment,
    2,
    row0,
    row,
    model,
    wetOD0=wetOD0,
    biomassFlux0=biomassFlux0,
    lgStream=logStrm
)
lGrowth_wet.append(wetOD / wetOD0)
lGrowth_silico.append(biomassFlux / biomassFlux0)

gL.toLog(logStrm, "\n>>> Wet time: 4h")
row = dfConsProd.iloc[2]
biomassFlux, wetOD = faL.doExperiment(
    o2Setting,
    "y8_pFBA_" + fileSuffixName + "_t4_" + timeStamp,
    "y8_FBA_" + fileSuffixName + "_t4_" + timeStamp,
    "y8_FVA_" + fileSuffixName + "_t4_" + timeStamp,
    condition,
    experiment,
    4,
    row0,
    row,
    model,
    wetOD0=wetOD0,
    biomassFlux0=biomassFlux0,
    lgStream=logStrm
)
lGrowth_wet.append(wetOD / wetOD0)
lGrowth_silico.append(biomassFlux / biomassFlux0)

gL.toLog(logStrm, "\n>>> Wet time: 6h")
row = dfConsProd.iloc[3]
biomassFlux, wetOD = faL.doExperiment(
    o2Setting,
    "y8_pFBA_" + fileSuffixName + "_t6_" + timeStamp,
    "y8_FBA_" + fileSuffixName + "_t6_" + timeStamp,
    "y8_FVA_" + fileSuffixName + "_t6_" + timeStamp,
    condition,
    experiment,
    6,
    row0,
    row,
    model,
    wetOD0=wetOD0,
    biomassFlux0=biomassFlux0,
    lgStream=logStrm
)
lGrowth_wet.append(wetOD / wetOD0)
lGrowth_silico.append(biomassFlux / biomassFlux0)

gL.toLog(logStrm, "\n>>> Wet time: 8h")
row = dfConsProd.iloc[4]
biomassFlux, wetOD = faL.doExperiment(
    o2Setting,
    "y8_pFBA_" + fileSuffixName + "_t8_" + timeStamp,
    "y8_FBA_" + fileSuffixName + "_t8_" + timeStamp,
    "y8_FVA_" + fileSuffixName + "_t8_" + timeStamp,
    condition,
    experiment,
    8,
    row0,
    row,
    model,
    wetOD0=wetOD0,
    biomassFlux0=biomassFlux0,
    lgStream=logStrm
)
lGrowth_wet.append(wetOD / wetOD0)
lGrowth_silico.append(biomassFlux / biomassFlux0)

fig = plt.figure(figsize=(13.3, 10))
ax = fig.add_subplot(111)

ax.scatter(
    [0, 2, 4, 6, 8], lGrowth_wet, alpha=0.7, marker="o", s=300, label="OD / OD t0"
)

ax.scatter(
    [0, 2, 4, 6, 8],
    lGrowth_silico,
    alpha=0.7,
    marker="x",
    s=300,
    label="Biomass flux / biomass flux t0",
)


plt.xlabel("Time (h)", size=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(fontsize=18)
plt.savefig(os.path.join(FIGUREDIR, "y8_" + fileSuffixName + "_" + timeStamp + ".png"))
plt.close()

logStrm.close()
