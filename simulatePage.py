from dash import html, callback, Input, Output, dcc, Dash, State, ctx
import pandas as pd
import cobra as cb
import os
import sys
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
import ast

MODULEDIR = os.path.split(os.path.dirname(__file__))[0]
sys.path.append(MODULEDIR)
import utils.genericLib as gL
import utils.configParamsLib as cpL
import utils.fluxAnalysisLib as faL

OUTDIR = gL.dDirs["out"]
MODELDIR = gL.dDirs["mods"]

configFileName = "parameters.toml"
dConfParams = cpL.loadConfigParams(confFileName=configFileName)

dfParamsDefs = pd.read_csv("parametersDefinition.tsv", sep = "\t")
dParams2Description = dfParamsDefs.set_index('Parameter').to_dict()["Definition"]

paramsFileName = "parametersSetting.toml"

def defineListModels2Simulate():
    if os.path.exists(paramsFileName) is True:
        dParams = cpL.loadConfigParams(confFileName=paramsFileName)
    else:
        dParams= {}
    # lModels2Simulate = ["backboneModel_postCuration_marine_ynb_0_2h.xml", "backboneModel_postCuration_marine_ynb_2_4h.xml",
    # "backboneModel_postCuration_marine_ynb_4_6h.xml", "backboneModel_postCuration_marine_ynb_6_8h.xml"]
    # dParams["models2Simulate"] = lModels2Simulate

    lModels2Simulate = ast.literal_eval(dParams["listConstrainedModels"])

    df = pd.DataFrame({'model': lModels2Simulate})
    df.to_csv(os.path.join(OUTDIR, "models2Simulate.tsv"), sep = "\t", index = False)

    return html.Div([ html.H5('Choose the model to simulate:'),
    dcc.Dropdown(
        id='dropdown_models2Simulate',
        options=[
            {'label':i, 'value':i} for i in df['model'].unique()
        ],
    ),
    dcc.Loading(
                id="loading-output_chosenModel",
                children=[html.Div([html.Div(id="output_chosenModel")])],
                type="circle",
            ),
    #dcc.Store(id='output_chosenModel'),
    html.Br(), html.Br()])


def simulateLayout():
    """
    :return: A Div containing dashboard content of simulate page loading.
    """
    layout = html.Div([
        html.H1('Simulate the constrained models'),
        html.Br(),

        html.Div([
        html.P(
            children=defineListModels2Simulate(),
            id='out-dropdownMenuModels'
        )]),

        html.Div(id="out-modelRxns"),
        dcc.Store(id='store-OFrxn'),
        html.Br(),
        html.Div(id="out-simulationType"),

        html.Div([
        html.H5('Input a string to label models simulations:'),
        dcc.Input(
            id="input_simsType",
            type="text",
            placeholder="",
        ),
        html.Br(),
        dcc.Loading(
                    id="loading-out-simsType",
                    children=[html.Div([html.Div(id="out-simsType")])],
                    type="circle",
                ),
        ]),
        html.Br(),

        dbc.Container([html.Button('Simulate models', id='btn-simulateModels', n_clicks=0),
        dcc.Loading(children=html.Div(id="out-simulationsMessage"), type="circle")]),

        html.Br(),

        html.Div(id="out-simulationsValue"),
        html.Br(),

        # dcc.Link('NEXT', href='/mapping'),
        html.Div(id="referenceMappingPage"),
        html.Div(id="simulation2Run")

    ])
    return layout


@callback(Output('output_chosenModel', 'children'),
      [Input('dropdown_models2Simulate', 'value')])
def model2Simulate(value):
    if value is not None:
        if os.path.exists(paramsFileName) is True:
            dParams = cpL.loadConfigParams(confFileName=paramsFileName)
        else:
            dParams= {}
        df = pd.read_csv(os.path.join(OUTDIR, "models2Simulate.tsv"), sep = "\t")
        dfChosenModel = df[df['model'] == value]
        model2Simulate = dfChosenModel.iloc[0]['model']
        dParams["model2Simulate"] = f"{model2Simulate}"
        paramsStream = open(paramsFileName, mode="w")
        for param in dParams:
            gL.toLog(paramsStream, dParams2Description[param])
            gL.toLog(paramsStream, f"{param}"  + ' = "' + f"{dParams[param]}" + '"')
            gL.toLog(paramsStream, "\n")
        paramsStream.close()
        return model2Simulate
    raise PreventUpdate

@callback(
    Output('out-modelRxns', 'children'),
    Input('output_chosenModel', 'children'),
)
def setOFreactionName(model2Simulate):
    dId2Name_rxns = {}
    model = cb.io.read_sbml_model(os.path.join(MODELDIR, model2Simulate))
    for rxn in model.reactions:
        dId2Name_rxns[rxn.name] = rxn.id

    return html.Div([html.H5('Choose the reaction to set as objective function:'),
        dcc.Dropdown(
           id='modelRxns-dropdown',
           options=[{'label': name, 'value':dId2Name_rxns[name]} for name in dId2Name_rxns]),
           html.Br(),html.Br()
           ])

@callback(Output('store-OFrxn', 'data'),
            Output("out-simulationType", "children"),
              Input('modelRxns-dropdown', 'value')
              )
def storeOFrxn(rxn):
    if os.path.exists(paramsFileName) is True:
        dParams = cpL.loadConfigParams(confFileName=paramsFileName)
    else:
        dParams= {}

    dParams["rxnOF"] = f"{rxn}"
    paramsStream = open(paramsFileName, mode="w")
    for param in dParams:
        gL.toLog(paramsStream, dParams2Description[param])
        gL.toLog(paramsStream, f"{param}"  + ' = "' + f"{dParams[param]}" + '"')
        gL.toLog(paramsStream, "\n")
    paramsStream.close()

    return rxn, html.Div([
       html.H5("Which type of simulation you want to run?"),
       dcc.RadioItems(
           options=['Flux Balance Analysis', 'parsimonious Flux Balance Analysis', 'Flux Variability Analysis'],
           value='Flux Balance Analysis',
           id='simulation2Run'),
           html.Br(),html.Br()
           ])


@callback(
    Output("out-simsType", "children"),
    Input("input_simsType", "value")
    )
def setSimsLabel(inputName):
    if inputName is not None:
        return f'You entered the following label:\t{inputName}'
    else:
        return f'ATTENTION: You have not entered any label.'

@callback(
    Output('out-simulationsValue', 'children'),
    Output('out-simulationsMessage', 'children'),
    Output('referenceMappingPage', 'children'),
    Input('output_chosenModel', 'children'),
    Input('simulation2Run', 'value'),
    State('store-OFrxn', 'data'),
    Input("out-simsType", "children"),
    Input('btn-simulateModels', 'n_clicks'),
)
def simulateConstrainedModels(modelFileName, simulationType, ofRxn, simLabel,simulateButton):
    if modelFileName is not None:
        model = cb.io.read_sbml_model(os.path.join(MODELDIR, modelFileName))
        # model.solver = "gurobi"
        model.solver = "glpk"
        if ofRxn is not None:
            model.reactions.get_by_id(ofRxn).lower_bound = 0
            model.reactions.get_by_id(ofRxn).upper_bound = 1000
            model.reactions.get_by_id(ofRxn).objective_coefficient = 1


            if os.path.exists(paramsFileName) is True:
                dParams = cpL.loadConfigParams(confFileName=paramsFileName)
            else:
                dParams= {}

            if simLabel == "ATTENTION: You have not entered any label.":
                simLabel = "_"
            else:
                simLabel = "_" + simLabel.split(":\t")[1].strip() + "_"

            if "btn-simulateModels" == ctx.triggered_id:
                timeStamp = gL.getTimeStamp()

                if simulationType == "Flux Balance Analysis":
                    try:
                        outputFBAname = "FBA" + simLabel + os.path.splitext(modelFileName)[0] + "_" + timeStamp
                        statusFBA = faL.computeFBAandSaveFlux(model,'FBA', outputFBAname)
                        fba = pd.read_csv(os.path.join(OUTDIR, outputFBAname + ".tsv"), sep = "\t")
                        ofFlux = round(fba[fba["Rxn"] == ofRxn]["Flux"].values[0], 3)
                        message_OfFlux = f"{modelFileName}. Objective function flux value:  {ofFlux} \n"
                        dParams["simsFileName"] = f"{outputFBAname}"+ ".tsv"
                    except:
                        message_OfFlux = f"{modelFileName}. Objective function flux value: Infeasible Solution\n"

                elif simulationType == "parsimonious Flux Balance Analysis":
                    try:
                        outputpFBAname = "pFBA" + simLabel + os.path.splitext(modelFileName)[0] + "_" + timeStamp
                        statuspFBA = faL.computeFBAandSaveFlux(model,'pFBA', outputpFBAname)
                        fba = pd.read_csv(os.path.join(OUTDIR, outputpFBAname + ".tsv"), sep = "\t")
                        ofFlux = round(fba[fba["Rxn"] == ofRxn]["Flux"].values[0], 3)
                        message_OfFlux = f"{modelFileName}. Objective function flux value: {ofFlux} \n"
                        dParams["simsFileName"] = f"{outputpFBAname}"+ ".tsv"
                    except:
                        message_OfFlux= f"{modelFileName}. Objective function flux value: Infeasible Solution\n"

                elif simulationType == "Flux Variability Analysis":
                    try:
                        outputFileNameFVA = "FVA" + simLabel + os.path.splitext(modelFileName)[0] + "_" + timeStamp
                        faL.computeFVAandSaveFlux(model, outputFileNameFVA)
                        fva = pd.read_csv(os.path.join(OUTDIR, outputFileNameFVA + ".tsv"), sep = "\t")
                        ofFluxMin = round(fva[fva["Rxn"] == ofRxn]["Min"].values[0], 3)
                        message_OfFlux = f"{modelFileName}. Objective function flux value:  {ofFluxMin}\n"
                        dParams["simsFileName"] = f"{outputFileNameFVA}" + ".tsv"
                    except:
                        message_OfFlux = f"{modelFileName}.Objective function flux value: Infeasible Solution\n"

                paramsStream = open(paramsFileName, mode="w")
                for param in dParams:
                    gL.toLog(paramsStream, dParams2Description[param])
                    gL.toLog(paramsStream, f"{param}"  + ' = "' + f"{dParams[param]}" + '"')
                    gL.toLog(paramsStream, "\n")
                paramsStream.close()

                return html.H6(message_OfFlux), html.H5(f'Simulation output has been saved in the {OUTDIR} directory. Proceed to the next panel to map fluxes.'), dcc.Link('NEXT', href='/mapping')
    raise PreventUpdate
