import dash
from dash import html, callback, Input, Output, dcc, Dash, State, ctx
import pandas as pd
import cobra as cb
import os
import sys
import ast
from dash.exceptions import PreventUpdate
MODULEDIR = os.path.split(os.path.dirname(__file__))[0]
sys.path.append(MODULEDIR)
import utils.genericLib as gL
import utils.configParamsLib as cpL
import utils.manipulateModelLib as mmL
import dash_bootstrap_components as dbc

MODELDIR = gL.dDirs["mods"]
RAWDIR = gL.dDirs["raw"]

configFileName = "parameters.toml"
dConfParams = cpL.loadConfigParams(confFileName=configFileName)

paramsFileName = "parametersSetting.toml"

dfParamsDefs = pd.read_csv("parametersDefinition.tsv", sep = "\t")
dParams2Description = dfParamsDefs.set_index('Parameter').to_dict()["Definition"]

def mergeLayout():
    """
    :return: A Div containing dashboard content of merge page loading.
    """
    layout = html.Div([
        html.H1('Create the integrated model'),
        html.Br(),

        html.H5('Click on the following button to merge the backbone and the new model'),

        dbc.Container([html.Button('Merge models', id='btn-mergeModels', n_clicks=0),
        dcc.Loading(children=html.Div(id="mergeModelSaving"), type="circle")]),

        # html.Div(id='mergeModelSaving'),
        html.Div(id='referenceConstrainPage'),

    ])
    return layout

@callback(Output('store-mergedModels', 'data'),
            Output('mergeModelSaving', 'children'),
            State('store-newModelName', 'data'),
            Input('btn-mergeModels', 'n_clicks'),
            )
def mergeModels(newModelName, mergeModelsButton):
    if os.path.exists(paramsFileName) is True:
        dParams = cpL.loadConfigParams(confFileName=paramsFileName)
    else:
        dParams= {}

    modelBackbone = cb.io.read_sbml_model(os.path.join(MODELDIR, dParams["finalCuratedModelName"]))

    newModel = cb.io.read_sbml_model(os.path.join(MODELDIR,newModelName + ".xml"))
    modelBackbone_plus_newModel = modelBackbone.merge(newModel)

    dSpecies2Annotation = mmL.getModelSpeciesAnnotation(os.path.join(MODELDIR, dParams["backboneModel"]))

    dfdSpecies2AnnotationNewModel = pd.read_csv(os.path.join(RAWDIR, "newModelSpeciesToAnnotation.tsv"), sep = "\t")
    dfdSpecies2AnnotationNewModel['dictAnnotation'] = dfdSpecies2AnnotationNewModel['dictAnnotation'].apply(ast.literal_eval)
    dSpecies2AnnotationNewModel = dfdSpecies2AnnotationNewModel.set_index('species').to_dict()
    dSpecies2Annotation.update(dSpecies2AnnotationNewModel["dictAnnotation"])

    dKeggId2ModelMet = mmL.getMet2KeggId(modelBackbone_plus_newModel, dSpecies2Annotation)
    finalMergedModelName = dConfParams["curatedBackboneModelName"].split(".")[0] + "_" + newModelName

    if "btn-mergeModels" == ctx.triggered_id:
        cb.io.write_sbml_model(modelBackbone_plus_newModel, os.path.join(MODELDIR, finalMergedModelName + ".xml"))

        dParams["mergedModelFileName"] = f"{finalMergedModelName}" + '.xml'
        paramsStream = open(paramsFileName, mode="w")
        for param in dParams:
            gL.toLog(paramsStream, dParams2Description[param])
            gL.toLog(paramsStream, f"{param}"  + ' = "' + f"{dParams[param]}" + '"')
            gL.toLog(paramsStream, "\n")
        paramsStream.close()

        return finalMergedModelName, html.Div([
            html.Br(), html.H5(f'The backbone and the new models have been merged. The resulting model has been saved in the {MODELDIR} directory. Proceed to the next panel.')])

@callback(Output('referenceConstrainPage', 'children'),
        Input('mergeModelSaving', 'children'),)
def goNextPage(lastCommand):
    if lastCommand is not None:
        if lastCommand["props"]["children"][1]["props"]["children"] == f'The backbone and the new models have been merged. The resulting model has been saved in the {MODELDIR} directory. Proceed to the next panel.':
            return dcc.Link('NEXT', href='/constrain')
