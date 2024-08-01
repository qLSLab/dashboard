import dash
from dash import html, callback, Input, Output, dcc, Dash, State, dash_table, ctx, no_update, Patch
import pandas as pd
import cobra as cb
import re
import os
import sys
import numpy as np
from dash.exceptions import PreventUpdate
import dash_ag_grid as dag
import dash_bootstrap_components as dbc
import ast


MODULEDIR = os.path.split(os.path.dirname(__file__))[0]
sys.path.append(MODULEDIR)
import utils.genericLib as gL
import utils.configParamsLib as cpL
import utils.experimentalSetupLib as expSL
import utils.manipulateModelLib as mmL

MODELDIR = gL.dDirs["mods"]
RAWDIR = gL.dDirs["raw"]


configFileName = "parameters.toml"
dConfParams = cpL.loadConfigParams(confFileName=configFileName)

paramsFileName = "parametersSetting.toml"

dfParamsDefs = pd.read_csv("parametersDefinition.tsv", sep = "\t")
dParams2Description = dfParamsDefs.set_index('Parameter').to_dict()["Definition"]


params_medium = ["Met", "Concentration"]
rowData_medium = gL.getRowdataList(params_medium)
columnDefs_medium = [
    {"headerName": "Met", "field": "Met", "checkboxSelection": True, "headerCheckboxSelection": True, "editable": True},
    {"headerName": "Concentration", "field": "Concentration", "editable": True}]

params_conditionMets = ["Met", "RxnID", "Concentration"]
rowData_conditionMets = gL.getRowdataList(params_conditionMets)
columnDefs_conditionMets = [
    {"headerName": "Met", "field": "Met", "checkboxSelection": True, "headerCheckboxSelection": True, "editable": True},
    {"headerName": "RxnID", "field": "RxnID", "editable": True},
    {"headerName": "Concentration", "field": "Concentration", "editable": True}]

params_extracFluxes = ["Time"]
rowData_extracFluxes = gL.getRowdataList(params_extracFluxes)
columnDefs_extracFluxes= [
    {"headerName": "Time", "field": "Time", "checkboxSelection": True, "headerCheckboxSelection": True, "editable": True}]

def constrainLayout():
    """
    :return: A Div containing dashboard content of constrain page loading.
    """
    layout = html.Div([
        html.H1('Constrain the integrated model with experimental data'),
        html.Br(),

        html.Div([
        html.H5('Input a label to define your condition:'),
        dcc.Input(
            id="input_condition",
            type="text",
            placeholder="",
        ),
        ]),

        html.Div(id="out_condition"),

        html.Br(), html.Br(),

        html.H5('Load the file including the wet medium composition:',
            # className="my-label"
            ),

        dcc.Upload(
        id='upload-medium',
        children=html.Div(['Drag and Drop or ', html.A('Select Files')]),
        #html.Button('Upload File')
        style={
        'width': '50%',
        'height': '60px',
        'lineHeight': '60px',
        'borderWidth': '1px',
        'borderStyle': 'dashed',
        'borderRadius': '5px',
        'textAlign': 'center',
        'margin': '10px'
        },
        # Allow multiple files to be uploaded
        multiple=True
        ),
        dcc.Loading(
                    id="loading-medium",
                    children=[html.Div([html.Div(id="out-mediumFile")])],
                    type="circle",
                ),
        html.Br(),

        html.Div([
        dag.AgGrid(
            id='table-editing-mediumFile',
            rowData=rowData_medium,
            columnDefs=columnDefs_medium,
            defaultColDef={"sortable": True, "filter": True},
            columnSize="sizeToFit",
            getRowId="params.data.Met",
            dashGridOptions={"rowSelection": "multiple"},
        ),

        html.Button('Add empty row', id='btn_add_medium'),
        html.Button("Remove selected rows", id="btn_remove_medium"),
        html.Button('Download table', id='btn-exportTablemedium', n_clicks=0),
        html.Div(id="mediumFile_output"),
        ]),

        html.Br(),
        html.Br(),

        html.H5('Load the file including the extracellular production and consumption rates:'),

        dcc.Upload(
        id='upload-extracFluxes',
        children=html.Div(['Drag and Drop or ', html.A('Select Files')]),
        #html.Button('Upload File')
        style={
        'width': '50%',
        'height': '60px',
        'lineHeight': '60px',
        'borderWidth': '1px',
        'borderStyle': 'dashed',
        'borderRadius': '5px',
        'textAlign': 'center',
        'margin': '10px'
        },
        multiple=True
        ),
        dcc.Loading(
                    id="loading-extracFluxes",
                    children=[html.Div([html.Div(id="out-extracFluxes")])],
                    type="circle",
                ),
        html.Br(),

        html.Div([
        dag.AgGrid(
            id='table-editing-extracFluxes',
            rowData=rowData_extracFluxes,
            columnDefs=columnDefs_extracFluxes,
            defaultColDef={"editable": True,"sortable": True, "filter": True},
            columnSize="sizeToFit",
            getRowId="params.data.Time",
            dashGridOptions={"rowSelection": "multiple"},
        ),

        html.Button('Add empty row', id='btn_add_extracFluxes'),
        html.Button("Remove selected rows", id="btn_remove_extracFluxes"),
        html.Button('Download table', id='btn-exportTableextracFluxes', n_clicks=0),

        html.Div([
            dcc.Input(
                id='editing-columns-name',
                placeholder='Enter a column name...',
                value='',
                style={'padding': 7}
            ),
            html.Button('Enter replica column name', id='editing-columns-button', n_clicks=0)
        ], style={'height': 50}),
        html.Div(id="extracFluxes_output"),
        ]),
        html.Br(), html.Br(),

        dbc.Container([html.Button('Constrain model', id='btn-constrainModel', n_clicks=0),
        dcc.Loading(children=html.Div(id="out-finalMessage"), type="circle")]),

        html.Div(id="referenceSimulatePage"),
        # dcc.Link('NEXT', href='/simulate')
        dcc.Store(id = "mediumFN", storage_type='session'),
        dcc.Store(id = "condMetsFN", storage_type='session'),
        dcc.Store(id = "extrFluxFN", storage_type='session'),
    ])
    return layout


# add or delete rows of table
@callback(
    Output("table-editing-mediumFile", "rowData", allow_duplicate=True),
    Input("btn_add_medium", "n_clicks"),
    Input("btn_remove_medium", "n_clicks"),
    State("table-editing-mediumFile", "rowData"),
    State("table-editing-mediumFile", "selectedRows"),
    prevent_initial_call=True,
)
def update_transaction(buttonAdd, buttonRemove, data, selection):
    if ctx.triggered_id == "btn_add_medium":
        new_row = {}
        for p in params_medium:
            new_row[p] = [None]
        df_new_row = pd.DataFrame(new_row)
        updated_table = pd.concat([pd.DataFrame(data), df_new_row])
        return updated_table.to_dict("records")
    elif selection is not None:
        if ctx.triggered_id == "btn_remove_medium":
            dfSelectedRows = pd.DataFrame(selection)
            dfComplete = pd.DataFrame(data)
            dfToKeep = pd.concat([dfComplete, dfSelectedRows]).drop_duplicates(keep = False).reset_index(drop = True)
            return dfToKeep.to_dict("records")

# add columns of table
@callback(
    Output("table-editing-extracFluxes", "columnDefs", allow_duplicate=True),
    State("editing-columns-name", "value"),
    Input("editing-columns-button", "n_clicks"),
    prevent_initial_call=True,
)
def addColumns(colName, buttonAdd):
    if ctx.triggered_id == "editing-columns-button":
        patched_grid = Patch()
        patched_grid.append({"headerName": colName, "field": colName, "editable": True, "checkboxSelection": True,})
        return patched_grid


# add or delete rows of table
@callback(
    Output("table-editing-extracFluxes", "rowData", allow_duplicate=True),
    Input("btn_add_extracFluxes", "n_clicks"),
    Input("btn_remove_extracFluxes", "n_clicks"),
    State("table-editing-extracFluxes", "rowData"),
    State("table-editing-extracFluxes", "selectedRows"),
    prevent_initial_call=True,
)
def update_transaction(buttonAdd, buttonRemove, data, selection):
    if ctx.triggered_id == "btn_add_extracFluxes":
        new_row = {}
        for p in params_extracFluxes:
            new_row[p] = [None]
        df_new_row = pd.DataFrame(new_row)
        updated_table = pd.concat([pd.DataFrame(data), df_new_row])
        return updated_table.to_dict("records")
    elif selection is not None:
        if ctx.triggered_id == "btn_remove_extracFluxes":
            dfSelectedRows = pd.DataFrame(selection)
            dfComplete = pd.DataFrame(data)
            dfToKeep = pd.concat([dfComplete, dfSelectedRows]).drop_duplicates(keep = False).reset_index(drop = True)
            return dfToKeep.to_dict("records")

@callback(
    Output('out-mediumFile', 'children'),
    Output('table-editing-mediumFile', 'rowData', allow_duplicate=True),
    Output("mediumFN", "data"),
    Input('upload-medium', 'contents'),
    State('upload-medium', 'filename'),
    prevent_initial_call=True,
)
def loadFile_newRxns(contents, fileName):
    if os.path.exists(paramsFileName) is True:
        dParams = cpL.loadConfigParams(confFileName=paramsFileName)
    else:
        dParams= {}
    if contents is not None:
        fileName = gL.parse_contents(contents, fileName)

        dParams["mediumFileName"] = f"{fileName}"
        paramsStream = open(paramsFileName, mode="w")
        for param in dParams:
            gL.toLog(paramsStream, dParams2Description[param])
            gL.toLog(paramsStream, f"{param}"  + ' = "' + f"{dParams[param]}" + '"')
            gL.toLog(paramsStream, "\n")
        paramsStream.close()
        df = pd.read_csv(os.path.join(RAWDIR, fileName), sep = "\t")

        fileLoadingMessage = f'You selected the following file:\t{fileName}'
        return fileLoadingMessage, df.to_dict("records"), fileName

    elif contents is None:
        fileLoadingMessage = 'ATTENTION: You have not selected any file. Please fill and then export the table below.'
        fileName = "wetMediumComposition.tsv"
        return fileLoadingMessage, no_update, fileName

@callback(
    Output("mediumFile_output", "children"),
    Input('btn-exportTablemedium', 'n_clicks'),
    Input("mediumFN", "data"),
    State("table-editing-mediumFile", "rowData"),
)
def saveFile_newRxns(exportButton, fileName, data):
    if os.path.exists(paramsFileName) is True:
        dParams = cpL.loadConfigParams(confFileName=paramsFileName)
    else:
        dParams= {}

    if "btn-exportTablemedium" == ctx.triggered_id:
        if fileName == "wetMediumComposition.tsv":
            dParams["conditionMetsFileName"] = f"{fileName}"
            paramsStream = open(paramsFileName, mode="w")
            for param in dParams:
                gL.toLog(paramsStream, dParams2Description[param])
                gL.toLog(paramsStream, f"{param}"  + ' = "' + f"{dParams[param]}" + '"')
                gL.toLog(paramsStream, "\n")
            paramsStream.close()
        pd.DataFrame(data).to_csv(os.path.join(RAWDIR, fileName), sep = "\t", index = False)
        return no_update
        # return True, {"fileName": fileName,
        #     "allColumns": True,
        #     "columnSeparator": "\t",
        #     "suppressCsvExport": True}

    raise PreventUpdate

@callback(
    Output('out-extracFluxes', 'children'),
    Output('table-editing-extracFluxes', 'columnDefs', allow_duplicate=True),
    Output('table-editing-extracFluxes', 'rowData', allow_duplicate=True),
    Output("extrFluxFN", "data"),
    Input('upload-extracFluxes', 'contents'),
    State('upload-extracFluxes', 'filename'),
    prevent_initial_call=True,
)
def loadFile_newRxns(contents, fileName):
    if os.path.exists(paramsFileName) is True:
        dParams = cpL.loadConfigParams(confFileName=paramsFileName)
    else:
        dParams= {}
    if contents is not None:
        fileName = gL.parse_contents(contents, fileName)

        dParams["extracFluxesFileName"] = f"{fileName}"
        paramsStream = open(paramsFileName, mode="w")
        for param in dParams:
            gL.toLog(paramsStream, dParams2Description[param])
            gL.toLog(paramsStream, f"{param}"  + ' = "' + f"{dParams[param]}" + '"')
            gL.toLog(paramsStream, "\n")
        paramsStream.close()
        df = pd.read_csv(os.path.join(RAWDIR, fileName), sep = "\t")

        fileLoadingMessage = f'You selected the following file:\t{fileName}'
        newCols = [{"headerName": i, "field": i} for i in list(df.columns)]
        return fileLoadingMessage, newCols, df.to_dict("records"), fileName

    elif contents is None:
        fileLoadingMessage = 'ATTENTION: You have not selected any file. Please fill and then export the table below.'
        fileName = "extracConsProd.tsv"
        return fileLoadingMessage, no_update, no_update, fileName

@callback(
    # Output("table-editing-extracFluxes", "exportDataAsCsv"),
    # Output("table-editing-extracFluxes", "csvExportParams"),
    Output("extracFluxes_output", "children"),
    Input('btn-exportTableextracFluxes', 'n_clicks'),
    Input("extrFluxFN", "data"),
    State("table-editing-extracFluxes", "rowData"),
)
def saveFile_newRxns(exportButton, fileName, data):
    if os.path.exists(paramsFileName) is True:
        dParams = cpL.loadConfigParams(confFileName=paramsFileName)
    else:
        dParams= {}
    if "btn-exportTableextracFluxes" == ctx.triggered_id:
        if fileName == "extracConsProd.tsv":
            dParams["extracFluxesFileName"] = f"{fileName}"
            paramsStream = open(paramsFileName, mode="w")
            for param in dParams:
                gL.toLog(paramsStream, dParams2Description[param])
                gL.toLog(paramsStream, f"{param}"  + ' = "' + f"{dParams[param]}" + '"')
                gL.toLog(paramsStream, "\n")
            paramsStream.close()

        pd.DataFrame(data).to_csv(os.path.join(RAWDIR, fileName), sep = "\t", index = False)
        return no_update
        # return True,{"fileName": fileName,
        #         "allColumns": True,
        #         "columnSeparator": "\t",
        #         "suppressCsvExport": True}

    raise PreventUpdate


@callback(
    Output("out_condition", "children"),
    Input('input_condition', 'value')
)
def setConditionLabel(label):
    if label is not None:
        return f'You entered the following label:\t{label}'
    else:
        return f'ATTENTION: You have not entered any label.'


@callback(
    Output('store-constrainedModels', 'data'),
    Output('out-finalMessage', 'children'),
    Output('referenceSimulatePage', 'children'),
    State('store-mergedModels', 'data'),
    Input('out-mediumFile', 'children'),
    Input('out-extracFluxes', 'children'),
    Input('btn-constrainModel', 'n_clicks'),
    Input("out_condition", "children"),
)
def constrainModel(mergedModelName, mediumFileName, extracFluxesFileName, constrainButton, modelLabel):
    lConstrainedModels = [mergedModelName + ".xml"]
    if os.path.exists(paramsFileName) is True:
        dParams = cpL.loadConfigParams(confFileName=paramsFileName)
    else:
        dParams= {}
    if "btn-constrainModel" == ctx.triggered_id:
        if modelLabel != "ATTENTION: You have not entered any label.":
            updatedModelMedium = modelLabel.split(":\t")[1] + ""
        else:
            updatedModelMedium = ""

        mergedModel = cb.io.read_sbml_model(os.path.join(MODELDIR, mergedModelName + ".xml"))
        dSpecies2Annotation = mmL.getModelSpeciesAnnotation(os.path.join(MODELDIR, dParams["backboneModel"]))

        dfdSpecies2AnnotationNewModel = pd.read_csv(os.path.join(RAWDIR, "newModelSpeciesToAnnotation.tsv"), sep = "\t")
        dfdSpecies2AnnotationNewModel['dictAnnotation'] = dfdSpecies2AnnotationNewModel['dictAnnotation'].apply(ast.literal_eval)
        dSpecies2AnnotationNewModel = dfdSpecies2AnnotationNewModel.set_index('species').to_dict()
        dSpecies2Annotation.update(dSpecies2AnnotationNewModel["dictAnnotation"])

        dKeggId2ModelMet = mmL.getMet2KeggId(mergedModel, dSpecies2Annotation)

        dfName2KeggId = pd.read_csv(os.path.join(RAWDIR,"name2KeggId.tsv"), sep = "\t")

        if mediumFileName is None:
            dfMedium = pd.DataFrame()
        elif mediumFileName != "ATTENTION: You have not selected any file. Please fill and then export the table below.":
            mediumFile = mediumFileName.split(":\t")[1].strip()
            if os.path.exists(os.path.join(RAWDIR,mediumFile)) is True:
                dfMedium = pd.read_csv(
                          os.path.join(RAWDIR,mediumFile),
                          sep="\t",
                          dtype={"Concentration": float})

        elif os.path.exists(os.path.join(RAWDIR,"wetMediumComposition.tsv")) is True:
            dfMedium = pd.read_csv(
                      os.path.join(RAWDIR,"wetMediumComposition.tsv"),
                      sep="\t",
                      dtype={"Concentration": float})


        if dfMedium.empty == False:
            lActiveMet = ["C00080", "C00001", "C00007", "C00009"]
            dfName2KeggId, mergedModel = expSL.addMediumData2Model(mergedModel, dfMedium, lActiveMet, dfName2KeggId, dKeggId2ModelMet)
            # dfName2KeggId, mergedModel = expSL.addMediumData2Model(mergedModel, dfMedium, dfConditionMets, lActiveMet, dfName2KeggId, dKeggId2ModelMet)
            if updatedModelMedium != "":
                updatedModelMedium = updatedModelMedium + "_" + os.path.splitext(mediumFile)[0]
            else:
                updatedModelMedium = os.path.splitext(mediumFile)[0] + ""

            cb.io.write_sbml_model(mergedModel, os.path.join(MODELDIR,  updatedModelMedium + ".xml"))
            lConstrainedModels = [updatedModelMedium + ".xml"]

        if extracFluxesFileName is None:
            dfConsProd = pd.DataFrame()
        elif extracFluxesFileName != "ATTENTION: You have not selected any file. Please fill and then export the table below.":

            extracFluxesFile = extracFluxesFileName.split(":\t")[1].strip()
            if os.path.exists(os.path.join(RAWDIR,extracFluxesFile)) is True:
                dfConsProd = pd.read_csv(
                          os.path.join(RAWDIR,extracFluxesFile),
                          sep="\t",
                          dtype= {"Time": str})


        elif os.path.exists(os.path.join(RAWDIR,"extracConsProd.tsv")) is True:
            dfConsProd = pd.read_csv(
                      os.path.join(RAWDIR,"extracConsProd.tsv"),
                      sep="\t",
                      dtype= {"Time": str})

        if dfConsProd.empty == False:
            dfConsProd.set_index("Time", inplace = True)
            ## define type columns
            dMeasuredMetsCols = {}
            for col in dfConsProd.columns:
                dMeasuredMetsCols[col] = float
            dfConsProd = dfConsProd.astype(dMeasuredMetsCols)
            lConstrainedModels = expSL.setExtracFluxes(dfConsProd, list(dMeasuredMetsCols.keys()), mergedModel, dKeggId2ModelMet, dfName2KeggId, updatedModelMedium)

        dfName2KeggId.to_csv(os.path.join(RAWDIR,"name2KeggId.tsv"), sep = "\t", index = False)

        if os.path.exists(paramsFileName) is True:
            dParams = cpL.loadConfigParams(confFileName=paramsFileName)
        else:
            dParams= {}
        dParams["listConstrainedModels"] = lConstrainedModels
        paramsStream = open(paramsFileName, mode="w")
        for param in dParams:
            gL.toLog(paramsStream, dParams2Description[param])
            gL.toLog(paramsStream, f"{param}"  + ' = "' + f"{dParams[param]}" + '"')
            gL.toLog(paramsStream, "\n")
        paramsStream.close()

        return lConstrainedModels, html.Div([html.Br(), html.H5(f'The integrated model has been constrained with the provided experimental data. The resulting models have been saved in the {MODELDIR} directory. Proceed to the next panel.')]), dcc.Link('NEXT', href='/simulate')
    raise PreventUpdate


def assignExch2NotUniqueMetMedium(mets):
    return html.Div([
       html.H5("Choose which metabolite needs to be associated to an exchange reaction:"),
       dcc.RadioItems(
           options=mets,
           value=mets[0],
           id=str(mets[0])
           )
       ]),
