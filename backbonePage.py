from dash import html, callback, Input, Output, dcc, Dash, State, dash_table, ctx, no_update
import pandas as pd
import os
import sys
import numpy as np
from cobra import Model
import cobra as cb
import dash_ag_grid as dag
from dash.exceptions import PreventUpdate
import time
import dash_bootstrap_components as dbc

MODULEDIR = os.path.split(os.path.dirname(__file__))[0]
sys.path.append(MODULEDIR)
import utils.keggLib as kL
import utils.genericLib as gL
import utils.configParamsLib as cpL
import utils.modelReconstructionLib as recL
import utils.manipulateModelLib as mmL

MODELDIR = gL.dDirs["mods"]
RAWDIR = gL.dDirs["raw"]

FILE_DIR = gL.CWDIR
configFileName = "parameters.toml"
dConfParams = cpL.loadConfigParams(confFileName=configFileName)

paramsFileName = "parametersSetting.toml"
dfParamsDefs = pd.read_csv("parametersDefinition.tsv", sep = "\t")
dParams2Description = dfParamsDefs.set_index('Parameter').to_dict()["Definition"]

## Load or create file including correspondence between metabolites name and their KEGG identifier
if os.path.exists(os.path.join(RAWDIR,"name2KeggId.tsv")) is True:
    dfName2KeggId = pd.read_csv(os.path.join(RAWDIR,"name2KeggId.tsv"), sep = "\t")
else:
    dfName2KeggId = pd.DataFrame()

## define all possible kegg codes
sOrgs = kL.getAllKEGGItemsInKEGGdb("organism")
dOrgs = {}
for org in sOrgs.strip().split("\n"):
    orgItemSplt = org.split("\t")
    dOrgs[orgItemSplt[2]] = orgItemSplt[1]

dOrgs["toymodel"] = "toymodel"
white_button_style = {'background-color': 'white',
                      'color': '#f2f2f2',
                      'border-color': '#f2f2f2',
                      'height': '40',
                      'width': '150px',
                      'margin-top': '5px',
                      'margin-left': '5px'}

active_button_style = {'background-color': '#d8d8d8',
                      'color': 'black',
                      'height': '40px',
                      'width': '150',
                      'margin-top': '5px',
                      'margin-left': '5px'}

def backboneModelLayout():
    """
    :return: A Div containing dashboard content of backbone model loading.
    """
    layout = html.Div([
        html.H1('The backbone model'),
        html.Br(),

        html.H5('Choose the organism name:'),
        dcc.Dropdown(
           id='keggCode-dropdown',
           options=[{'label': organismIds, 'value':dOrgs[organismIds]} for organismIds in dOrgs],
           # value='Homo sapiens (human)'
           value='hsa'),
        html.Br(),html.Br(),

        html.Div([
           html.H5("Do you want to upload a backbone model?"),
           dcc.RadioItems(
               options=['Yes', 'No'],
               value='No',
               id='uploadBackboneModelChoice'
           )
       ]),

       html.Div(id='modelLoading'),html.Br(),html.Br(),
       # html.Div([
       # html.H5('Input path of your downloads folder in order to automatically retrieve exported files:'),
       # dcc.Input(
       #     id="input_downloadFolder",
       #     type="text",
       #     placeholder="",
       # ),
       # ]),
       # dcc.Store(id='store-downloadFolder'),
       # html.Button('Move files', id='btn_moveFiles_B', n_clicks=0),
       # html.Div(id='downloadFolderPosition_output_b'), html.Br(),
       dbc.Container([
       html.Button('Proceed', id='btn-proceedButton', n_clicks=0),
       dcc.Loading(children=html.Div(id="backboneModelSavingNo"), type="circle")]),
       dbc.Container([
       html.Button('Save and Proceed', id='btn-saveBackboneModelButton', n_clicks=0,  style = dict(display='none')),
       dcc.Loading(children=html.Div(id="backboneModelSavingYes"), type="circle")]),

       html.Div(id='referenceCreateModelPage'),

       html.Div(id='newRxnsToAddFileName'),
       html.Div(id='modelPath'),
       html.Div(id='boundariesToChangeFileName'),
       dcc.Store(id='gprRefinementFileName'),
       html.Div(id="out-simsType"),
       html.Div(id="curateChoice"),
       html.Div(id="btn-simulateModels"),
       dcc.Store(id = "newRxnsFN"),
       dcc.Store(id = "newBoundsFN"),
       dcc.Store(id = "newGprFN"),

    ])
    return layout

finalCuratedModelName = dConfParams["curatedBackboneModelName"]



@callback(Output('modelLoading', 'children'),
        Output('backboneModelSavingNo', 'children'),
        Output('btn-proceedButton', 'style'),
        Output('btn-saveBackboneModelButton', 'style'),
        # Output('btn-saveBackboneModelButton', 'style', allow_duplicate=True),
        # Output('btn-saveBackboneModelButton', 'style'),
        Input('uploadBackboneModelChoice', 'value'),
        Input('btn-proceedButton', 'n_clicks'),
        )


def choose2LoadModel(loadingChoice, proceedButton):
    if os.path.exists(paramsFileName) is True:
        dParams = cpL.loadConfigParams(confFileName=paramsFileName)
    else:
        dParams= {}


    if loadingChoice == "Yes":
        dParams["backboneModelLoadingChoice"] = f"{loadingChoice}"
        paramsStream = open(paramsFileName, mode="w")
        for param in dParams:
            gL.toLog(paramsStream, dParams2Description[param])
            gL.toLog(paramsStream, f"{param}"  + ' = "' + f"{dParams[param]}" + '"')
            gL.toLog(paramsStream, "\n")
        paramsStream.close()

        return html.Div([html.Br(),html.Br(), html.H5('Load the backbone model:'),
        dcc.Upload(
         id='upload-baseModel',
         children=html.Div(['Drag and Drop or ', html.A('Select Files')]),
         style={
         'width': '50%',
         'height': '60px',
         'lineHeight': '60px',
         'borderWidth': '1px',
         'borderStyle': 'dashed',
         'borderRadius': '5px',
         'textAlign': 'center',
         'margin': '10px'
         }),
         # html.Div(id='modelPath'),
         dcc.Loading(
                     id="loading-backboneModel",
                     children=[html.Div([html.Div(id="modelPath")])],
                     type="circle",
                 ),
         html.Br(),html.Br(),
         html.H5("Do you want to curate the backbone model?"),
         dcc.RadioItems(
                options=['Yes', 'No'],
                value='No',
                id='curateChoice'),

        html.Br(),html.Br(),
        html.Div(id='addCurationFilesSlot')
        ]), "",  dict(display='none'), active_button_style
    else:

        dParams["backboneModelLoadingChoice"] = f"{loadingChoice}"
        paramsStream = open(paramsFileName, mode="w")
        for param in dParams:
            gL.toLog(paramsStream, dParams2Description[param])
            gL.toLog(paramsStream, f"{param}"  + ' = "' + f"{dParams[param]}" + '"')
            gL.toLog(paramsStream, "\n")
        paramsStream.close()

        modelGW = Model()
        if "btn-proceedButton" == ctx.triggered_id:
            cb.io.write_sbml_model(modelGW, os.path.join(MODELDIR, finalCuratedModelName))

            return None, html.Div([
                html.Br(),
                html.H5('Proceed to the next panel.')]), dict(display='none'), dict(display='none')
    raise PreventUpdate



@callback(
    Output('modelPath', 'children'),
    Input('upload-baseModel', 'contents'),
    State('upload-baseModel', 'filename'),
    )
def getPath(list_of_contents, list_of_names):
    if os.path.exists(paramsFileName) is True:
        dParams = cpL.loadConfigParams(confFileName=paramsFileName)
    else:
        dParams= {}


    if list_of_contents is not None:
        modelPath = gL.parse_contents(list_of_contents, list_of_names)

        dParams["backboneModel"] = f"{modelPath}"
        # print("\n\dParams\n", dParams, "\n\n")

        paramsStream = open(paramsFileName, mode="w")
        for param in dParams:
            gL.toLog(paramsStream, dParams2Description[param])
            gL.toLog(paramsStream, f"{param}"  + ' = "' + f"{dParams[param]}" + '"')
            gL.toLog(paramsStream, "\n")
        paramsStream.close()

        fileLoadingMessage = f"The selected backbone model is:\t{modelPath}"
    else:
        fileLoadingMessage = 'ATTENTION: You have not selected any backbone model.'

        dParams["backboneModel"] = "Any backbone model has been selected."
        paramsStream = open(paramsFileName, mode="w")
        for param in dParams:
            gL.toLog(paramsStream, dParams2Description[param])
            gL.toLog(paramsStream, f"{param}"  + ' = "' + f"{dParams[param]}" + '"')
            gL.toLog(paramsStream, "\n")
        paramsStream.close()

    return fileLoadingMessage


@callback(Output('addCurationFilesSlot', 'children'),
              Input('curateChoice', 'value'),
              )
def choose2AddCurationFiles(curateChoice):
    if os.path.exists(paramsFileName) is True:
        dParams = cpL.loadConfigParams(confFileName=paramsFileName)
    else:
        dParams= {}


    if curateChoice == "Yes":

        dParams["curateChoice"] = f"{curateChoice}"
        paramsStream = open(paramsFileName, mode="w")
        for param in dParams:
            gL.toLog(paramsStream, dParams2Description[param])
            gL.toLog(paramsStream, f"{param}"  + ' = "' + f"{dParams[param]}" + '"')
            gL.toLog(paramsStream, "\n")
        paramsStream.close()

        return html.Div([html.Br(),
            html.H5("Do you want to add external files to curate backbone model?"),
            dcc.RadioItems(
                options=['Yes', 'No'],
                value='No',
                id = "addCurationFiles"
            ),
            html.Div(id='newRxnsToAddFile'),
            # html.Div(id='newRxnsToAddFileName'),
            dcc.Loading(
                        id="loading-newRxns",
                        children=[html.Div([html.Div(id="newRxnsToAddFileName")])],
                        type="circle",
                    ),
            html.Br(),
            html.Div(id='newRxnsToAddTable'),
            html.Div(id='exportTableRxnsButton'),
            html.Br(),html.Br(),

            html.Div(id='boundariesToChangeFile'),
            dcc.Loading(
                        id="loading-bounds",
                        children=[html.Div([html.Div(id="boundariesToChangeFileName")])],
                        type="circle",
                    ),
            html.Br(),
            html.Div(id='boundariesToChangeTable'),
            html.Div(id='exportTableboundariesToChangeButton'),
            html.Br(), html.Br(),

            html.Div(id='gprRefinementFile'),
            dcc.Loading(
                        id="loading-gpr",
                        children=[html.Div([html.Div(id="gprRefinementFileName")])],
                        type="circle",
                    ),
            html.Br(),
            html.Div(id='gprRefinementTable'),
            html.Div(id='exportTablegprRefinementButton'),
            html.Br(), html.Br(),

            html.Div(id='export'),

            ])

params_rxns2Add = ["Rxn", "Equation", "Compartment", "Gpr"]
rowData_rxns2Add = gL.getRowdataList(params_rxns2Add)
columnDefs_rxns2Add = [
    {"headerName": "Rxn", "field": "Rxn", "checkboxSelection": True, "headerCheckboxSelection": True, "editable": True},
    {"headerName": "Equation", "field": "Equation", "editable": True},
    {"headerName": "Compartment", "field": "Compartment", "editable": True},
    {"headerName": "Gpr", "field": "Gpr", "editable": True}]

params_boundary = ['Id', 'Lb', 'Ub']
rowData_boundary = gL.getRowdataList(params_boundary)
columnDefs_boundary = [
    {"headerName": "Id", "field": "Id", "checkboxSelection": True, "headerCheckboxSelection": True, "editable": True},
    {"headerName": "Lb", "field": "Lb", "editable": True},
    {"headerName": "Ub", "field": "Ub", "editable": True}]

params_gpr = ['Id', 'Gpr']
rowData_gpr = gL.getRowdataList(params_gpr)
columnDefs_gpr = [
    {"headerName": "Id", "field": "Id", "checkboxSelection": True, "headerCheckboxSelection": True, "editable": True},
    {"headerName": "Gpr", "field": "Gpr", "editable": True}]


@callback([Output('newRxnsToAddFile', 'children'),
        Output('newRxnsToAddTable', 'children'),
        Output('boundariesToChangeFile', 'children'),
        Output('boundariesToChangeTable', 'children'),
        Output('gprRefinementFile', 'children'),
        Output('gprRefinementTable', 'children'),],
        Input('addCurationFiles', 'value')
        )
def choose2AddCurationFiles(addCurationFileChoice):
    if os.path.exists(paramsFileName) is True:
        dParams = cpL.loadConfigParams(confFileName=paramsFileName)
    else:
        dParams= {}


    if addCurationFileChoice == "Yes":

        dParams["addCurationFileChoice"] = f"{addCurationFileChoice}"
        paramsStream = open(paramsFileName, mode="w")
        for param in dParams:
            gL.toLog(paramsStream, dParams2Description[param])
            gL.toLog(paramsStream, f"{param}"  + ' = "' + f"{dParams[param]}" + '"')
            gL.toLog(paramsStream, "\n")
        paramsStream.close()

        return html.Div([
            html.Br(),
            html.H6('Load the file including new reactions to add to the background model or fill and then export the table below. Any change here introduced in the uploaded file needs table to be exported.',
                ),
            dcc.Upload(
            id='upload-newRxnsToAddFileName',
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
            }),
            ]), html.Div([
            dag.AgGrid(
                id='table-editing-rxns',
                rowData=rowData_rxns2Add,
                columnDefs=columnDefs_rxns2Add,
                defaultColDef={"sortable": True, "filter": True},
                columnSize="sizeToFit",
                getRowId="params.data.Rxn",
                dashGridOptions={"rowSelection": "multiple"},
            ),

            html.Button('Add empty row', id='btn_add_newRxns'),
            html.Button("Remove selected rows", id="btn_remove_newRxns"),
            html.Button('Download table', id='btn-exportTableRxns', n_clicks=0),
            html.Div(id='newRxns_output'),
            ]), html.Div([
            html.H6('Load the file including boundaries to change in the background model. Any change here introduced in the uploaded file needs table to be exported.',
                ),
            dcc.Upload(
            id='upload-boundariesToChangeFileName',
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
            }),
            ]), html.Div([
            dag.AgGrid(
                id='table-editing-boundaries',
                rowData=rowData_boundary,
                columnDefs=columnDefs_boundary,
                defaultColDef={"sortable": True, "filter": True},
                columnSize="sizeToFit",
                getRowId="params.data.Id",
                dashGridOptions={"rowSelection": "multiple"},
            ),
            html.Button('Add empty row', id='btn_add_boundary'),
            html.Button("Remove selected rows", id="btn_remove_boundary"),
            html.Button('Download table', id='btn-exportTableboundary', n_clicks=0),
            html.Div(id='boundaries_output'),
            ]), html.Div([
            html.H6('Load the file including GPR to fix in the background model. Any change here introduced in the uploaded file needs table to be exported.',
                ),
            dcc.Upload(
            id='upload-gprRefinementFileName',
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
            }),
            ]), html.Div([
            dag.AgGrid(
                id='table-editing-gpr',
                rowData=rowData_gpr,
                columnDefs=columnDefs_gpr,
                defaultColDef={"sortable": True, "filter": True},
                columnSize="sizeToFit",
                getRowId="params.data.Id",
                dashGridOptions={"rowSelection": "multiple"},
            ),
            html.Button('Add empty row', id='btn_add_gpr'),
            html.Button("Remove selected rows", id="btn_remove_gpr"),
            html.Button('Download table', id='btn-exportTablegpr', n_clicks=0),
            html.Div(id='gpr_output'),
            ])
    raise PreventUpdate


# add or delete rows of table
@callback(
    Output("table-editing-rxns", "rowData", allow_duplicate=True),
    Input("btn_add_newRxns", "n_clicks"),
    Input("btn_remove_newRxns", "n_clicks"),
    State("table-editing-rxns", "rowData"),
    State("table-editing-rxns", "selectedRows"),
    prevent_initial_call=True,
)
def update_transaction(buttonAdd, buttonRemove, data, selection):
    if ctx.triggered_id == "btn_add_newRxns":
        new_row = {}
        for p in params_rxns2Add:
            new_row[p] = [None]
        df_new_row = pd.DataFrame(new_row)
        updated_table = pd.concat([pd.DataFrame(data), df_new_row])
        return updated_table.to_dict("records")
    elif selection is not None:
        if ctx.triggered_id == "btn_remove_newRxns":
            dfSelectedRows = pd.DataFrame(selection)
            dfComplete = pd.DataFrame(data)
            dfToKeep = pd.concat([dfComplete, dfSelectedRows]).drop_duplicates(keep = False).reset_index(drop = True)
            return dfToKeep.to_dict("records")


# add or delete rows of table
@callback(
    Output("table-editing-boundaries", "rowData", allow_duplicate=True),
    Input("btn_add_boundary", "n_clicks"),
    Input("btn_remove_boundary", "n_clicks"),
    State("table-editing-boundaries", "rowData"),
    State("table-editing-boundaries", "selectedRows"),
    prevent_initial_call=True,
)
def update_transaction(buttonAdd, buttonRemove, data, selection):
    if ctx.triggered_id == "btn_add_boundary":
        new_row = {}
        for p in params_boundary:
            new_row[p] = [None]
        df_new_row = pd.DataFrame(new_row)
        updated_table = pd.concat([pd.DataFrame(data), df_new_row])
        return updated_table.to_dict("records")
    elif selection is not None:
        if ctx.triggered_id == "btn_remove_boundary":
            dfSelectedRows = pd.DataFrame(selection)
            dfComplete = pd.DataFrame(data)
            dfToKeep = pd.concat([dfComplete, dfSelectedRows]).drop_duplicates(keep = False).reset_index(drop = True)
            return dfToKeep.to_dict("records")

# add or delete rows of table
@callback(
    Output("table-editing-gpr", "rowData", allow_duplicate=True),
    Input("btn_add_gpr", "n_clicks"),
    Input("btn_remove_gpr", "n_clicks"),
    State("table-editing-gpr", "rowData"),
    State("table-editing-gpr", "selectedRows"),
    prevent_initial_call=True,
)
def update_transaction(buttonAdd, buttonRemove, data, selection):
    if ctx.triggered_id == "btn_add_gpr":
        new_row = {}
        for p in params_gpr:
            new_row[p] = [None]
        df_new_row = pd.DataFrame(new_row)
        updated_table = pd.concat([pd.DataFrame(data), df_new_row])
        return updated_table.to_dict("records")
    elif selection is not None:
        if ctx.triggered_id == "btn_remove_gpr":
            dfSelectedRows = pd.DataFrame(selection)
            dfComplete = pd.DataFrame(data)
            dfToKeep = pd.concat([dfComplete, dfSelectedRows]).drop_duplicates(keep = False).reset_index(drop = True)
            return dfToKeep.to_dict("records")


@callback(
    Output('newRxnsToAddFileName', 'children'),
    Output('table-editing-rxns', 'rowData', allow_duplicate=True),
    Output("newRxnsFN", "data"),
    Input('upload-newRxnsToAddFileName', 'contents'),
    State('upload-newRxnsToAddFileName', 'filename'),
    prevent_initial_call=True,
)
def loadFile_newRxns(contents, fileName):
    if os.path.exists(paramsFileName) is True:
        dParams = cpL.loadConfigParams(confFileName=paramsFileName)
    else:
        dParams= {}


    if contents is not None:
        fileName = gL.parse_contents(contents, fileName)

        dParams["refineModel_rxns"] = f"{fileName}"
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
        fileName = "modelRefinement_newRxns2Add.tsv"
        return fileLoadingMessage, no_update, fileName


@callback(
    Output("newRxns_output", "children"),
    Input('btn-exportTableRxns', 'n_clicks'),
    Input("newRxnsFN", "data"),
    State("table-editing-rxns", "rowData"),
)
def saveFile_newRxns(exportButton, fileName, data):
    if os.path.exists(paramsFileName) is True:
        dParams = cpL.loadConfigParams(confFileName=paramsFileName)
    else:
        dParams= {}
    if "btn-exportTableRxns" == ctx.triggered_id:
        if fileName == "modelRefinement_newRxns2Add.tsv":
            dParams["refineModel_rxns"] = fileName
            paramsStream = open(paramsFileName, mode="w")
            for param in dParams:
                gL.toLog(paramsStream, dParams2Description[param])
                gL.toLog(paramsStream, f"{param}"  + ' = "' + f"{dParams[param]}" + '"')
                gL.toLog(paramsStream, "\n")
            paramsStream.close()

        pd.DataFrame(data).to_csv(os.path.join(RAWDIR, fileName), sep = "\t", index = False)
        return no_update
    raise PreventUpdate

@callback(
    Output('boundariesToChangeFileName', 'children'),
    Output('table-editing-boundaries', 'rowData', allow_duplicate=True),
    Output("newBoundsFN", "data"),
    Input('upload-boundariesToChangeFileName', 'contents'),
    State('upload-boundariesToChangeFileName', 'filename'),
    prevent_initial_call=True,
)
def loadFile_newBounds(contents, fileName):
    if os.path.exists(paramsFileName) is True:
        dParams = cpL.loadConfigParams(confFileName=paramsFileName)
    else:
        dParams= {}
    if contents is not None:
        fileName = gL.parse_contents(contents, fileName)

        dParams["refineModel_bounds"] = f"{fileName}"
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
        fileName = "modelRefinement_bounds2Change.tsv"
        return fileLoadingMessage, no_update, fileName

@callback(
    Output("boundaries_output", "children"),
    Input('btn-exportTableboundary', 'n_clicks'),
    Input("newBoundsFN", "data"),
    State("table-editing-boundaries", "rowData"),
)
def saveFile_newRxns(exportButton, fileName, data):
    if os.path.exists(paramsFileName) is True:
        dParams = cpL.loadConfigParams(confFileName=paramsFileName)
    else:
        dParams= {}
    if "btn-exportTableboundary" == ctx.triggered_id:
        if fileName == "modelRefinement_bounds2Change.tsv":
            dParams["refineModel_bounds"] = f"{fileName}"
            paramsStream = open(paramsFileName, mode="w")
            for param in dParams:
                gL.toLog(paramsStream, dParams2Description[param])
                gL.toLog(paramsStream, f"{param}"  + ' = "' + f"{dParams[param]}" + '"')
                gL.toLog(paramsStream, "\n")
            paramsStream.close()
        pd.DataFrame(data).to_csv(os.path.join(RAWDIR, fileName), sep = "\t", index = False)
        return no_update

    raise PreventUpdate


@callback(
    Output('gprRefinementFileName', 'children'),
    Output('table-editing-gpr', 'rowData', allow_duplicate=True),
    Output("newGprFN", "data"),
    Input('upload-gprRefinementFileName', 'contents'),
    State('upload-gprRefinementFileName', 'filename'),
    prevent_initial_call=True,
)
def loadFile_newGpr(contents, fileName):
    if os.path.exists(paramsFileName) is True:
        dParams = cpL.loadConfigParams(confFileName=paramsFileName)
    else:
        dParams= {}
    if contents is not None:
        fileName = gL.parse_contents(contents, fileName)

        dParams["refineModel_gpr"] = f"{fileName}"
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
        fileName = "modelRefinement_newGpr.tsv"
        return fileLoadingMessage, no_update, fileName

@callback(
    Output("gpr_output", "children"),
    Input('btn-exportTablegpr', 'n_clicks'),
    Input("newGprFN", "data"),
    State("table-editing-gpr", "rowData"),
)
def saveFile_newRxns(exportButton, fileName, data):
    if os.path.exists(paramsFileName) is True:
        dParams = cpL.loadConfigParams(confFileName=paramsFileName)
    else:
        dParams= {}
    if "btn-exportTablegpr" == ctx.triggered_id:
        if fileName == "modelRefinement_newGpr.tsv":
            dParams["refineModel_gpr"] = f"{fileName}"
            paramsStream = open(paramsFileName, mode="w")
            for param in dParams:
                gL.toLog(paramsStream, dParams2Description[param])
                gL.toLog(paramsStream, f"{param}"  + ' = "' + f"{dParams[param]}" + '"')
                gL.toLog(paramsStream, "\n")
            paramsStream.close()
        pd.DataFrame(data).to_csv(os.path.join(RAWDIR, fileName), sep = "\t", index = False)
        return no_update
    raise PreventUpdate


@callback(Output('store_keggOrgCode', 'data'),
              Input('keggCode-dropdown', 'value'),
              )
def storeKeggOrgCode(orgCode):
    if os.path.exists(paramsFileName) is True:
        dParams = cpL.loadConfigParams(confFileName=paramsFileName)
    else:
        dParams= {}

    dParams["orgCode"] = f"{orgCode}"
    paramsStream = open(paramsFileName, mode="w")
    for param in dParams:
        gL.toLog(paramsStream, dParams2Description[param])
        gL.toLog(paramsStream, f"{param}"  + ' = "' + f"{dParams[param]}" + '"')
        gL.toLog(paramsStream, "\n")
    paramsStream.close()

    return orgCode

@callback(Output('backboneModelSavingYes', 'children'),
            Output('btn-saveBackboneModelButton', 'style', allow_duplicate = True),
              Input('modelPath', 'children'),
              Input('keggCode-dropdown', 'value'),
              Input('curateChoice', 'value'),
              Input('newRxnsToAddFileName', 'children'),
              Input('boundariesToChangeFileName', 'children'),
              Input('gprRefinementFileName', 'children'),
              Input('btn-saveBackboneModelButton', 'n_clicks'),
              prevent_initial_call=True,
              )
def loadBackgroundModel(model, orgCode, curateChoice, curateFile_newRxns, curateFile_bounds2Change, curateFile_grp2Refine, saveButton):
    if os.path.exists(paramsFileName) is True:
        dParams = cpL.loadConfigParams(confFileName=paramsFileName)
    else:
        dParams= {}
    finalCuratedModelName = dConfParams["curatedBackboneModelName"]

    dParams["finalCuratedModelName"] = f"{finalCuratedModelName}"
    paramsStream = open(paramsFileName, mode="w")
    for param in dParams:
        gL.toLog(paramsStream, dParams2Description[param])
        gL.toLog(paramsStream, f"{param}"  + ' = "' + f"{dParams[param]}" + '"')
        gL.toLog(paramsStream, "\n")
    paramsStream.close()


    if model is not None and model != 'ATTENTION: You have not selected any backbone model.':
        modelName = model.split(":\t")[1].strip()

        dSpecies2Annotation = mmL.getModelSpeciesAnnotation(os.path.join(MODELDIR, modelName))

        modelGW = cb.io.read_sbml_model(os.path.join(MODELDIR, modelName))
        backModelPath = os.path.split(model.split(":\t")[1])[1]
        message = 'Backbone model has been saved in the ' + f"{MODELDIR}" + ' directory as ' + f"{finalCuratedModelName}" + '. Proceed to the next panel.'

        if curateChoice == "Yes":
            if "btn-saveBackboneModelButton" == ctx.triggered_id:
                modelGW = recL.curateModel(modelGW,orgCode, dfName2KeggId, modelName, curateFile_newRxns, curateFile_bounds2Change, curateFile_grp2Refine, dSpecies2Annotation)
                cb.io.write_sbml_model(modelGW, os.path.join(MODELDIR, finalCuratedModelName))
                time.sleep(1)
                return  html.Div([
                    html.Br(),
                    html.H5(message)]), dict(display='none')
        else:
            if "btn-saveBackboneModelButton" == ctx.triggered_id:
                cb.io.write_sbml_model(modelGW, os.path.join(MODELDIR, finalCuratedModelName))
                time.sleep(1)
                return  html.Div([
                    html.Br(),
                    html.H5(message)]), dict(display='none')
    raise PreventUpdate


@callback(Output('referenceCreateModelPage', 'children'),
        Input('backboneModelSavingYes', 'children'),
        Input('backboneModelSavingNo', 'children'),)
def goNextPage(saveYes, saveNo):
    if saveYes is not None:
        if os.path.exists(paramsFileName) is True:
            dParams = cpL.loadConfigParams(confFileName=paramsFileName)
        else:
            dParams= {}
        finalCuratedModelName = dParams["finalCuratedModelName"]
        message = 'Backbone model has been saved in the ' + f"{MODELDIR}" + ' directory as ' + f"{finalCuratedModelName}" + '. Proceed to the next panel.'

        if saveYes["props"]["children"][1]["props"]["children"] == message:
            return dcc.Link('NEXT', href='/createModel')
    elif saveNo is not None:
        return dcc.Link('NEXT', href='/createModel')
