from dash import html, callback, Input, Output, dcc, Dash, State, dash_table, ctx, no_update
import pandas as pd
import cobra as cb
import os
import sys
import numpy as np
import dash_ag_grid as dag
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc

MODULEDIR = os.path.split(os.path.dirname(__file__))[0]
sys.path.append(MODULEDIR)
import utils.genericLib as gL
import utils.configParamsLib as cpL
import utils.manipulateModelLib as mmL
import utils.modelReconstructionLib as recL

FILE_DIR = gL.CWDIR
MODELDIR = gL.dDirs["mods"]
RAWDIR = gL.dDirs["raw"]

configFileName = "parameters.toml"
dConfParams = cpL.loadConfigParams(confFileName=configFileName)



paramsFileName = "parametersSetting.toml"

dfParamsDefs = pd.read_csv("parametersDefinition.tsv", sep = "\t")
dParams2Description = dfParamsDefs.set_index('Parameter').to_dict()["Definition"]

####### queste due variabiiqui di seguito sono da sistemare non appenna avrò accesso alla struttura dei file

# params_newrxns = ["Id", "Rxns"]
params_newrxns = ["Rxn", "Equation", "Compartment", "Gpr"]
rowData_newrxns = gL.getRowdataList(params_newrxns)
# columnDefs_newrxns = [
#     {"headerName": "Id", "field": "Id", "checkboxSelection": True, "headerCheckboxSelection": True, "editable": True},
#     {"headerName": "Rxns", "field": "Rxns", "editable": True}]
columnDefs_newrxns = [
        {"headerName": "Rxn", "field": "Rxn", "checkboxSelection": True, "headerCheckboxSelection": True, "editable": True},
        {"headerName": "Equation", "field": "Equation", "editable": True},
        {"headerName": "Compartment", "field": "Compartment", "editable": True},
        {"headerName": "Gpr", "field": "Gpr", "editable": True}]


# params_newExrxns = ["Id", "Rxns"]
params_newExrxns = ["Rxn", "Equation", "Compartment", "Gpr"]
rowData_newExrxns = gL.getRowdataList(params_newExrxns)
columnDefs_newExrxns = [
    {"headerName": "Id", "field": "Id", "checkboxSelection": True, "headerCheckboxSelection": True, "editable": True},
    {"headerName": "Rxns", "field": "Rxns", "editable": True}]

def newModelLayout():
    """

    :return: A Div containing dashboard content of homepage loading.
    """
    layout = html.Div([
        html.H1('Generate the new model'),
        html.Br(),

        html.Div([
        html.H5('Input the name of the new model:'),
        dcc.Input(
            id="input_newModelName",
            type="text",
            placeholder="",
        ),
        html.Br(),
        dcc.Loading(
                    id="loading-modelName",
                    children=[html.Div([html.Div(id="out-newModelName")])],
                    type="circle",
                ),
        ]),
        html.Br(), html.Div([
        html.H5('Load the file including the list of reactions of the new model. Any change here introduced in the uploaded file needs table to be exported.',
            # className="my-label"
            ),
        dcc.Upload(
        id='upload-rxnsNewModel',
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
        dcc.Loading(
                    id="loading-newModelRxns",
                    children=[html.Div([html.Div(id="rxnsNewModelPath")])],
                    type="circle",
                ),
        ]),
        html.Br(),
        html.Div([
        dag.AgGrid(
            id='table-editing-rxnsNewModel',
            rowData=rowData_newrxns,
            columnDefs=columnDefs_newrxns,
            defaultColDef={"sortable": True, "filter": True},
            columnSize="sizeToFit",
            getRowId="params.data.Rxn",
            dashGridOptions={"rowSelection": "multiple"},
        ),

        html.Button('Add empty row', id='btn_add_internalRxns'),
        html.Button("Remove selected rows", id="btn_remove_internalRxns"),
        html.Button('Download table', id='btn-exportTableinternalRxns', n_clicks=0),
        html.Div(id='rxnsNewModel_output'),
        ]),
        html.Br(),html.Br(),
        dbc.Container([html.Button('Generate the model', id='btn-createNewModel', n_clicks=0),
        dcc.Loading(children=html.Div(id="newModelSaving"), type="circle")]),
        html.Div(id='tableProposedExchRxns'), html.Br(), html.Br(),
        html.Br(),html.Br(),
        html.Button('Close the newly generated model', id='btn-closeNewModel', style = dict(display='none')),
        html.Div(id='exchNewModelPath'), html.Br(),
        html.Div(id='referenceMergePage'),
        dcc.Store(id = "intRxnsFN"),
        dcc.Store(id = "exchRxnsFN"),
    ])
    return layout

@callback(
    Output("out-newModelName", "children"),
    Input("input_newModelName", "value")
    )
def setModelName(inputName):
    if os.path.exists(paramsFileName) is True:
        dParams = cpL.loadConfigParams(confFileName=paramsFileName)
    else:
        dParams= {}
    if inputName is not None:
        dParams["newModelName"] = f"{inputName}"
        paramsStream = open(paramsFileName, mode="w")
        for param in dParams:
            gL.toLog(paramsStream, dParams2Description[param])
            gL.toLog(paramsStream, f"{param}"  + ' = "' + f"{dParams[param]}" + '"')
            gL.toLog(paramsStream, "\n")
        paramsStream.close()
        return f'You entered the following name:\t{inputName}'
    else:
        return f'ATTENTION: You have not entered any name.'


# add or delete rows of table
@callback(
    Output("table-editing-rxnsNewModel", "rowData", allow_duplicate=True),
    Input("btn_add_internalRxns", "n_clicks"),
    Input("btn_remove_internalRxns", "n_clicks"),
    State("table-editing-rxnsNewModel", "rowData"),
    State("table-editing-rxnsNewModel", "selectedRows"),
    prevent_initial_call=True,
)
def update_transaction(buttonAdd, buttonRemove, data, selection):
    if ctx.triggered_id == "btn_add_internalRxns":
        new_row = {}
        for p in params_newrxns:
            new_row[p] = [None]
        df_new_row = pd.DataFrame(new_row)
        updated_table = pd.concat([pd.DataFrame(data), df_new_row])
        return updated_table.to_dict("records")
    elif selection is not None:
        if ctx.triggered_id == "btn_remove_internalRxns":
            dfSelectedRows = pd.DataFrame(selection)
            dfComplete = pd.DataFrame(data)
            dfToKeep = pd.concat([dfComplete, dfSelectedRows]).drop_duplicates(keep = False).reset_index(drop = True)
            return dfToKeep.to_dict("records")


@callback(
    Output('rxnsNewModelPath', 'children'),
    Output('table-editing-rxnsNewModel', 'rowData', allow_duplicate=True),
    Output("intRxnsFN", "data"),
    Input('upload-rxnsNewModel', 'contents'),
    State('upload-rxnsNewModel', 'filename'),
    Input("out-newModelName", "children"),
    prevent_initial_call=True,
)
def loadFile_newRxns(contents, fileName, experimentModelName):
    if os.path.exists(paramsFileName) is True:
        dParams = cpL.loadConfigParams(confFileName=paramsFileName)
    else:
        dParams= {}
    # print("\ndParams\n", dParams, "\n")
    if contents is not None:
        fileName = gL.parse_contents(contents, fileName)

        dParams["newModelInternalRxns"] = f"{fileName}"
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
        if experimentModelName != "ATTENTION: You have not entered any name.":
            experimentModel = experimentModelName.split(":\t")[1].strip()
            fileName = experimentModel + "_internalRxns.tsv"
            return fileLoadingMessage, no_update, fileName
    raise PreventUpdate

@callback(
    Output("rxnsNewModel_output", "children"),
    Input('btn-exportTableinternalRxns', 'n_clicks'),
    Input("out-newModelName", "children"),
    Input("intRxnsFN", "data"),
    State("table-editing-rxnsNewModel", "rowData"),
)
def saveFile_newRxns(exportButton, experimentModelName, fileName, data):
    if os.path.exists(paramsFileName) is True:
        dParams = cpL.loadConfigParams(confFileName=paramsFileName)
    else:
        dParams= {}
    if "btn-exportTableinternalRxns" == ctx.triggered_id:
        if experimentModelName != "ATTENTION: You have not entered any name.":
            experimentModel = experimentModelName.split(":\t")[1].strip()
        if fileName == experimentModel + "_internalRxns.tsv":
            dParams["newModelInternalRxns"] = f"{fileName}"
            paramsStream = open(paramsFileName, mode="w")
            for param in dParams:
                gL.toLog(paramsStream, dParams2Description[param])
                gL.toLog(paramsStream, f"{param}"  + ' = "' + f"{dParams[param]}" + '"')
                gL.toLog(paramsStream, "\n")
            paramsStream.close()
        pd.DataFrame(data).to_csv(os.path.join(RAWDIR, fileName), sep = "\t", index = False)
        return no_update
    raise PreventUpdate


lColumns_exch = ["MetId", "MetName", "MetCompartment", "ExchRxnId", "Lb", "Ub"]
columnDefs_exch = [
    {"headerName": "MetId", "field": "MetId", "checkboxSelection": True, "headerCheckboxSelection": True, "editable": True},
    {"headerName": "MetName", "field": "MetName", "editable": True},
    {"headerName": "MetCompartment", "field": "MetCompartment", "editable": True},
    {"headerName": "ExchRxnId", "field": "ExchRxnId", "editable": True},
    {"headerName": "Lb", "field": "Lb", "editable": True},
    {"headerName": "Ub", "field": "Ub", "editable": True}]

active_button_style = {'background-color': '#d8d8d8',
                      'color': 'black',
                      'width': '300px',
                      'height': '40px',
                      'margin-top': '5px',
                      'margin-left': '5px'}


@callback(Output('store-newModelName', 'data'),
            Output('newModelSaving', 'children'),
            Output('tableProposedExchRxns', 'children'),
            Output("btn-closeNewModel", "style"),
            Input('rxnsNewModelPath', 'children'),
            Input("out-newModelName", "children"),
            State('store_keggOrgCode', 'data'),
            Input('btn-createNewModel', 'n_clicks'),
            )
def generateModel_internalRxns(rxnsToAddFileName, experimentModelName, keggOrgCode, createModelButton):
    if os.path.exists(paramsFileName) is True:
        dParams = cpL.loadConfigParams(confFileName=paramsFileName)
    else:
        dParams= {}
    if "btn-createNewModel" == ctx.triggered_id:
        experimentModel = experimentModelName.split(":")[1].strip()
        dfName2KeggId = pd.read_csv(os.path.join(RAWDIR,"name2KeggId.tsv"), sep = "\t")

        dSpecies2Annotation = mmL.getModelSpeciesAnnotation(os.path.join(MODELDIR, dParams["backboneModel"]))

        modelGW = cb.io.read_sbml_model(os.path.join(MODELDIR, dConfParams["curatedBackboneModelName"]))
        dKeggId2ModelMet = mmL.getMet2KeggId(modelGW, dSpecies2Annotation)

        newModel, dKeggId2ModelMet, dfName2KeggId = recL.createNewModel_part1_addInternalRxns(dfName2KeggId, modelGW, rxnsToAddFileName, experimentModel, keggOrgCode, dSpecies2Annotation, dizKeggCompound2ModelCompound = dKeggId2ModelMet)

        dRxns2AnnotationNewModel = {}
        for s in newModel.reactions:
            if s.annotation != {}:
                dRxns2AnnotationNewModel[s.id] = s.annotation

        dfdRxns2AnnotationNewModel = pd.DataFrame(dRxns2AnnotationNewModel.items(), columns=['reaction', 'dictAnnotation'])
        dfdRxns2AnnotationNewModel.to_csv(os.path.join(RAWDIR, "newModelReactionsToAnnotation.tsv"), sep = "\t", index = False)

        dSpecies2AnnotationNewModel = {}
        for s in newModel.metabolites:
            if s.annotation != {}:
                dSpecies2AnnotationNewModel[s.id] = s.annotation

        dfdSpecies2AnnotationNewModel = pd.DataFrame(dSpecies2AnnotationNewModel.items(), columns=['species', 'dictAnnotation'])
        dfdSpecies2AnnotationNewModel.to_csv(os.path.join(RAWDIR, "newModelSpeciesToAnnotation.tsv"), sep = "\t", index = False)

        cb.io.write_sbml_model(newModel,os.path.join(MODELDIR,experimentModel + ".xml"))
        dfName2KeggId.to_csv(os.path.join(RAWDIR,"name2KeggId.tsv"), sep = "\t", index = False)

        newModel, dfProposedExchangeRxns = recL.createNewModel_part2_addExchangeRxns(newModel)
        cb.io.write_sbml_model(newModel,os.path.join(MODELDIR,experimentModel + ".xml"))

        if os.path.exists(paramsFileName) is True:
            dParams = cpL.loadConfigParams(confFileName=paramsFileName)
        else:
            dParams= {}

        dParams["newModelFileName"] = f"{experimentModel}" + ".xml"
        dParams["exchangeNewModel"] = "exchangeReactions_newModel.tsv"
        paramsStream = open(paramsFileName, mode="w")
        for param in dParams:
            gL.toLog(paramsStream, dParams2Description[param])
            gL.toLog(paramsStream, f"{param}"  + ' = "' + f"{dParams[param]}" + '"')
            gL.toLog(paramsStream, "\n")
        paramsStream.close()

        return experimentModel, html.Div([
            html.Br(),
            html.H5('New model has been generated. Revise the proposed exchange reactions and then click on the "Close the newly generated model" button to finalize the model.')]), html.Div([
            dag.AgGrid(
                id='table-putativeExchRxns',
                rowData=dfProposedExchangeRxns.to_dict("records"),
                columnDefs=columnDefs_exch,
                defaultColDef={"sortable": True, "filter": True},
                columnSize="sizeToFit",
                getRowId="params.data.MetId",
                dashGridOptions={"rowSelection": "multiple"},
            ),
            html.Button('Add empty row', id='btn_add_exch'),
            html.Button("Remove selected rows", id="btn_remove_exch"),
            html.Button('Download table', id='btn_exportTable_exch', n_clicks=0),
            html.Div(id='putativeExch_output')
            ]), active_button_style

    raise PreventUpdate

# add or delete rows of table
@callback(
    Output("table-putativeExchRxns", "rowData", allow_duplicate=True),
    Input("btn_add_exch", "n_clicks"),
    Input("btn_remove_exch", "n_clicks"),
    State("table-putativeExchRxns", "rowData"),
    State("table-putativeExchRxns", "selectedRows"),
    prevent_initial_call=True,
)
def update_transaction(buttonAdd, buttonRemove, data, selection):
    if ctx.triggered_id == "btn_add_exch":
        new_row = {}
        for p in lColumns_exch:
            new_row[p] = [None]
        df_new_row = pd.DataFrame(new_row)
        updated_table = pd.concat([pd.DataFrame(data), df_new_row])
        return updated_table.to_dict("records")
    elif selection is not None:
        if ctx.triggered_id == "btn_remove_exch":
            dfSelectedRows = pd.DataFrame(selection)
            dfComplete = pd.DataFrame(data)
            dfToKeep = pd.concat([dfComplete, dfSelectedRows]).drop_duplicates(keep = False).reset_index(drop = True)
            return dfToKeep.to_dict("records")

@callback(
    # Output("table-putativeExchRxns", "exportDataAsCsv"),
    # Output("table-putativeExchRxns", "csvExportParams"),
    Output("putativeExch_output", "children"),
    Input('btn_exportTable_exch', 'n_clicks'),
    State("table-putativeExchRxns", "rowData"),
)
def saveFile_putativeExch(exportButton, data):
    if "btn_exportTable_exch" == ctx.triggered_id:

        pd.DataFrame(data).to_csv(os.path.join(RAWDIR, "exchangeReactions_newModel.tsv"), sep = "\t", index = False)
        return no_update

        # return True,{"fileName": "exchangeReactions_newModel.tsv",
        #         "allColumns": True,
        #         "columnSeparator": "\t",
        #         "suppressCsvExport": True}
    raise PreventUpdate

@callback(Output('exchNewModelPath', 'children'),
        Output('referenceMergePage', 'children'),
            Input("out-newModelName", "children"),
            Input('btn-closeNewModel', 'n_clicks'),
            )
def generateModel(experimentModelName, closeModelButton):
    if os.path.exists(paramsFileName) is True:
        dParams = cpL.loadConfigParams(confFileName=paramsFileName)
    else:
        dParams= {}

    if "btn-closeNewModel" == ctx.triggered_id:
        dfProposedExchangeRxns = pd.read_csv(os.path.join(RAWDIR, dParams["exchangeNewModel"]), sep ="\t")
        newModel = cb.io.read_sbml_model(os.path.join(MODELDIR, dParams["newModelFileName"]))

        newModel = recL.createNewModel_part3_addRevisedExchangeRxns(newModel, dfProposedExchangeRxns)
        cb.io.write_sbml_model(newModel,os.path.join(MODELDIR,dParams["newModelFileName"]))

        dRxns2AnnotationExchNewModel = {}
        for row in dfProposedExchangeRxns.itertuples():
            dRxns2AnnotationExchNewModel["EX_" + row.MetId] = {}

        dfdRxns2AnnotationExchNewModel = pd.DataFrame(dRxns2AnnotationExchNewModel.items(), columns=['reaction', 'dictAnnotation'])
        dfdRxns2AnnotationExchNewModel.to_csv(os.path.join(RAWDIR, "newModelExchangeReactionsToAnnotation.tsv"), sep = "\t", index = False)

        return html.Div([
            html.Br(),
            html.H5(f'New model has been saved in the {MODELDIR} directory. Proceed to the next panel.')]), dcc.Link('NEXT', href='/merge')
    raise PreventUpdate
