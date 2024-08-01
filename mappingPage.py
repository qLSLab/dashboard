import dash
from dash import html, callback, Input, Output, dcc, Dash, State, ctx, no_update
import pandas as pd
import cobra as cb
import re
import os
import sys
import matplotlib as mpl
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
import ast

MODULEDIR = os.path.split(os.path.dirname(__file__))[0]
sys.path.append(MODULEDIR)
import utils.keggLib as kL
import utils.genericLib as gL
import utils.configParamsLib as cpL
import utils.fluxAnalysisLib as faL
from Bio.Graphics.KGML_vis import KGMLCanvas
from IPython.display import Image, HTML
import utils.manipulateModelLib as mmL


# Setting the parameters
OUTDIR = gL.dDirs["out"]
MODELDIR = gL.dDirs["mods"]
FIGUREDIR = gL.dDirs["figs"]
UTILSDIR =  gL.dDirs["utils"]
MAPDIR =  gL.dDirs["map"]
SCRIPTDIR = gL.dDirs["scripts"]
RAWDIR = gL.dDirs["raw"]

# setting create model params
configFileName = "parameters.toml"
dConfParams = cpL.loadConfigParams(confFileName=configFileName)

paramsFileName = "parametersSetting.toml"

dfParamsDefs = pd.read_csv("parametersDefinition.tsv", sep = "\t")
dParams2Description = dfParamsDefs.set_index('Parameter').to_dict()["Definition"]

def mappingLayout():
    """
    :return: A Div containing dashboard content of mapping page loading.
    """
    layout = html.Div([
        html.H1('Map the flux distribution'),
        html.Br(),

        html.Div([
             html.H5("Where do you want to map fluxes?"),
             dcc.RadioItems(
                 options=['KEGG database maps', 'Own SVG map file'],
                 value='KEGG database maps',
                 id='mapChoice'
             )
         ]),

        html.Div(id='mapRxnsOnKeggMaps'),
        html.Div(id='mapRxnsOnSvg'),
        html.Div(id='fvaFlux2Map'),

        dbc.Container([html.Button('Map fluxes on SVG map', id='btn-mapFluxesOnSvg', style = dict(display='none')),
        dcc.Loading(children=html.Div(id="out-svgMessage"), type="circle")]),
    ])
    return layout

active_button_style = {'background-color': '#d8d8d8',
                      'color': 'black',
                      'width': '150px',
                      'height': '40px',
                      'margin-top': '5px',
                      'margin-left': '5px'}

@callback(Output('mapRxnsOnKeggMaps', 'children'),
            Output("btn-mapFluxesOnSvg", "style"),
              Input('mapChoice', 'value'),
              Input('store_keggOrgCode', 'data')
              )
def chooseWhere2MapFluxes(whereMapChoice,keggOrgCode):
    if os.path.exists(paramsFileName) is True:
        dParams = cpL.loadConfigParams(confFileName=paramsFileName)
    else:
        dParams= {}

    simFileName = dParams["simsFileName"]
    model2Simulate = dParams["model2Simulate"]

    lVisitedRxns=[] #to avoid ripetute call
    if whereMapChoice == "Own SVG map file":
        ## carica SVG
        return html.Div([html.Br(),
        html.Br(),
        html.H5('Load your own SVG file map:',),
        dcc.Upload(
        id='upload-svg',
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
        html.Br(),

        html.Div(id="out-svg"),
        html.Div(id="out-svgMessage"),

        html.Br(),
        ]), active_button_style

    elif whereMapChoice == "KEGG database maps":
        return html.Div([html.Br(), html.Div([
           html.H5("If multiple reactions identifiers map on the same reaction arrow in the map, you want to map: "),
           dcc.RadioItems(
               options=['Sum of their fluxes', 'Mean of their fluxes', 'Minimum of their fluxes', 'Maximum of their fluxes'],
               value='Sum of their fluxes',
               id='strategy2Map'),
               html.Br(),html.Br()
               ]),
               dbc.Container([html.Button('Map fluxes on KEGG maps', id='btn-mapFluxesOnKeggMaps', n_clicks=0),
               dcc.Loading(children=html.Div(id="out-mapOnKegg"), type="circle")]), html.Br()]), dict(display='none')


@callback(
    Output('out-svgMessage', 'children'),
    Output('out-svg', 'children'),
    Output('btn-mapFluxesOnSvg', 'style', allow_duplicate=True),
    Input('upload-svg', 'contents'),
    State('upload-svg', 'filename'),
    prevent_initial_call=True,
    )
def loadSvgFile(contents, fileName):
    if os.path.exists(paramsFileName) is True:
        dParams = cpL.loadConfigParams(confFileName=paramsFileName)
    else:
        dParams= {}

    simFileName = dParams["simsFileName"]
    if contents is not None:
        filePath = gL.parse_contents(contents, fileName)
        dParams["svgMapFileName"] = f"{filePath}"
        # dParams.update(dParams_initial)
        paramsStream = open(paramsFileName, mode="w")
        for param in dParams:
            gL.toLog(paramsStream, dParams2Description[param])
            gL.toLog(paramsStream, f"{param}"  + ' = "' + f"{dParams[param]}" + '"')
            gL.toLog(paramsStream, "\n")
        paramsStream.close()

        fileLoadingMessage = f'You selected the following SVG file:\t{filePath}'

        if simFileName.startswith("FVA_"):
            buttonMapFluxOnSvg = html.Div([html.Br(),
            html.H5("Since you run a Flux Variability Analysis, you want to map:"),
            dcc.RadioItems(
                options=['Minimum fluxes', 'Maximum fluxes', 'Flux variations as arrows width', 'Minimum fluxes with flux variation as arrows width', 'Maximum fluxes with flux variation as arrows width'],
                value='Minimum fluxes',
                id='fvaFlux2Map'),
                html.Br()])
            return fileLoadingMessage, buttonMapFluxOnSvg, active_button_style
        else:
            buttonMapFluxOnSvg = html.Div([html.Br()])
            return fileLoadingMessage, buttonMapFluxOnSvg, active_button_style
    else:
        fileLoadingMessage = 'ATTENTION: You have not selected any file.'
        buttonMapFluxOnSvg = ""
        return fileLoadingMessage, buttonMapFluxOnSvg, no_update

@callback(
    Output('mapRxnsOnSvg', 'children'),
    Input('btn-mapFluxesOnSvg', 'n_clicks'),
    Input('out-svgMessage', 'children'),
    Input('fvaFlux2Map', 'value'),
)
def mapFluxesOnSvg(mapButton, mapFile, fvaFlux2Map):
    if os.path.exists(paramsFileName) is True:
        dParams = cpL.loadConfigParams(confFileName=paramsFileName)
    else:
        dParams= {}

    simFileName = dParams["simsFileName"]

    if "btn-mapFluxesOnSvg" == ctx.triggered_id:
        fluxFile = os.path.join(OUTDIR, simFileName)
        if simFileName.startswith("FBA_") or simFileName.startswith("pFBA_"):
            fbaFlux = str(2)
            fvaFlux = str(2)
        elif simFileName.startswith("FVA_"):
            if fvaFlux2Map == "Minimum fluxes":
                fbaFlux = str(1)
                fvaFlux = str(1)
            elif fvaFlux2Map == "Maximum fluxes":
                fbaFlux = str(2)
                fvaFlux = str(2)
            elif fvaFlux2Map == "Flux variations as arrows width":
                fbaFlux = str(3)
                fvaFlux = str(3)
            elif fvaFlux2Map == "Minimum fluxes with flux variation as arrows width":
                fbaFlux = str(1)
                fvaFlux = str(3)
            elif fvaFlux2Map == "Maximum fluxes with flux variation as arrows width":
                fbaFlux = str(2)
                fvaFlux = str(3)

        scriptPath = os.path.join(SCRIPTDIR, "colorFluxes_py3.py")

        if mapFile != 'ATTENTION: You have not selected any file.':
            mapPath = mapFile.split(":\t")[1].strip()
            mapColouredPath = os.path.join(MAPDIR, os.path.splitext(simFileName)[0] + ".svg")
            command = "python " + scriptPath + " " +  mapPath + " " + fba_fva_flux + " " + fbaFlux + " " + fvaFlux + " " + mapColouredPath + " " + str('gist_rainbow_r')
            os.system(str(command))

        return html.H5(f'The resulting map has been saved in the {MAPDIR} directory.')


@callback(
    Output('out-mapOnKegg', 'children'),
    Input('btn-mapFluxesOnKeggMaps', 'n_clicks'),
    Input('strategy2Map', 'value'),
    State('store_keggOrgCode', 'data'),
)
def mapFluxesOnKeggMaps(mapButton, flux2Map,keggOrgCode):
    if os.path.exists(paramsFileName) is True:
        dParams = cpL.loadConfigParams(confFileName=paramsFileName)
    else:
        dParams= {}

    simFileName = dParams["simsFileName"]
    model2Simulate = dParams["model2Simulate"]

    dRxns2Annotation = mmL.getModelRxnsAnnotation(os.path.join(MODELDIR, dParams["backboneModel"])) ## uso il modello salvato con questa chiave perche' poi il salvataggio con cobrapy annulla le info contenute nel tag annotation

    dfdRxns2AnnotationNewModel = pd.read_csv(os.path.join(RAWDIR, "newModelReactionsToAnnotation.tsv"), sep = "\t")
    dfdRxns2AnnotationNewModel['dictAnnotation'] = dfdRxns2AnnotationNewModel['dictAnnotation'].apply(ast.literal_eval)
    dRxns2AnnotationNewModel = dfdRxns2AnnotationNewModel.set_index('reaction').to_dict()
    dRxns2Annotation.update(dRxns2AnnotationNewModel["dictAnnotation"])

    dfdRxns2AnnotationExchNewModel = pd.read_csv(os.path.join(RAWDIR, "newModelExchangeReactionsToAnnotation.tsv"), sep = "\t")
    dfdRxns2AnnotationExchNewModel['dictAnnotation'] = dfdRxns2AnnotationExchNewModel['dictAnnotation'].apply(ast.literal_eval)
    dRxns2AnnotationExchNewModel = dfdRxns2AnnotationExchNewModel.set_index('reaction').to_dict()
    dRxns2Annotation.update(dRxns2AnnotationExchNewModel["dictAnnotation"])

    if "btn-mapFluxesOnKeggMaps" == ctx.triggered_id:
        model = cb.io.read_sbml_model(os.path.join(MODELDIR, model2Simulate))
        lPaths = []
        dRxnId2KeggId = {}
        for rxn in model.reactions:
            if rxn.id in dRxns2Annotation and "kegg.reaction" in dRxns2Annotation[rxn.id]:
                kRxn = dRxns2Annotation[rxn.id]["kegg.reaction"]
                dInfoRxn = kL.getKeggInfo("rn:" + kRxn)
                if "PATHWAY" in dInfoRxn:
                    for path in list(dInfoRxn["PATHWAY"].keys()):
                        new= re.sub('[a-z]+', keggOrgCode,path)
                        lPaths.append(new)
                dRxnId2KeggId[rxn.id] = "rn:" + dRxns2Annotation[rxn.id]["kegg.reaction"]

            if rxn.id in dRxns2Annotation and "kegg.pathway" in dRxns2Annotation[rxn.id]:
                lPaths += [dRxns2Annotation[rxn.id]["kegg.pathway"]]

            lPaths = gL.unique(lPaths)

        lPaths = gL.unique(lPaths)

        jet = cm = plt.get_cmap('viridis_r')

        dfFlux = pd.read_csv(os.path.join(OUTDIR, simFileName), sep = "\t")
        dfFlux_wKeggRxnId = dfFlux.replace({"Rxn": dRxnId2KeggId})

        if simFileName.startswith("FBA_") or simFileName.startswith("pFBA_"):
            highestFlux = dfFlux.Flux.max()
        elif simFileName.startswith("FVA_"):
            highestFlux = max(dfFlux.Min.max() + dfFlux.Max.max())

        norm = mpl.colors.LogNorm(vmin = 0.1, vmax = (abs(highestFlux)), clip = False)
        m = cmx.ScalarMappable(norm=norm, cmap=jet)

        for path in lPaths:
            try:
                grayMap = kL.downloadPdfMaps(path, "#cccccc", MAPDIR)
                gene = grayMap.genes
                for g in gene:
                    lRxns = g.reaction.split()
                    if any(r in list(dfFlux_wKeggRxnId["Rxn"]) for r in lRxns) is True:
                        ## identificare il loro flusso
                        dfFlux_currentRxns = dfFlux_wKeggRxnId[dfFlux_wKeggRxnId["Rxn"].isin(lRxns)]
                        if simFileName.startswith("FBA_") or simFileName.startswith("pFBA_"):
                            if flux2Map == 'Sum of their fluxes':
                                resultingFlux = dfFlux_currentRxns.Flux.sum()
                            elif flux2Map == 'Mean of their fluxes':
                                resultingFlux = dfFlux_currentRxns.Flux.mean()
                            elif flux2Map == 'Minimum of their fluxes':
                                resultingFlux = dfFlux_currentRxns.Flux.min()
                            elif flux2Map == 'Maximum of their fluxes':
                                resultingFlux = dfFlux_currentRxns.Flux.max()
                        elif simFileName.startswith("FVA_"):
                            if flux2Map == 'Sum of their fluxes':
                                resultingFlux = dfFlux_currentRxns.Min.sum() + dfFlux_currentRxns.Max.sum()
                            elif flux2Map == 'Mean of their fluxes':
                                resultingFlux = (dfFlux_currentRxns.Min + dfFlux_currentRxns.Max).mean()
                            elif flux2Map == 'Minimum of their fluxes':
                                resultingFlux = dfFlux_currentRxns.Min.min()
                            elif flux2Map == 'Maximum of their fluxes':
                                resultingFlux = dfFlux_currentRxns.Max.max()

                        color = mpl.colors.rgb2hex(m.to_rgba(resultingFlux))

                        for gg in g.graphics:
                            gg.fgcolor = "#ffffff"
                            gg.bgcolor = color

                canvas = KGMLCanvas(grayMap)
                if path != "ko01100":
                    canvas.import_imagemap = True
                else:
                    canvas.import_imagemap = False
                canvas.draw(os.path.join(MAPDIR,  os.path.splitext(simFileName)[0] + "_" + path + ".pdf"))
                kL.PDF(os.path.join(MAPDIR, os.path.splitext(simFileName)[0] + "_" + path + ".pdf"))
                for g in gene:
                    for gg in g.graphics:
                        gg.fgcolor = "#cccccc"
                        gg.bgcolor = "#cccccc"

            except:
                print(path, " does not exist!\n")
        return html.H5(f'The resulting maps have been saved in the {MAPDIR} directory.')
