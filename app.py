import dash
from dash import html, register_page, callback, Input, Output, dcc, Dash, State, dash_table, ctx
import dash_bootstrap_components as dbc
import pandas as pd
import cobra as cb
import re
import bioservices.kegg as kegg
import os
import sys
import numpy as np
import dash_uploader as du

MODULEDIR = os.path.split(os.path.dirname(__file__))[0]
sys.path.append(MODULEDIR)
import utils.keggLib as kL
import utils.genericLib as gL
import utils.configParamsLib as cpL
import utils.experimentalSetupLib as expSL
import utils.manipulateModelLib as mmL
import utils.fluxAnalysisLib as faL
import utils.figuresLib as figL
import utils.modelReconstructionLib as recL
from Bio.Graphics.KGML_vis import KGMLCanvas
from IPython.display import Image, HTML

import homePage as homeL
import backbonePage as backboneL
import newModelPage as newModelL
import mergePage as mergeL
import constrainPage as constrainL
import simulatePage as simulateL
import mappingPage as mappingL
import chatbotPage as chatbotL

# Setting the parameters
RAWDIR = gL.dDirs["raw"]
OUTDIR = gL.dDirs["out"]
MODELDIR = gL.dDirs["mods"]
FIGUREDIR = gL.dDirs["figs"]
UTILSDIR =  gL.dDirs["utils"]
MAPDIR =  gL.dDirs["map"]

timeStamp = gL.getTimeStamp()

# setting create model params
configFileName = "parameters.toml"
dConfParams = cpL.loadConfigParams(confFileName=configFileName)

# logStrm = gL.logFileOpen(gL.dDirs["logs"], timeStamp, gL.getBaseName(sys.argv[0]))
# gL.toLog(logStrm, f"\ntimeStamp: {timeStamp}")
# cpL.toLogConfigParams(confFileName=configFileName, logStream=logStrm)

## Load or create file including correspondence between metabolites name and their KEGG identifier
dfName2KeggId = pd.read_csv("name2KeggId.tsv", sep = "\t")

lActiveMet = {"C00080": "H+", "C00001":"H2O", "C00007": "Oxygen", "C00009": "Phosphate"}
for activeMet in lActiveMet:
    if lActiveMet[activeMet].lower() not in dfName2KeggId["Name"].values:
        dfName2KeggId.loc[len(dfName2KeggId)] = [activeMet, lActiveMet[activeMet]]
dfName2KeggId.to_csv("name2KeggId.tsv", sep = "\t", index = False)

## manca:
## 1. come scegliere quali rxns mappare sulle mappe di kegg quando ci sono rxns nel modello che corrispondono alla stessa rxn
##    in kegg essendo kegg stesso non compartimentalizzato?
## 3. sistemare file parameters
## 4. nel new model: quando genero il modello e il modello esiste già lo strumento mi carica la versione già esistente. Chiedere in tal caso all utente se caricare versione esistente o se cmq rigeerare modello.

############################################################################
############################################################################
# ## creation log file of parameters
# paramsStream = open("parametersSetting_" + timeStamp + ".toml", mode="w")
# dParams = {}
# for param in dParams:
#     gL.toLog(paramsStream, param)
#     gL.toLog(paramsStream, dParams[param])
#     gL.toLog(paramsStream, "\n")
# paramsStream.close()
############################################################################
############################################################################


# Code source: https://medium.com/@mcmanus_data_works/how-to-create-a-multipage-dash-app-261a8699ac3f

# Toggle the themes at [dbc.themes.LUX]
# The full list of available themes is:
# BOOTSTRAP, CERULEAN, COSMO, CYBORG, DARKLY, FLATLY, JOURNAL, LITERA, LUMEN,
# LUX, MATERIA, MINTY, PULSE, SANDSTONE, SIMPLEX, SKETCHY, SLATE, SOLAR,
# SPACELAB, SUPERHERO, UNITED, YETI, ZEPHYR.
# To see all themes in action visit:
# https://dash-bootstrap-components.opensource.faculty.ai/docs/themes/explorer/

#NAVBAR = create_navbar()
# To use Font Awesome Icons
FA621 = "https://use.fontawesome.com/releases/v6.2.1/css/all.css"
APP_TITLE = "Pipeline"

app = dash.Dash(
    __name__,
    suppress_callback_exceptions=True,
    external_stylesheets=[
        dbc.themes.LUX,  # Dash Themes CSS
        FA621,  # Font Awesome Icons CSS
    ],
    title=APP_TITLE,
    #use_pages=True,  # New in Dash 2.7 - Allows us to register pages
)

du.configure_upload(app, gL.CWDIR)


# the style arguments for the sidebar. We use position:fixed and a fixed width
SIDEBAR_STYLE = {
    "position": "fixed",
    "top": 0,
    "left": 0,
    "bottom": 0,
    "width": "16rem",
    "padding": "2rem 1rem",
    "background-color": "#f8f9fa",
}

# the styles for the main content position it to the right of the sidebar and
# add some padding.
CONTENT_STYLE = {
    "margin-left": "18rem",
    "margin-right": "2rem",
    "padding": "2rem 1rem",
}

sidebar = html.Div(
    [
        # html.H2("Pipeline navigation menu", className="display-4"),
        html.H2("Navigation menu"),
        html.Hr(),
        # html.P(
        #     "A simple sidebar layout with navigation links", className="lead"
        # ),
        dbc.Nav(
            [
                dbc.NavLink("Backbone model", href="/backboneModel", active="exact"),
                dbc.NavLink("New model", href="/createModel", active="exact"),
                dbc.NavLink("Create integrated model", href="/merge", active="exact"),
                dbc.NavLink("Constrain integrated model", href="/constrain", active="exact"),
                dbc.NavLink("Simulate model", href="/simulate", active="exact"),
                dbc.NavLink("Mapping simulated fluxes", href="/mapping", active="exact"),
                dbc.NavLink("Chatbot-based exploratory analysis", href="/chatbot", active="exact"),
            ],
            vertical=True,
            pills=True,
            className="bg-light",
        ),
    ],
    style=SIDEBAR_STYLE,
)

content = html.Div(id="page-content", style=CONTENT_STYLE)


app.layout = html.Div([dcc.Location(id="url"), sidebar, content, dcc.Store(id="store_keggOrgCode", storage_type='session'),
        dcc.Store(id='store-df', storage_type='session'), dcc.Store(id='store-dict', storage_type='session'),
        dcc.Store(id='store-newModelName', storage_type='session'), dcc.Store(id='store-mergedModels', storage_type='session'),
        dcc.Store(id='store-constrainedModels', storage_type='session'), dcc.Store(id='store-downloadFolder', storage_type='session')])


@app.callback(Output("page-content", "children"), [Input("url", "pathname")])
def render_page_content(pathname):
    if pathname == "/":
        return homeL.homeLayout()
    elif pathname == "/backboneModel":
        return backboneL.backboneModelLayout()
    elif pathname == "/createModel":
        return newModelL.newModelLayout()
    elif pathname == "/merge":
        return mergeL.mergeLayout()
    elif pathname == "/constrain":
        return constrainL.constrainLayout()
    elif pathname == "/simulate":
        return simulateL.simulateLayout()
    elif pathname == "/mapping":
        return mappingL.mappingLayout()
    elif pathname == "/chatbot":
        return chatbotL.chatbotLayout()
    # If the user tries to reach a different page, return a 404 message
    return html.Div(
        [
            html.H1("404: Not found", className="text-danger"),
            html.Hr(),
            html.P(f"The pathname {pathname} was not recognised..."),
        ],
        className="p-3 bg-light rounded-3",
    )

server = app.server

if __name__ == '__main__':
    app.run_server(debug=True)
