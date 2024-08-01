from dash import html, callback, Input, Output, dcc, Dash, State, dash_table, ctx, no_update
import dash_table as dt
from dash.dash_table.Format import Group
import base64
import subprocess
import threading
import time
from uploadModel import upload

# Run the Streamlit app
def run_streamlit():
    try:
        process = subprocess.Popen(["streamlit", "run", "app_chatbot.py", "--server.port", "8501", "--server.headless", "true"],
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # Wait for streamlit to start
        time.sleep(10)
        return process
    except Exception as e:
        print(f"Error al iniciar Streamlit: {e}")
        return None

# Define the URL for streamlit
def get_streamlit_url():
    return "http://localhost:8501"

def chatbotLayout():
    layout = html.Div([
        html.H1('Chatbot'),

        html.H5('Load a new model.'),
        
        dcc.Upload(
            id='upload-rxnsNewModel',
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
            },
            multiple=False
        ),
        
        dcc.Loading(
            id="loading-newModelRxns",
            children=[html.Div([html.Div(id="ModelPathForChatbot")])],
            type="circle",
        ),

        html.Div(id='error-message', style={'color': 'red'}),

        html.Iframe(
            id='xml-viewer',
            src=get_streamlit_url(),  # URL de Streamlit
            style={
                'width': '90%',
                'height': '500px',
                'border': '1px solid black',
                'display': 'none'  # Initially hidden
            },
        ),
    ])
    return layout

@callback(
    Output('ModelPathForChatbot', 'children'),
    Output('xml-viewer', 'style'),
    Output('error-message', 'children'),
    Output('upload-rxnsNewModel', 'disabled'),
    Output('upload-rxnsNewModel', 'style'),  
    Input('upload-rxnsNewModel', 'contents'),
    State('upload-rxnsNewModel', 'filename')
)
def update_output(content, filename):
    if content is not None:
        if filename.endswith('.xml'):
            # Decode file content
            content_type, content_string = content.split(',')
            decoded = base64.b64decode(content_string).decode('utf-8')

            # Start streamlit in separate thread
            streamlit_thread = threading.Thread(target=run_streamlit)
            streamlit_thread.start()

            print(filename)
            upload(pathModel= filename, uri='bolt://localhost:11007', authUser='neo4j', authPass='1')


            return (filename, 
                    {'width': '90%', 'height': '500px', 'border': '1px solid black', 'display': 'block'}, 
                    '', 
                    True,  # Deactivate component of upload
                    {'width': '50%', 'height': '60px', 'lineHeight': '60px', 'borderWidth': '1px', 
                     'borderStyle': 'solid', 'borderRadius': '5px', 'textAlign': 'center', 
                     'margin': '10px', 'backgroundColor': '#e0e0e0', 'color': '#a0a0a0'}  # Estilo de desactivado
                   )
        else:
            # If the file is no XML, show error messagge
            return '', {'display': 'none'}, 'El archivo subido no tiene la extensión .xml. Por favor, sube un archivo con la extensión correcta.', False, {}
    else:
        return '', {'display': 'none'}, '', False, {}
