from dash import html, dcc

def homeLayout():
    """

    :return: A Div containing dashboard content of homepage loading.
    """
    layout = html.Div([
        html.H1('GPRuler'),
        html.Br(),
        html.H5('Abstract'),
        html.H6('Genome-scale metabolic models are increasingly used tools for providing full functional state readout of given organisms. Metabolic models are valued allies in supporting wet experiments for the design of new engineered organism platforms able to perform novel functions and produce new desired compounds. Nevertheless, tools able to reconstruct the corresponding digital twins by merging content of new pathways according to the knowledge stored into already available metabolic models are still missing. In this regard, we propose an open-source python-based framework that can be used to load, curate, reconstruct and integrate metabolic networks, yielding SBML integrated models that are ready for further analyses.'),
        html.Br(), html.Br(),
        html.H5('Documentation'),
        html.Div([html.H6('The documentation to install and use GPRuler is included in the '), html.A("Wiki", href="https://git-ricerca.unimib.it/myCoding/replay"), html.H6('.')]),
        html.Br(), html.Br(),
        html.H5('Getting Help'),
        html.H6('For support, please contact: '),
        html.A([html.H6('Dario Pescini')], title ='email_dario', href='mailto:dario.pescini@unimib.it'),
        html.A([html.H6('Marzia Di Filippo')], title ='email_dario', href='mailto:marzia.difilippo@unimib.it'),
        html.A([html.H6('Ennio Guzmán')], title ='email_dario', href='mailto:e.guzmancolmenares@campus.unimib.it'),
        html.Br(), html.Br(),
        dcc.Link('START', href='/backboneModel'),
    ])
    return layout
