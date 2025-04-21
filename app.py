import dash
from dash import dcc, html, dash_table
from dash.dependencies import Input, Output
import pandas as pd
import sqlite3
import re
import requests


# --- Load Data ---

# Adjust paths to point to files in the same directory as app.py or a cloud storage URL
metabolite_db_path = "https://drive.google.com/uc?export=download&id=1rH9mEdVT31x-oDvORNMUrHKoVAi0svmw"

# Orthology mapping
orth_conv = pd.read_csv("https://drive.google.com/uc?export=download&id=14KseMAFevC4iSChtyiSBZFDpo8tLw2KG", sep="\t", compression="gzip", dtype=str)
m2h = dict(zip(orth_conv["mouse_symbol"], orth_conv["human_symbol"]))

# LR Pairs
mouse_LRs = pd.read_csv("https://drive.google.com/uc?export=download&id=1WdxNVUZ8KIm_goVGmJgRXj7RcUa710cw")
mouse_LRs.rename(columns={"LIGAND_1": "Ligand", "RECEPTOR_1": "Receptor"}, inplace=True)
human_LRs = mouse_LRs[["LRI", "Ligand", "Receptor", "DataBase", "is_hormone", "In_Plasma", "age", "Exosome", "Secretomics"]].copy()
human_LRs['Ligand'] = human_LRs['Ligand'].map(m2h)
human_LRs['Receptor'] = human_LRs['Receptor'].map(m2h)
human_LRs = human_LRs.dropna(subset=["Ligand", "Receptor"])
human_LRs['LRI'] = human_LRs['Ligand'] + ":" + human_LRs['Receptor']

mouse_LR_between = pd.read_csv(resources_path + "Long_dist_LRs.csv")
mouse_LR_between.rename(columns={"LIGAND_1": "Ligand", "RECEPTOR_1": "Receptor"}, inplace=True)
human_LR_between = mouse_LR_between[["LRI", "Ligand", "Receptor", "DataBase", "is_hormone", "In_Plasma", "age", "Exosome", "Secretomics"]].copy()
human_LR_between['Ligand'] = human_LR_between['Ligand'].map(m2h)
human_LR_between['Receptor'] = human_LR_between['Receptor'].map(m2h)
human_LR_between = human_LR_between.dropna(subset=["Ligand", "Receptor"])
human_LR_between['LRI'] = human_LR_between['Ligand'] + ":" + human_LR_between['Receptor']

# Pathways
pathways = pd.read_csv("https://drive.google.com/uc?export=download&id=1eTUZS32OxTMqkV_O0lmG8AZd471Oxv-q")
pathways["Pathway"] = pathways["Pathway"].str.replace(r' - Mus musculus.*', '', regex=True)

# Metabolites Database
url = "https://drive.google.com/uc?export=download&id=1rH9mEdVT31x-oDvORNMUrHKoVAi0svmw"
db_path = "metabolite_db.sqlite"
with requests.get(url, stream=True) as r:
    r.raise_for_status()
    with open(db_path, 'wb') as f:
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)
conn = sqlite3.connect(db_path)

query = "SELECT name FROM sqlite_master WHERE type='table';"
tables = pd.read_sql_query(query, conn)
dataframes = {}
for table in tables['name']:
    dataframes[table] = pd.read_sql_query(f"SELECT * FROM {table};", conn)
conn.close()
metabolite_gene_pathway = pd.read_csv("https://drive.google.com/uc?export=download&id=1nTZoK9tHV5OS0V3bi-Jr8tMKAAuWpSQs")
metabolite_gene_pathway.drop(columns="pathway", inplace=True)
metabolite_gene_pathway.rename(columns={"description": "pathway", "HMDB_ID": "hmdb"}, inplace=True)

# Hormones
hormone_organs = pd.read_csv("https://drive.google.com/uc?export=download&id=1k3BDs2tYrpynPDvfJ90s4Lh5GIbkYlZK")
hormones_base = pd.read_json("https://drive.google.com/uc?export=download&id=1Wqg2nmmQI2KREcXuwEo18dtEdPhYbul3").T

# Get unique values for dropdowns
all_tissues = dataframes['tissue_location']['tissue_location'].unique().tolist()
all_biospecimens = dataframes['biospecimen_location']['biospecimen_location'].unique().tolist()
all_diseases = dataframes['disease']['disease'].unique().tolist()
all_organs = sorted(set([org.strip() for organs in hormone_organs["Organs_Involved"].dropna() for org in organs.split(',')]))

# --- Helper Functions ---

def get_hormones_by_organs(hormone_organs, must_have, equivalent_organs):
    pattern = "".join(f"(?=.*{re.escape(org)})" for org in must_have)
    for group in equivalent_organs:
        group_pattern = "|".join(map(re.escape, group))
        pattern += f"(?=.*({group_pattern}))"
    mask = hormone_organs["Organs_Involved"].str.contains(pattern, case=False, na=False, regex=True)
    return list(hormone_organs.loc[mask, "Hormone"])

def hormone_signalling(hormones, prior_hormone_signalling, LR_between):
    hormones = hormones[hormones.index.str.contains('|'.join(prior_hormone_signalling), na=False, case=False)]
    source_hormones = set().union(*[hormones['source'].get(key, []) for key in prior_hormone_signalling])
    target_hormones = set().union(*[hormones['target'].get(key, []) for key in prior_hormone_signalling])
    hormones_list = list(source_hormones) + list(target_hormones)
    hormone_net = LR_between[(LR_between['is_hormone'] == 1) & ((LR_between['Ligand'].isin(hormones_list)) | (LR_between['Receptor'].isin(hormones_list)))]
    return hormones, hormone_net

def proper_title_case(lst):
    return [word[0].upper() + word[1:].lower() if word else word for word in lst]

# --- Initialize App ---

app = dash.Dash(__name__)

# --- Layout ---

app.layout = html.Div([
    dcc.Tabs(id='tabs', value='lr-pairs', children=[
        dcc.Tab(label='LR Pairs', value='lr-pairs'),
        dcc.Tab(label='Pathways', value='pathways'),
        dcc.Tab(label='Metabolites', value='metabolites'),
        dcc.Tab(label='Hormones', value='hormones'),
    ]),
    html.Div(id='tabs-content')
])

# --- Callbacks ---

@app.callback(
    Output('tabs-content', 'children'),
    Input('tabs', 'value')
)
def render_content(tab):
    if tab == 'lr-pairs':
        return html.Div([
            html.H3("Ligand-Receptor Pairs"),
            dcc.Dropdown(id='species-dropdown', options=[{'label': 'Mouse', 'value': 'mouse'}, {'label': 'Human', 'value': 'human'}], value='mouse'),
            html.Div([html.Label("Filter by is_hormone:"), dcc.Dropdown(id='is-hormone-filter', options=[{'label': '0', 'value': 0}, {'label': '1', 'value': 1}])]),
            html.Div([html.Label("Filter by In_Plasma:"), dcc.Dropdown(id='in-plasma-filter', options=[{'label': str(v), 'value': v} for v in mouse_LRs['In_Plasma'].unique()])]),
            dash_table.DataTable(id='table-lr', page_size=20)
        ])
    elif tab == 'pathways':
        return html.Div([
            html.H3("Pathways and Genes"),
            dcc.Input(id='pathway-search', type='text', placeholder='Search Pathway'),
            dash_table.DataTable(id='table-pathways', page_size=20, columns=[{"name": i, "id": i} for i in pathways.columns])
        ])
    elif tab == 'metabolites':
        return html.Div([
            html.H3("Metabolites"),
            html.Div([html.Label("Tissues:"), dcc.Dropdown(id='tissues-dropdown', options=[{'label': t, 'value': t} for t in all_tissues], multi=True)]),
            html.Div([html.Label("Biospecimen Locations:"), dcc.Dropdown(id='biospecimens-dropdown', options=[{'label': b, 'value': b} for b in all_biospecimens], multi=True)]),
            html.Div([html.Label("Diseases:"), dcc.Dropdown(id='diseases-dropdown', options=[{'label': d, 'value': d} for d in all_diseases], multi=True)]),
            html.Div([html.Label("Display:"), dcc.Dropdown(id='met-display-dropdown', options=[{'label': 'Metabolites', 'value': 'metabolites'}, {'label': 'Metabolite Net', 'value': 'met_net'}], value='metabolites')]),
            dash_table.DataTable(id='table-metabolites', page_size=20)
        ])
    elif tab == 'hormones':
        return html.Div([
            html.H3("Hormones"),
            html.Div([html.Label("Must-Have Organs:"), dcc.Dropdown(id='must-have-dropdown', options=[{'label': o, 'value': o} for o in all_organs], multi=True)]),
            html.Div([html.Label("Equivalent Organs (group per line, comma-separated):"), dcc.Textarea(id='equivalent-organs-input', placeholder='e.g., brain,hypothalamus\nliver,kidney', style={'width': '100%', 'height': 100})]),
            html.Div([html.Label("Display:"), dcc.Dropdown(id='hor-display-dropdown', options=[{'label': 'Hormones', 'value': 'hormones'}, {'label': 'Hormone Net', 'value': 'hor_net'}], value='hormones')]),
            dash_table.DataTable(id='table-hormones', page_size=20)
        ])

# LR Pairs Tab
@app.callback(
    Output('table-lr', 'data'),
    [Input('species-dropdown', 'value'),
     Input('is-hormone-filter', 'value'),
     Input('in-plasma-filter', 'value')]
)
def update_lr_table(species, is_hormone, in_plasma):
    df = mouse_LRs if species == 'mouse' else human_LRs
    if is_hormone is not None:
        df = df[df['is_hormone'] == is_hormone]
    if in_plasma is not None:
        df = df[df['In_Plasma'] == in_plasma]
    return df.to_dict('records')

@app.callback(
    Output('table-lr', 'columns'),
    Input('species-dropdown', 'value')
)
def update_lr_columns(species):
    df = mouse_LRs if species == 'mouse' else human_LRs
    return [{"name": i, "id": i} for i in df.columns]

# Pathways Tab
@app.callback(
    Output('table-pathways', 'data'),
    Input('pathway-search', 'value')
)
def update_pathways_table(search):
    df = pathways
    if search:
        df = df[df['Pathway'].str.contains(search, case=False, na=False)]
    return df.to_dict('records')

# Metabolites Tab
@app.callback(
    [Output('table-metabolites', 'data'),
     Output('table-metabolites', 'columns')],
    [Input('tissues-dropdown', 'value'),
     Input('biospecimens-dropdown', 'value'),
     Input('diseases-dropdown', 'value'),
     Input('met-display-dropdown', 'value')]
)
def update_metabolites_table(tissues, biospecimens, diseases, display):
    df_tissue = dataframes['tissue_location']
    df_biospecimen = dataframes['biospecimen_location']
    df_disease = dataframes['disease']
    if tissues:
        pattern = '|'.join(tissues)
        df_tissue = df_tissue[df_tissue['tissue_location'].str.contains(pattern, case=False, na=False)]
    if biospecimens:
        df_biospecimen = df_biospecimen[df_biospecimen['biospecimen_location'].isin(biospecimens)]
    if diseases:
        pattern = '|'.join(diseases)
        df_disease = df_disease[df_disease['disease'].str.contains(pattern, case=False, na=False)]
    metabolites = pd.merge(pd.merge(df_tissue, df_biospecimen, on="hmdb"),
                           pd.merge(df_disease, dataframes['pathway'], on="hmdb"), on="hmdb")
    metabolites = pd.merge(metabolites, metabolite_gene_pathway[["pathway", "hmdb", "symbol", "organism"]], on=["pathway", "hmdb"])
    if display == 'metabolites':
        df = metabolites
    else:  # 'met_net'
        df = metabolites[(metabolites['symbol'].isin(set(human_LR_between["Ligand"]))) & (metabolites['symbol'].isin(set(human_LR_between["Receptor"])))]
    columns = [{"name": i, "id": i} for i in df.columns]
    return df.to_dict('records'), columns

# Hormones Tab
@app.callback(
    [Output('table-hormones', 'data'),
     Output('table-hormones', 'columns')],
    [Input('must-have-dropdown', 'value'),
     Input('equivalent-organs-input', 'value'),
     Input('hor-display-dropdown', 'value')]
)
def update_hormones_table(must_have, equivalent_organs_text, display):
    must_have = must_have or []
    equivalent_organs = [line.split(',') for line in equivalent_organs_text.split('\n') if line.strip()] if equivalent_organs_text else []
    prior_hormone_signalling = get_hormones_by_organs(hormone_organs, must_have, equivalent_organs)
    hormones = hormones_base.applymap(lambda x: proper_title_case(x) if isinstance(x, list) else x)
    hormones_df, hormone_net = hormone_signalling(hormones, prior_hormone_signalling, human_LR_between)
    df = hormones_df if display == 'hormones' else hormone_net
    columns = [{"name": i, "id": i} for i in df.columns]
    return df.to_dict('records'), columns

# --- Run App ---

if __name__ == '__main__':
    app.run_server(host='0.0.0.0', port=8050, debug=False)