#To write to html use fig.write_html('file_name.html')

import plotly.graph_objects as go
import json
import os
import sys
import requests
import io
import pandas as pd
import numpy as np
import rdflib
from SPARQLWrapper import SPARQLWrapper, JSON

def heatmap(csv):
    
    #Store necessary info
    
    color_scale = [[0,'#566573'],[0.5,'#566573'],[0.5,'#1abc9c'],[1,'#1abc9c']]
    ticktext = ['Absent','Present']
    tickvals = [0,1]

    #Store in DataFrame
    
    df = pd.read_csv(csv,header=None)
    df_new = pd.DataFrame(list(map(str,df.values.tolist())))

    #Draw figure
    
    trace = go.Heatmap(
        x = df.loc[0],
        y = df[0].tolist(),
        z = df, 
        colorscale = color_scale,
        colorbar = dict(thickness=25, 
                        tickvals=tickvals, 
                        ticktext=ticktext))
    
    #hoverlabel not made yet
    
    fig = go.Figure(data = [trace])
    return fig

#credit: Doug Fils

def get_sparql_dataframe(service, query):
    """
    Helper function to convert SPARQL results into a Pandas data frame.
    """
    sparql = SPARQLWrapper(service)
    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    result = sparql.query()

    processed_results = json.load(result.response)
    cols = processed_results['head']['vars']
    
    out = []
    for row in processed_results['results']['bindings']:
        item = []
        for c in cols:
            item.append(row.get(c, {}).get('value'))
        out.append(item)
    return pd.DataFrame(out, columns=cols)

def get_taxa(protein_list,url_query_3):
    values_string = '"^^xsd:token "'.join(protein_list)
    url_query_4 = url_query_3.replace('list_goes_here',values_string)
    #print(url_query_4)
    return(get_sparql_dataframe(OPP_SERVE,namespaces + url_query_4))

def depth_profile(protein_list):
    #Define the SPARQL endpoints and main namespace

    global OPP_SERVE
    OPP_SERVE = "https://kg.oceanproteinportal.org/sparql"  #OPP SPARQL Endpoint

    #Define expected namespaces from knowledge graph

    global namespaces
    namespaces = """
    PREFIX view: <http://schema.oceanproteinportal.org/v2/views/>
    PREFIX opp: <http://schema.oceanproteinportal.org/v2/>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    PREFIX owl: <http://www.w3.org/2002/07/owl#>
    PREFIX dc: <http://purl.org/dc/elements/1.1/>
    PREFIX dcterms: <http://purl.org/dc/terms/>
    PREFIX geo: <http://www.w3.org/2003/01/geo/wgs84_pos#>
    PREFIX geosparql: <http://www.opengis.net/ont/geosparql#>
    PREFIX obo: <http://purl.obolibrary.org/obo/>
    PREFIX odo: <http://ocean-data.org/schema/>
    PREFIX odo_dt: <http://ocean-data.org/schema/data-type/v1.0/>
    PREFIX prov: <http://www.w3.org/ns/prov#>
    PREFIX schema: <http://schema.org/>
    PREFIX sf: <http://www.opengis.net/ont/sf#>
    PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
    PREFIX sosa: <http://www.w3.org/ns/sosa/>
    PREFIX ssn: <http://www.w3.org/ns/ssn/>
    PREFIX ssn_system: <http://www.w3.org/ns/ssn/systems/>
    PREFIX time: <http://www.w3.org/2006/time#>
    PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
    """

    url_query = """SELECT DISTINCT ?label ?sample ?feature STR(?depth) as ?depth ?station ?lat ?lon STR(?spectralCount) as ?spectralCount ?orfName ?proteinDescription
    WHERE 
      { <urn:opp:dataset:metzyme-0.2> opp:storesResultsForSample ?sample .
        ?sample rdfs:label ?label . 
        ?sample sosa:isSampleOf ?feature .
        ?feature opp:depth ?depth .
        ?feature opp:inVicinityOfStation ?stationInfo .
        ?stationInfo opp:stationName ?station .	
    OPTIONAL { 
        ?feature geo:lat ?lat .
        ?feature geo:lon ?lon .
             }
        ?observation opp:observationInSample ?sample .
        ?observation schema:value ?spectralCount . 
        ?observation rdf:type opp:ProteinTotalSpectralCount .
        ?protein opp:sampleProteinProperty ?observation .
        ?protein opp:proteinIdentifier ?orfName . 
        ?protein opp:proteinName ?proteinDescription .
    FILTER REGEX(?proteinDescription, "ftsz", "i") . 
    }
    ORDER BY ?label"""

    url_query_3 = """SELECT DISTINCT ?tn ?pn
    WHERE {
      VALUES ?pn {"list_goes_here"^^xsd:token}
      ?p opp:proteinIdentifier ?pn .
      ?p opp:identifier ?t .
      ?t schema:name ?tn .
      ?t rdf:type opp:NCBITaxonIdentifier .
      }
    ORDER BY ?tn"""

    # submit query and return as df
    df = get_sparql_dataframe(OPP_SERVE, namespaces + url_query)
    
    df = df.astype({'station':'str'})
    df = df.astype({'depth':'int'})
    df = df.astype({'spectralCount':'int'})
    
    options_list = df['orfName'].unique()
    protein_groups = df.groupby('orfName',sort=False)
    
    df_taxa = get_taxa(protein_list,url_query_3)
    
    #Create figure
    
    fig = go.Figure()
    
    for ptn in protein_list:
        df_new = protein_groups.get_group(ptn)
        groups = df_new.groupby('station', sort=False)
        #print(list(df_taxa['pn']))
        for name,group in groups:
            fig.add_trace(go.Scatter(
                x=group['spectralCount'],
                y=group['depth'],
                mode='lines+markers',
                name = 'St. ' + name + ' ' + ptn,
                line = dict(width=2.5),
                marker = dict(size=8),
                hovertemplate='Station: ' + name + '<br>Protein: ' + ptn + '<br>Taxon: ' + df_taxa['tn'][list(df_taxa['pn']).index(ptn)].replace('"','') + '<br>Spectral Counts' + ': %{x}<br>Depth: %{y}<extra></extra>'))

    #Invert y-axis and add labels

    fig['layout']['yaxis']['autorange'] = 'reversed'
    fig.update_layout(legend_title_text='Station and Protein')
    fig.update_yaxes(title_text='Depth (m)',tickfont=dict(size=18),title_font_size=18)
    fig.update_xaxes(title_text='Spectral Counts',tickfont=dict(size=18),title_font_size=18)
    fig.update_layout(title='Depth Profile',title_font_size=20)

    return fig

def getphylocounts(level,mtwithtaxon,list1,list2,list3):
        global allTaxaList
        allTaxaList = list1
        global xlist
        xlist = list2
        global sankeyList
        sankeyList = list3
        df2 = mtwithtaxon.groupby(level).size()
        df2_list = list(zip(df2,df2.index))
        for entry in df2_list:
            count = entry[0]
            phylo1 = entry[1][-1]
            phylo2 = entry[1][-2]
            if not phylo1 in allTaxaList:
                allTaxaList.append(phylo1)
                xlist.append(1/len(level))
            if not phylo2 in allTaxaList:
                allTaxaList.append(phylo2)
                xlist.append(1/len(level))
            #global sankeyList
            sankeyList['source'].append(allTaxaList.index(phylo2))
            sankeyList['target'].append(allTaxaList.index(phylo1))
            sankeyList['value'].append(count)

def sankey_plot(taxons_file,mtoutput_file):
    '''for the Sankey to work we need to return the taxon information for all the genome hits in metatryp
    this is essentially the information needed to return the LCA
    I assume this is relatively easy to output from metatryp API, but I am recreating here for testing the Sankey'''

    #this is the taxon table we gave to David for LCA analysis in metatryp
    taxons = pd.read_csv(taxons_file)

    #example metatryp output for peptide LSHQAIAEAIGSTR
    mtoutput = pd.read_csv(mtoutput_file)

    #add taxon information to mtoutput 
    mtwithtaxon = pd.merge(mtoutput, taxons, left_on='Taxon Id', right_on='taxon_oid')

    #IMPORTANT: to avoid circularity issues, need to replace replace lower level "unclassified":
    mtwithtaxon['Genus'] = mtwithtaxon['Genus'].replace('unclassified', 'unclassified Genus')
    mtwithtaxon['Family'] = mtwithtaxon['Family'].replace('unclassified', 'unclassified Family')
    mtwithtaxon['Order'] = mtwithtaxon['Order'].replace('unclassified', 'unclassified Order')
    mtwithtaxon['Class'] = mtwithtaxon['Class'].replace('unclassified', 'unclassified Class')
    mtwithtaxon['Phylum'] = mtwithtaxon['Phylum'].replace('unclassified', 'unclassified Phylum')
    mtwithtaxon['Domain'] = mtwithtaxon['Domain'].replace('unclassified', 'unclassified Domain')

    level0 = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    level1 = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus']
    level2 = ['Domain', 'Phylum', 'Class', 'Family', 'Order']
    level3 = ['Domain', 'Phylum', 'Class', 'Family']
    level4 = ['Domain', 'Phylum', 'Class']
    level5 = ['Domain', 'Phylum']

    xlist = []
    allTaxaList = []
    sankeyList = {'source':[],'target':[],'value':[]}
              
    levels = [level0, level1, level2, level3, level4 ,level5]
    for level in levels:
        getphylocounts(level,mtwithtaxon,allTaxaList,xlist,sankeyList)

    fig = go.Figure(data=[go.Sankey(
        node = {
          'line': dict(color = "black", width = 0.5),
          'label': allTaxaList,
          'color': "blue",
            'x': xlist,
        },
        link = sankeyList,
    )])
    return fig

#Potentially change these to classes?

hm = heatmap('peptide_heatmap_example.csv')
hm.write_html('graphics_heatmap_output.html')

sp = sankey_plot('taxontable.csv','examplepeptide.csv')
sp.write_html('graphics_sankey_output.html')

dp = depth_profile(['NODE_162056_length_786_cov_1.9316_1_786_+','NODE_84267_length_1274_cov_3.33388_125_1252_+'])
dp.write_html('graphics_depth_profile_output.html')
