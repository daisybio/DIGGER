#!/usr/bin/env python
# coding: utf-8


import csv
import networkx as nx
import numpy as np
import pickle
import pandas as pd
import os

from django.conf import settings
from domain.Process import process_data as pr
from domain.Process.load_data import DomainG_all as Graphs
from sqlalchemy import text
from django.urls import reverse

# --- Get database connection aka 'SQLAlchemie engine'
engine = settings.DATABASE_ENGINE


# load the network
def load_obj(name):
    with open('domain/data/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)


# PPI= pd.read_csv( "domain/data/PPI_interface_mapped_to_exon.csv")

# function for visualization  DomainView
def vis_node_(node, organism, missing=False):
    p_id = '"' + node.split(".")[1] + '"'
    ppp = node.split(".")[1]
    node = node.split(".")[0] + "/" + node.split(".")[1]
    DomainG = Graphs[organism]
    G = nx.Graph()
    if DomainG.has_node(node):
        # copy over the edges with their attributes
        for edge in DomainG.edges(node, data=True):
            G.add_edge(edge[0], edge[1], **edge[2])

    # get the amount of edges that are predicted and how many are original
    predicted = 0
    original = 0
    for edge in G.edges(data=True):
        if edge[2]['confidence'] == 'original':
            original += 1
        else:
            predicted += 1
    print(f"Node {node} has {predicted} predicted edges and {original} original edges")

    g = nx.Graph()
    g.add_node(node)
    for n in G.nodes():
        domain = n.split("/")[1]
        if n == node:
            continue
        # add confidence attribute if it is not original in G, this is done seperately as to not overwrite the original
        if G.has_edge(n, node) and G[node][n]['confidence'] != 'original':
            if not g.has_edge(n, domain):
                g.add_edge(n, domain, confidence=G[n][node]['confidence'])
            if not g.has_edge(node, domain):
                g.add_edge(node, domain, confidence=G[n][node]['confidence'])
        else:
            g.add_edge(n, domain, confidence='original')
            g.add_edge(node, domain, confidence='original')


    # nodes
    N = {'original': [], 'high': [], 'mid': [], 'low': []}
    # edges
    E = {'original': [], 'high': [], 'mid': [], 'low': []}
    for n in g.nodes:
        # find if there is an edge connected to this node that is not predicted
        confidence = set()
        for e in g.edges(n, data=True):
            confidence.add(e[2]['confidence'])

        for confid in confidence:
            # Domain node
            if len(n.split("/")) == 1:
                # print(n)
                label_name = pr.Domain_name(n)[0] + ' (' + n + ')'
                N[confid].append(f'{{id: "{n + ppp}", label: "{label_name}", group: "Domain", physics: true, source: {p_id}, value: "4"}},')
                continue

            gene = n.split("/")[0]
            # Main Domain of interest
            if n == node:
                try:
                    domain_name = n.split("/")[1]
                    domain_name = pr.Domain_name(domain_name)[0] + ' (' + domain_name + ')'
                    # if the node is missing in the isoform, make it red
                    color = 'RED' if missing else 'BLUE'
                    N[confid].append(f'{{id: "{n + ppp}", label: "{pr.entrez_to_name(gene, organism)} - {domain_name}", group: "MDomain", color: {color}, physics: false, source: {p_id}, value: "5"}},')
                except KeyError:
                    print("KeyError with", node)
            # a gene node
            else:
                try:
                    N[confid].append(f'{{id: "{n + ppp}", label: "{pr.entrez_to_name(gene, organism)}", group: "protein", physics: true, source: {p_id}, value: "2"}},')
                except KeyError:
                    print("KeyError with", node)

    for e in g.edges(data=True):

        # edge to the main domain
        if any(x == node for x in e):
            colour = 'GREEN' if e[2]['confidence'] == 'original' else 'LIGHTGREEN'
            dashes = 'true' if missing else 'false'
            E[e[2]['confidence']].append(f'{{from: "{e[0] + ppp}", to: "{e[1] + ppp}", length:  L1, dashes: {dashes}, '
                                         f'width: WIDTH_SCALE * 2, color:  {colour} }},')

        else:
            colour = 'YELLOW' if e[2]['confidence'] == 'original' else 'YELLOW'

            E[e[2]['confidence']].append(f'{{from: "{e[0] + ppp}", to: "{e[1] + ppp}", length:  L2, width: '
                                         f'WIDTH_SCALE, color:  {colour} }},')

    return N, E, pr.entrez_to_name(node.split("/")[0], organism), node.split("/")[1]


# internal fucntion to get info about node
def vis_node(graph, node):
    # add name of gene later
    G = nx.ego_graph(graph, node)
    N = []
    E = []
    for n in G.nodes():
        # visualize node with entres ID
        N.append("{id: \"" + n + "\", label:  \"" + n.split("/")[0] + "-" + n.split("/")[1] + "\"},")

        # OR visulaize with gene names (it takes more time to convert)
        # N.append("{id: \""+n+"\", label:  \""+pr.entrez_to_name(n.split("/")[0])+"-"+n.split("/")[1]+"\"},")

    for e in G.edges(data=True):
        E.append("{from: \"" + e[0] + "\", to: \"" + e[1] + "\"},")
    return N, E, len(G) - 1


def expand_table(table, unique_domains, entrezID, organism):
    Graph = Graphs[organism]

    # ADD information about Pfam domain later here:

    # number of domains with known interactions
    table["Interactions mediated by the domain"] = np.zeros(table.shape[0])
    Text_nodes = {}
    text_edges = {}
    for domain in unique_domains:
        print("Unique domains:", entrezID, domain)
        node = entrezID + "/" + domain
        if Graph.has_node(node):
            # add a function that convert from entrez ID to 
            nodes, edges, degree = vis_node(Graph, node)

            filters = table['Pfam ID'].isin([domain])
            table.at[filters, "Interactions mediated by the domain"] = degree
            Text_nodes[domain] = nodes
            text_edges[domain] = edges

        else:
            pass

    return table, Text_nodes, text_edges


def exon_3D(exon_IDs, Ensemble_transID, organism):
    # NUMBER OF INTERACTION IN EVERY EXON INTERFACE
    N = []
    exons_in_interface = []
    # p1=PPI[ PPI['Transcript stable ID_x']==Ensemble_transID]
    # p2=PPI[ PPI['Transcript stable ID_y']==Ensemble_transID]
    text_edges = []
    partners = []

    query = """
                  SELECT * 
                  FROM ppi_data_""" + organism + """  
                  WHERE  "Transcript stable ID_x"=:ensemble_trans_id
                  """
    tr_1 = pd.read_sql_query(sql=text(query), con=engine, params={'ensemble_trans_id': Ensemble_transID})

    partners = tr_1['Transcript stable ID_y'].unique().tolist()

    query = """   
                  SELECT * 
                  FROM ppi_data_""" + organism + """  
                  WHERE "Transcript stable ID_y"=:ensemble_trans_id
                  """
    tr_2 = pd.read_sql_query(sql=text(query), con=engine, params={'ensemble_trans_id': Ensemble_transID})

    co_partners = []
    partner_dict = {Ensemble_transID: pr.tranID_convert(Ensemble_transID, organism)[3]}
    if len(partners) != 0:
        partners = set(partners + tr_2['Transcript stable ID_x'].unique().tolist())

        if len(partners) != 0:
            for x in partners:
                try:
                    partner_dict[x] = pr.tranID_convert(x, organism)[3]
                    co_partners.append(partner_dict[x])
                except TypeError:
                    pass

        # get unique exons where we extract the domains
        unqiue_exons = set(tr_1['Exon stable ID_y'].unique().tolist() + tr_2['Exon stable ID_x'].unique().tolist() +
                           tr_1['Exon stable ID_x'].unique().tolist())
        query = f"""
        SELECT "Exon stable ID", "Pfam ID"
        FROM exons_to_domains_data_{organism}
        WHERE "Exon stable ID" IN ({','.join([f"'{x}'" for x in unqiue_exons])})
        """
        domains = pd.read_sql_query(sql=text(query), con=engine)
        exon_to_domain = domains.set_index('Exon stable ID').to_dict()['Pfam ID']

        edge_tuples = set()
        # make tuples of the partners in fashion of ((transcript_x, exon_x), (transcript_y, exon_y))
        for idx, row in pd.concat([tr_1, tr_2]).iterrows():
            edge_tuples.add(((partner_dict.get(row['Transcript stable ID_x'], None),
                              exon_to_domain.get(row['Exon stable ID_x'], None)),
                             (partner_dict.get(row['Transcript stable ID_y'], None),
                              exon_to_domain.get(row['Exon stable ID_y'], None))))

        # convert the edge tuples to strings such that we can use it in the graph
        text_edges = [(f"{x[0][0]}/{x[0][1]}", f"{x[1][0]}/{x[1][1]}")
                      for x in edge_tuples if all([x[0][0], x[0][1], x[1][0], x[1][1]])]

    # avoid duplicates as well as False values
    co_partners = list(set([x for x in co_partners if type(x) == str]))
    print(f"Co partners: {co_partners}")

    for exon_ID in exon_IDs:
        # print(exon_ID)
        # --- Get tables from database

        '''
        query = """
                
                FROM ppi_data 
                WHERE "Exon stable ID_x"=:exon_id AND "Transcript stable ID_x"=:ensemble_trans_id
                """
        p1 = pd.read_sql_query(sql=text(query), con=engine, params={'exon_id': exon_ID,
                                                                    'ensemble_trans_id': Ensemble_transID})

        query = """
                SELECT DISTINCT "Transcript stable ID_x", "u_ac_1", "Transcript stable ID_y", "u_ac_2"
                FROM ppi_data 
                WHERE "Exon stable ID_y"=:exon_id AND "Transcript stable ID_y"=:ensemble_trans_id
                """
        p2 = pd.read_sql_query(sql=text(query), con=engine, params={'exon_id': exon_ID,
                                                                    'ensemble_trans_id': Ensemble_transID})
        '''

        # Compare the new and old dataframes
        # p1_old=PPI[ (PPI['Exon stable ID_x']==exon_ID)  &  (PPI['Transcript stable ID_x']==Ensemble_transID)].drop(columns=['Exon stable ID_x','Exon stable ID_y']).drop_duplicates()
        # p2_old=PPI[ (PPI['Exon stable ID_y']==exon_ID)  &  (PPI['Transcript stable ID_y']==Ensemble_transID)].drop(columns=['Exon stable ID_y','Exon stable ID_x']).drop_duplicates()
        # assert (np.array_equal(p1.values, p1_old.values))
        # assert (np.array_equal(p2.values, p2_old.values))

        p1 = tr_1[tr_1['Exon stable ID_x'] == exon_ID].drop(
            columns=['Exon stable ID_x', 'Exon stable ID_y']).drop_duplicates()
        p2 = tr_2[tr_2['Exon stable ID_y'] == exon_ID].drop(
            columns=['Exon stable ID_x', 'Exon stable ID_y']).drop_duplicates()

        p2 = p2[['Transcript stable ID_y', 'u_ac_2', 'Transcript stable ID_x', 'u_ac_1']]
        p2 = p2.rename(columns={
            'Transcript stable ID_y': 'Transcript stable ID_x',
            'u_ac_2': 'u_ac_1',
            'u_ac_1': 'u_ac_2',
            'Transcript stable ID_x': 'Transcript stable ID_y',
        })

        p1 = p1.append(p2, ignore_index=True).drop_duplicates()
        p1 = p1['u_ac_2'].unique()

        n = len(p1)
        N.append(n)
        if n != 0: exons_in_interface.append(exon_ID)
    return N, exons_in_interface, co_partners, text_edges


# if the input is a transcript ID:

def input_transcript(Ensemble_transID, organism):
    # check if the transcript have known domains:
    # ********** If not just visualize the exons and proteins interection

    # if the transcript have known domains:
    exons, domains, unique_domains = pr.transcript(Ensemble_transID, organism)

    output = pr.tranID_convert(Ensemble_transID, organism)
    if output == 0: return 0
    if any(x == False for x in output):
        return 0
    else:
        tran_name, gene_name, Ensemble_geneID, entrezID, gene_description = output

    text1 = "The transcript " + tran_name + ' has ' + str(exons.shape[0]) + ' exons and ' + str(
        len(unique_domains)) + " unique protein domains."
    # HTML code to visualize the table

    domains, Text_nodes, text_edges = expand_table(domains, unique_domains, entrezID, organism)
    domains = domains.sort_values(['Exon rank in transcript', 'Pfam start'], ascending=[True, True])

    # Link to visualize the network
    domains["Link to other databases"] = ""
    h = reverse('home') + "graph/"
    h2 = 'target="'
    h3 = '_blank"'
    # df_filter =domains['Interactions mediated by the domain']!=0

    domains["Link to other databases"] = ' <a href="https://www.ebi.ac.uk/interpro/entry/pfam/' + domains[
        'Pfam ID'] + '  "target="_blank">Pfam  </a>   &nbsp;&nbsp;&nbsp; <a href="https://3did.irbbarcelona.org/dispatch.php?type=domain&value=' + \
                                         domains['Pfam ID'] + '"target="_blank">3did  </a>      </h5 class> '

    # domains.at[df_filter,"Visualization of the domain interactions"]='<a target="'+'_blank"href="'+h+entrezID+"."+domains['Pfam ID']+'">'+gene_name+'-'+domains['Pfam ID']+'</a>'

    pd.set_option('display.max_colwidth', 1000)

    exon_info = exons[['Exon stable ID', 'Exon rank in transcript']]
    n, exons_in_interface, co_partners, co_partner_edges = exon_3D(exon_info['Exon stable ID'].tolist(),
                                                                   Ensemble_transID, organism)
    exon_info.loc[:, 'Number of interaction interface mapped to the exon'] = n

    # print(exons_in_interface)
    droped1 = exon_info.merge(domains, how='left', left_on='Exon stable ID', right_on='Exon stable ID')
    droped1 = droped1.rename(columns={"Pfam ID": "Corresponding domain ID", "Exon stable ID": "Exon  ID",
                                      "Exon rank in transcript_x": 'Exon rank in transcript'})

    droped1 = droped1.drop(columns=['CDS start', 'CDS end', 'Pfam start', 'Pfam end', "Link to other databases",
                                    "Interactions mediated by the domain", 'Transcript stable ID',
                                    'Chromosome/scaffold name', "Strand", "Genomic coding start", "Genomic coding end",
                                    'Exon rank in transcript_y'])

    droped2 = domains.drop(
        columns=['Exon rank in transcript', 'Exon stable ID', 'CDS start', 'CDS end', 'Pfam start', 'Pfam end',
                 'Transcript stable ID', 'Chromosome/scaffold name', "Strand", "Genomic coding start",
                 "Genomic coding end"])

    droped2 = droped2.drop_duplicates()

    # prepare the html table
    droped1["Corresponding domain ID"] = droped1["Corresponding domain ID"].fillna('-')

    h4 = reverse('home') + "ID/exon/" + organism + "/"

    droped1["Protein features encoded by the exon"] = '<a target="' + '_blank"href="' + h4 + droped1[
        "Exon  ID"] + '">Exon Page</a>'

    droped1.loc[droped1[
                    "Corresponding domain ID"] != '-', "Corresponding domain ID"] = '<a target="' + '_blank"href="' + 'https://www.ebi.ac.uk/interpro/entry/pfam/' + \
                                                                                    droped1[
                                                                                        "Corresponding domain ID"] + '">' + \
                                                                                    droped1[
                                                                                        "Corresponding domain ID"] + '</a>'

    droped1["Number of interaction interface mapped to the exon"] = droped1[
        "Number of interaction interface mapped to the exon"].astype(str)

    # add domain name column here

    droped2["Interactions mediated by the domain"] = droped2["Interactions mediated by the domain"].astype(int).astype(
        str)

    # s=[pr.Domain_name(x) for x  in droped1["Corresponding domain ID"].unique()]

    droped2["Symbol"], droped2["Summary"] = zip(*droped2['Pfam ID'].map(pr.Domain_name))

    droped2 = droped2[
        ['Pfam ID', 'Interactions mediated by the domain', 'Symbol', 'Summary', 'Link to other databases']]

    return (domains, unique_domains, exons, text1, domains.to_html(escape=False), Text_nodes, text_edges, tran_name,
            gene_name, Ensemble_geneID, entrezID, gene_description, exons, droped1.to_html(
        **settings.TO_HTML_PARAMETERS), droped2.to_html(**settings.TO_HTML_PARAMETERS), exons_in_interface, co_partners,
            co_partner_edges)
