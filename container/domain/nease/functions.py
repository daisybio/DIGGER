import ast
import timeit

from .process import *
import scipy.stats as stats
import plotly.graph_objects as go
from .annotated_graph import filter_ppi_graph, ppi_interactions, filter_by_ddi, pathway_node_degree

import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
import numpy as np


def create_digger_link(row, organism):
    if row['Interacting domain']:
        organism_name = organism.lower()
        return f"<a href='{DIGGER}{organism_name}/{row['Exon stable ID']}' target='_blank'>{row['Exon stable ID']}</a>"
    else:
        return ''


def create_elm_link(row):
    return f"<a href='http://elm.eu.org/elms/{row['ELMIdentifier']}' target='_blank'>{row['ELMIdentifier']}</a>"


# main functions for nease
def exons_to_edges(mapped, G, elm_interactions, organism):
    # check if domains have known interactions/binding:
    mapped['domain'] = mapped['NCBI gene ID'] + '/' + mapped['Pfam ID']

    # domain in ddi
    mapped['DDI'] = mapped['domain'].apply(lambda x: G.has_node(x))

    # domaimn in elm
    mapped['elm'] = mapped['domain'].apply(lambda x: x in list(elm_interactions['interactor 2'].unique()))

    mapped['Interacting domain'] = mapped.apply(lambda x: bool(x['elm'] + x['DDI']), axis=1)

    def interactiontypye(ddi, elm):
        if (ddi & elm):
            return "DDI and DMI"
        elif ddi:
            return "DDI"
        elif elm:
            return "DMI"

    mapped['Interaction type'] = mapped.apply(lambda x: interactiontypye(x['DDI'], x['elm']), axis=1)

    # mapped['Interacting domain']=mapped['domain'].apply(lambda x: G.has_node(x))

    mapped = mapped.rename(columns={"max_change": "dPSI",
                                    "domain": "Domain ID"}).reset_index(drop=True)
    mapped['Visualization link'] = ''
    mapped['Visualization link'] = mapped.apply(lambda x: create_digger_link(x, organism), axis=1)
    return mapped


def affected_edges(nease_data, Join, only_DDIs):
    data = nease_data.data
    mapping = nease_data.mapping
    elm_interactions = nease_data.elm_interactions

    # get domains with DDIs
    interacting_domains = data[data['Interacting domain']]

    # Identify binding of affected domains = Edges in the PPI

    # helper function to search in ddi
    #retun nothing in case there is no node in joint graph
    def get_neighbors(graph, node):
        try:
            return graph.neighbors(node)
        except:
            return []

    def get_elm(node):

        interactors = list(elm_interactions[elm_interactions['interactor 2'] == node]['Interator gene 1'].unique())
        return [str(x) for x in interactors]

    # search in DDI and ELM
    t = lambda node: list(set([x.split('/')[0] for x in [n for n in get_neighbors(Join, node)]]
                              + get_elm(node)
                              ))

    interacting_domains['Affected binding (NCBI)'] = interacting_domains['Domain ID'].apply(t)
    interacting_domains = interacting_domains.rename(columns={"Pfam ID": "Identifier"}).reset_index(drop=True)

    # get edges from elm nad pdb
    if not only_DDIs:

        elm_affected = nease_data.elm_affected
        pdb_affected = nease_data.pdb

        # get elm edges
        t = lambda x: [str(x) for x in
                       list(elm_interactions[elm_interactions['interactor 1'] == x]['Interator gene 2'].unique())]
        elm_affected['Affected binding (NCBI)'] = elm_affected['ID'].apply(t)

        # only elm with interactions
        elm_affected = elm_affected[elm_affected['Affected binding (NCBI)'].map(lambda d: len(d)) > 0]

        # get elm edges

        if ~elm_affected.empty:
            elm_affected = elm_affected.rename(columns={"ELMIdentifier": "Identifier",
                                                        "entrezgene": "NCBI gene ID"}).drop(
                columns=['Gene stable ID', 'ID']).reset_index(drop=True)

            interacting_domains = pd.concat([interacting_domains, elm_affected], ignore_index=True)

    return interacting_domains


def gene_to_edges(data, pdb, only_DDIs):
    #For every gene get all edges
    gene_edges = {}
    for gene in data['NCBI gene ID'].unique():
        edges = data[data['NCBI gene ID'] == gene]['Affected binding (NCBI)']
        edges = [item for sublist in edges for item in sublist]

        gene_edges[gene] = list(set(edges))

    # get pdb interactions
    if not only_DDIs:

        for gene in pdb['NCBI gene ID'].unique():
            edges = pdb[pdb['NCBI gene ID'] == gene]['entrezgene']
            edges = [item for sublist in edges for item in sublist]

            if gene in gene_edges:

                gene_edges[gene] = list(set(gene_edges[gene] + edges))
            else:
                gene_edges[gene] = list(set(edges))

    return gene_edges


def pathway_enrichment(g2edges, paths, mapping, ppi_ddi_graph, organism, p_value_cutoff, only_DDIs):
    # General enrichment analysis
    pathway_genes = []
    # Totat degree of structural network for human (pre-computer)
    # For statistical test: edge enrichment

    #  every pathway degree two structural PPIs and
    # 'Degree in the structural PPI' : degree in ppi annotated with DDI.DMI/PDB
    # 'Degree in the PPI/DDI' : degree in ppi annotated with DDI only
    ppi_set = ppi_interactions(PPI[organism])

    if only_DDIs:
        filtered_ppis = filter_by_ddi(ppi_ddi_graph, ppi_set)
        n = len(filtered_ppis)
        ppi_type = 'Degree in the PPI/DDI'

    else:
        filtered_ppis = filter_ppi_graph(ppi_set, ppi_ddi_graph, elm_interactions[organism], pdb[organism], mapping)
        n = len(filtered_ppis)
        ppi_type = 'Degree in the structural PPI'

    # convert filtered ppi to networkx graph
    annotated_graph = nx.Graph()
    annotated_graph.add_edges_from(filtered_ppis)

    # number of effected edges
    affected_edges = len([item for sublist in g2edges.values() for item in sublist])

    # make a dictionary for faster lookup:
    mapping_dict = dict(zip(mapping['NCBI gene ID'], mapping['Gene name']))

    # for every path :
    path_name = []
    path_id = []
    source = []
    genes = []
    p_values = []
    score = []
    c = 0

    for path in list(paths['pathway']):
        # count of affected edges connected to the pathway
        # specific to that pathway list p
        # initialise the variable for every path 
        connected = 0
        genes_tmp = []
        gene_count = 0
        c += 1
        #if c == 10:
        #    break

        try:
            # get path total degree "p" and gene list
            path_genes = list(paths[paths['pathway'] == path]['entrez_gene_ids'])[0]
            p = pathway_node_degree(annotated_graph, path_genes)

        except Exception as e:
            print(e)
            path_genes = []
            pass

        for gene in g2edges:
            # for every affected gene
            # count affected gene edges connected to the
            # specific to the gene and to the pathway list
            tmp = len([x for x in g2edges[gene] if x in path_genes])

            if tmp > 0:

                # increment for path edges
                connected = connected + tmp

                # add gene to the gene list of the pathway
                # Entrez_to_name(gene, mapping_dict=mapping_dict) + " (" + str(tmp) + ")"
                genes_tmp.append(f"{Entrez_to_name(gene, mapping_dict=mapping_dict)} ({tmp})")

                # gene specific test
                _, p_gene = edge_enrich(tmp, len(g2edges[gene]) - tmp, p, n)

                if p_gene <= 0.05:
                    # gene with edges siginifically connected to the pathway

                    gene_count = gene_count + 1

        #  affected edges not connected to that pathway
        not_connected = affected_edges - connected

        # Join function is slow can be optimized later 
        if genes_tmp == []:
            genes_tmp = ''
        else:
            genes_tmp = (", ").join(genes_tmp)

        # Run hypergeometric test on affected edges
        _, p_value_tmp = edge_enrich(connected, not_connected, p, n)

        p_values.append(p_value_tmp)

        # compute combined score
        if p_value_tmp < p_value_cutoff:
            s = -(np.sqrt(gene_count) * np.log10(p_value_tmp))
        else:
            s = 0

        score.append(s)
        path_name.append(path)
        source.append(list(paths[paths['pathway'] == path]['source'])[0])
        path_id.append(list(paths[paths['pathway'] == path]['external_id'])[0])
        genes.append(genes_tmp)

    # save results
    Enrichment = pd.DataFrame(list(zip(path_id, path_name, source, genes, p_values, score)),
                              columns=['Pathway ID', 'Pathway name', 'Source',
                                       'Spliced genes (number of interactions affecting the pathway)', "p_value",
                                       'Nease score'])

    return Enrichment.sort_values(['p_value'], ascending=True)


def single_path_enrich(path_id, Pathways, g2edges, mapping, organism, only_DDIs, entrez_name_map, ppi_ddi_network):
    # Totat degree of structural network for human (pre-computer)
    # For statistical test: edge enrichment
    ppi_set = ppi_interactions(PPI[organism])

    if only_DDIs:
        filtered_ppis = filter_by_ddi(ppi_ddi_network, ppi_set)
        n = len(filtered_ppis)
        ppi_type = 'Degree in the PPI/DDI'

    else:
        filtered_ppis = filter_ppi_graph(ppi_set, ppi_ddi_network,
                                         elm_interactions[organism], pdb[organism], mapping)
        n = len(filtered_ppis)
        ppi_type = 'Degree in the structural PPI'

    # convert filtered ppi to networkx graph
    annotated_graph = nx.Graph()
    annotated_graph.add_edges_from(filtered_ppis)

    path_genes = list(Pathways[Pathways['external_id'] == path_id]['entrez_gene_ids'])[0]
    p = pathway_node_degree(annotated_graph, path_genes)
    # convert path_genes to str
    if isinstance(path_genes, str):
        path_genes_str = ast.literal_eval(path_genes)
    else:
        path_genes_str = path_genes

    #collect:
    spliced_genes = []
    spliced_genes_entrez = []
    gene_association = []
    num = []
    affected_edges = []
    affected_edges_entrez = []
    p_val = []

    # graph to save affected edges
    G = nx.Graph()

    for g in g2edges:

        # affected edges of the gene g
        affected = g2edges[g]

        # edges connected to the pathway
        edges = [x for x in affected if x in path_genes_str]
        a = len(edges)

        # Not connected
        b = len(affected) - a

        if a != 0:
            # calculate gene specific p_value:
            _, p_value = edge_enrich(a, b, p, n)

            # Save results
            spliced_genes.append(Entrez_to_name(g, mapping_dict=entrez_name_map))
            spliced_genes_entrez.append(g)
            gene_association.append(g in path_genes_str)
            num.append(f"{a}/{a + b} ({a / (a + b) * 100:.2f}%)")
            affected_edges.append(', '.join([Entrez_to_name(x, mapping_dict=entrez_name_map) for x in edges]))
            affected_edges_entrez.append(', '.join(edges))
            p_val.append(p_value)

            # save affected edges
            G.add_edges_from([(g, x) for x in edges])

    Enrichment = pd.DataFrame(list(
        zip(spliced_genes, spliced_genes_entrez, gene_association, num, p_val, affected_edges, affected_edges_entrez)),
        columns=['Spliced genes', 'NCBI gene ID', 'Gene is known to be in the pathway',
                 'Edges associated with the pathway (%)', 'p_value',
                 'Affected binding (edges)', 'Affected binding (NCBI)'])

    return Enrichment, G


# stats functions
# Statistical test
# Edge enrichment test

def edge_enrich(a, b, p, n):
    # function to calculate P value
    # test if affected edges by AS are significally enriched in a pathway
    # fisher exact test
    # a+b affected edges, with a the one linked to pathway p
    # p total degree of pathway p
    # n total edges in the ppi used.

    # linked to pathway but not affected
    c = p - a
    if c < 0:
        return np.nan, 1.0
    # background of test: not linked to p and not affected edges
    d = (2 * n) - p - b

    # retun oddsratio and pvalue from fisher exact test
    try:
        return stats.fisher_exact([[a, b], [c, d]], alternative='greater')
    except ValueError:
        print("Wrong numbers in fisher exact test")
        return np.nan, 1.0


# Visualization Function

def extract_subnetwork(path_genes,
                       ppi,
                       affected_genes,
                       all_spliced_genes,
                       k,
                       mapping,
                       affected_graph,
                       significant,
                       entrez_name_map,
                       organism):
    # Affected_genes: genes with lost/gained interaction in the pathway
    # all_spliced_genes: all genes affected with splicing
    ensembl_map = mapping.set_index('Gene stable ID')['NCBI gene ID'].to_dict()
    translated = set()
    missing = set()
    missing_flag = False
    for gene in all_spliced_genes:
        if gene in ensembl_map:
            translated.add(str(ensembl_map.get(gene, gene)))
        else:
            missing.add(gene)
    # look online once for missing genes
    online_result = Ensemb_to_entrez(missing, organism)
    print(f"Fetched {len(online_result)} of {len(missing)} missing genes online")
    for gene in missing:
        if gene in online_result:
            translated.add(str(ensembl_map.get(gene, gene)))
        else:
            missing_flag = True

    all_spliced_genes = translated

    # Extract the pathway module for the complete PPI
    # We would like to visualize the pathway with affected edges:
    G = ppi.subgraph(path_genes + affected_genes)
    G = nx.Graph(G)
    # add connections from affected_graph
    G.add_edges_from(affected_graph.edges())
    G.remove_nodes_from(list(nx.isolates(G)))

    # Position nodes using Fruchterman-Reingold force-directed algorithm.
    pos = nx.spring_layout(G, k=k, iterations=100)
    # pos = nx.kamada_kawai_layout(G)

    # Prepare the visualization
    for n, p in pos.items():
        G.nodes[n]['pos'] = p

    # Define node and edges in plott
    node_trace = go.Scatter(x=[],
                            y=[],
                            text=[],
                            mode='markers+text',
                            hoverinfo='text',
                            textposition='top center',
                            #textfont_size = 22,
                            marker=dict(
                                reversescale=True,
                                color=[],
                                size=[],
                                line=dict(width=0)))

    colored_node_trace = go.Scatter(x=[],
                                    y=[],
                                    text=[],
                                    mode='markers+text',
                                    hoverinfo='text',
                                    textposition='top center',
                                    #textfont_size = 22,
                                    marker=dict(
                                        reversescale=True,
                                        color=[],
                                        size=[],
                                        line=dict(width=0)))

    edge_trace = go.Scatter(x=[],
                            y=[],
                            mode='lines',
                            line=dict(width=0.7,
                                      color='#888'),
                            hoverinfo='none')

    colored_edge_trace = go.Scatter(x=[],
                                    y=[],
                                    mode='lines',
                                    line=dict(width=4,
                                              color='red'),
                                    hoverinfo='none')

    nodes_with_red_edges = set()
    for edge in G.edges():
        x0, y0 = G.nodes[edge[0]]['pos']
        x1, y1 = G.nodes[edge[1]]['pos']

        # Check if edge is affected
        if affected_graph.has_edge(*edge):
            colored_edge_trace['x'] += tuple([x0, x1, None])
            colored_edge_trace['y'] += tuple([y0, y1, None])
            nodes_with_red_edges.add(edge[0])
            nodes_with_red_edges.add(edge[1])

        else:
            edge_trace['x'] += tuple([x0, x1, None])
            edge_trace['y'] += tuple([y0, y1, None])

    for node in G.nodes():
        x, y = G.nodes[node]['pos']
        node_info = Entrez_to_name(node, mapping_dict=entrez_name_map)

        if node not in path_genes:
            # Node not in pathway
            color = 'orange'
            if node in significant:
                size = 50
            else:
                size = 30

        elif str(node) in all_spliced_genes:
            # spliced and part of the pathway
            color = 'red'
            if node in significant:
                size = 50
            else:
                size = 30

        else:
            # other pathway nodes
            color = '#888'
            if node in nodes_with_red_edges:
                color = '#888887'
            if node in significant:
                size = 50
            else:
                size = 20

        trace_type = colored_node_trace if color != '#888' or size > 20 else node_trace

        trace_type['x'] += tuple([x])
        trace_type['y'] += tuple([y])
        trace_type['text'] += tuple([node_info])

        trace_type['marker']['color'] += tuple([color])
        trace_type['marker']['size'] += tuple([size])

    return [colored_edge_trace, node_trace, edge_trace, colored_node_trace], missing_flag


# plot domains stats

def stats_domains(affecting_percentage,
                  number_of_features,
                  domain_number,
                  elm_number,
                  pdb_number,
                  file_path):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 7))
    fig.subplots_adjust(wspace=0)

    # pie chart parameters
    ratios_pie = [affecting_percentage, 1 - affecting_percentage]
    labels_pie = ['Affecting protein features', 'Not affecting any feature']
    pie_colors = ['#377DB0', '#E68A00']

    ax1.set_facecolor('white')
    ax1.pie(ratios_pie, labels=labels_pie, autopct='%1.1f%%', startangle=0, wedgeprops={'edgecolor': 'white'},
            colors=pie_colors, textprops={'fontsize': 12, 'color': 'black'})  # Set font size and color for percentages
    ax1.set_title("Genes with AS affecting protein features", fontsize=14)

    ratios_bar = [round(elm_number / number_of_features, 2) * 100, round(pdb_number / number_of_features, 2) * 100,
                  round(domain_number / number_of_features, 2) * 100]

    labels_bar = ['Linear motifs', 'Residues', 'Domains']

    colors = ['#344552', '#377DB0', '#92C0DD'][::-1]

    ax2.set_facecolor('white')
    ax2.barh(labels_bar, ratios_bar, color=colors)
    ax2.set_xlim(0, 100)
    ax2.set_title('Affected features', fontsize=14)

    # Remove unnecessary spines and ticks, while keeping the clean look
    ax2.tick_params(bottom=False)
    ax2.tick_params(left=False)
    ax2.spines[['left', 'bottom']].set_color('#D2D2D2')
    ax2.spines[['left', 'bottom']].set_linewidth(2)
    ax2.spines['right'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.grid(False)
    ax2.yaxis.set_ticks_position('none')
    ax2.xaxis.set_visible(False)

    for i, v in enumerate(ratios_bar):
        ax2.text(v + 1, i, f"{v:.1f}%", va='center', color='black', fontweight='bold')

    plt.tight_layout()

    plt.savefig(f"{file_path}_thumb.jpg", format='jpg', bbox_inches='tight')
    plt.savefig(f"{file_path}.jpg", format='jpg', bbox_inches='tight', dpi=1200)
    plt.clf()
    plt.close()
    # Save the figure

    return


def get_node_depth(G, node, n=0):
    pred = list(G.predecessors(node))
    if not pred:  # No predecessors
        return n
    return get_node_depth(G, pred[0], n + 1)


def interpolate_line(x0, y0, x1, y1, num_points=20):
    """Interpolate points between two coordinates for smoother hover behavior."""
    x_vals = np.linspace(x0, x1, num_points)
    y_vals = np.linspace(y0, y1, num_points)
    return x_vals, y_vals


# create network for the pathways and their connected genes
def all_pathway_network(enrichment_table: pd.DataFrame, pathways: nx.DiGraph, k=0.8, db_name=None):
    if len(enrichment_table) < 2:
        return None

    pathway_map = {}
    edge_weights = []
    node_source = {}
    for i, row in enrichment_table.iterrows():
        id = row['Pathway ID']
        name = row['Pathway name']
        genes = row['Spliced genes (number of interactions affecting the pathway)'].split(',')
        genes = {x.split('(')[0].strip() for x in genes}
        node_source[id] = row['Source']
        pathway_map[id] = name, genes

    G = nx.Graph()
    for pathway_s in pathway_map:
        for pathway_d in pathway_map:
            if (pathway_s == pathway_d or G.has_edge(pathway_s, pathway_d) or
                    not pathways.has_node(pathway_s) or not pathways.has_node(pathway_d)):
                continue
            # check if pathway_s is a predecessor of pathway_d in the pathways
            parent_pathway_s = list(pathways.in_edges(pathway_s))
            parent_pathway_d = list(pathways.in_edges(pathway_d))
            if (not nx.has_path(pathways, pathway_s, pathway_d) and
                    ((len(parent_pathway_s) == 0 or len(parent_pathway_d) == 0) or
                     parent_pathway_s[0][0] != parent_pathway_d[0][0])):
                continue
            _, genes_s = pathway_map[pathway_s]
            _, genes_d = pathway_map[pathway_d]
            overlap = genes_s.intersection(genes_d)
            if len(overlap) > 0:
                G.add_edge(pathway_s, pathway_d, weight=len(overlap))
                edge_weights.append(len(overlap))

    print(f"Graph has {len(G.nodes)} nodes and {len(G.edges)} edges")
    if len(G.nodes) == 0:
        return None

    # prepare for visualization
    top_pathways = [p for p in enrichment_table.head(8)['Pathway ID'].tolist() if p in G.nodes]
    edge_weights.sort()
    depths = {n: get_node_depth(pathways, n) for n in G.nodes}
    max_depth = max(depths.values(), default=0)
    num_dbs = len(db_name.split(','))
    iterations = 80 if num_dbs == 1 else 80 * num_dbs

    # position nodes using Fruchterman-Reingold force-directed algorithm.
    if db_name == 'KEGG':
        pos = nx.spring_layout(G, k=k, iterations=100, fixed=[top_pathways[0]], pos={top_pathways[0]: (0, 0)})
    else:
        pos = nx.spring_layout(G, k=k, iterations=iterations)
    for n, p in pos.items():
        G.nodes[n]['pos'] = p

    # create nodes and edges. Edges are weighted by the number of shared genes
    node_trace = go.Scatter(x=[],
                            y=[],
                            text=[],
                            mode='markers+text',
                            hoverinfo='text',
                            textposition='top center',
                            opacity=1,
                            marker=dict(
                                reversescale=True,
                                color=[],
                                size=[],
                                opacity=1,
                                line=dict(width=0)))

    # convert to dict for individual line widths
    edge_traces = {}
    for weight in edge_weights:
        edge_traces[weight] = go.Scatter(x=[], y=[], text=f"Shared genes: {weight}",
                                         mode='lines',
                                         hoverinfo='text',
                                         opacity=0.5,
                                         line=dict(width=min(((weight * 0.5) ** 2) + 1, 50), color='#B0C4DE'))

    for node in G.nodes():
        x, y = G.nodes[node]['pos']
        node_info = pathway_map[node][0]
        if 'Homo sapiens (human)' in node_info:
            node_info = node_info.replace(' - Homo sapiens (human)', '')
        elif 'Mus musculus (mouse)' in node_info:
            node_info = node_info.replace(' - Mus musculus (mouse)', '')

        if num_dbs > 1:
            node_info = f"{node_info} ({node_source[node]})"

        if node in top_pathways:
            node_info = f"<b>{node_info}</b>"

        node_trace['x'] += tuple([x])
        node_trace['y'] += tuple([y])
        node_trace['text'] += tuple([node_info])
        level = depths[node]
        base_size = 30 if max_depth < 5 else 10

        node_trace['marker']['color'] += tuple(['#98bdd7' if level != 0 else '#4298d6'])
        node_trace['marker']['size'] += tuple([base_size + (max_depth - level) * 5])

    for edge in G.edges(data=True):
        weight = edge[2]['weight']
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]

        # this is for a better hover experience
        x_vals, y_vals = interpolate_line(x0, y0, x1, y1, num_points=20)

        edge_traces[weight]['x'] += tuple(x_vals) + (None,)
        edge_traces[weight]['y'] += tuple(y_vals) + (None,)

    return [*edge_traces.values(), node_trace]
