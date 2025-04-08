import os
import pickle
import random
import re
import time
import timeit
import traceback

import pandas as pd

from django.db import connection
from django.conf import settings
from django.http import HttpResponse, JsonResponse
from django.shortcuts import render, redirect
from django.urls import reverse
from django.utils.html import escape
from io import StringIO

from pandas.errors import ParserError

from .Process import exonstodomain as exd
from .Process import exon as ex
from .Process import process_data as pr
from .Process import transcript as tr
from .Process import gene as g
from .Process import network_analysis as nt
from .Process import mutliple_query as mq
from .Process import process_data as proc_data
from .Process import nease_output as no
from .models import NeaseSavedRun

# --- Create folder
# Global jobs path
jobs_path = os.path.join(settings.MEDIA_ROOT, 'jobs')
if not os.path.exists(jobs_path):
    os.makedirs(jobs_path)


# Display transcripts of a gene
def gene(request, gene_ID, organism):
    print("Currently in gene view")
    transcript_table, gene_name = g.input_gene(gene_ID, organism)

    if transcript_table == []:
        return render(request, 'visualization/no_info.html', {'name': gene_name,
                                                              'info': 'Without any known Pfam domains'})

    context = {
        'tb': transcript_table,
        'name': gene_name,
    }
    return render(request, 'visualization/gene.html', context)


# if the user inputs multiple queries, this function will handle them
def multiple_queries(request, inputs, organism):
    print("Currently in multiple queries view")
    # retrieve get data from the form
    input_names = [x.strip() for x in inputs.split(",")]
    print(input_names)
    transcript_table = {}

    for query in input_names:
        try:
            if re.match(r'ENS\w*E\d+$', query):
                print("Checking exon")
                transcript_table[query] = ex.exon_table(query, organism)
            elif re.match(r'ENS\w*T\d+$', query):
                transcript_table[query] = tr.transcript_table(query, organism)
            elif re.match(r'ENS\w*P\d+$', query):
                transcript_table[query] = tr.transcript_table(query, organism, True)
            else:
                with connection.cursor() as cursor:
                    cursor.execute(
                        """SELECT ensembl_id FROM domain_gene_""" + organism + """ WHERE gene_symbol ILIKE %s""",
                        [query])
                    row = cursor.fetchone()
                    query_id = row[0] if row else query
                    transcript_table[query] = g.input_gene(query_id, organism)[0]
        except Exception as e:
            print(e)
            continue

    matches = 0
    not_found = []
    filtered_table = {}
    for query in transcript_table.keys():
        if len(transcript_table[query]) > 0:
            matches += 1
            filtered_table[query] = transcript_table[query]
        else:
            # remove the query from the transcript table
            not_found.append(query.strip())

    transcript_table = filtered_table

    if matches == 0:
        origin = request.GET.get('origin', None)
        original_input = "~".join(inputs.split(","))
        return redirect(f"{reverse(origin)}?search={original_input}&organism={organism}")

    # iterate through all query ids
    id_list = []
    for query in transcript_table.keys():
        if len(transcript_table[query]) == 0:
            continue
        # search for the id in the transcript table string. This is suboptimal but it works for now
        co_partners = []
        find_co_partners = True
        try:
            transcript_ids = re.findall(r'ENS\w*[T,P]\d+', transcript_table[query])
            for t_id in transcript_ids:
                tran_name, _, _, entrez_id, _ = proc_data.tranID_convert(t_id, organism)
                exons, _, unique_domains = proc_data.transcript(t_id, organism)
                if find_co_partners:
                    _, _, co_partners, _ = exd.exon_3D(exons['Exon stable ID'].tolist(), t_id, organism)
                    find_co_partners = False
                id_list.append([tran_name, entrez_id, unique_domains, co_partners])
        except AttributeError:
            print("Attribute error for: ", query)
            pass
    # for each id, get the graph data
    combined_nodes = {}
    combined_edges = {}
    try:
        # convert id_list ids to entrez ids
        subgraph_g, confirmed_proteins, missing_domains = mq.create_subgraph(id_list, organism)

        nodes, edges = mq.vis_nodes_many(subgraph_g, id_list, confirmed_proteins, missing_domains)
        combined_nodes = {k: combined_nodes.get(k, []) + v for k, v in nodes.items()}
        combined_edges = {k: combined_edges.get(k, []) + v for k, v in edges.items()}
    except Exception as e:
        print(e)
        pass

    context = {
        'show_network': len(combined_nodes['original']) > 0,
        'trans_table': transcript_table,
        'combined_nodes': combined_nodes,
        'combined_edges': combined_edges,
        'switch': [],
        'not_found': not_found,
    }

    return render(request, 'visualization/multiple_queries.html', context)


# Display information of an exon (Exon-Level Analysis)
def exon(request, organism, exon_ID):
    print("Currently in exon view")
    v = ex.input_exon(exon_ID, organism)

    if v is None:
        return render(request, 'visualization/no_info.html', {'name': exon_ID,
                                                              'info': 'Wrong entry or exon in a gene without any known '
                                                                      'Pfam domains'})
    else:
        _, domains, gene_name, Ensemble_geneID, entrezID, tb_transc, table_domains, number = v

    # only if the exon code for domains with known interactions
    nodes_domainV = {}
    edges_domainV = {}
    switcher = []
    switcher_js = []
    first = []
    maxx = 0

    # Interactionview
    Interactiveview_select = {}
    first_victim = []
    first_name = ""

    if number > 0:
        # ProteinView
        nodes, edges, pd_interaction = ex.vis_exon(domains, entrezID, gene_name, exon_ID, organism)

        # DomainView
        first = domains[0]
        for pfams in domains:
            n, e, _, _ = exd.vis_node_(entrezID + "." + pfams, organism)
            if len(e['original']) > maxx:
                maxx = len(e)
                first = pfams
            if len([val for sublist in e.values() for val in sublist]) != 0:
                # initialise keys with empty list if they don't exist
                for key in n.keys():
                    if key not in nodes_domainV:
                        nodes_domainV[key] = []
                for key in e.keys():
                    if key not in edges_domainV:
                        edges_domainV[key] = []

                nodes_domainV = {k: nodes_domainV[k] + v for k, v in n.items()}
                edges_domainV = {k: edges_domainV[k] + v for k, v in e.items()}
                switcher.append('<option value="' + pfams + '"> ' + pfams + '</option>')
                switcher_js.append('case "' + pfams + '": return node.source === "' + pfams + '";')

    else:
        nodes, edges, pd_interaction = {}, {}, []

    # PPI res interfaces on the exon:
    # table: HTML table with all PPIs that have res interface mapped to the exon
    # number_of_PPI: number of interactions

    table, number_of_PPI = ex.PPI_inter(exon_ID, gene_name, organism)

    if number > 0 and len(pd_interaction) > 0:
        # added to combine evidence of DDI and Residue in one final table
        if number_of_PPI > 0:
            ppi_from_res = table['Partner Protein'].unique()
            f = pd_interaction['Partner Protein'].isin(ppi_from_res)
            pd_interaction.loc[f, 'Residue evidence'] = '<center>&#9989;</center>'

        # InteractionView
        pd_interaction['_'] = pd_interaction["Partner Protein"]

        # This also affects the pd_interaction dataframe, watch out
        Interactiveview_select = pr.interactive_select(pd_interaction)

        # the first protein to show, select a protein with confidence original
        try:
            first_victim = pd_interaction[pd_interaction['Confidence'] == 'Original']['NCBI gene ID'].tolist()[0]
            first_name = pd_interaction[pd_interaction['Confidence'] == 'Original']['_'].tolist()[0]
        except IndexError:
            # if there is no original confidence, select the first protein
            first_victim = pd_interaction['NCBI gene ID'].tolist()[0]
            first_name = pd_interaction['_'].tolist()[0]

        pd_interaction = pd_interaction[
            ["Affected Protein", 'Partner Protein', 'NCBI gene ID', 'Retained DDIs', 'Lost DDIs',
             'Percentage of lost domain-domain interactions', 'Residue evidence', "Protein-protein interaction",
             'Confidence']]

        pd_interaction = pd_interaction.rename(columns={

            "Percentage of lost domain-domain interactions": "% of missing DDIs",
            "Retained DDIs": "Retained domain-domain interactions",
            "Lost DDIs": "Missing domain-domain interactions",
            "Protein-protein interaction": "Protein-protein interaction",
            'Residue evidence': 'Residue-level evidence*',
        })

        pd_interaction = pd_interaction.to_html(table_id='Interaction_table', **settings.TO_HTML_RESPONSIVE_PARAMETERS)

    table = table.to_html(**settings.TO_HTML_PARAMETERS)

    context = {

        'tb1': tb_transc,
        'tb2': table_domains,
        'tb3': pd_interaction,
        'tb4': table,
        'name': gene_name,
        'exon_ID': exon_ID,
        'entrezID': entrezID,
        'gID': Ensemble_geneID,
        'dis': number > 0,

        "dis2": number == -1,
        'dis3': number_of_PPI != 0,

        # only a self loop for the domains>> no interactionView
        'dis4': number > 0 and len(pd_interaction) == 0,
        'long_table': number_of_PPI > 25,
        'pv_nodes': nodes,
        'pv_edges': edges,

        'first_domain': first,
        'switch1': switcher,
        'switch2': switcher_js,
        'Domainview_edges': edges_domainV,
        'Domainview_nodes': nodes_domainV,

        'Interactiveview_select': Interactiveview_select,
        'first_vict': first_victim,
        'first_name': first_name,

        'enable_Proteinview': len(edges_domainV.get('original')) < 70 if edges_domainV.get('original') else True,
    }
    return render(request, 'visualization/exon.html', context)


# Display information of a transcript or a protein (Isoform-Level Analysis)
def transcript(request, P_id, organism):
    print("Currently in transcript view")

    out = tr.Protein_view(P_id, organism)

    if out == 0: return render(request, 'visualization/no_info.html',
                               {'name': P_id, 'info': 'Wrong entry or protein without any known Pfam domains'})
    if out == 1: return render(request, 'visualization/no_info.html',
                               {'name': P_id, 'info': 'The selected protein does not have any interaction in the '
                                                      'current PPI database'})

    nodes, edges, unique, text1, tran_name, Ensemble_geneID, entrezID, gene_description, droped1, droped2, trID, p, \
        missed, pd_interaction, isoforms, co_partners = out

    # Interactionview
    Interactiveview_select = {}
    first_victim = []
    first_name = []

    if len(pd_interaction) != 0:
        pd_interaction['Residue evidence'] = ''

        pd_interaction.loc[
            pd_interaction["NCBI gene ID"].isin(co_partners), 'Residue evidence'] = '★'

        pd_interaction = pd_interaction.sort_values('Protein name')
        pd_interaction['_'] = pd_interaction['Protein name'] + " " + pd_interaction['Residue evidence']

        Interactiveview_select = pr.interactive_select(pd_interaction)

        # the first protein to show, select a protein with confidence original, alternatively the first protein
        try:
            first_victim = pd_interaction[pd_interaction['Confidence'] == 'Original']['NCBI gene ID'].tolist()[0]
            first_name = pd_interaction[pd_interaction['Confidence'] == 'Original']['_'].tolist()[0]
        except IndexError:
            first_victim = pd_interaction['NCBI gene ID'].tolist()[0]
            first_name = pd_interaction['_'].tolist()[0]

        pd_interaction = pd_interaction.rename(columns={

            "Protein name": "Partner Protein",
            'Residue evidence': 'Residue-level evidence',
            "Percentage of lost domain-domain interactions": "% of missing DDIs",
            "Retained DDIs": "Retained domain-domain interactions",
            "Lost DDIs": "Missing domain-domain interactions",
            "Protein-protein interaction": "Protein-protein interaction",
        })

        pd_interaction = pd_interaction[
            ["Selected Protein variant", 'Partner Protein', 'NCBI gene ID', 'Retained domain-domain interactions',
             'Missing domain-domain interactions', '% of missing DDIs', 'Residue-level evidence',
             "Protein-protein interaction", 'Confidence']]
        pd_interaction = pd_interaction.to_html(table_id='Interaction_table', **settings.TO_HTML_RESPONSIVE_PARAMETERS)

    # Get ID of missing domains with interactions
    if len(missed) != 0:
        missing_domains = missed['Pfam ID'].unique()
        missed = missed.to_html(**settings.TO_HTML_PARAMETERS)

    nodes_domainV = {}
    edges_domainV = {}
    switcher = []
    switcher_js = []
    first = unique[0]
    maxx = 0

    # DomainView for retained domains
    for pfams in unique:
        n, e, _, _ = exd.vis_node_(entrezID + "." + pfams, organism)
        if len(e['original']) > maxx:
            maxx = len(e)
            first = pfams
        if len([val for sublist in e.values() for val in sublist]) != 0:
            # initialise keys with empty list if they don't exist
            for key in n.keys():
                if key not in nodes_domainV:
                    nodes_domainV[key] = []
            for key in e.keys():
                if key not in edges_domainV:
                    edges_domainV[key] = []

            nodes_domainV = {k: nodes_domainV[k] + v for k, v in n.items()}
            edges_domainV = {k: edges_domainV[k] + v for k, v in e.items()}

            switcher.append('<option value="' + pfams + '"> ' + pfams + '</option>')
            switcher_js.append('case "' + pfams + '": return node.source === "' + pfams + '";')

    # DomainView for missing domains
    switcher_m = []
    if len(missed) != 0:
        for pfams in missing_domains:
            n, e, _, _ = exd.vis_node_(entrezID + "." + pfams, organism, missing=True)
            if len(e['original']) > maxx:
                maxx = len(e)
                first = pfams
            if len([val for sublist in e.values() for val in sublist]) != 0:
                # initialise keys with empty list if they don't exist
                for key in n.keys():
                    if key not in nodes_domainV:
                        nodes_domainV[key] = []
                for key in e.keys():
                    if key not in edges_domainV:
                        edges_domainV[key] = []

                nodes_domainV = {k: nodes_domainV[k] + v for k, v in n.items()}
                edges_domainV = {k: edges_domainV[k] + v for k, v in e.items()}
                switcher_m.append('<option value="' + pfams + '"> ' + pfams + ' (missing in the isoform) </option>')
                switcher_js.append('case "' + pfams + '": return node.source === "' + pfams + '";')

    # get total length of edges_domainV to determine disable_Proteinview
    len_edges_domainV = len([val for sublist in edges_domainV.values() for val in sublist])

    context = {
        'dt': droped1,
        'text1': text1,
        'tran_name': tran_name,
        'gene_description': gene_description,
        'trID': trID,
        'gID': Ensemble_geneID,
        'path': p,
        'pv_nodes': nodes,
        'pv_edges': edges,
        'entrezID': entrezID,
        'dt2': droped2,
        'dt3': missed,
        'dt4': pd_interaction,
        "dt5": isoforms,

        'dis1': missed != [],
        'dis2': pd_interaction != [],

        'dis3': isoforms != [],

        'Interactiveview_select': Interactiveview_select,
        'first_vict': first_victim,
        'first_name': first_name,

        'first_domain': first,
        'switch1': switcher,
        'switch1_missing': switcher_m,
        'switch2': switcher_js,
        'Domainview_edges': edges_domainV,
        'Domainview_nodes': nodes_domainV,

        # make it make sense
        'enable_Proteinview': (len_edges_domainV < 200 or len(unique) + len(missed) == 1),
    }

    return render(request, 'visualization/transcript.html', context)


# isoform_level search box
def isoform_level(request):
    if "search" in request.GET:  # If the form is submitted
        # Get and sanitize the search_query
        search_query = request.GET['search'].strip()

        # Get organism (human, mouse)
        organism = request.GET.get('organism', None)

        # allow user to input multiple queries separated by commas
        multiple = search_query.split(",")
        if len(multiple) > 1:
            url = reverse('multiple-queries', kwargs={'inputs': search_query, 'organism': organism})
            return redirect(f"{url}?origin=isoform-level")

        # Try and parse the search_query as gene name from the database
        with connection.cursor() as cursor:
            cursor.execute("""SELECT ensembl_id FROM domain_gene_""" + organism + """ WHERE gene_symbol ILIKE %s""",
                           [search_query])
            row = cursor.fetchone()

        if row:
            search_query = row[0]

        search_query = search_query.split("+")[0]
        search_query = search_query.split("%")[0]
        search_query = search_query.split(".")[0]

        # regex to check if search query starts with ENS and ends with T or P
        if re.match(r'ENS\w*[T,P]\d+$', search_query):
            print("User input is a protein")
            return redirect(transcript, P_id=search_query, organism=organism)

        elif re.match(r'ENS\w*G\d+$', search_query):
            print("User input is a gene")
            return redirect(gene, gene_ID=search_query, organism=organism)

        search = request.GET['search'].strip()
        if "~" in search:
            search = ", ".join(search.split("~"))
        msg = {"error": f"No database entries contain {search}."}

    else:
        msg = {}

    return render(request, 'setup/isoform_level.html', msg)


# exon_level search box
def exon_level(request):
    msg = {"error": "No entries found in our database for your search query."}
    if "search" in request.GET:  # If the form is submitted
        # Input and Exon ID
        print('-----------------------------------------------------------')
        search_query = request.GET['search']
        organism = request.GET.get('organism', None)
        search_query = search_query.replace(" ", "")
        search_query = search_query.split("+")[0]
        search_query = search_query.split("%")[0]
        search_query = search_query.split(".")[0]
        print(search_query)

        # allow user to input multiple queries separated by commas
        multiple = search_query.split(",")
        if len(multiple) > 1:
            url = reverse('multiple-queries', kwargs={'inputs': search_query, 'organism': organism})
            return redirect(f"{url}?origin=exon-level")

        if re.match(r'ENS\w*[E]\d+$', search_query):
            return redirect(exon, organism=organism, exon_ID=search_query)
    elif "search 2" in request.GET:  # If the form is submitted
        # Input coordinate of the exon
        # Check if coordinate are correct
        #  Example   ' ENSG00000266028  206437964 206437042 '

        print('-----------------------------------------------------------')
        search_query = request.GET['search 2']
        organism = request.GET.get('organism', None)

        multiple = search_query.split(",")
        if len(multiple) > 1:
            search_query = multiple[0]

        search_query = search_query.split(" ")
        search_query = [x for x in search_query if x != '']
        # search_query[0]=search_query[0].split(".")[0]
        if len(search_query) == 3 and re.match(r'ENS\w*G\d+$', search_query[0]) and \
                search_query[1].isdigit() and search_query[2].isdigit():

            gene_ID = search_query[0]
            s1 = int(search_query[1])
            e1 = int(search_query[2])

            # Correct for very big inputs
            if abs(s1 - e1) < 3000:

                # map coordinates to exon
                exonID = pr.coordinate_to_exonID(gene_ID, s1, e1, organism)

                if exonID != []:
                    return redirect(exon, organism=organism, exon_ID=exonID)
                    # return exon(request,exonID)
                else:
                    return render(request, 'visualization/no_info.html', {'name': gene_ID,
                                                                          'info': 'The exon is not in the database'})
    elif "search 3" in request.GET:  # If option 3 is selected

        # Get and sanitize the search_query
        search_query = request.GET['search 3'].strip()
        organism = request.GET.get('organism', None)

        # allow user to input multiple queries separated by commas
        multiple = search_query.split(",")
        if len(multiple) > 1:
            url = reverse('multiple-queries', kwargs={'inputs': search_query, 'organism': organism})
            return redirect(f"{url}?origin=exon-level")

        with connection.cursor() as cursor:
            cursor.execute("""SELECT ensembl_id FROM domain_gene_""" + organism + """ WHERE gene_symbol ILIKE %s""",
                           [search_query])
            row = cursor.fetchone()

        if row:
            search_query = row[0]

        search_query = search_query.split("+")[0]
        search_query = search_query.split("%")[0]
        search_query = search_query.split(".")[0]

        # Input search is a protein:
        if re.match(r'ENS\w*[TP]\d+$', search_query):
            return redirect(transcript, organism=organism, P_id=search_query)

        # Input search is a gene:
        elif re.match(r'ENS\w*G\d+$', search_query):
            return redirect(gene, organism=organism, gene_ID=search_query)

    else:
        msg = {}

    return render(request, 'setup/exon_level.html', msg)


# PPI network analysis
def network(request):
    error_message = ""
    jump_div = ""

    # Option 1: List of Ensembl IDs
    if "option1" in request.POST:
        organism = request.POST.get('organism', None)
        input_query = []
        for element in request.POST['input'].split('\n'):
            element = element.strip()
            if element:
                input_query.append(element)

        input_query = list(set(input_query))

        # max input IDs
        if 2000 > len(input_query) > 1:
            if re.match(r'^ENS\w*[GTP]', input_query[0]):
                job_num = str(random.randrange(500))
                with open(f'{jobs_path}/{job_num}.txt', "wb") as fp:
                    pickle.dump(input_query, fp)
                return redirect(Multi_proteins, organism=organism, job=job_num)
        else:
            error_message = f"The input list must contain between 2 and 2000 Ensembl IDs."
            jump_div = 'option1'

    # Option 2: Upload file
    if "option2" in request.POST and 'gene-count-file' in request.FILES:
        error_message_suffix = ""
        organism = request.POST.get('organism', None)
        try:
            # --- Check input file for correct format
            # Try to decode as UTF-8, sanitize and parse as table
            try:
                file_string = escape(request.FILES['gene-count-file'].read().decode('UTF-8'))
                file_buffer = StringIO(file_string)
                # Parse as pandas dataframe
                transcript_count_df = pd.read_table(file_buffer)
            except UnicodeDecodeError:
                error_message_suffix = "could not be parsed as an text file"
                raise RuntimeError

            except ParserError:
                error_message_suffix = f"could not be parsed as an table file (CSV or TSV)"
                raise RuntimeError

            # Check input shape
            if transcript_count_df.shape[0] < 2 or transcript_count_df.shape[1] < 2:
                error_message_suffix = f"could not be parsed as table or has less than two rows and columns"
                raise RuntimeError

            # Kevin: Zakaria please insert the magic down below:)
            # Zaka: And this is where the magic happens :p
            # Elias: lets hope this magic still works after I changed the format checking

            # Check if the first row corresponds to transcript Ensembl IDs
            if not (re.match(r'^ENS\w*T', transcript_count_df.iloc[0, 0]) or re.match(r'^ENS\w*T',
                                                                                      transcript_count_df.iloc[1, 0])):
                error_message_suffix = f"must have Ensembl transcript IDs in the first column starting with \"ENST\""
                raise RuntimeError

            # --- Try parsing counts for the different options (search for FPKM, tpm or counts)
            # max_isoforms: the max number of isoforms to consider:
            max_isoforms = int(request.POST['transcript-count-max'])

            column_names = transcript_count_df.columns

            # Cufflinks file (or a similar thing)
            if "FPKM" in column_names:
                transcript_count_df = transcript_count_df.sort_values(by=['FPKM'], ascending=False)
                cut_rows = transcript_count_df.iloc[:, 0].unique()[:max_isoforms]
                print('Input matches cufflinks output ')

            # Kallisto output counts in tpm
            elif "tpm" in column_names:
                transcript_count_df = transcript_count_df.sort_values(by=['tpm'], ascending=False)
                cut_rows = transcript_count_df.iloc[:, 0].unique()[:max_isoforms]
                print('Input matches kallisto output ')

            # Generic count matrix
            elif "counts" in column_names:
                transcript_count_df = transcript_count_df.sort_values(by=['counts'], ascending=False)
                cut_rows = transcript_count_df.iloc[:, 0].unique()[:max_isoforms]
                print('Input with counts column ')

            # Could not find the row
            else:
                error_message_suffix = "does not contain a column with the counts. The column must be named either \"FPKM\", \"tpm\" or \"counts\""
                raise RuntimeError

            # and let DIGGER do the magic ;)
            job_num = str(random.randrange(500))
            with open(f'{jobs_path}/{job_num}.txt', "wb") as fp:
                pickle.dump(cut_rows, fp)  # Pickling
                print(f"Starting network analysis with {len(cut_rows)} rows")
            return redirect(Multi_proteins, organism=organism, job=job_num)

        except RuntimeError:
            print("Could not parse uploaded file correctly")
            error_message = f"The uploaded file \"{request.FILES['gene-count-file']}\" {error_message_suffix}."
            jump_div = 'option2'

    return render(request, 'setup/network.html', context={
        'error_message': error_message,
        'jump_div': jump_div
    })


# network level analysis
def Multi_proteins(request, organism, job='0'):
    with open(f'{jobs_path}/{job}.txt', "rb") as fp:  # Unpickling
        inputs = pickle.load(fp)

    if re.match(r'^ENS\w*G', inputs[0]):
        info = nt.analysis_input_genes(inputs, organism)

    elif re.match(r'^ENS\w*[TP]', inputs[0]):
        info = nt.analysis_input_isoforms(inputs, organism)
    else:
        return render(request, 'visualization/no_info.html', {'name': inputs[0],
                                                              'info': 'Wrong entry or protein without any known Pfam '
                                                                      'domains'})
    if not info:
        return render(request, 'visualization/no_info.html', {'name': 'your job',
                                                              'info': 'Too many inputs (max=2000 genes)'})

    else:
        print("info exists")
        genes, missing, num_isoforms = info

        Net = nt.Construct_network(genes, missing, job, organism)

        if Net == 0:
            return render(request, 'visualization/no_info.html', {'name': 'your job',
                                                                  'info': 'There is no known interaction between '
                                                                          'the provided proteins'})

        else:
            nodes, edges, tab, tb_html = Net

    context = {

        'pv_nodes': nodes,
        'pv_edges': edges,
        "tab": tb_html,
        'ID': job,
        "genes_number": len(missing),
        "isoforms_num": num_isoforms,
        'interacted_nodes': len(nodes),

    }

    return render(request, 'visualization/network.html', context)


def set_previous_analysis(request, post_request=True):
    if post_request:
        run_id = request.POST.get('previousAnalysis')
    else:
        run_id = request.GET.get('runId')
    print("got previous analysis with run ID:", run_id)
    try:
        events, info_tables = no.get_nease_events(run_id)

    except FileNotFoundError:
        context = {'error_msg': "Could not find this analysis, please run it again."}
        return render(request, 'setup/nease_setup.html', context)

    except Exception as e:
        context = {'error_msg': str(e)}
        return render(request, 'setup/nease_setup.html', context)

    for key, value in info_tables.items():
        info_tables[key] = value.to_html(table_id=f"{key}_table", **settings.TO_HTML_RESPONSIVE_PARAMETERS)

    save_info = NeaseSavedRun.objects.get(run_id=run_id)
    current_duration = save_info.saved_for_days
    time_left = save_info.days_left()

    context = {
        'input_name': save_info.file_name,
        'custom_name': save_info.get_custom_name(),
        **events.summary,
        **info_tables,
        'stats': f"{run_id}_thumb.jpg",
        'stats_full': f"{run_id}.jpg",
        'shareable_link': request.build_absolute_uri(reverse('nease-analysis')) + "?runId=" + run_id,
        'run_id': run_id,
        'current_duration': current_duration,
        'time_left': time_left,
        'organism': save_info.organism,
        'input_format': save_info.input_format,
        'predicted_ddis': save_info.predicted_DDIs,
        'p_value_cutoff': save_info.p_value_cutoff,
        'min_delta': save_info.min_delta,
        'majiq_confidence': save_info.majiq_confidence,
        'only_ddis': save_info.only_ddis,
        'remove_not_in_frame': save_info.remove_not_in_frame,
        'only_divisible_by_3': save_info.only_divisible_by_three,
        **events.get_databases()
    }
    try:
        if request.POST.get('example', None) == 'true':
            return render(request, 'visualization/nease_result_example.html', context)

        return render(request, 'visualization/nease_result.html', context)
    except Exception as e:
        print(e)
        traceback.print_exc()
        return render(request, 'setup/nease_setup.html', context)


# this does the initial nease run or loads a previous analysis
def setup_nease(request):
    # handle previous analysis
    if request.POST.get('previousAnalysis', None):
        return set_previous_analysis(request, True)

    # handle previous analysis as a get request
    if request.GET.get('runId', None):
        return set_previous_analysis(request, False)

    # otherwise continue with new analysis
    if not request.FILES:
        return render(request, 'setup/nease_setup.html')

    # reject files that are too big
    MAX_FILE_SIZE_MB = 100
    if request.FILES['splicing-events-file'].size > MAX_FILE_SIZE_MB * 1_000_000:
        return HttpResponse("File is too big. Please upload a file smaller than 100 MB.", status=400)

    # Get the input file and post data
    input_data = request.FILES
    if 'splicing-events-file' not in input_data:
        return HttpResponse("No input file provided", status=400)

    organism = request.POST.get('organism', 'human')
    database_type = request.POST.get('inputType', 'Standard')
    custom_name = request.POST.get('analysis_name', None)

    confidences = []
    for confidence in ["high", "mid", "low"]:
        if request.POST.get(f"predicted-checkbox-{confidence}", False):
            confidences.append(confidence)

    # get form data in the correct format
    enrich_dbs = request.POST.getlist('databases_to_enrich')
    p_value = float(request.POST.get('p_value_cutoff', 0.05))
    min_delta = float(request.POST.get('min_delta', 0.05))
    majiq_confidence = float(request.POST.get('Majiq_confidence', 0.95))
    # Check if values are in POST Request (They are set to True by default), if not, set them False
    # We are just checking if the value is in the request, not the actual value (this is the same as done for the
    # confidences). So our default values are in the user interface, but not here.
    only_ddis = request.POST.get('only_DDIs', False)
    rm_not_in_frame = request.POST.get('remove_non_in_frame', False)
    divisible_by_3 = request.POST.get('only_divisible_by_3', False)

    # Run the NEASE job
    print(f"{time.asctime(time.localtime(time.time()))} Submitted NEASE job with params: {organism}, {database_type}, "
          f"{p_value}, {rm_not_in_frame}, {divisible_by_3}, {min_delta}, {majiq_confidence}, {only_ddis}, {confidences}")
    context = {
        'error_msg': None
    }

    try:
        if input_data['splicing-events-file'].name.endswith('.json'):
            table = pd.read_json(input_data['splicing-events-file'])
        else:
            if not no.file_needs_cleaning(input_data['splicing-events-file']):
                table = pd.read_table(input_data['splicing-events-file'])
            else:
                table = no.read_extra_spaces(input_data['splicing-events-file'])
    except Exception as e:
        print(e)
        table = None

    try:
        if not isinstance(table, pd.DataFrame) or table.empty:
            raise ValueError("Could not parse the input file. Please make sure the file is in the correct format.")

        events, info_tables, run_id = no.run_nease(table, organism, {'db_type': database_type,
                                                                     'enrich_dbs': enrich_dbs,
                                                                     'p_value': p_value,
                                                                     'rm_not_in_frame': rm_not_in_frame,
                                                                     'divisible_by_3': divisible_by_3,
                                                                     'min_delta': min_delta,
                                                                     'majiq_confidence': majiq_confidence,
                                                                     'only_ddis': only_ddis,
                                                                     'confidences': confidences},
                                                   input_data['splicing-events-file'].name,
                                                   custom_name)

        for key, value in info_tables.items():
            info_tables[key] = value.to_html(table_id=f"{key}_table", **settings.TO_HTML_RESPONSIVE_PARAMETERS)

        save_info = NeaseSavedRun.objects.get(run_id=run_id)
        current_duration = save_info.saved_for_days
        time_left = save_info.days_left()

        context = {
            'input_name': input_data['splicing-events-file'].name,
            'custom_name': custom_name,
            **events.summary,
            **info_tables,
            'stats': run_id + "_thumb.jpg",
            'stats_full': run_id + ".jpg",
            'shareable_link': request.build_absolute_uri(reverse('nease-analysis')) + "?runId=" + run_id,
            'run_id': run_id,
            'current_duration': current_duration,
            'time_left': time_left,
            'organism': save_info.organism,
            'input_format': save_info.input_format,
            'predicted_ddis': save_info.predicted_DDIs,
            'p_value_cutoff': save_info.p_value_cutoff,
            'min_delta': save_info.min_delta,
            'majiq_confidence': save_info.majiq_confidence,
            'only_ddis': save_info.only_ddis,
            'remove_not_in_frame': save_info.remove_not_in_frame,
            'only_divisible_by_3': save_info.only_divisible_by_three,
            **events.get_databases(),
        }
        # return the example page that includes explanations for the user
        if request.POST.get('example', None) == 'true':
            return render(request, 'visualization/nease_result_example.html', context)

        return render(request, 'visualization/nease_result.html', context)

    except Exception as e:
        traceback.print_exc()
        # check for the filename and see if there is an indication that the input format is wrong
        possible_true_format = no.match_name_with_format(input_data['splicing-events-file'].name)
        error_out = str(e)
        if "format" in error_out and possible_true_format:
            if database_type == 'Standard':
                error_out = "Invalid standard format. Are you sure you selected the correct input format?"
            error_out += f" Try to use the {possible_true_format}!"
        elif "format" in error_out and not possible_true_format:
            error_out += f" Try to use the Standard input!"
        context['error_msg'] = error_out
    return render(request, 'setup/nease_setup.html', context)


# extra functions for the NEASE output once the analysis is done
def nease_extra_functions(request):
    function_name = request.GET.get('func', None)
    if not function_name:
        return HttpResponse("No function provided", status=400)
    run_id = request.GET.get('runId', None)
    if not run_id:
        return HttpResponse("No run ID provided", status=400)
    databases = request.GET.get('databases', None)
    pathway = request.GET.get('pathway', None)
    k = request.GET.get('k', None)

    if databases:
        databases = databases.split(",")
        databases = list({db.strip() for db in databases})

    # get the enrichment table
    try:
        print("Got function: ", function_name)
        table_name = "enrich"
        if function_name == 'classic':
            out_table = no.nease_classic_enrich(no.get_nease_events(run_id), databases, run_id)
        elif function_name == 'nease':
            out_table = no.nease_enrichment(no.get_nease_events(run_id), databases, run_id)
        elif function_name == 'pathway':
            out_table = no.pathway_info(no.get_nease_events(run_id), pathway, run_id)
            table_name = "path"
        elif function_name == 'visualise':
            out_table = no.visualise_path(no.get_nease_events(run_id), pathway, k)
        elif function_name == 'save':
            duration = request.GET.get('duration', None)
            if not duration:
                return HttpResponse("No duration provided", status=400)
            out_table = no.change_save_timing(run_id, duration)
        else:
            return HttpResponse(f"Unknown function: {function_name}", status=400)
    except Exception as e:
        traceback.print_exc()
        return HttpResponse(f"Error: {str(e)}", status=500)

    if isinstance(out_table, pd.DataFrame):
        return HttpResponse(out_table.to_html(table_id=f"{function_name}_{table_name}",
                                              **settings.TO_HTML_RESPONSIVE_PARAMETERS))
    elif isinstance(out_table, str):
        return HttpResponse(out_table)
    elif isinstance(out_table, tuple):
        return JsonResponse({'table': out_table[0].to_html(table_id=f"{function_name}_{table_name}",
                                                           **settings.TO_HTML_RESPONSIVE_PARAMETERS),
                             'plot': out_table[1]}, status=200)
    else:
        return JsonResponse(out_table, status=400)


# retrieve the list of organisms for the organism selection dropdown based on the available data folders
def get_organisms(request):
    directory_path = os.path.join(settings.PROJECT_ROOT, 'domain/data')

    if os.path.exists(directory_path) and os.path.isdir(directory_path):
        organisms = os.listdir(directory_path)
        true_organisms = []
        trivial_names = []
        for orga in organisms:
            if not os.path.isdir(os.path.join(directory_path, orga)):
                continue
            true_organisms.append(orga.split("[")[0])
            trivial_names.append(orga.split("[")[1][:-1])
        return JsonResponse({'organisms': true_organisms, 'trivial_names': trivial_names})
    else:
        return JsonResponse({'organisms': [], 'trivial_names': []})
