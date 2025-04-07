from domain.Process import process_data as pr
from domain.Process import exonstodomain as exd
from domain.Process import proteininfo as info
import pandas as pd
import os
from django.urls import reverse
from sqlalchemy import text

from django.conf import settings

from domain.Process import network_analysis as nt
from domain.Process.load_data import DomainG_all

# --- Get database connection aka 'SQLAlchemie engine'
engine = settings.DATABASE_ENGINE



def TranscriptsID_to_table(transcripts, organism, entrez='0'):
    if len(transcripts) >= 1:
        # print('1111111111')
        ID = []
        name = []
        pfams = []
        missing_PPI = []

        DomainG = DomainG_all[organism]
        g2d = nt.g2d_all[organism]

        all_pfams = g2d[entrez]
        # print(transcripts)
        for tr in transcripts:
            query = """
                          SELECT * 
                          FROM exons_to_domains_data_""" + organism + """ 
                          WHERE "Transcript stable ID"=:transcript_id 
                          """
            tdata = pd.read_sql_query(sql=text(query), con=engine, params={'transcript_id': tr})

            # tdata=tdata.drop(columns=["Unnamed: 0"]).drop_duplicates()

            # df_filter = pr.data['Transcript stable ID'].isin([tr])
            # tdata=pr.data[df_filter]

            # print(tdata)
            if len(tdata) != 0:

                tmp = pr.tranID_convert(tr, organism)
                if tmp == 0: continue
                n = tmp[0]
                name.append(n)
                ID.append(tr)
                p = tdata["Pfam ID"].unique()
                p = p[~pd.isnull(p)]
                p = sorted(p)

                # look for interesting isoforms with missing PPI using the join graph
                missing = []
                missing = [x for x in all_pfams if x not in p]
                missing = [entrez + '/' + x for x in missing]
                missing_count = 0
                not_missing_count = len(all_pfams) - len(missing)
                for x in missing:
                    if DomainG.has_node(x):
                        # remove the second last element of the string
                        missing_count += 1
                for x in [f"{entrez}/{x}" for x in all_pfams if x in p]:
                    if DomainG.has_node(x):
                        # remove the second last element of the string
                        not_missing_count += 1
                total_count = missing_count + not_missing_count
                percent_retained = (not_missing_count / total_count)
                green_boxes = int(10 * percent_retained)
                missing_interaction = '<span class="text-success">█</span>' * green_boxes + ('<span '
                                                                    'class="text-danger">█</span>') * (10 - green_boxes)

                missing_PPI.append(missing_interaction)

                # add hyperlink
                p = [nt.link(x) for x in p]
                pfams.append(', '.join(p))

        if ID != []:
            pd_isoforms = pd.DataFrame(list(zip(name, ID, pfams, missing_PPI)),
                                       columns=['Transcript name', 'Transcript ID', 'Pfam domains',
                                                '<span class="text-success">Present</span> / <span class="text-danger">'
                                                'Missing</span> interacting domains in the isoform'])
            pd_isoforms['length'] = pd_isoforms['Pfam domains'].str.len()
            pd_isoforms.sort_values('length', ascending=False, inplace=True)
            pd_isoforms = pd_isoforms.drop(columns=['length'])

            h = reverse('home') + "ID/" + organism + "/"
            pd_isoforms["Link"] = '<a class="visualize" href="' + h + pd_isoforms["Transcript ID"] + '">' + " Visualize " + '</a>'
            pd_isoforms["Effect on Interactions"] = '<a class="visualize" href="' + h + pd_isoforms["Transcript ID"] + '#v-pills-interaction">' + " Show Interactions " + '</a>'

            pd.set_option('display.max_colwidth', 1000)

            pd_isoforms = pd_isoforms.to_html(**settings.TO_HTML_PARAMETERS)

    return pd_isoforms, n.split('-')[0]

    # changed


def input_gene(gene_ID, organism):
    # get a list of all transcripts of the selected gene
    pd_isoforms = []
    n = ''

    transcripts = pr.gene_to_all_transcripts(gene_ID, organism)
    entrez = str(int(nt.ensembl_to_entrez(gene_ID, organism)))

    if len(transcripts) == 0:
        return [], []

    pd_isoforms, n = TranscriptsID_to_table(transcripts, organism, entrez)

    return pd_isoforms, n
