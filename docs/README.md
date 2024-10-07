<p align="center">
  <img src="https://exbio.wzw.tum.de/digger/static/image/DIGGER.png" height="130">
</p>

Domain Interaction Graph Guided ExploreR (DIGGER) helps to investigate the impact of alternative splicing on protein-protein interactions (PPI) on the level of an exon, event, transcript, or the whole network. DIGGER is ideally suited to investigate the difference in interactions
between the isoforms, analyze the effect of isoform switch, or to explore how alternative splicing events such as exon skipping lead to altered
interactions of protein isoforms.   

# The concept of DIGGER

Experimental evidence suggests that the majority of protein isoforms have different interaction partners. However, conducting experiments for all possible isoforms will be time- and resource-consuming. Domain Interaction Graph Guided ExploreR (DIGGER) integrates protein-protein and domain-domain interactions (DDI) into a joint graph and maps interacting residues to exons  (Figure 1, **A**). The concept is the following: if an alternative splicing event leads to splicing out an interacting domain or residues, the interaction between proteins might be lost or impaired. 

<p align="center">
<img alt="DIGGER workflow" src="https://exbio.wzw.tum.de/digger/static/image/figure%201.png" width="800"/>
</p>
<p align="center">
  <em>Figure 1: The overview of DIGGER</em>
</p>

DIGGER uses the following sources for PPI, DDIs, and interacting residues:

- BioGRID (PPIs)
- 3did (DDIs)
- DOMINE (DDIs)
- PDB (interacting residues)

Domain-domain interaction predictions are based on the predictions of [PPIDM](10.1371/journal.pcbi.1008844).

DIGGER provides PPI visualizations in a structural context using a dynamic graph visualization that can be toggled between a protein
isoform and a domain-centric view (Figure 1, **B**). Furthermore, the tool maps the protein features encoded by a selected exon,
to judge the functional role of individual exons in the PPI (Figure 1, **C**).  DIGGER’s joint PPI and domain-domain interaction network can also be used for subnetwork extraction, providing a basis for network analysis (Figure 1, **D**).

Finally, DIGGER provides a web interface for **NEASE** (Network Enrichment method for Alternative Splicing Events) - a tool for the functional enrichment of alternative splicing exons/events.

# Step-by-step functional analysis of alternative splicing

## Isoform-level analysis

The isoform-level analysis investigates the difference in interactions between isoforms or the consequences of isoform switches. The documentation will use an example of BAG1 and NCK2 genes from the [DIGGER publication](https://doi.org/10.1093/nar/gkaa768). 

As input (Figure 2), isoform-level analysis accepts the name of a gene, an Ensebml identifier of a gene, and an Ensembl identifier of a transcript. A user can also put multiple identifiers separated by a comma. Then a user should select an organism. Currently, DIGGER contains the data from Homo sapiens and Mus musculus.

<img src="./Figure2.png">
<p align="center">
  <em>Figure 2: Input genes are BAG1 and NCK2; the chosen organism is Homo sapiens.</em>
</p>

If there is information about interactions for the isoforms of interest, the resulting window will show a summary of information about them in a network ('Network') and tabular views ('Queries'). A user can switch between them in the menu in the upper left corner.

A network view (Figure 3) provides information about known/predicted interactions of query isoforms. Predicted DDIs toggle bars allow the addition of predicted DDIs of different levels of confidence. Interactions lost or impaired due to splicing are marked as dashed lines. A user could also adjust the view of the network by enabling/disabling physics and adjusting the height of a canvas.
<img src="./Figure3.png">
<p align="center">
  <em>Figure 3: A network view for chosen isoforms of BAG1 and NCK2.</em>
</p>

The tabular view (Figure 4) presents the information about query isoforms or all annotated protein-coding isoforms if a query consists of a gene name/identifier. The Pfam domains column contains information about known PFAM domains of a resulting protein. "Present / Missing interacting domains in the isoform" column demonstrates the percentage of interactions possibly lost or impaired for this particular isoform due to alternative splicing. To further investigate the impact of a splicing event, a user can click 'Visualize'.

<img src="./Figure4.png">
<p align="center">
  <em>Figure 4: a tabular view on isoforms of BAG1 and NCK2.
  A transcript BAG-207 possibly lost a PPI due to alternative splicing that leads to removing the part of the protein that encodes a domain PF02179. A transcript NCK-202 also lost a PPI due to splicing of domains PF00017 and PF14604.</em>
</p>

A "Visualize" window, provides five different visualization panels: Overview, Protein View, Interaction View, Domain View, and Other Isoforms.

The Overview window header (Figure 5) shows the general information about the transcript.

<img src="./Figure5.png">
<p align="center">
  <em>Figure 5: General information about a transcript</em>
</p>

This page also contains a scheme of transcript exons and mapped Pfam domains, and the table of all exons and Pfam domains mapped on exons (Figures 6, 7). 

<img src="./Figure6A.png">
<img src="./Figure6B.png">
<p align="center">
  <em>Figure 6: Comparing BAG1-210 and BAG-207 isoforms, one can detect that these isoforms differ in exons that encode a domain PF02179 as well as exons that contain interacting residues</em>
</p>

<img src="./Figure7A.png">
<img src="./Figure7B.png">
<p align="center">
  <em>Figure 7: Comparing NCK-201 and NCK-202 isoforms, one can detect that these isoforms differ in exons that encode domains PF14604, the second PF00018, PF00017 as well as exons that contain interacting residues</em>
</p>

The Protein view shows all known direct interactors of an isoform. View mode can be changed from PPI to PPI-DDI view (Figure 8). 

<img src="./Figure8.png">
<p align="center">
  <em>Figure 8: PPI-DDI view for BAG1-207 depicts possibly missing interactions with BAG3, HSPA1A, HSPA4, and HSPA8</em>
</p>

Interaction View (Figure 9) allows investigation into details of the PPI between a query isoform and all known and predicted interactors. The interactors can be chosen from the **Select Interaction Partner Menu** (just tick the box with the interaction of interest there).

<p align="center">
  <img src="./Figure9.png" width="400" align="center">
</p>
<p align="center">
  <em>Figure 9: Interaction View of BAG1-207 isoform. An interaction with HSPA8 might be lost due to alternative splicing of PF02179</em>
</p>

Domain View (Figure 10) shows all interactors known for a particular domain.

<p align="center">
  <img src="./Figure10.png" width="400" align="center">
</p>
<p align="center">
  <em>Figure 10: Splicing out a domain PF02179 in the isoform BAG1-207 leads to missing PPIs with BAG3, HSPA1A, HSPA4, and HSPA8</em>
</p>

Other isoforms View (Figure 11) again depicts the table of all other known protein-coding isoforms of a gene.
<p align="center">
  <img src="./Figure11.png">
</p>
<p align="center">
  <em>Figure 11: All transcripts of BAG1</em>
</p>

## Exon-level analysis

The Exon-level analysis investigates the consequences of exon skipping events to interactions between isoforms.

In addition to input (Figure 12) as in isoform-level analysis, exon-level analysis accepts an Ensebml identifier of an exon, an Ensembl identifier of a gene with the coordinates of an exon. A user can also put multiple identifiers separated by a comma. Then a user should select an organism.

<p align="center">
  <img src="./Figure12.png">
</p>
<p align="center">
  <em>Figure 12: An input exon is exon ENSE00003569638 from BAG1; the organism is Homo sapiens</em>
</p>

Results of the exon analysis consist of DomainView, ProteinView, InteractionView, and the table of Proteins that use this exon. DomainView, ProteinView, and InteractionView are similar to Isoform-level analysis but for an exon of interest. DomainView additionally shows the information about the corresponding Pfam domain with a link to the Pfam and 3did database (Figure 13).

<p align="center">
  <img src="./Figure13.png">
</p>
<p align="center">
  <em>Figure 13: A DomainView for and exon ENSE00003569638 </em>
</p>

A table "Proteins that use this exon" (Figure 14) reports the list of all protein-coding isoforms that encodes this exon:

<p align="center">
  <img src="./Figure14.png">
</p>
<p align="center">
  <em>Figure 14: A table "Proteins that use this exon" for an exon ENSE00003569638 </em>
</p>


## Network-level analysis
Network analysis investigates the systematic impact of alternative splicing events on the PPIs. As input (Figure 15), Network-level analysis accepts the list of genes, transcripts, proteins, or a transcript count file and constructs a subnetwork from this list.
<p align="center">
  <img src="./Figure15.png">
</p>
<p align="center">
  <em>Figure 15: An example input for Network-level analysis </em>
</p>

If interactions between isoforms are known, the isoforms are connected to the subnetworks, and the interactions are marked as lost or present as well as based on the interaction evidence (Figure 16). For further exploration, the network can be downloaded in SIF or GraphML formats.
<p align="center">
  <img src="./Figure16.png">
</p>
<p align="center">
  <em>Figure 16: An example output for Network-level analysis </em>
</p>

# Running NEASE

**NEASE** (Network Enrichment method for Alternative Splicing Events) is a splicing-based functional network enrichment tool. Currently, NEASE supports the following formats:
- Standard (List of genes with exon coordinates)
- MAJIQ *deltapsi.tsv
- Whippet *.diff
- rMATS SE.MATS.JC.txt

The Standard input format contains three columns: gene ID, exon start coordinate, and exon end coordinate. E.g.:

**Gene ID** | **Start** | **End**  
--- | --- | --- 
ENSG00000033627 | 42508572	| 42508589 | 290 | 286 | 289 | 285 | 287 | 287 | 272 | 276 | 269

If you have a lot of files to run with NEASE, consider using the Python package: https://github.com/louadi/NEASE.

To start the analysis, upload your file (**A**), define an organism (**B**) (currently, Homo sapiens or Mus musculus), the input file type (**C**), and adjust the parameters.

For all input files except MAJIQ:
- Usage of predicted DDIs (**D**)
- Analysis name (**E**)
- P-value cutoff (**F**)
- Minimum delta PSI (Min delta, **G**)

For MAJIQ:
- Minimum delta PSI will start by default from 0.2 (Min delta, **G**)
- Confidence interval (0.95 corresponds to P-value 0.05, **H**)
<p align="center">
</p>

<p align="center">
  <img src="./FigureNEASE.png">
</p>
<p align="center">
  <em>Figure 17: NEASE starting page </em>
</p>

The results of the analysis will be stored for **7 days** and visible from the NEASE starting page (Previous Analyses). For large files, the running time might take up to several minutes.

The resulting 'Overview' page contains: 
- the general information about the number of genes with AS that affect any known protein feature;
- the table of affected protein domains;
- the table with affected interactions;
- the table with affected linear features;
- the table with affected interacting residues.

All tables are downloadable in a CSV format.

NEASE contains three further types of analysis: enrichment analysis, pathway analysis, and visualization.

NEASE edge enrichment supports the following ontologies:

#### Homo sapiens: 
- PharmGKB
- HumanCyc
- Wikipathways
- Reactome
- KEGG
  
#### Mus musculus
- KEGG
- Reactome
- MouseCyc

Type in the name of the ontology of interest and click the "Run enrichment" button. The resulting table will show the list of pathways, the list of genes that lost/impaired interactions with these pathways due to alternative splicing, and the significance assessment with p-value and adjusted p-value. NEASE score re-weights the pathways according to the hub bias. For more details, check the [NEASE publication](https://doi.org/10.1186/s13059-021-02538-1).

To investigate the enriched pathways in detail, use Pathway analysis. Type in the pathway ID from NEASE edge enrichment analysis, e.g., "path:mmu04010". The resulting table contains the information about each gene with list/impaired interactions with this pathway due to alternative splicing: whether it belongs to the pathway, how many interactions were lost/impaired, and the p-value.

Finally, "Visualize pathways" option provides the visualization of the resulting network with lost/impaired interactions. Be aware, the visualization might take up to several minutes. 

# Hosting your own instance

If you want, you can also host your own instance of DIGGER. To do so, follow the instructions below.

## Deploying DIGGER
To install DIGGER 2.0 (running it for the first time) follow these steps:
```shell script
# First, follow this link, if you want to install docker-compose: 
# https://docs.docker.com/compose/install/

#  Clone this repository and change into the created directory
git clone https://github.com/louadi/DIGGER.git && cd DIGGER

# Download a copy of all the data files into domain/data/Homo sapiens[human]/
# If you have more organisms, add the files into respective foldes, e.g. domain/data/Mus musculus[mouse]/
wget https://zenodo.org/record/3973368/files/data.zip    # extract the files


# Create a copy of the .env.sample file and edit the .env file
cp .env.sample .env   # now edit the .env file 

# Deploy and build the containers
docker-compose up -d --build

# Apply migrations (make sure the containers are up an running, else you will get an error)
docker-compose exec web python manage.py migrate --noinput 

# Import all the datasets into the database
docker-compose exec web python manage.py import_datasets

# Collect all the static files
docker-compose exec web python manage.py collectstatic --no-input

# Enjoy your instance of DIGGER

```


## Extending DIGGER by DDIs
If you downloaded the data provided by us, the extended DDIs are already included and you will **not** have to do 
this step. Only extend DIGGER like this if you want to add new organisms and have sufficient data to do so.  
To extend DIGGER using generated Domain-Domain Interactions (DDIs) follow these steps:
````bash
# First, create a conda environment to ensure you have all dependencies installed
conda env create -f DIGGER_env.yml

# Activate the environment
conda activate DIGGER

# Make sure you have the necessary files in the sourcedata folder
# more info about what these files should look like can be found in the sourcedata README.md

# copy the example database_sources to have a backup and list of all options. 
# Edit the database_sources.yml file to your needs
cp preprocess/sourcedata/example.database_sources.yml preprocess/sourcedata/database_sources.yml

# Run the prediction script
cd preprocess
python main_ddi_extend.py

# Once this has run, depending on your settings, restarting the DIGGER container will show the new data
cd ..
docker-compose up -d --force-recreate
````
# Frequently asked questions

**Question 1:** Could I input multiple isoforms?

*Answer 1:* Yes, you can input multiple isoform IDs, gene names, and gene IDs into the query field and even mix them. The first query page will then start by showing the interactions between input genes. To investigate each gene individually, use the "Queries" option in the upper left corner of the page.

**Question 2:** What are high, medium, and low confidence levels for DDI?

*Answer 2:* The confidence level is defined based on the prediction score from the PPIDM algorithm.

**Question 3:** I would like to use the output images in my manuscript. How could I do that?

*Answer 3:* You can save figures with the right click of the mouse -> 'Save Image as' option. Also, for NEASE plots, you can use a download button to download high-resolution figures.

**Question 4:** What is 'min delta' option in the NEASE analysis?

*Answer 4:* 'Min delta' refers to the deltaPSI value that defines the level of differential splicing between conditions. Here PSI value means 'Percentage Spliced In' and reflects the ratio between isoforms that use and splice out alternative exon. The calculation of PSI values varies from tool to tool, thus for more details consult the manuscript about the tools of interest. deltaPSI is simply the difference between the PSI value of an event in one condition and the PSI value of the same event in the other condition.
  The choice of the threshold depends on the study design. If only highly differentially spliced events should be considered, the threshold could be set to 0.1 or even 0.2. E.g., MAJIQ uses 0.2 as a default threshold.

**Question 5:** What is "MAJIQ confidence" option?

*Answer 5:* MAJIQ does not output p-values. Instead, for each alternative junction, the tool reports "Confidence Level" and 0.95 is suggested as a default threshold. For more details, check the MAJIQ publication.

**Question 6:** How long will my analysis be available?

*Answer 6:* By default, previous analyses will be available for 7 days. This period can be extended.

# Cite

If you use DIGGER, please cite:

Zakaria Louadi, Kevin Yuan, Alexander Gress, Olga Tsoy, Olga Kalinina, Jan Baumbach, Tim Kacprowski*, Markus List*. DIGGER: exploring the functional role of alternative splicing in protein interactions, Nucleic Acids Research, https://doi.org/10.1093/nar/gkaa768  (*joint last authors)

If you use NEASE, please cite:

Zakaria Louadi, Maria L. Elkjaer, Melissa Klug, Chit Tong Lio, Amit Fenn, Zsolt Illes, Dario Bongiovanni, Jan Baumbach, Tim Kacprowski, Markus List*, Olga Tsoy*. Functional enrichment of alternative splicing events with NEASE reveals insights into tissue identity and diseases, Genome Biology, https://doi.org/10.1186/s13059-021-02538-1 (*joint last authors)

# Contact us
Elias Albrecht: elias.albrecht@in.tum.de  
Olga Tsoy: olga.tsoy@uni-hamburg.de
