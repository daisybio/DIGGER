# DIGGER
<p align="center">
  <img src="https://exbio.wzw.tum.de/digger/static/image/DIGGER.png" height="130">
</p>



Protein-protein interaction (PPI) networks are a key resource for systems biology. However, they do not consider the influence of alternative splicing, even though experimental evidence suggests that interaction partners are different for isoforms of the same protein. Domain Interaction Graph Guided ExploreR (DIGGER) integrates protein-protein interactions and domain-domain interactions into a joint graph and maps interacting residues to exons. DIGGER allows the users to query exons or isoforms individually or as a set to visually explore their interactions and it is available at: https://exbio.wzw.tum.de/digger




## Deploying DIGGER
To install DIGGER 2.0 (running it for the first time) follow these steps:
```shell script
# First, follow this link, if you want to install docker-compose: 
# https://docs.docker.com/compose/install/

#  Clone this repository and change into the created directory
git clone https://github.com/daisybio/DIGGER.git && cd DIGGER

# Download a copy of all the data files into domain/data/
# If you have more organisms, add the files into respective foldes, e.g. domain/data/Mus musculus[mouse]/
wget https://zenodo.org/records/14770005/files/data.zip

# unzip this data to the appropriate directory
unzip data.zip -d container/domain/

# repeat this with the nease data:
wget https://zenodo.org/records/14770005/files/nease_data.zip
unzip nease_data.zip -d container/domain/nease/


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

# Create the superuser for the admin page
# This can be found under https://WEBSITE/secure-admin and is used to manage the django models
docker-compose exec web python manage.py createsuperuser

# Enjoy your instance of DIGGER

```

## Extending DIGGER by DDIs for new organisms
>[!Warning]
>This is only necessary if you add a new organism and want to predict DDIs for it

To extend DIGGER using generated Domain-Domain Interactions (DDIs) for a new organism follow these steps:
````bash
# First, create a conda environment to ensure you have all dependencies installed
conda env create -f DIGGER_env.yml

# Activate the environment
conda activate DIGGER

# Make sure you have the necessary files for your new organism in the sourcedata folder
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


## Cite

If you use DIGGER, please cite:


Zakaria Louadi, Kevin Yuan, Alexander Gress, Olga Tsoy, Olga Kalinina, Jan Baumbach, Tim Kacprowski*, Markus List*. DIGGER: exploring the functional role of alternative splicing in protein interactions, Nucleic Acids Research, https://doi.org/10.1093/nar/gkaa768  (*joint last author)


## Development

First released by: Zakaria Louadi

Maintainers: Elias Albrecht / Konstantin Pelz


## Contact us
Elias Albrecht: elias.albrecht@tum.de
Konstantin Pelz: konstantin.pelz@tum.de
