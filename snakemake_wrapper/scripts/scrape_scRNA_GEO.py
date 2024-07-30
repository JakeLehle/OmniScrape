#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 09:35:27 2024

@author: jlehle
"""
# %%


import pandas as pd
from bs4 import BeautifulSoup
import xml.etree.ElementTree as ET
from random import randint
from random import seed
import time
import json
import csv
import re




import tarfile

#%%
import os, sys
#os.chdir(os.environ['HOME']+'/SC/fastq')
os.chdir(os.environ['HOME']+'/OmniScrape/tmp')
#%%
#######
# THIS IS WHERE THINGs STARTED WORKING #
#######
#%%
import pandas as pd
pd.options.display.max_colwidth = 400
#X = pd.read_fwf('/master/jlehle/OmniScrape/tmp/gds_result.txt', header=None)
X = pd.read_fwf('/master/jlehle/OmniScrape/tmp/Smart-seq-GEO.txt', header=None)
ftp_links = X[X[0].str.startswith('FTP download')]
ftp_links = ftp_links.reset_index(drop=True)
ftp_links = ftp_links.loc[0:,0]
ftp_links
# %% scrape the ftp from the gds file
ftp_list = []
# loop through the rows using iterrows()
for index, row in ftp_links.items():
    try:
        tmp = row.split(') ', -1)[1]
        ftp_list.append(tmp)
    except:
        pass


# %% check the ftp links
# import module
import requests
# create a function
# pass the url

def url_ok(foo_url):
    foo_url = 'https' + foo_url[3:] + 'soft/'
    # pass the url into
    # request.hear
    response = requests.head(foo_url)
    # check the status code
    if response.status_code == 200:
        return True
    else:
        return False


# driven code
status = []
for test_url in ftp_list:
    status.append(url_ok(test_url))


def countX(lst, x):
    count = 0
    for ele in lst:
        if (ele == x):
            count = count + 1
    return count


#%%
print(f'''Of the initiaimport wgetl {len(ftp_links)} datasets provided from the GEO search,
      {len(ftp_list)} point to GEO FTP pages.
      Of these {len(ftp_list)} FTP pages {countX(status, True)} returned 200 GET statuses for having GEO *.soft files.''')

#%%
from itertools import compress
ftp_list = list(compress(ftp_list, status))
len(ftp_list)
#%%
import wget
import GEOparse

gse_dict = {}
for url in ftp_list:
    wget_this = url + 'soft/' + url[:-1].split('/', -1)[-1] + '_family.soft.gz'
    filename = wget.download(wget_this)
    gse = GEOparse.get_GEO(filepath=filename)
    gse_dict[url[:-1].split('/', -1)[-1]] = gse.phenotype_data


#%%
num_rows = []
for key in gse_dict.keys():
    num_rows.append(len(gse_dict[key]))

#%%
print(f'''From the {countX(status, True)} FTP pages that returned 200 GET 
      statuses for having GEO *.soft files. 
      There are an available {sum(num_rows)} samples that are available for download.''')

#%%
#Start filtering the 25,202 samples

#Figure out some way to write everything you have collected so fart into table this will be the first target file of the script. 
I will need to have GSE the 

#%%

sample = []
organism = []
library_se = []
library_st = []
gse_dict_copy = gse_dict
for key in gse_dict_copy.keys():
    gse = GEOparse.get_GEO(filepath= key+"_family.soft.gz")
    for index, gsm in gse.phenotype_data['geo_accession'].items():
        sample.append(gse.phenotype_data.loc[gsm]['geo_accession'])
        organism.append(gse.phenotype_data.loc[gsm]['organism_ch1'])
        library_se.append(gse.phenotype_data.loc[gsm]['library_selection'])
        library_st.append(gse.phenotype_data.loc[gsm]['library_strategy'])
    
#%%  
# function to get unique values
def unique(list1):
    unique_list = pd.Series(list1).drop_duplicates().tolist()
    for x in unique_list:
        print(x)
#%%
organism
unique(organism)
unique(library_se)
unique(library_st)
type(organims)
#%%
#######
# THIS IS WHERE I PUT THE CODE THAT I'M WORKING ON #
#######

#%%
#Mohadeseh needs more samples so I'm gonna set that up really fast

wget_this = 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE160nnn/GSE160269/soft/GSE160269_family.soft.gz'
filename = wget.download(wget_this)
gse = GEOparse.get_GEO(filepath=filename)
gse.phenotype_data
SRA_file = pd.read_csv("test.tsv", sep="\t")
SRA_file.loc[0]
SRA_file.run_accession
print(SRA_file.experiment_title)
SRA_file.biosample

data = {'SRR': SRA_file.run_accession,
        'GSM': SRA_file.experiment_title}


#%% Bladder cancer datasets
import wget
import GEOparse

wget_this = 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE222nnn/GSE222315/soft/GSE222315_family.soft.gz'
filename = wget.download(wget_this)
gse = GEOparse.get_GEO(filepath=filename)

data = {'Sample': ["P1_BCa_scRNAseq", "P1_BCa_scRNAseq", "P1_BCa_scRNAseq", "P1_BCa_scRNAseq", \
"P2_NAT_scRNAseq", "P2_NAT_scRNAseq", "P2_NAT_scRNAseq", "P2_NAT_scRNAseq", \
"P2_BCa_scRNAseq", "P2_BCa_scRNAseq", "P2_BCa_scRNAseq", "P2_BCa_scRNAseq", \
"P3_NAT_scRNAseq", "P3_NAT_scRNAseq", "P3_NAT_scRNAseq", "P3_NAT_scRNAseq", \
"P3_BCa_scRNAseq", "P3_BCa_scRNAseq", "P3_BCa_scRNAseq", "P3_BCa_scRNAseq", \
"P4_NAT_scRNAseq", "P4_NAT_scRNAseq", "P4_NAT_scRNAseq", "P4_NAT_scRNAseq", \
"P4_BCa_scRNAseq", "P4_BCa_scRNAseq", "P4_BCa_scRNAseq", "P4_BCa_scRNAseq", \
"P5_NAT_scRNAseq", "P5_NAT_scRNAseq", "P5_NAT_scRNAseq", "P5_NAT_scRNAseq", \
"P5_BCa_scRNAseq", "P5_BCa_scRNAseq", "P5_BCa_scRNAseq", "P5_BCa_scRNAseq", \
"P6_BCa_scRNAseq", \
"P7_BCa_scRNAseq", \
"P8_BCa_scRNAseq", \
"P9_BCa_scRNAseq"],
        'Condition': np.repeat(gse.phenotype_data.source_name_ch1, [4, 4, 4, 4, 4, 4, 4, 4, 4, 1, 1, 1, 1]).tolist(),
        'Source': np.repeat('Bladder', [40]).tolist(),
        'Cancer Stage': np.repeat(gse.phenotype_data['characteristics_ch1.3.tumor stage'],  [4, 4, 4, 4, 4, 4, 4, 4, 4, 1, 1, 1, 1]).tolist(),
        'Sex': np.repeat(gse.phenotype_data['characteristics_ch1.1.Sex'], [4, 4, 4, 4, 4, 4, 4, 4, 4, 1, 1, 1, 1]).tolist(),
        'Age': np.repeat(gse.phenotype_data['characteristics_ch1.2.age'], [4, 4, 4, 4, 4, 4, 4, 4, 4, 1, 1, 1, 1]).tolist()       
        }
df = pd.DataFrame(data)
#%%
len(data['Sample'])
len(data['Source'])
df = pd.DataFrame(data)
df.Sample.tolist()
#%%    
    
    
SRA_file = pd.read_csv("test.tsv", sep="\t")
SRA_file.loc[0]
SRA_file.run_accession
print(SRA_file.experiment_title)
SRA_file.biosample

data = {'SRR': SRA_file.run_accession,
        'GSM': SRA_file.experiment_title}
#%%
new_string_final = []
for index, string in data['GSM'].items():
    new_string = string.split(' 10X', 1)[0]
    new_string_final.append(new_string)
#%%    
data['Sample'] = new_string_final    
data['SRR'].to_string()
    
gse_dict[url[:-1].split('/', -1)[-1]] = gse.phenotype_data

#%%
#the sample
len(gse.phenotype_data['geo_accession'])
gse.phenotype_data.loc['GSM3327702']['organism_ch1']
#organism is
gse.phenotype_data['organism_ch1']
library info
gse.phenotype_data['library_selection']
gse.phenotype_data['library_strategy']


# %%
gse_dict.keys()
gse = GEOparse.get_GEO(filepath="GSE118389_family.soft.gz")
type(gse.phenotype_data)
gse.phenotype_data
gse.gpls
gse.gsms
gse.gsms['GSM4909314'].show_metadata()
gse.gsms['GSM6129415'].get_metadata_attribute('Sample_title')
gse.gsms['GSM4909314'].download_SRA(
    'jlehle@txbiomed.org', directory=GSE_UID, nproc=1,)
gse.gsms['GSM4909314'].download_supplementary_files()

pd.set_option('display.max_columns', None)
print(gse.phenotype_data)
print(gse.phenotype_data['geo_accession'])

for accession in gse.phenotype_data['geo_accession']; :
    gse.gsms[accession].download_SRA(
        'jlehle@txbiomed.org', directory=gse.phenotype_data['series_id'], nproc=1)
print(gse.phenotype_data['series_id'])


# %%
#######
# THIS IS WHERE I PUT THE CODE THAT FAILED ME#
#######

# %%
def extract(tar_url, extract_path='.'):
    print(tar_url)
    tar = tarfile.open(tar_url, 'r')
    for item in tar:
        tar.extract(item, extract_path)
        if item.name.find(".tgz") != -1 or item.name.find(".tar") != -1:
            extract(item.name, "./" + item.name[:item.name.rfind('/')])


try:

    extract(sys.argv[1] + '.tgz')
    print('Done.')
except:
    name = os.path.basename(sys.argv[0])
    print(name[:name.rfind('.')], '<filename>')
# %%
extract('GSE252723_family.xml.tgz')

# %%
# %%
# Reading the data inside the xml
# file to a variable under the name
# data
# %%
with open('GSE252723_family.xml', 'r') as f:
    data = f.read()

# %%
# Passing the stored data inside
# the beautifulsoup parser, storing
# the returned object
Bs_data = BeautifulSoup(data, "xml")
print(Bs_data)
# %%
xmlTree = ET.parse('GSE252723_family.xml')
print(xmlTree)
tags = {elem.tag for elem in xmlTree.iter()}
tags
# %%
# Finding all instances of tag
# `unique`
b_unique = Bs_data.find_all('Data-Table')
b_unique
# %%
headers = {
    'authority': 'screen-beta-api.wenglab.org',
    'accept': 'application/json',
    'accept-language': 'en-US,en;q=0.9',
    'content-type': 'application/json',
    'origin': 'https://screen.wenglab.org',
    'referer': 'https://screen.wenglab.org/',
    'sec-ch-ua': '\".Not/A)Brand\";v=\"99\", \"Google Chrome\";v=\"103\", \"Chromium\";v=\"103\"',
    'sec-ch-ua-mobile': '?1',
    'sec-ch-ua-platform': '\"Android\"',
    'sec-fetch-dest': 'empty',
    'sec-fetch-mode': 'cors',
    'sec-fetch-site': 'same-site',
    'user-agent': 'Mozilla/5.0 (Linux; Android 6.0; Nexus 5 Build/MRA58N) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/103.0.0.0 Mobile Safari/537.36'
}
# %%
url = 'https://www.ncbi.nlm.nih.gov/gds/?term=200264681'
# %%
url_output = requests.get(url)
url_output.headers['Content-Type']
url_output.text
geo_json = url_output.json()

pd.readhttps: // www.ncbi.nlm.nih.gov/gds /?term = 200264681
type(X)
X
X.loc[0, 0].startswith('1.')


fname = 'guppy-0.1.10.tar.gz'
url = 'https://pypi.python.org/packages/source/g/guppy/' + fname
r = requests.get(url)
open(fname, 'wb').write(r.content)
