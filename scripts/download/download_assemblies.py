"""
This script downloads the assemblies from NCBI 
using python's NCBI-dataset-API for sequence 
download.

python Download_Assemblies.py --assembly_file=<example.txt> 
                              --Entrez_email=foo@bar.ca 
                              --get_links
"""

from genericpath import isdir, isfile
import argparse
import os
import sys
from tqdm import tqdm

import urllib.request
from Bio import Entrez
import time
import subprocess

def get_links(args):
    if os.path.isfile('../data/assembly_links.txt'):
        os.remove('../data/assembly_links.txt')

    sys.stdout.write('....Opening assemblies file.....\n')
    sys.stdout.flush()
    with open(args["assembly_file"], 'r') as f:
        sys.stdout.write('....Collecting assembly links.....\n')
        sys.stdout.flush()

        lines = f.readlines()
        for ass in tqdm(lines):
            try:
                ass = ass[:-1]
                handle = Entrez.esearch(db="assembly", term=ass)
                record = Entrez.read(handle)
                for id in record['IdList']:
                    #print(id)
                    # Get Assembly Summary
                    esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
                    esummary_record = Entrez.read(esummary_handle)

                    # Parse Accession Name from Summary
                    link = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
                    if link == "":
                        link = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
                    
                    with open(os.path.join('../data/assembly_links.txt'),'a') as f:
                        f.write(link)
                        f.write('\n')
                           
                time.sleep(1)
        
            except: "urllib.error.HTTPError"

def run(args):
    
    Entrez.email = args["Entrez_email"]
    if args['get_links']:
            get_links(args)
    
    if not os.path.isdir('../data/Assemblies'):
        os.mkdir('Assemblies')

    with open('../data/assembly_links.txt','r') as f:
        links = f.readlines()
    
    sys.stdout.write('....Downloading assemblies, this might take a while .....\n')
    sys.stdout.flush()
    
    for link in tqdm(links):
        link = link.split('\n')[0]
        name = "_".join((link.split('/')[-1]).split('_')[:2])
        
        if not os.path.isfile('../data/Assemblies/'+name+'.fna'):
            output = subprocess.call(['wget','-q', '-r', link])
            time.sleep(4.25)
            stream = os.popen('find ftp.ncbi.nlm.nih.gov -name "*fna.gz"')
            output = stream.read()
            output = subprocess.call(['mv', output.split('\n')[-2], '../data/Assemblies/'+name+'.fna.gz'])
            time.sleep(1)
            output = subprocess.call(['gzip', '-d', '../data/Assemblies/'+name+'.fna.gz'])
            time.sleep(1)
            output = subprocess.call(['rm', '-r', 'ftp.ncbi.nlm.nih.gov'])
            time.sleep(1)
    
    sys.stdout.write('....All assemblies have been downloaded......\n')
    sys.stdout.flush()
    

    print('Check the following list for missing records:')

    with open(args["assembly_file"], 'r') as f:

        lines = f.readlines()
        for ass in lines:
            name = ass[:-1]
            name_1 = 'GCA'+name[3:]
            name_2 = 'GCF'+name[3:]

            if (not os.path.isfile('../data/Assemblies/'+name_1+'.fna')) and (not os.path.isfile('../data/Assemblies/'+name_2+'.fna')):
                print('-',name)

    print('Try running the scripts again without the "get_links" flag to download missing records. \n If the problem persist, you will have to download them manually')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--assembly_file', action='store', type=str, required=True)
    parser.add_argument('--Entrez_email', action='store', type=str, required=True)
    parser.add_argument('--get_links', action='store_true')
    args= parser.parse_args()
    args=vars(args)
    run(args)


if __name__ =='__main__':
    main()