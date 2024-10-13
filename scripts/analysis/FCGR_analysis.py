"""Generating avoided pattern plot for the paper"""
import json
import matplotlib.pyplot as plt
import scripts.src.CGR_utils as CGR_utils
from Bio import SeqIO
import os
import pandas as pd
import re

ENVS = ["Temperature", "pH"]
NUM_CLUSTERS = {"Temperature": 4, "pH": 2}
k = 8

def draw_fcgrs(sequence, id, len, domain, env_label, env, species, results_folder, mode):

    fcgr = CGR_utils.FCGR(k=k, bits=8)
    chaos1 = fcgr(sequence[0])
    chaos2 = fcgr(sequence[1])
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

    # Display the first image
    ax1.imshow(fcgr.plot(chaos1), cmap='gray')
    ax1.set_title(f'{domain[0]}_\n{env_label[0]}\n{species[0]}')
    ax1.axis('off')  # Hide the axes

    # Display the second image
    ax2.imshow(fcgr.plot(chaos2), cmap='gray')
    ax2.set_title(f'{domain[1]}_\n{env_label[1]}\n{species[1]}')
    ax2.axis('off')  # Hide the axes
    plt.tight_layout(pad=3.0)
    if mode == 'all':
        plt.savefig(f"{results_folder}/FCGRs/all_candidates/k{k}/{id[0]}_{id[1]}.png", dpi=300)
        plt.xticks([])
        plt.yticks([])
        plt.close()
    if mode == 'selected':
        plt.savefig(f"{results_folder}/FCGRs/accepted_candidates/k{k}/{id[0]}_{id[1]}.png", dpi=300)
        plt.xticks([])
        plt.yticks([])
        plt.close()

def read_fasta(file_path):

    id_2_sequences = {}
    # Open the FASTA file
    with open(file_path, 'r') as file:
        # Iterate over each record
        for record in SeqIO.parse(file, 'fasta'):
            id_2_sequences[str(record.id)] = str(record.seq)

    return id_2_sequences

def clean_sequence(sequence):
    # Replace any character not A, C, G, T, or N with N
    return re.sub(r'[^ACGTN]', 'N', sequence)

def run(candidates, env, fragment_length, results_folder, data_folder, mode = "all"):


    data_folder = f"{data_folder}/fragments_{fragment_length}"
    fasta_file = os.path.join(data_folder, env, f'Extremophiles_{env}.fas')
    summary_path = os.path.join(data_folder, env, f'Extremophiles_{env}_Summary.tsv')
    summary_data = pd.read_csv(summary_path, sep='\t')
    id_2_sequences = read_fasta(fasta_file)
    for key in candidates.keys():
        try:
            if mode == 'selected':
                id1, id2 = key.split('_')[0]+"_"+key.split('_')[1], key.split('_')[2]+"_"+key.split('_')[3]
            else:
                id1, id2 = candidates[key][0], candidates[key][1]
            domain1 = str(list(summary_data[summary_data['Assembly'] == id1]['Domain'])[0])
            domain2 = str(list(summary_data[summary_data['Assembly'] == id2]['Domain'])[0])
            species1 = str(list(summary_data[summary_data['Assembly'] == id1]['species'])[0])
            species2 = str(list(summary_data[summary_data['Assembly'] == id2]['species'])[0])
            env_label1 = str(list(summary_data[summary_data['Assembly'] == id1][env])[0])
            env_label2 = str(list(summary_data[summary_data['Assembly'] == id2][env])[0])
            sequence_id1 = clean_sequence(id_2_sequences[id1])
            sequence_id2 = clean_sequence(id_2_sequences[id2])
            draw_fcgrs((sequence_id1, sequence_id2), (id1, id2), fragment_length, (domain1, domain2), (env_label1, env_label2), env, (species1, species2), results_folder, mode)
        except:
            continue

    return id_2_sequences


##########################################################################################

# results_folder = os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..', 'results'))
# data_folder = os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..', 'data'))
# fragment_length = 100000
# for env in ["Temperature", "pH"]:
#     with open(f"{results_folder}/candidates/pairs_ids_candidates_{env}.json", 'r') as f:
#         candidates = json.load(f)
#     id_2_sequences = run(candidates, env, fragment_length, results_folder, data_folder)


