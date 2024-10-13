# Description: This script is used to identify the candidate pairs of species that are in different domains and have been clustered together by at least 3 clustering algorithms.
import os
import json
from Bio import SeqIO
import re
from scripts.src.utils import clean_sequence


GOOD_ALGO = {"Temperature":
                 ["VAE+IM", "VAE+affinity_propagation", "VAE+HDBSCAN", "UMAP+HDBSCAN", "iDeLUCS"],
             "pH":
                 ["VAE+IM", "VAE+affinity_propagation", "VAE+HDBSCAN", "iDeLUCS"]}
k = 8
ENVS = ["Temperature", "pH"]
ids_2_dist = {}

def candidates_identification(summary_dataset, env):
    pairs = {}
    file_path = os.path.join(f'/content/candidates_{env}.json')
    if os.path.isfile(file_path):
        with open(file_path, 'r') as file:
            pairs = json.load(file)

    for name in GOOD_ALGO[env]:
        assignment_algo = name

        for group in summary_dataset.groupby(assignment_algo):

            if (group[0] != -1):  # (group[0] in GOOD_CLUSTERS[name])
                species = group[1]["Species"].values
                ids = group[1]["Assembly"].values
                for i in range(len(species) - 1):
                    for j in range(i + 1, len(species)):
                        if group[1]["Domain"].values[i] != group[1]["Domain"].values[j]:
                            pair = tuple(sorted((species[i], species[j])))
                            pair_species = str(pair[0]) + "_" + str(pair[1])
                            if pair_species in pairs:
                                pairs[pair_species]["algo"].append(name)
                                pairs[pair_species]["pair_id_first"].append(ids[i])
                                pairs[pair_species]["pair_id_second"].append(ids[j])
                                pairs[pair_species]["counts"] += 1
                            else:
                                pairs[pair_species] = {}
                                pairs[pair_species]["pair_species_first"] = [species[i]]
                                pairs[pair_species]["pair_species_second"] = [species[j]]
                                pairs[pair_species]["algo"] = [name]
                                pairs[pair_species]["pair_id_first"] = [ids[i]]
                                pairs[pair_species]["pair_id_second"] = [ids[j]]
                                pairs[pair_species]["counts"] = 1

    with open(file_path, 'w') as file:
        json.dump(pairs, file, indent=2)


def analyse_candidates(file_dir, env):
    file_path = os.path.join(f'{file_dir}/candidates_{env}.json')

    with open(file_path, 'r') as file:
        pairs = json.load(file)

    final_pairs = {}
    just_pairs_ids = {}

    for i in pairs:
        clustering_models = list(set(pairs[i]["algo"]))
        count = pairs[i]["counts"]

        if len(clustering_models) >= 3 and count >= 5:
            final_pairs[i] = pairs[i]
            final_pairs[i]["algo"] = list(set(pairs[i]["algo"]))
            final_pairs[i]["pair_id_first"] = list(set(pairs[i]["pair_id_first"]))
            final_pairs[i]["pair_id_second"] = list(set(pairs[i]["pair_id_second"]))

            just_pairs_ids[final_pairs[i]["pair_species_first"][0]+","+final_pairs[i]["pair_species_second"][0]] = (final_pairs[i]["pair_id_first"][0], final_pairs[i]["pair_id_second"][0])

    print(len(just_pairs_ids))
    final_file_path = os.path.join(f'{file_dir}/final_candidates_{env}.json')
    with open(final_file_path, 'w') as file:
        json.dump(final_pairs, file, indent=2)

    final_file_path = os.path.join(f'{file_dir}/pairs_ids_candidates_{env}.json')
    with open(final_file_path, 'w') as file:
        json.dump(just_pairs_ids, file, indent=2)

    return just_pairs_ids

def read_fasta(file_path):

    id_2_sequences = {}
    # Open the FASTA file
    with open(file_path, 'r') as file:
        # Iterate over each record
        for record in SeqIO.parse(file, 'fasta'):
            id_2_sequences[str(record.id)] = str(record.seq)

    return id_2_sequences

for env in ENVS:
    ids_2_dist = {}
    just_pairs_ids = analyse_candidates("../candidates",env)
    print(ids_2_dist)
    max_key = max(ids_2_dist, key=ids_2_dist.get)
    max_value = ids_2_dist[max_key]
    print(max_key, max_value)

    min_key = min(ids_2_dist, key=ids_2_dist.get)
    min_value = ids_2_dist[min_key]
    print(min_key, min_value)

    json_file_path = f'../Distances/candidates_{env}.json'
    # Write the dictionary to the JSON file
    with open(json_file_path, 'w') as json_file:
        json.dump(ids_2_dist, json_file, indent=4)