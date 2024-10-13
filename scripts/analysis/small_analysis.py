# # # # # import pandas as pd
# # # # # import os
# # # # # import json
# # # # # from Bio import SeqIO
# # # # #
# # # def read_fasta(file_path):
# # #
# # #     id_2_sequences = {}
# # #     # Open the FASTA file
# # #     with open(file_path, 'r') as file:
# # #         # Iterate over each record
# # #         for record in SeqIO.parse(file, 'fasta'):
# # #             id_2_sequences[str(record.id)] = str(record.seq)
# # #
# # #     return id_2_sequences
# # # # #
# # # # #
# # # # #
# # # # # results_folder = os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..', 'results'))
# # # # # filtered_folder = os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..', 'Filtered_FCGRs'))
# # # # #
# # # # #
# # # # # fragment_length = 100000
# # # # # s = []
# # # # # g = []
# # # # # for env in ["Temperature", "pH"]:
# # # # #     data_folder = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'data'))
# # # # #     # candidates = json.load(open(f"{filtered_folder}/candidates/majority_candidates.json", 'r'))
# # # # #     candidates = os.listdir(filtered_folder)
# # # # #     print(candidates)
# # # # #     data_folder = f"{data_folder}/fragments_{fragment_length}"
# # # # #     fasta_file = os.path.join(data_folder, env, f'Extremophiles_{env}.fas')
# # # # #     summary_path = os.path.join(data_folder, env, f'Extremophiles_{env}_Summary.tsv')
# # # # #     summary_data = pd.read_csv(summary_path, sep='\t')
# # # # #     id_2_sequences = read_fasta(fasta_file)
# # # # #
# # # # #     for key in candidates:
# # # # #         try:
# # # # #             print(key)
# # # # #             id1, id2 = key.split('_')[0]+"_"+key.split('_')[1], key.split('_')[2]+"_"+key.split('_')[3].replace('.png', '')
# # # # #             domain1 = str(list(summary_data[summary_data['Assembly'] == id1]['Domain'])[0])
# # # # #             domain2 = str(list(summary_data[summary_data['Assembly'] == id2]['Domain'])[0])
# # # # #             genus1 = str(list(summary_data[summary_data['Assembly'] == id1]['genus'])[0])
# # # # #             genus2 = str(list(summary_data[summary_data['Assembly'] == id2]['genus'])[0])
# # # # #             species1 = str(list(summary_data[summary_data['Assembly'] == id1]['species'])[0])
# # # # #             species2 = str(list(summary_data[summary_data['Assembly'] == id2]['species'])[0])
# # # # #             env_label1 = str(list(summary_data[summary_data['Assembly'] == id1][env])[0])
# # # # #             env_label2 = str(list(summary_data[summary_data['Assembly'] == id2][env])[0])
# # # # #             s.append(species1)
# # # # #             s.append(species2)
# # # # #             g.append(genus1)
# # # # #             g.append(genus2)
# # # # #
# # # # #
# # # # #         except:
# # # # #             continue
# # # # #
# # # # #     print(len(candidates))
# # # # #     print(len(list(set(s))))
# # # # #     print(len(list(set(g))))
# # # #
# # # #
# # # #
# # # # import matplotlib.pyplot as plt
# # # # import numpy as np
# # # #
# # # # # Example data (replace with actual data extraction)
# # # # data = {
# # # #     'GCA_000211475.1': {'k-mers': ['AAA', 'AAC', 'AAG', 'AAT'], 'deviation': [0.02, 0.01, 0.00, -0.01]},
# # # #     'GCA_000246985.3': {'k-mers': ['AAA', 'AAC', 'AAG', 'AAT'], 'deviation': [0.015, 0.010, 0.005, 0.000]},
# # # #     'GCA_000446015.1': {'k-mers': ['AAA', 'AAC', 'AAG', 'AAT'], 'deviation': [0.02, 0.01, 0.00, -0.01]},
# # # #     'GCA_000512735.1': {'k-mers': ['AAA', 'AAC', 'AAG', 'AAT'], 'deviation': [0.020, 0.015, 0.010, 0.005]},
# # # #     'GCA_000725425.1': {'k-mers': ['AAA', 'AAC', 'AAG', 'AAT'], 'deviation': [0.020, 0.015, 0.010, 0.005]},
# # # #     'GCA_002214605.1': {'k-mers': ['AAA', 'AAC', 'AAG', 'AAT'], 'deviation': [0.020, 0.015, 0.010, 0.005]}
# # # # }
# # # #
# # # # # Create a combined dataframe
# # # # all_kmers = sorted(set(kmer for values in data.values() for kmer in values['k-mers']))
# # # # combined_data = {kmer: [] for kmer in all_kmers}
# # # #
# # # # # Fill the combined_data dictionary with deviations for each genome
# # # # for genome, values in data.items():
# # # #     for kmer in all_kmers:
# # # #         if kmer in values['k-mers']:
# # # #             combined_data[kmer].append(values['deviation'][values['k-mers'].index(kmer)])
# # # #         else:
# # # #             combined_data[kmer].append(0)  # Use 0 if k-mer is not present in this genome
# # # #
# # # # # Convert combined_data to a format suitable for plotting
# # # # kmers = list(combined_data.keys())
# # # # deviations = np.array(list(combined_data.values())).T
# # # # genomes = list(data.keys())
# # # #
# # # # # Plotting
# # # # fig, ax = plt.subplots(figsize=(15, 8))
# # # # bar_width = 0.15
# # # # index = np.arange(len(kmers))
# # # #
# # # # # Grouped bar plot
# # # # for i, genome in enumerate(genomes):
# # # #     ax.bar(index + i * bar_width, deviations[i], bar_width, label=genome)
# # # #
# # # # ax.set_xlabel('K-mers')
# # # # ax.set_ylabel('Deviation')
# # # # ax.set_title('K-mer Deviation Across Genomes')
# # # # ax.set_xticks(index + bar_width * (len(genomes) - 1) / 2)
# # # # ax.set_xticklabels(kmers)
# # # # ax.legend()
# # # #
# # # # plt.xticks(rotation=45)
# # # # plt.tight_layout()
# # # # plt.savefig('grouped_bar_chart.png')
# # # # plt.show()
# # # #
# # #
# # #
# # # import os
# # # import pandas as pd
# # #
# # # def list_files_in_folder(folder_path):
# # #     # List to store file names
# # #     file_names = []
# # #
# # #     # Iterate over all files in the folder
# # #     for file_name in os.listdir(folder_path):
# # #         # Append the file name to the list
# # #         file_names.append(file_name.split("_")[0]+"_"+file_name.split("_")[1].replace(".png", ""))
# # #         file_names.append(file_name.split("_")[2] + "_" + file_name.split("_")[3].replace(".png", ""))
# # #
# # #
# # #     return file_names
# # #
# # # results_folder = os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..', 'Filtered_FCGRs'))
# # # file_names = list(set(list_files_in_folder(results_folder)))
# # #
# # # fragment_length = 100000
# # # b = []
# # # a = []
# # # g_b = []
# # # g_a = []
# # # for env in ["Temperature", "pH"]:
# # #     data_folder = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'data'))
# # #     data_folder = f"{data_folder}/fragments_{fragment_length}"
# # #     fasta_file = os.path.join(data_folder, env, f'Extremophiles_{env}.fas')
# # #     summary_path = os.path.join(data_folder, env, f'Extremophiles_{env}_Summary.tsv')
# # #     summary_data = pd.read_csv(summary_path, sep='\t')
# # #     for id in file_names:
# # #         try:
# # #             domain = str(list(summary_data[summary_data['Assembly'] == id]['Domain'])[0])
# # #             genus = str(list(summary_data[summary_data['Assembly'] == id]['genus'])[0])
# # #             if domain == "Bacteria":
# # #                 b.append(id)
# # #                 g_b.append(genus)
# # #             if domain == "Archaea":
# # #                 a.append(id)
# # #                 g_a.append(genus)
# # #         except:
# # #             continue
# # #
# # #
# # # print(len(set(a)))
# # # print(len(set(b)))
# # #
# # # print(len(list(set(file_names))))
# # # print(len(list(set(g_a))))
# # # print(len(list(set(g_b))))
# #
# # from scripts.pipelines.build_signature_dataset import *
# #
# # def save_fasta(whole_genome, file, cluster, name):
# #     assembly_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'clusters', 'concat'))
# #     assembly_path2 = os.path.join(assembly_path, cluster)
# #     with open(os.path.join(assembly_path2 ,f'{file}.fas'), 'w') as f:
# #         f.write(f'>{name}\n')
# #         f.write(f'{whole_genome}\n')
# #
# # c = ["cluster 6", "cluster 12"]
# #
# # for cluster in c:
# #     assembly_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'clusters'))
# #     assembly_path2 = os.path.join(assembly_path, cluster)
# #     for file in os.listdir(assembly_path2):
# #         assembly_path3 = os.path.join(assembly_path2, file)
# #         fasta_file = os.listdir(assembly_path3) if os.listdir(assembly_path3) else None
# #         _min = 0
# #         fasta_path = os.path.join(assembly_path3, fasta_file[0])
# #         seq_names, seqs, plasmid_names = summary_fasta(fasta_path, _min)
# #         print(len(seq_names))
# #         whole_genome = produce_fragment(seq_names, seqs, 100000, is_whole_genome=True)[:-1]
# #         print(len(whole_genome))
# #         save_fasta(whole_genome, file, cluster, seq_names[0])
#
#
#
# import os
# from Bio import SeqIO
#
# def combine_fasta_files(input_folder, output_file):
#     with open(output_file, 'w') as outfile:
#         for filename in os.listdir(input_folder):
#             print(filename)
#             if filename.endswith(".fasta") or filename.endswith(".fas"):
#                 filepath = os.path.join(input_folder, filename)
#                 with open(filepath, 'r') as infile:
#                     for record in SeqIO.parse(infile, 'fasta'):
#                         SeqIO.write(record, outfile, 'fasta')
#
# assembly_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'clusters'))
# # Define the input folder containing the FASTA files and the output file name
# input_folder = f"{assembly_path}/concat/cluster 12"
# output_file = f'{assembly_path}/combined_sequences_cluster12.fasta'
#
# # Combine the FASTA files
# combine_fasta_files(input_folder, output_file)
#

import os
import json
import numpy as np

# results_folder = os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..', 'results'))
# dist = json.load(open(f"{results_folder}/distances/Domain_Temperature_level_dist_MLDSP.json", 'r'))
# avg = {}
# for i in dist:
#     print(i)
#     avg[i] = np.mean(list(dist[i].values()))
#
# print(avg)
# print(max(avg.values()))
# print(min(avg.values()))
# print(np.mean(list(avg.values())))

def trim_outliers(data, percent):
    # Calculate the lower and upper percentiles
    lower_percentile = percent / 2
    upper_percentile = 100 - lower_percentile

    # Calculate the lower and upper threshold values
    lower_threshold = np.percentile(data, lower_percentile)
    upper_threshold = np.percentile(data, upper_percentile)

    # Filter out the outliers
    trimmed_data = [x for x in data if lower_threshold <= x <= upper_threshold]

    return trimmed_data


results_folder = os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..', 'results'))
dist = json.load(open(f"{results_folder}/distances/genus_pH_level_dist_Desc.json", 'r'))
avg = {}
for i in dist:
    # print(i)
    avg[i] = np.mean(list(dist[i].values()))

print(avg)
trimmed_avg = trim_outliers(list(avg.values()), 5)
print(max(trimmed_avg))
# print(max(avg.values()))
print(min(trimmed_avg))
print(np.mean(list(trimmed_avg)))



# import json
# import os
# import shutil
#
# def read_json(file_path):
#     with open(file_path, 'r') as file:
#         data = json.load(file)
#
#     return data
#
# results_folder = os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..', 'results'))
# dist = read_json(f"{results_folder}/candidates/majority_candidates_new_threshold.json")
# fcgr_folder = f"{results_folder}/FCGRs"
# d = {}
# for m, k in enumerate(dist):
#     print(k)
#     new_k = k.replace("'","")
#     shutil.copy(f"{fcgr_folder}/all_candidates/k8/{new_k}.png", f"{fcgr_folder}/accepted_candidates_new/{k}.png")
#     d[m] = dist[k]
#
# with open(f"{results_folder}/candidates/final_candidates.json", 'w') as file:
#     json.dump(d, file, indent=4)

