import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
import shap
from extprime.utils.utils import kmersFasta
import json
import matplotlib.pyplot as plt
import os
import matplotlib.colors as mcolors
from scipy.stats import spearmanr
import seaborn as sns
from adjustText import adjust_text  # For better label placement
from pathlib import Path

k = 3
kmers_label = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'ATA', 'ATC', 'ATG', 'CAA',
               'CAC', 'CAG', 'CCA', 'CCC', 'CCG', 'CGA', 'CGC', 'CTA', 'CTC', 'GAA', 'GAC', 'GCA', 'GCC', 'GGA', 'GTA',
               'TAA', 'TCA']



temperature = [0.05315811, 0.03050008, 0.03524833, 0.0402312, 0.025798, 0.02891119,
        0.02633873, 0.02382808, 0.03098142, 0.03082318, 0.02975528, 0.03463928,
        0.03377405, 0.03096574, 0.03451577, 0.02486449, 0.03002391, 0.03156688,
        0.02686612, 0.03266274, 0.03263328, 0.03274237, 0.0198588, 0.0302579,
        0.03705325, 0.0261334, 0.02926867, 0.03520646, 0.02986665, 0.02337804,
        0.03441187, 0.03373675]

pH = [0.04398287, 0.02889037, 0.03174242, 0.03457252, 0.02425236, 0.02925751,
      0.03019546, 0.02211669, 0.02953973, 0.03104185, 0.02966024, 0.02946116,
      0.03500242, 0.03042681, 0.03153062, 0.02559669, 0.03128845, 0.03253462,
      0.0277128, 0.04004679, 0.04029354, 0.04148019, 0.0174389, 0.03188939,
      0.03699027, 0.03024128, 0.03045014, 0.04029907, 0.03268185, 0.02109142,
      0.02668636, 0.03160523]


def load_data(env, path, fragment_length=100000):
    data_path = f'{path}/fragments_{fragment_length}/{env}/Extremophiles_{env}.fas'
    label_path = f'{path}/fragments_{fragment_length}/{env}/Extremophiles_{env}_GT_{mode}.tsv'
    summary_data = pd.read_csv(f'{path}/fragments_{fragment_length}/{env}/Extremophiles_{env}_Summary.tsv', sep='\t')
    labels = pd.read_csv(label_path, sep='\t')['cluster_id']
    _, kmers = kmersFasta(data_path, k=k, transform=None, reduce=True)
    kmers_normalized = np.transpose((np.transpose(kmers) / np.linalg.norm(kmers, axis=1)))
    return summary_data, labels, kmers_normalized

def encode_labels(labels):
    unique_labels = list(labels.unique())
    labels_encoder = {k: i for i, k in enumerate(unique_labels)}
    y = [labels_encoder[i] for i in labels]
    return y, labels_encoder

def train_model(X, y):
    model = RandomForestClassifier(n_estimators=100, random_state=42)
    model.fit(X, y)
    return model

def compute_shap_values(model, X):
    explainer = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(X)
    return shap_values

def extract_features(summary_data, labels_encoder, shap_values, clusters, kmers_label, env):
    features_all = {}
    for c in clusters:
        for ids in clusters[c]:
            features = {"pos": []}
            # features = {}
            wanted_index = summary_data[summary_data["Assembly"] == ids].index[0]
            l = list(summary_data[summary_data["Assembly"] == ids][env])
            num_label = labels_encoder[l[0]]
            sorted_shap_values = list(enumerate([np.abs(i) for i in shap_values[wanted_index, :, num_label]]))
            sorted_shap_values = sorted(sorted_shap_values, key=lambda x: x[1], reverse=True)
            print(sorted_shap_values)
            kmer_indices = [x[0] for x in sorted_shap_values[:15]]
            features["pos"] = [kmers_label[i] for i in kmer_indices]
            # pos = [kmers_label[i] for i in kmer_indices]

            # for i in sorted_shap_values:
            #     if kmers_label[i[0]] in pos:
            #         features[kmers_label[i[0]]] = i[1]
            #     else:
            #         features[kmers_label[i[0]]] = -1

            features_all[ids] = features


    return features_all


def correlation_groups(clusters, k_mer_dict, domain_dict, genus_dict, species_dict, plot_title_size=20, plot_axis_numbers_size=15,
                      plot_axis_title_size=20):
    """
    Create correlation plots for k-mer profiles within clusters with min-max line.
    """
    comparison_label = "3-mer normalized counts"
    for c in clusters:
        print(f"Cluster {c}")
        comparison_label = f" {c}"

        bacteria = [i for i in clusters[c] if domain_dict[i] == "Bacteria"]
        archaea = [i for i in clusters[c] if domain_dict[i] == "Archaea"]


        for id_1 in bacteria:
            k_mer_ids = k_mer_dict[id_1]
            for id_2 in archaea:
                if id_1 != id_2:
                    k_mer_ids_2 = k_mer_dict[id_2]
                    cor_test_rho, cor_test_pvalue = spearmanr(k_mer_ids, k_mer_ids_2)

                    # Create plot title with correlation results
                    if cor_test_pvalue < 1e-05:
                        plot_title = f"{comparison_label}, rho={cor_test_rho:.2f}, p < 1e-05"
                    else:
                        plot_title = f"{comparison_label}, rho={cor_test_rho:.2f}, p={cor_test_pvalue:.5f}"


                    plot_df = pd.DataFrame({
                        'mean_counts_g1': k_mer_ids,
                        'mean_counts_g2': k_mer_ids_2,
                        'sequence': kmers_label
                    })

                    # Create the plot
                    plt.figure(figsize=(10, 10))

                    # Plot all points
                    plt.scatter(plot_df['mean_counts_g1'], plot_df['mean_counts_g2'],
                              alpha=0.6, color='#45b39d', s=100)

                    # Add line connecting min and max points
                    plot_df_sorted = plot_df.sort_values('mean_counts_g1')
                    plt.plot([plot_df_sorted['mean_counts_g1'].iloc[0],
                             plot_df_sorted['mean_counts_g1'].iloc[-1]],
                            [plot_df_sorted['mean_counts_g2'].iloc[0],
                             plot_df_sorted['mean_counts_g2'].iloc[-1]],
                            'r--', alpha=1, color = "#F08080")


                    texts = []
                    for _, row in plot_df.iterrows():
                        texts.append(plt.text(row['mean_counts_g1'], row['mean_counts_g2'],
                                              row['sequence'], fontsize=10))
                    adjust_text(texts, arrowprops=dict(arrowstyle='->', color='black', lw=1))
                    species_1 = species_dict[id_1]
                    species_2 = species_dict[id_2]

                    if species_2 == "Methanoculleus_A taiwanensis": species_2 = "Methanoculleus taiwanensis"
                    if species_2 == "Methanolinea_B mesophila": species_2 = "Methanolinea mesophila"
                    # if species_1 == "Pseudothermotoga_B elfii": species_1 = "Pseudothermotoga elfii"

                    x_label = f"{species_1}"
                    y_label = f"{species_2}"
                    # Customize plot appearance
                    plt.title(plot_title, fontsize=plot_title_size)
                    plt.xlabel(x_label, fontsize=plot_axis_title_size, fontstyle='italic')
                    plt.ylabel(y_label, fontsize=plot_axis_title_size, fontstyle='italic')
                    plt.xticks(fontsize=plot_axis_numbers_size)
                    plt.yticks(fontsize=plot_axis_numbers_size)
                    plt.grid(True, alpha=0.3)






                    # Make plot square and set equal scales
                    plt.axis('square')

                    plt.tight_layout()
                    print(comparison_label)
                    plt.savefig(f"{results_folder}/statistics/{comparison_label}_{id_1}_{id_2}.png", format='png')
                    plt.show()
                    plt.close()


def significance_test(summary_data, clusters, env, k_mer_dict):
    for c in clusters:
        print(f"Cluster {c}")
        num_microbes = len(clusters[c])
        correlation_matrix = np.zeros((num_microbes, num_microbes))
        for id_1 in clusters[c]:
            k_mer_ids = k_mer_dict[id_1]
            for id_2 in clusters[c]:
                # if id_1 != id_2:
                k_mer_ids_2 = k_mer_dict[id_2]
                coef, p_value = spearmanr(k_mer_ids, k_mer_ids_2)
                correlation_matrix[clusters[c].index(id_1), clusters[c].index(id_2)] = coef
                print(f"Microbe {id_1} vs Microbe {id_2}: Spearman Coefficient = {coef}, p-value = {p_value}")

        # Use microbe IDs directly as labels instead of numeric indices
        labels = clusters[c]  # Use the actual microbe IDs

        # Plot heatmap with microbe IDs
        fig, ax = plt.subplots(figsize=(12, 10))  # Increased figure size to accommodate longer labels
        ax = sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm',
                         xticklabels=labels, yticklabels=labels,
                         vmin=-1, vmax=1)

        plt.title('Spearman Correlation Heatmap of 3-mer Counts - Cluster ' + c)
        ax.set_xlabel("Microbe ID")
        ax.set_ylabel("Microbe ID")

        # Rotate x-axis labels for better readability
        plt.xticks(rotation=45, ha='right')
        plt.yticks(rotation=0)

        # Adjust layout to prevent label cutoff
        plt.tight_layout()
        plt.savefig(f"{results_folder}/statistics/{c}_correlation.png", format='png',
                    bbox_inches='tight')  # Added bbox_inches='tight' to prevent label cutoff
        plt.show()


def save_features(features_all, file_path):
    with open(file_path, 'w') as json_file:
        json.dump(features_all, json_file, indent=4)




import matplotlib.pyplot as plt

import matplotlib.pyplot as plt
import numpy as np

import matplotlib.pyplot as plt
import numpy as np

import matplotlib.pyplot as plt
import numpy as np

def plot_deviations(clusters, temp, pH, kmers_label, features_important, labels, names, kmers, bacteria):
    for c in clusters:
        ds = {}
        max_list = []
        min_list = []
        b_key = bacteria[c]  # Bacterial reference key for this cluster
        b_rel = features_important[b_key]['pos']
        print(f"--- Cluster {c} ---")
        print("Bacterial relevant k-mers:", b_rel)

        # Build deviation values
        for n, k in zip(names, kmers):
            if n in clusters[c]:
                deviations = ([elem2 - elem1 for elem1, elem2 in zip(temp, k)]
                              if c != 'c10'
                              else [elem2 - elem1 for elem1, elem2 in zip(pH, k)])
                ds[n] = deviations
                max_list.append(max(deviations))
                min_list.append(min(deviations))

        max_all = max(max_list)
        min_all = min(min_list)

        num_plots = len(ds)
        num_cols = 2
        num_rows = (num_plots + num_cols - 1) // num_cols

        fig, axs = plt.subplots(num_rows, num_cols, figsize=(32, num_rows * 5))
        axs = axs.flatten()

        # Step 1: identify bold-worthy kmers from all Archaea
        ref_devs = ds[b_key]
        bold_kmers_from_archaea = set()

        for a in clusters[c]:
            domain = list(labels[labels["Assembly"] == a]["Domain"])
            if domain and domain[0] != "Bacteria":
                a_rel = features_important[a]['pos']
                cur_devs = ds[a]

                common_kmers = set(b_rel).intersection(a_rel)
                same_sign_kmers = []

                for i, kmer in enumerate(kmers_label):
                    if kmer in common_kmers:
                        if np.sign(ref_devs[i]) == np.sign(cur_devs[i]) and np.sign(ref_devs[i]) != 0:
                            bold_kmers_from_archaea.add(kmer)
                            same_sign_kmers.append(kmer)

                # ✅ PRINT THE RATIO
                species_name = labels[labels["Assembly"] == a]["species"].values[0]
                if common_kmers:
                    ratio = len(same_sign_kmers) / len(common_kmers)
                    print(f"{a} ({species_name}) | Same-sign/Intersection: {len(same_sign_kmers)}/{len(common_kmers)} = {ratio:.2f}")
                else:
                    print(f"{a} ({species_name}) | No common relevant k-mers.")

        for index, (key, ax) in enumerate(zip(ds, axs)):
            intersect = set(b_rel).intersection(set(features_important[key]['pos']))
            l_domain = list(labels[labels["Assembly"] == key]["Domain"])
            l_env = list(labels[labels["Assembly"] == key]["species"])

            bars = ax.bar(kmers_label, ds[key], color='blue')

            for bar, l in zip(bars, kmers_label):
                bar.set_color('#45b39d' if l in features_important[key]['pos'] else '#F08080')

            cur_devs = ds[key]

            new_labels = []
            if l_domain and l_domain[0] == "Bacteria":
                for i, l in enumerate(kmers_label):
                    if l in bold_kmers_from_archaea:
                        new_labels.append(r"$\bf{" + l + "}$")
                    else:
                        new_labels.append(l)
            else:
                for i, l in enumerate(kmers_label):
                    if l in intersect:
                        if np.sign(ref_devs[i]) == np.sign(cur_devs[i]) and np.sign(ref_devs[i]) != 0:
                            new_labels.append(r"$\bf{" + l + "}$")
                        else:
                            new_labels.append(l)
                    else:
                        new_labels.append(l)

            ax.set_xticklabels(new_labels, rotation=45, fontsize=13)
            if l_env[0] == "Methanoculleus_A taiwanensis": l_env[0] = "Methanoculleus taiwanensis"
            if l_env[0] == "Methanolinea_B mesophila": l_env[0] = "Methanolinea mesophila"
            if l_env[0] == "Pyrococcus sp000211475": l_env[0] = "Pyrococcus furiosus"
            if l_env[0] == "Thermofilum_B adornatum": l_env[0] = "Thermofilum adornatum"
            if l_env[0] == "Thermococcus_A litoralis": l_env[0] = "Thermococcus litoralis"
            if l_env[0] == "Methanobacterium_C paludis": l_env[0] = "Methanobacterium paludis"
            if  l_env[0] == "Pseudothermotoga_B elfii":  l_env[0] = "Pseudothermotoga elfii"

            ax.set_title(f'{l_env[0]}', fontsize=18, fontstyle='italic')
            ax.set_xlabel('k-mers', fontsize=16)
            ax.set_ylabel('Deviation', fontsize=16)
            ax.tick_params(axis='x', rotation=45)
            ax.set_ylim(min_all - 0.01, max_all + 0.01)
            ax.set_facecolor('#e0f7fa' if l_domain and l_domain[0] == "Bacteria" else '#f6ddcc')

            legend_labels = ['Relevant 3-mers', 'Not Relevant 3-mers']
            custom_lines = [
                plt.Line2D([0], [0], color='#45b39d', lw=4),
                plt.Line2D([0], [0], color='#F08080', lw=4)
            ]
            ax.legend(custom_lines, legend_labels, loc='upper left', fontsize='large', prop={'style': 'normal'})

        for ax in axs[num_plots:]:
            ax.remove()

        custom_lines = [
            plt.Line2D([0], [0], color='#e0f7fa', lw=8),
            plt.Line2D([0], [0], color='#f6ddcc', lw=8)
        ]
        fig.legend(custom_lines, ['Bacteria', 'Archaea'], loc='upper right', bbox_to_anchor=(0.95, 0.94),
                   ncol=1, fontsize=12, title='Domains', title_fontsize=12, borderpad=0.8,
                   edgecolor='black', prop={'style': 'normal'})

        plt.subplots_adjust(left=0.12, right=0.85, top=0.9, bottom=0.1)
        plt.tight_layout(rect=[0.12, 0.03, 0.85, 0.95])
        plt.savefig(f"{results_folder}/kmer_importance/{c}_new.pdf", format='pdf')
        plt.show()



def save_deviation_data_to_json(clusters, temp, pH, kmers_label, features_important, labels, names, kmers, output_file):
    deviation_data = {}

    for c in clusters:
        cluster_data = {}
        for n, k in zip(names, kmers):
            if n in clusters[c]:
                sample_data = {}
                deviations = [elem2 - elem1 for elem1, elem2 in zip(temp, k)] if c != 'c10' else [elem2 - elem1 for
                                                                                                  elem1, elem2 in
                                                                                                  zip(pH, k)]
                for kmer, deviation in zip(kmers_label, deviations):
                    sample_data[kmer] = {
                        "Deviation": deviation,
                        "Relevant": kmer in features_important[n]['pos']
                    }
                cluster_data[n] = sample_data
        deviation_data[c] = cluster_data

    with open(output_file, 'w') as json_file:
        json.dump(deviation_data, json_file, indent=4)
    print(f"Deviation data saved to {output_file}")


def plot_deviations_with_heatmap(clusters, temp, pH, kmers_label, features_important, labels, names, kmers):
    cmap = plt.get_cmap('Greens')
    for c in clusters:
        ds = {}
        max_list = []
        min_list = []
        for n, k in zip(names, kmers):
            if n in clusters[c]:
                deviations = [elem2 - elem1 for elem1, elem2 in zip(temp, k)] if c != 'c10' else [elem2 - elem1 for elem1, elem2 in zip(pH, k)]
                ds[n] = deviations
                max_list.append(max(deviations))
                min_list.append(min(deviations))

        max_all = max(max_list)
        min_all = min(min_list)

        fig, axs = plt.subplots(len(ds), 1, figsize=(15, len(ds) * 5))

        for index, key in enumerate(ds):
            l_domain = list(labels[labels["Assembly"] == key]["Domain"])
            l_env = list(labels[labels["Assembly"] == key]["Temperature"])

            bars = axs[index].bar(kmers_label, ds[key], color='blue')

            importance_values = [features_important[key][kmer] for kmer in kmers_label if features_important[key][kmer] != -1]

            # Set bar color based on feature importance with a heatmap effect

            norm = mcolors.Normalize(vmin=0, vmax=max(importance_values))

            for bar, kmer in zip(bars, kmers_label):
                if features_important[key][kmer] == -1:
                    bar.set_color('#F08080')
                else:
                    importance_value = norm(features_important[key][kmer])   # Scale value
                    bar.set_color(cmap(importance_value)) # Apply heatmap coloring based on importance

            axs[index].set_title(f'{key} - {l_domain[0]} - {l_env[0]}')
            axs[index].set_xlabel('k-mers', fontsize=14)
            axs[index].set_ylabel('Deviation', fontsize=14)
            axs[index].tick_params(axis='x', rotation=45)
            axs[index].set_ylim(min_all-0.01, max_all+0.01)
            axs[index].set_facecolor('#e0f7fa')

            legend_labels = ['High Importance', 'Low Importance']
            custom_lines = [plt.Line2D([0], [0], color=cmap(norm(1)), lw=4),
                            plt.Line2D([0], [0], color=cmap(norm(0)), lw=4)]
            axs[index].legend(custom_lines, legend_labels, loc='lower left', fontsize='large')

        plt.tight_layout()
        plt.savefig(f"{results_folder}/kmer_importance/{c}_heatmap_bars.pdf", format='pdf')
        plt.show()



results_folder = "outputs/statistics"
data_folder = "outputs/exp2"

data_root = Path(data_folder).expanduser().resolve()
output_root = Path(results_folder).expanduser().resolve()

mode = "Env"
fragment_length = 100000

for env in ["Temperature"]:
    if env == "Temperature":
        clusters = {
            # "test1":["GCA_000023965.1", "GCA_000026225.1"],
            # # "test2":["GCA_000759775.1", "GCA_900112975.1", "GCA_000217815.1", "GCA_004402405.1", "GCA_000425505.1", "GCA_017874375.1", "GCA_000018365.1", "GCA_000336995.1", "GCA_000017945.1", "GCA_002214545.1"],
            # "Confirmed Pairs Group 3": ["GCA_000512735.1", "GCA_000211475.1", "GCA_000446015.1", "GCA_000725425.1", "GCA_002214605.1",
            #            "GCA_000246985.3"],
            "Confirmed Pairs Group 5": ["GCA_003568865.1", "GCA_900095385.1", "GCA_000304355.2", "GCA_001571405.1", "GCA_001602375.1", "GCA_004102725.1", "GCA_017873855.1"]
            #"Confirmed Pairs Group 2": ["GCA_000016785.1", "GCA_000789255.1"],
            #"Confirmed Pairs Group 4": ["GCA_000504085.1", "GCA_000214725.1", "GCA_000969905.1"]
            #"Confirmed Pairs Group 6": ["GCA_000147695.3", "GCA_000166095.1", "GCA_000212395.1", "GCA_003722315.1"]
        }
    elif env == "pH":
          clusters = {"Confirmed Pairs Group 1": ["GCA_000145615.1", "GCA_000317795.1"]}

    # Main script
    summary_data, labels, kmers_normalized = load_data(env, data_folder, fragment_length)

    # labels = pd.read_csv(f"{data_folder}/fragments_{fragment_length}/Temperature/Extremophiles_Temperature_Summary.tsv", sep="\t")
    names, kmers = kmersFasta(f"{data_folder}/fragments_{fragment_length}/{env}/Extremophiles_{env}.fas", k=k, transform=None, reduce=True)
    domin_dict = {names[i]: list(summary_data[summary_data["Assembly"] == names[i]]["Domain"])[0] for i in range(len(names))}
    genus_dict = {names[i]: list(summary_data[summary_data["Assembly"] == names[i]]["genus"])[0] for i in range(len(names))}
    species_dict = {names[i]: list(summary_data[summary_data["Assembly"] == names[i]]["species"])[0] for i in range(len(names))}
    k_mer_dict = {names[i]: kmers[i] for i in range(len(names))}
    correlation_groups(clusters, k_mer_dict, domin_dict, genus_dict, species_dict)


    # significance_test(summary_data, clusters, env, k_mer_dict)


    y, labels_encoder = encode_labels(labels)

    model = train_model(kmers_normalized, y)
    shap_values = compute_shap_values(model, kmers_normalized)

    features_all = extract_features(summary_data, labels_encoder, shap_values, clusters, kmers_label, env)
    save_features(features_all, f"{results_folder}/kmer_importance/features_data_{env}.json")

    with open(f"{results_folder}/kmer_importance/features_data_{env}.json", 'r') as json_file:
        features_important = json.load(json_file)

    labels = pd.read_csv(f"{data_folder}/fragments_{fragment_length}/Temperature/Extremophiles_Temperature_Summary.tsv", sep="\t")
    names, kmers = kmersFasta(f"{data_folder}/fragments_{fragment_length}/{env}/Extremophiles_{env}.fas", k=k, transform=None, reduce=True)
    output_file = f"{results_folder}/kmer_importance/k_deviation_data_{env}.json"
    save_deviation_data_to_json(clusters, temperature, pH, kmers_label, features_important, labels, names, kmers, output_file)
    plot_deviations(clusters, temperature, pH, kmers_label, features_important, labels, names, kmers, {"Confirmed Pairs Group 5": "GCA_003568865.1"})
