import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
import shap
from scripts.src.utils import kmersFasta
import json
import matplotlib.pyplot as plt
import os
import matplotlib.colors as mcolors

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


def save_features(features_all, file_path):
    with open(file_path, 'w') as json_file:
        json.dump(features_all, json_file, indent=4)

def plot_deviations(clusters, temp, pH, kmers_label, features_important, labels, names, kmers):
    for c in clusters:
        ds = {}
        max_list = []
        min_list = []
        for n, k in zip(names, kmers):
            if n in clusters[c]:
                deviations = [elem2 - elem1 for elem1, elem2 in zip(temp, k)] if c != 'c10' else [elem2 - elem1 for
                                                                                                  elem1, elem2 in
                                                                                                  zip(pH, k)]
                ds[n] = deviations
                max_list.append(max(deviations))
                min_list.append(min(deviations))

        max_all = max(max_list)
        min_all = min(min_list)

        num_plots = len(ds)
        num_cols = 2
        num_rows = (num_plots + num_cols - 1) // num_cols  # Calculate number of rows needed

        fig, axs = plt.subplots(num_rows, num_cols, figsize=(25, num_rows * 5))
        axs = axs.flatten()  # Flatten the grid for easy iteration

        for index, (key, ax) in enumerate(zip(ds, axs)):
            l_domain = list(labels[labels["Assembly"] == key]["Domain"])
            l_env = list(labels[labels["Assembly"] == key]["Temperature"])
            l_env = list(labels[labels["Assembly"] == key]["species"])
            bars = ax.bar(kmers_label, ds[key], color='blue')
            for bar, l in zip(bars, kmers_label):
                bar.set_color('#45b39d' if l in features_important[key]['pos'] else '#F08080')
            ax.set_title(f'{l_env[0]}', fontsize=18, fontstyle='italic')
            ax.set_xlabel('k-mers', fontsize=16)
            ax.set_ylabel('Deviation', fontsize=16)
            ax.tick_params(axis='x', rotation=45)
            ax.set_ylim(min_all - 0.01, max_all + 0.01)
            if l_domain[0] == "Bacteria":
                ax.set_facecolor('#e0f7fa')
            else:
                ax.set_facecolor('#f6ddcc')

            legend_labels = ['Relevant 3-mers', 'Not Relevant 3-mers']
            custom_lines = [plt.Line2D([0], [0], color='#45b39d', lw=4),
                            plt.Line2D([0], [0], color='#F08080', lw=4)]
            ax.legend(custom_lines, legend_labels, loc='lower left', fontsize='large', prop={'style': 'normal'})

        # Hide unused subplots if any
        for ax in axs[num_plots:]:
            ax.remove()  # Remove any extra subplots that are not used

        # Create a single legend for the domains
        legend_labels_all = ['Bacteria', 'Archaea']
        custom_lines = [plt.Line2D([0], [0], color='#e0f7fa', lw=8),
                        plt.Line2D([0], [0], color='#f6ddcc', lw=8)]

        fig.legend(custom_lines, legend_labels_all, loc='upper right', bbox_to_anchor=(0.95, 0.94), ncol=1,
                   fontsize=12,
                   title='Domains', title_fontsize=12, borderpad=0.8, edgecolor='black', prop={'style': 'normal'})

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



results_folder = os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..', 'results'))
data_folder = os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..', 'data'))
mode = "Env"
fragment_length = 100000

for env in ["Temperature", "pH"]:
    if env == "Temperature":
        clusters = {
            "c3": ["GCA_000512735.1", "GCA_000211475.1", "GCA_000446015.1", "GCA_000725425.1", "GCA_002214605.1",
                   "GCA_000246985.3"],
            "c6": ["GCA_003568865.1", "GCA_000015825.1", "GCA_000275865.1", "GCA_000304355.2", "GCA_001017125.1",
                   "GCA_001571405.1", "GCA_001602375.1", "GCA_002503885.1", "GCA_004102725.1", "GCA_017873855.1",
                   "GCA_900095385.1"],
            "c11": ["GCA_000016785.1", "GCA_015163485.1", "GCA_000789255.1"],
            "c12": ["GCA_000504085.1", "GCA_002287235.1", "GCA_000214725.1", "GCA_000969905.1", "GCA_017873625.1"]
        }
    elif env == "pH":
        clusters = {"c10": ["GCA_000145615.1", "GCA_000317795.1"]}

    # Main script
    summary_data, labels, kmers_normalized = load_data(env, data_folder, fragment_length)
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
    plot_deviations(clusters, temperature, pH, kmers_label, features_important, labels, names, kmers)
