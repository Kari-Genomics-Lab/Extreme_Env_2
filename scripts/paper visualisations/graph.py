import argparse

import pandas as pd
import json
import matplotlib.pyplot as plt
import numpy as np
import os

PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'results/exp2'))

font = {'family': 'serif',
        'weight': 'bold',
        'size': 14}
# Function to load and transform JSON data into a DataFrame
def load_json(file_paths):
    all_data = pd.DataFrame()
    for file_label, file_path in file_paths.items():
        with open(f"{PATH}/{file_path}", 'r') as file:
            data = json.load(file)

        results = []
        for k, models in data.items():
            for model, scores in models.items():
                print(scores)
                results.append({
                    'file': file_label,
                    'k': int(k),
                    'model': model,
                    'env_score': float('{0:.2f}'.format(scores['env'])),
                    'taxa_score': float('{0:.2f}'.format(scores['tax']))
                })

        temp_df = pd.DataFrame(results)
        all_data = pd.concat([all_data, temp_df], ignore_index=True)

    return all_data


def run(env, scenario="normal"):
    if scenario == "normal":
        file_paths = {
            '10000': f'fragments_10000/all_results_{env}.json',
            '50000': f'fragments_50000/all_results_{env}.json',
            '100000': f'fragments_100000/all_results_{env}.json',
            '250000': f'fragments_250000/all_results_{env}.json',
            '500000': f'fragments_500000/all_results_{env}.json',
            '1000000': f'fragments_1000000/all_results_{env}.json'
        }
    if scenario == "challenging":
        file_paths = {
            '10000': f'fragments_10000/Challenging_all_results_{env}.json',
            '50000': f'fragments_50000/Challenging_all_results_{env}.json',
            '100000': f'fragments_100000/Challenging_all_results_{env}.json',
            '250000': f'fragments_250000/Challenging_all_results_{env}.json',
            '500000': f'fragments_500000/Challenging_all_results_{env}.json',
            '1000000': f'fragments_1000000/Challenging_all_results_{env}.json'
        }

    # Load all data
    df = load_json(file_paths)
    # Initialize a list to store best k values and scores for each model and file

    # Assuming 'df' is already loaded and prepared from previous steps
    # models = df['model'].unique()
    models = df['model'].unique()
    files = df['file'].unique()
    find_best_k_for_each_fragment(df, env, scenario)
    find_best_fragment_for_each_k(df, env, scenario)
    create_plots('env', models, files, df, env, scenario)
    create_plots('taxa', models, files, df, env, scenario)


def find_best_k_for_each_fragment(df, env, scenario):
    best_ks = []
    # print(df)
    # Group data by file and model to compute best k values
    for file_label, file_group in df.groupby('file'):

        for model, model_group in file_group.groupby('model'):
            max_acc_env = max(list(model_group['env_score']))
            best_k_env = list(model_group[model_group["env_score"] == max_acc_env]['k'])
            max_acc_taxa = max(list(model_group['taxa_score']))
            best_k_taxa = list(model_group[model_group["taxa_score"] == max_acc_taxa]['k'])

            best_ks.append({
                'length': file_label,
                'model': model,
                'best_k_env': best_k_env,
                'env_score_at_best_k': max_acc_env,
                'best_k_taxa': best_k_taxa,
                'taxa_score_at_best_k': max_acc_taxa
            })

    # Convert list of best ks to DataFrame
    best_ks_df = pd.DataFrame(best_ks)

    # Save to CSV
    best_ks_df.to_csv(f'./results/best_k_values_and_scores_{env}_{scenario}.csv', index=False)
    return best_ks_df


def find_best_fragment_for_each_k(df, env, scenario):
    best_lengths = []
    # print(df)
    # Group data by file and model to compute best k values
    for k_label, k_group in df.groupby('k'):

        for model, model_group in k_group.groupby('model'):
            max_acc_env = max(list(model_group['env_score']))
            best_k_env = list(model_group[model_group["env_score"] == max_acc_env]['file'])
            max_acc_taxa = max(list(model_group['taxa_score']))
            best_k_taxa = list(model_group[model_group["taxa_score"] == max_acc_taxa]['file'])

            best_lengths.append({
                'k': k_label,
                'model': model,
                'best_length_env': best_k_env,
                'env_score_at_best_length': max_acc_env,
                'best_length_taxa': best_k_taxa,
                'taxa_score_at_best_length': max_acc_taxa
            })

    # Convert list of best ks to DataFrame
    best_lengths_df = pd.DataFrame(best_lengths)

    # Save to CSV
    best_lengths_df.to_csv(f'./results/best_lengths_values_and_scores_{env}_{scenario}.csv', index=False)
    return best_lengths_df


# Function to create plots for each score type
def create_plots(score_type, models, files, df, env, senario):
    print(models)
    # Create a figure with a subplot for each model
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(18, 12))
    # Adjusted for better aspect ratio
    # fig.suptitle(f'Performance Across k Values for {score_type.capitalize()} Scores', fontsize=16, y=1.02)  # Adjust for better title placement

    # Enhanced aesthetics
    markers = ['o', '^', 's', 'p', '*', 'h', 'x']  # Different markers for each file
    colors = plt.cm.viridis(np.linspace(0, 1, len(files)))  # Use a perceptually uniform colormap

    # Determine the global y-axis limits
    all_scores = df[score_type + '_score'].dropna().astype(float)
    min_score, max_score = all_scores.min(), all_scores.max()
    margin = (max_score - min_score) * 0.05  # Tighter margin for publication
    y_min, y_max = min_score - margin, max_score + margin
    # Flatten the axes array for easier iteration
    models = ['SVM', 'ANN', 'RF']

    axes = axes.flatten()
    for ax, model in zip(axes, models):
        model_group = df[df['model'] == model]
        for file_idx, file_label in enumerate(files):
            file_group = model_group[model_group['file'] == file_label]
            line, = ax.plot(file_group['k'], file_group[score_type + '_score'].astype(float),
                    marker=markers[file_idx % len(markers)],
                    linestyle='-', label=f'{file_label}', color=colors[file_idx])

        ax.set_title(f'{model} Model Classification', fontsize=14)
        ax.set_xlabel('k-mer Length', fontsize=14)
        ax.set_ylabel('Accuracy (%)', fontsize=14)
        ax.set_ylim(y_min, y_max)  # Set uniform y-limits
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)  # Finer grid
        ax.legend(title='Fragment Length', fontsize=10)
    plt.subplots_adjust(left=0.06, right=0.97, top=0.85, bottom=0.15)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(f"./graph/{score_type}_{env}_{senario}.pdf", format='pdf', dpi=300)
    plt.tight_layout()
    plt.show()


def create_plots_one_model():
    font = {'family': 'sans-serif',
            'sans-serif': 'Helvetica',
            'weight': 'bold',
            'size': 14}

    # plt.rc('font', **font)
    plt.rc('axes', titlesize=16)
    plt.rc('axes', titlepad=20)  # Increased title padding for subplot titles

    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(15, 12))

    markers = ['o', '^', 's', 'p', '*', 'h', 'x']
    colors = plt.cm.viridis(np.linspace(0, 1, 6))

    model = 'SVM'
    axes = axes.flatten()

    ENVS = [("Temperature", "taxa"), ("Temperature", "env"), ("pH", "taxa"), ("pH", "env")]
    handles_labels = []
    index = 0
    letters = ['a)', 'b)', 'c)', 'd)']

    for ax, e in zip(axes, ENVS):
        handles_labels = []
        file_paths = {
            '10000': f'fragments_10000/Challenging_all_results_{e[0]}.json',
            '50000': f'fragments_50000/Challenging_all_results_{e[0]}.json',
            '100000': f'fragments_100000/Challenging_all_results_{e[0]}.json',
            '250000': f'fragments_250000/Challenging_all_results_{e[0]}.json',
            '500000': f'fragments_500000/Challenging_all_results_{e[0]}.json',
            '1000000': f'fragments_1000000/Challenging_all_results_{e[0]}.json'
        }

        df = load_json(file_paths)
        files = df['file'].unique()
        model_group = df[df['model'] == model]

        for file_idx, file_label in enumerate(files):
            file_group = model_group[model_group['file'] == file_label]
            line, = ax.plot(file_group['k'], file_group[e[1] + '_score'].astype(float),
                            marker=markers[file_idx % len(markers)],
                            linestyle='-', label=f'{file_label}', color=colors[file_idx])
            handles_labels.append((line, file_label))

        if e[1] == "env":
            ax.set_title(f'{letters[index]} {e[0]} Dataset - Environment Label', fontsize=17, loc='center')
        elif e[1] == "taxa":
            ax.set_title(f'{letters[index]} {e[0]} Dataset - Taxonomy Label', fontsize=17, loc='left')

        ax.set_xlabel('k-mer Length', fontsize=14)
        ax.set_ylabel('Accuracy (%)', fontsize=14)
        ax.set_ylim(30, 105)
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        # ax.annotate(letters[index], xy=(0.05, 0.95), xycoords='axes fraction', fontsize=14, fontweight='bold', va='top',
        #             ha='right')
        index += 1

    handles, labels = zip(*handles_labels)
    fig.legend(handles, labels, loc='upper left', bbox_to_anchor=(0.01, 0.92), ncol=1, fontsize=12,
               title='Fragment Length', title_fontsize=14, edgecolor = 'black')

    plt.subplots_adjust(left=0.12, right=0.97, top=0.9, bottom=0.1)
    plt.tight_layout(rect=[0.12, 0.03, 1, 0.95])
    plt.savefig("./scripts/results/graph/svm.pdf", format='pdf', dpi=300)

    plt.show()


# Create plots for 'env' and 'taxa' scores

create_plots_one_model()
ENVS = ["Temperature", "pH"]
# ENVS = ["Radio_label"]
for env in ENVS:
    run(env, "normal")
    run(env, "challenging")
    # create_plots_one_model()


