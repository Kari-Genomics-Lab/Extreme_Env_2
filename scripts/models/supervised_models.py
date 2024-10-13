import json
import os
import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold
from scripts.src.utils import kmersFasta

import warnings
from multiprocessing import Lock

warnings.filterwarnings("ignore")

lock = Lock()

def save_results(result_folder, result, dataset, run):

    os.makedirs(result_folder, exist_ok=True)
    file_path = os.path.join(result_folder, f'Supervised_Results_{dataset}.json')
    lock.acquire()  # Acquire lock before accessing the file
    try:
        existing_data = {}
        if os.path.isfile(file_path):
            with open(file_path, 'r') as file:
                existing_data = json.load(file)
        existing_data[run] = result
        with open(file_path, 'w') as file:
            json.dump(existing_data, file, indent=2)
    finally:
        lock.release()  # Ensure the lock is always released


def load_json_results(path, continue_flag):
    if os.path.isfile(path) and continue_flag:
        with open(path, 'r') as file:
            data = json.load(file)
        return {int(k): v for k, v in data.items()}
    return {}


def run_supervised_classification(data_folder, result_folder, max_k, env, exp, classifiers):
    results_json = {}
    fasta_file =  os.path.join(data_folder, f'Extremophiles_{env}.fas')
    for k in range(1, max_k + 1):
        results_json[k] = {}
        _, kmers = kmersFasta(fasta_file, k=k, transform=None, reduce=True)
        kmers_normalized = np.transpose((np.transpose(kmers) / np.linalg.norm(kmers, axis=1)))
        results_json = perform_classification(data_folder, kmers_normalized, k, results_json, env, classifiers)
        print(f"Finished processing k = {k}", flush=True)
        del kmers_normalized
    save_results(result_folder, results_json, env, exp)


def perform_classification(data_folder, kmers, k, results_json, env, classifiers):

    env_file = os.path.join(data_folder, env, f'Extremophiles_{env}_GT_Env.tsv')
    tax_file = os.path.join(data_folder, env, f'Extremophiles_{env}_GT_Tax.tsv')

    for name, (algorithm, params) in classifiers.items():
        results_json[k][name] = [0, 0]
        for index, label_file in enumerate([env_file, tax_file]):
            label_data = pd.read_csv(label_file, sep='\t')
            results_json[k][name][index] = cross_validate_model(kmers, label_data, algorithm, params)

    return results_json


def cross_validate_model(kmers, label_data, algorithm, params):
    label_data.reset_index(drop=True, inplace=True)  # Reset index if not already aligned
    kmers_df = pd.DataFrame(kmers)
    Dataset = pd.concat([label_data['cluster_id'], kmers_df], axis=1).dropna()

    # Prepare data for cross-validation
    skf = StratifiedKFold(n_splits=10, shuffle=True)
    unique_labels = list(label_data['cluster_id'].unique())
    label_index_map = {label: idx for idx, label in enumerate(unique_labels)}

    scores = []

    for train_idx, test_idx in skf.split(Dataset, Dataset['cluster_id']):
        model = algorithm(**params)
        x_train, y_train = Dataset.iloc[train_idx].drop('cluster_id', axis=1), Dataset.iloc[train_idx]['cluster_id']
        x_test, y_test = Dataset.iloc[test_idx].drop('cluster_id', axis=1), Dataset.iloc[test_idx]['cluster_id']

        # Convert labels to indices
        y_train = [label_index_map[label] for label in y_train]
        y_test = [label_index_map[label] for label in y_test]

        model.fit(x_train, y_train)
        scores.append(model.score(x_test, y_test))

        del model

    average_score = sum(scores) / len(scores)
    print(average_score)
    return average_score
