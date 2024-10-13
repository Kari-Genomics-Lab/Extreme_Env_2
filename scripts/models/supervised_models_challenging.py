import json
import os
import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedGroupKFold
from scripts.src.utils import SummaryFasta, kmersFasta
import warnings
from multiprocessing import Lock

warnings.filterwarnings("ignore")
lock = Lock()

def save_results(result, dataset, result_folder, run):
    file_path = os.path.join(result_folder, f'Challenging_Supervised_Results_{dataset}.json')
    lock.acquire()  # Acquire lock before accessing the file
    try:
        existing_data = {}
        if os.path.isfile(file_path):
            with open(file_path, 'r') as file:
                existing_data = json.load(file)
        if run in existing_data.keys():
            existing_data = {}
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


def preprocess_data(fasta_file, summary_file):
    names, _, _, _ = SummaryFasta(fasta_file)
    summary_dataset = pd.read_csv(summary_file, sep='\t',
                                  usecols=["Domain", "Radio_label", "Assembly", "Genus", "Species"])
    assembly_dict = {name: summary_dataset[summary_dataset["Assembly"] == name].iloc[0]
                     for name in names if summary_dataset["Assembly"].eq(name).any()}
    data = pd.DataFrame({
        "Domain": [info["Domain"] for info in assembly_dict.values()],
        "Genus": [info["Genus"] for info in assembly_dict.values()],
        "Species": [info["Species"] for info in assembly_dict.values()],
        "Assembly": [info["Assembly"] for info in assembly_dict.values()],
        "sequence_id": list(assembly_dict.keys())
    })
    return data


def run_supervised_classification_challenging(path, fasta_file, max_k, result_folder, env, exp, classifiers):

    data = preprocess_data(fasta_file, f"{path}/Extremophiles_GTDB_Radio.tsv")
    results_json = {}
    for k in range(1, max_k + 1):
        results_json[k] = {}
        _, kmers = kmersFasta(fasta_file, k=k, transform=None, reduce=True)
        kmers_normalized = np.transpose((np.transpose(kmers) / np.linalg.norm(kmers, axis=1)))
        results_json = perform_classification(kmers_normalized, k, results_json, result_folder, env, data, classifiers)
        print(f"Finished processing k = {k}", flush=True)
        del kmers_normalized
    save_results(results_json, env, result_folder, exp)


def perform_classification(kmers, k, results_json, result_folder, env, data, classifiers):

    env_file = os.path.join(result_folder, env, f'Extremophiles_{env}_GT_Env.tsv')
    tax_file = os.path.join(result_folder, env, f'Extremophiles_{env}_GT_Tax.tsv')

    for name, (algorithm, params) in classifiers.items():
        results_json[k][name] = [0, 0]
        for index, label_file in enumerate([env_file, tax_file]):
            label_data = pd.read_csv(label_file, sep='\t')
            results_json[k][name][index] = cross_validate_model(kmers, label_data, algorithm, params, data)

    return results_json


def cross_validate_model(kmers, label_data, algorithm, params, data):
    label_data.reset_index(drop=True, inplace=True)  # Reset index if not already aligned
    kmers_df = pd.DataFrame(kmers)
    Dataset = pd.concat([label_data['cluster_id'], kmers_df], axis=1).dropna()

    # Prepare data for cross-validation
    skf = StratifiedGroupKFold(n_splits=10)
    unique_labels = list(label_data['cluster_id'].unique())
    label_index_map = {label: idx for idx, label in enumerate(unique_labels)}

    scores = []

    for train_idx, test_idx in skf.split(Dataset, Dataset['cluster_id'], groups=data["Genus"].values):
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
    return average_score
