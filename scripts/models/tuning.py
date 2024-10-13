from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score


import json
import os
import pandas as pd
import numpy as np
from src.utils import kmersFasta

import warnings
from multiprocessing import Lock

warnings.filterwarnings("ignore")
lock = Lock()




def run_supervised_classification_tuning(fasta_file, max_k, result_folder, env, exp):
    results_json = {}
    for k in range(1, max_k + 1):
        results_json[k] = {}
        _, kmers = kmersFasta(fasta_file, k=k, transform=None, reduce=True)
        kmers_normalized = np.transpose((np.transpose(kmers) / np.linalg.norm(kmers, axis=1)))
        results_json = perform_classification(kmers_normalized, k, results_json, result_folder, env)
        print(f"Finished processing k = {k}", flush=True)
        del kmers_normalized



def perform_classification(kmers, k, results_json, result_folder, env):

    env_file = os.path.join(result_folder, env, f'Extremophiles_{env}_GT_Env.tsv')
    tax_file = os.path.join(result_folder, env, f'Extremophiles_{env}_GT_Tax.tsv')

    for index, label_file in enumerate([env_file, tax_file]):
        label_data = pd.read_csv(label_file, sep='\t')
        run_hyper_parameter_tuning(kmers, label_data['cluster_id'], k, env, index, result_folder)


    return results_json



def run_hyper_parameter_tuning(X, y, k, env, index, result_folder):
    #
    # unique_labels = list(label_data['cluster_id'].unique())
    # label_index_map = {label: idx for idx, label in enumerate(unique_labels)}

    # Split the data
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # SVM Hyperparameter tuning
    svm_params = {
        'C': [0.1, 1, 10, 100],
        'gamma': [1, 0.1, 0.01, 0.001],
        'kernel': ['rbf', 'linear']
    }
    svm = SVC()
    svm_grid = GridSearchCV(svm, svm_params, cv=5, scoring='accuracy', n_jobs=-1)
    svm_grid.fit(X_train, y_train)

    print(f"Best parameters for SVM: {svm_grid.best_params_}")
    print(f"Best SVM accuracy: {svm_grid.best_score_}")

    # ANN Hyperparameter tuning
    ann_params = {
        'hidden_layer_sizes': [(10,), (50,), (100,)],
        'activation': ['tanh', 'relu'],
        'solver': ['sgd', 'adam'],
        'alpha': [0.0001, 0.05],
        'learning_rate': ['constant', 'adaptive']
    }
    ann = MLPClassifier(max_iter=200)
    ann_grid = GridSearchCV(ann, ann_params, cv=5, scoring='accuracy', n_jobs=-1)
    ann_grid.fit(X_train, y_train)

    print(f"Best parameters for ANN: {ann_grid.best_params_}")
    print(f"Best ANN accuracy: {ann_grid.best_score_}")

    # Random Forest Hyperparameter tuning
    rf_params = {
        'n_estimators': [10, 50, 100],
        'max_features': ['auto', 'sqrt'],
        'max_depth': [10, 50, 100, None],
        'criterion': ['gini', 'entropy']
    }
    rf = RandomForestClassifier()
    rf_grid = GridSearchCV(rf, rf_params, cv=5, scoring='accuracy', n_jobs=-1)
    rf_grid.fit(X_train, y_train)

    print(f"Best parameters for Random Forest: {rf_grid.best_params_}")
    print(f"Best Random Forest accuracy: {rf_grid.best_score_}")

    # Evaluating on the test set
    best_svm = svm_grid.best_estimator_
    best_ann = ann_grid.best_estimator_
    best_rf = rf_grid.best_estimator_

    svm_pred = best_svm.predict(X_test)
    ann_pred = best_ann.predict(X_test)
    rf_pred = best_rf.predict(X_test)

    print(f"Test Accuracy for SVM: {accuracy_score(y_test, svm_pred)}")
    print(f"Test Accuracy for ANN: {accuracy_score(y_test, ann_pred)}")
    print(f"Test Accuracy for Random Forest: {accuracy_score(y_test, rf_pred)}")

    best_params = {
        'SVM': svm_grid.best_params_,
        'ANN': ann_grid.best_params_,
        'RandomForest': rf_grid.best_params_
    }
    file_path = os.path.join(result_folder, f'best_params_{k}_{env}_{index}.json')
    with open(file_path, 'w') as f:
        json.dump(best_params, f, indent=4)


