import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import umap
import hdbscan
from sklearn.metrics import homogeneity_completeness_v_measure as cluster_quality
import os
import vamb
from src.idelucs.cluster import iDeLUCS_cluster
from src.utils import kmersFasta
from sklearn.cluster import MeanShift
from sklearn.cluster import AffinityPropagation

# Configuration and data preparation

TAX_LEVEL = 'Genus'
# Dimensionality reduction functions
def VAE(sequence_file):
    with vamb.vambtools.Reader(sequence_file) as filehandle:
        composition = vamb.parsecontigs.Composition.from_file(filehandle)

    rpkms = np.ones((composition.matrix.shape[0], 1), dtype=np.float32)
    dataloader = vamb.encode.make_dataloader(
        rpkms,
        composition.matrix,
        composition.metadata.lengths,
    )
    vae = vamb.encode.VAE(nsamples=1, nlatent=32)
    vae.trainmodel(dataloader)
    latent = vae.encode(dataloader)
    names = composition.metadata.identifiers

    return names, latent

def UMAP(sequence_file=None):
    names, kmers = kmersFasta(sequence_file, k=6)
    return names, umap.UMAP(random_state=42, n_neighbors=30, min_dist=0.0, n_components=64).fit_transform(kmers)


# Clustering functions
def IM(latent, names, sequence_file):
    # clusters: Iterable[tuple[str, Iterable[str]]], separator: str
    latent = np.float32(latent)
    # Cluster and output clusters
    with vamb.vambtools.Reader(sequence_file) as filehandle:
        composition = vamb.parsecontigs.Composition.from_file(filehandle)
    clusters = {}
    labels = [-1] * len(names)
    for (n, cluster) in enumerate(vamb.cluster.ClusterGenerator(latent, composition.metadata.lengths)):
        members = list(cluster.members)
        for i in members:
            labels[i] = n

    return labels


def meanshift(latent, names):
    mean_shift = MeanShift(bandwidth=1).fit(latent)
    labels = mean_shift.labels_

    return labels


def affinity_propagation(latent, names):
    clustering = AffinityPropagation(random_state=5).fit(latent)
    labels = clustering.labels_
    return labels


def DBSCAN(latent, names):
    clusterer = hdbscan.HDBSCAN(min_cluster_size=3, gen_min_span_tree=True).fit(latent)
    return clusterer.labels_


def iDeLUCS(sequence_file=None, params=None):
    model = iDeLUCS_cluster(**params)
    return model.fit_predict(sequence_file)


def insert_assignment(summary_dataset, assignment_algo, GT_file, labels):
    """
    Updates the summary dataset with new cluster assignments and calculates cluster quality metrics using the specified taxonomic level.

    Parameters:
        summary_dataset (DataFrame): The dataset to update.
        assignment_algo (str): The name of the algorithm used for the assignment.
        names (list): List of sequence names or identifiers.
        labels (list): List of labels corresponding to the names.
    """

    summary_dataset[assignment_algo] = pd.Series(dtype='int')
    df = pd.read_csv(GT_file, sep='\t')

    for i, acc in enumerate(df["Assembly"].values):
        summary_dataset.loc[summary_dataset["Assembly"] == acc, assignment_algo] = labels[i]

    # Calculate and print the clustering quality metrics.
    if TAX_LEVEL in summary_dataset:
        GT = summary_dataset[TAX_LEVEL].dropna().to_numpy()
        unique_labels = np.unique(GT)
        label_mapping = {label: idx for idx, label in enumerate(unique_labels)}
        y = np.array([label_mapping[genus] for genus in GT])
        y_pred = summary_dataset[assignment_algo].to_numpy()
        results = cluster_quality(y, y_pred)
        summary_dataset.to_csv("assigned_summary.csv")

        print(f'Clustering Quality for {assignment_algo}:', results)
        return summary_dataset
    else:
        print(f'Column "{TAX_LEVEL}" not found in summary_dataset for clustering quality calculation.')


def run_models(env, path, fragment_length, k):
    result_folder = f"{path}/exp2/0/fragments_{fragment_length}"

    summary_file = f'{path}/Extremophiles_GTDB.tsv'
    summary_dataset = pd.read_csv(summary_file, sep='\t',
                                  usecols=["Assembly", "Domain", env,
                                           "pH", "Phylum", "Class", "Order",
                                           "Family", "Genus", "Species"])
    summary_dataset.dropna(subset=[env], inplace=True)
    GT_file = os.path.join(result_folder, env, f'Extremophiles_{env}_GT_Tax.tsv')

    sequence_file = os.path.join(result_folder, env, f'Extremophiles_{env}.fas')
    idelucs_params = {'sequence_file': sequence_file, 'n_clusters': 200, 'n_epochs': 50, 'n_mimics': 3,
                      'batch_sz': 512, 'k': 6, 'weight': 0.25, 'n_voters': 5}

    models = {
        "VAE": VAE,
        "UMAP": UMAP,
        "iDeLUCS": ""
    }

    clust_algorithms = {"IM": IM, "HDBSCAN": DBSCAN, "affinity_propagation": affinity_propagation,
                        "meanshift": meanshift}

    for model_name, model_func in models.items():
        print(f"......... Processing {model_name} ...............")
        if model_name == "iDeLUCS":
            model = iDeLUCS_cluster(**idelucs_params)
            labels, latent = model.fit_predict(sequence_file)

        else:
            names, latent = model_func(sequence_file)
            for clust_name, clust_func in clust_algorithms.items():
                if clust_name == "IM":
                    labels = clust_func(latent, names, sequence_file)
                    assignment_algo = f'{model_name}+{clust_name}'
                    insert_assignment(summary_dataset, assignment_algo, GT_file, labels)
                else:
                    labels = clust_func(latent, names)
                    assignment_algo = f'{model_name}+{clust_name}'
                    insert_assignment(summary_dataset, assignment_algo, GT_file, labels)
        if model_name == "iDeLUCS":
            insert_assignment(summary_dataset, model_name, GT_file, labels)

    algo_names = [f'{dim_name}+{clust_name}' for dim_name in models for clust_name in clust_algorithms if
                  dim_name != "iDeLUCS"]
    algo_names.append('iDeLUCS')  # Add the standalone iDeLUCS algorithm

    return summary_dataset, algo_names


def analyze_result(summary_dataset, dataset_column):
    """
    Analyzes the dataset to count significant taxa for each environmental condition
    and generates a list of algorithm names based on dimensionality and clustering methods.

    Parameters:
        summary_dataset (pd.DataFrame): The DataFrame containing the dataset.
        dataset_column (str): The name of the column in the DataFrame to analyze.
        dim_algorithms (list): A list of tuples representing dimensionality reduction algorithms.
        clust_algorithms (list): A list of tuples representing clustering algorithms.
    """
    env_taxonomic_counts = {}

    for env_value in summary_dataset[dataset_column].unique():
        env_subset = summary_dataset[summary_dataset[dataset_column] == env_value]
        taxonomic_counts = env_subset[TAX_LEVEL].value_counts()
        significant_taxa = taxonomic_counts[taxonomic_counts >= 2]
        env_taxonomic_counts[env_value] = len(significant_taxa)

    df_env_taxa = pd.DataFrame(list(env_taxonomic_counts.items()), columns=[dataset_column, 'Significant Taxa Count'])

    return df_env_taxa


def analyze_clustering(algo_names, summary_dataset, env):

    GOOD_CLUSTERS = {}
    df = analyze_result(summary_dataset, env)

    print(algo_names)
    for name in algo_names:
        results_df = pd.DataFrame({
            name: summary_dataset[name],
            'true_genus': summary_dataset[TAX_LEVEL],
            'mode_cluster': summary_dataset.groupby(name)[TAX_LEVEL].transform(
                lambda x: x.value_counts().idxmax()),
            'cluster_size': summary_dataset.groupby(name)[TAX_LEVEL].transform('count')
        })

        # Calculate the true size for each mode cluster
        genus_size = summary_dataset.groupby(TAX_LEVEL)[TAX_LEVEL].count().to_dict()
        results_df['true_size'] = results_df['mode_cluster'].map(genus_size)

        # Calculate the mode size for each cluster
        mode_size = summary_dataset[results_df['mode_cluster'] == results_df['true_genus']].groupby(name)[
            name].count().to_dict()

        results_df['mode_size'] = results_df[name].map(mode_size)

        # Calculate completeness and contamination
        results_df['completeness'] = np.minimum(results_df['mode_size'], results_df['true_size']) / results_df[
            'true_size']
        results_df['contamination'] = (results_df['cluster_size'] - results_df['mode_size']) / results_df[
            'cluster_size']

        # Update results_df with dataset info
        unique_dataset = summary_dataset.drop_duplicates(subset=[TAX_LEVEL]).set_index(TAX_LEVEL)[env]
        results_df[env] = results_df['mode_cluster'].map(unique_dataset)

        # Filtering results for visualization
        clustering_results = results_df[results_df['cluster_size'] >= 2]
        HQ = clustering_results[
            (clustering_results['completeness'] >= 0.50) & (clustering_results['contamination'] <= 0.50)]
        GOOD_CLUSTERS[name] = HQ[name].values

        # Count and display results
        gt_counts = HQ.groupby(env)['mode_cluster'].nunique()
        df2 = pd.DataFrame({env: gt_counts.index, 'Count': gt_counts.values})
        print(df2)
        print(df2['Count'].sum())

        # Update df for plotting
        df[name] = np.zeros(len(df))
        look = df2.set_index(env)['Count'].to_dict()
        for key in look:
            df.loc[df[env] == key, name] = look[key]

    # Plot the merged dataframe
    plot_results(df, env)


def plot_results(df, env):
    df.plot.bar(x=env, rot=0, figsize=(8, 6.8), legend=True, width=0.4, fontsize=14)
    plt.title("No. True Genera vs. No. Recovered Genera by each Algorithm", fontsize=16)
    plt.ylabel("No. Genera", fontsize=14)
    plt.ylim(0, 37)
    plt.xlabel(env, fontsize=14)
    plt.savefig(f"Unsupervised_{env}.png", format="png")



