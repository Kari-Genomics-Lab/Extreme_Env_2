import matplotlib.pyplot as plt
import json
import os
from FCGR_analysis import run
# This script is to get the candidates that are similar in the majority distance metrics

distance_metrics = ["descriptor", "dssim", "lpips"]
thresholds = [0.20, 0.50, 0.50]
majority_num = 3
RES_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..', 'results'))
DATA_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..', 'data'))


def plot_dist(dist_list, dist_metric):
    x_values = list(range(1, len(dist_list) + 1))
    # Create bar chart
    plt.bar(x_values, dist_list)
    plt.title(f'candidates bar chart plot/ temp / {dist_metric}')
    plt.xlabel('Index')
    plt.ylabel('distances')
    plt.grid(True)
    plt.show()

def majority_candidates():
    dist_data = {}
    limited_dist = {}
    for index, d in enumerate(distance_metrics):
        limited_dist[d] = []
        with open(f"{RES_PATH}/distances/candidates_{d}_all.json", 'r') as json_file:

            dist_data[d] = json.load(json_file)

            for id in dist_data[d]:
                n1, n2 = id.split("_")[1], id.split("_")[3]
                if n1 == "000145615.1" and n2 == "000317795.1":
                    print(dist_data[d][id])
                # dist = float("{:.2f}".format(dist_data[d][id]))
                if dist_data[d][id] <= thresholds[index]:
                    limited_dist[d].append(id)

            print(d, len(limited_dist[d]))

    majority = {}
    for i in limited_dist:

        limited_dist[i] = set(limited_dist[i])
    for i in dist_data[d]:
        count = 0
        if i in limited_dist[distance_metrics[0]]:
            count += 1
        if i in limited_dist[distance_metrics[1]]:
            count += 1
        if i in limited_dist[distance_metrics[2]]:
            count += 1
        if count >= majority_num:
            majority[i] = {d: dist_data[d][i] for d in distance_metrics if i in limited_dist[d]}

    print(len(majority))
    json.dump(majority, open(f"{RES_PATH}/candidates/majority_candidates_new_threshold.json", "w"), indent=4)
    return majority

fragment_length = 100000
candidates = majority_candidates(RES_PATH)
for env in ["Temperature", "pH"]:
    id_2_sequences = run(candidates, env, fragment_length, results_folder, data_folder, 'selected')