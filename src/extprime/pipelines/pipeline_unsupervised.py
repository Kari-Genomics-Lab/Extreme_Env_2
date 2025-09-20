import argparse
from extprime.models.unsupervised_non_parametric import run_models, analyze_clustering
from extprime.analysis.candidates_identification import candidates_identification

ENVS = ["Temperature", "pH"]
NUM_CLUSTERS = {"Temperature": 4, "pH": 2}
EXP_NUM = 10


def run_pipeline(args):
    fragement_length = args["fragement_length"]
    k = args["k_mer"]
    path = args["fragment_path"]
    data_path = args["data_root"]


    if args["exp_type"] == "parametric":
        n_clusters = args["n_clusters"]
        results_folder = args["outputs_root"]
        Env = args["Env"]
        for env in ENVS:
            if Env:
                run_models(n_clusters, results_folder, env, k, "Env")
            else:
                run_models(n_clusters, results_folder, env, k, "Tax")

    elif args["exp_type"] == "non-parametric":
        results_folder = args["outputs_root"]
        for env in ENVS:
          for exp in range(EXP_NUM):
            summary_dataset, algo_names = run_models(env, path, data_path, fragement_length, k, result_folder)
            # for algo in algo_names:
            analyze_clustering(algo_names, summary_dataset, env, results_folder)
            candidates_identification(summary_dataset, env, results_folder)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--exp_type', action='store', type=str)
    parser.add_argument('--k_mer', action='store', type=int)
    parser.add_argument('--fragement_length', action='store', type=int)
    parser.add_argument('--n_clusters', action='store', type=int, default=4)
    parser.add_argument('--outputs_root', action='store', type=str)
    parser.add_argument('--env', action='store', type=str)
    parser.add_argument('--fragment_path', type=str)
    parser.add_argument('--data_root', type=str)

    args = vars(parser.parse_args())

    run_pipeline(args)

if __name__ == '__main__':
    main()
