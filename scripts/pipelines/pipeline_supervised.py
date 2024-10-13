import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.svm import SVC
import os
from build_signature_dataset import run_fragment_builder
from scripts.models.supervised_models import run_supervised_classification
from scripts.models.supervised_models_challenging import run_supervised_classification_challenging
from scripts.models.tuning import run_supervised_classification_tuning

ENVS = ["Temperature", "pH"]
NUM_CLUSTERS = {"Temperature": 4, "pH": 2}
FRAGMENT_LENGTHS = [10000, 50000, 100000, 250000, 500000, 1000000]
DATA_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'data'))
RES_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'results'))


def experiment_task(args, env, exp, fragment_length):
    print("\n Running the pipeline is started:")

    # Building the fragments
    fragment_file = f"{args['exp_type']}/{exp}/fragments_{fragment_length}"
    print(f"\n Building fragment with length {fragment_length} is started.")
    run_fragment_builder(DATA_PATH, fragment_file, fragment_length, args['whole_genome'], env)
    print(f"\n Fragment with length {fragment_length} has been created.", flush=True)

    # Run the supervised classification under the first scenario (not challenging)

    data_folder = os.path.join(f"{DATA_PATH}/{args['exp_type']}/{exp}/fragments_{fragment_length}", env)
    result_folder = f"{RES_PATH}/{args['exp_type']}/{exp}/fragments_{fragment_length}"

    print(f"\n Classification is started (scenario 1).")
    run_supervised_classification(data_folder, args['max_k'], result_folder, env, exp, args['classifiers'])
    print(f"\n Classification ended (scenario 1).", flush=True)

    # Run the supervised classification under the 2nd scenario (challenging)
    print(f"\n Classification is started (scenario 2).")
    run_supervised_classification_challenging(data_folder, result_folder, args['max_k'], env, exp, args['classifiers'])
    print(f"\n Classification ended (scenario 2).", flush=True)

def hyper_prameter_tuning(args, env, exp, fragment_length):
    print("\n Running the pipeline is started:")
    # Building the fragments
    fragment_file = f"{args['exp_type']}/{exp}/fragments_{fragment_length}"
    print(f"\n Building fragment with length {fragment_length} is started.")
    run_fragment_builder(DATA_PATH, fragment_file, fragment_length, args['whole_genome'], env)
    print(f"\n Fragment with length {fragment_length} has been created.", flush=True)
    result_folder = f"{DATA_PATH}/{args['exp_type']}/{exp}/fragments_{fragment_length}"
    fasta_file = os.path.join(result_folder, env, f'Extremophiles_{env}.fas')
    print(f"\n Classification is started (scenario 1).")
    run_supervised_classification_tuning(fasta_file, result_folder, args['max_k'], result_folder, env, exp)
    print(f"\n Classification ended (scenario 1).", flush=True)
    
    
def run_pipeline(args):

    if args["exp_type"] == "tuning":
        num_exp = 1
        print(f"\n number of excuters: {os.cpu_count()}")
        with ProcessPoolExecutor() as executor:
            futures = []
            for env in ENVS:
                for exp in range(num_exp):
                    for fragment_length in FRAGMENT_LENGTHS:
                        # Submit tasks to the process pool
                        future = executor.submit(hyper_prameter_tuning, args, env, exp, fragment_length)
                        futures.append(future)

            # Waiting for all futures to complete
            for future in as_completed(futures):
                future.result()  # You can handle exceptions here if needed

        
    elif args["exp_type"] == "exp1":
        classifiers = {"SVM": (SVC, {'kernel': 'rbf', 'class_weight': 'balanced', 'C': 10})}
        args['classifiers'] = classifiers
        start = args["continue_exp"]
        num_exp = args["num_exp"]

    elif args["exp_type"] == "exp2":
        classifiers = {
            "SVM": (SVC, {'kernel': 'rbf', 'class_weight': 'balanced', 'C': 10}),
            "Random Forest": (RandomForestClassifier, {}),
            "ANN": (MLPClassifier, {'hidden_layer_sizes': (256, 64), 'solver': 'adam', 'activation': 'relu', 'alpha': 1, 'learning_rate_init': 0.001, 'max_iter': 300, 'n_iter_no_change': 10})
        }
        args['classifiers'] = classifiers
        start = 0
        num_exp = args["num_exp"]

    print(f"\n number of excuters: {os.cpu_count()}")
    with ProcessPoolExecutor() as executor:
        futures = []
        for env in ENVS:
            for exp in range(start, num_exp, 1):
                for fragment_length in FRAGMENT_LENGTHS:
                    # Submit tasks to the process pool
                    future = executor.submit(experiment_task, args, env, exp, fragment_length)
                    futures.append(future)

        # Waiting for all futures to complete
        for future in as_completed(futures):
            future.result()  # You can handle exceptions here if needed

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--exp_type', action='store', type=str)
    parser.add_argument('--max_k', action='store', type=int)
    parser.add_argument('--whole_genome', action='store_true')
    parser.add_argument('--num_exp', action='store', type=int)
    parser.add_argument('--continue_exp', action='store', type=int)

    args = vars(parser.parse_args())

    run_pipeline(args)

if __name__ == '__main__':
    main()
