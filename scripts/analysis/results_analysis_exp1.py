import os
import json
import numpy as np

# Base directory where your folders are extracted
base_dir = '../results/supervised_classification/exp1/'

# Experiment folders
experiment_folders = ['0', '1', '2', '3', '4']

# Fragment sizes
fragment_sizes = ['10000', '50000', '100000', '250000', '500000', '1000000']

def statistical_resutls(env):
    # Initialize a dictionary to store all values for averaging and variance calculation
    data_storage = {fragment: {str(k): {'env_values': [], 'tax_values': []} for k in range(1, 10)} for fragment in
                    fragment_sizes}

    # Iterate through each experiment folder
    for exp_folder in experiment_folders:
        exp_folder_path = os.path.join(base_dir, exp_folder)

        # Iterate through each fragment size
        for fragment in fragment_sizes:
            fragment_path = os.path.join(exp_folder_path, f'fragments_{fragment}')
            temperature_file = f'Challenging_Supervised_Results_{env}.json'

            # Construct full path to the JSON file
            json_path = os.path.join(fragment_path, temperature_file)

            # Check if the file exists and read data
            if os.path.exists(json_path):
                with open(json_path, 'r') as file:
                    data = json.load(file)[exp_folder]  # Assuming '0' is the fixed experiment identifier

                # Store values for each k-mer size
                for k in range(1, 10):
                    k_str = str(k)
                    if k_str in data:
                        data_storage[fragment][k_str]['env_values'].append(data[k_str]['SVM'][0])
                        data_storage[fragment][k_str]['tax_values'].append(data[k_str]['SVM'][1])

    # Compute averages and variances for each fragment and k-mer size
    results = {fragment: {str(k): {} for k in range(1, 10)} for fragment in fragment_sizes}
    for fragment in fragment_sizes:
        for k in range(1, 10):
            k_str = str(k)
            env_values = data_storage[fragment][k_str]['env_values']
            tax_values = data_storage[fragment][k_str]['tax_values']
            if env_values and tax_values:  # Ensure there are values to compute stats
                results[fragment][k_str]['environment_avg'] = np.mean(env_values) * 100
                results[fragment][k_str]['environment_var'] = np.var(env_values) * 100
                results[fragment][k_str]['taxonomy_avg'] = np.mean(tax_values) * 100
                results[fragment][k_str]['taxonomy_var'] = np.var(tax_values) * 100

    # Output the final average results
    for fragment in fragment_sizes:
        print(f"Fragment Size: {fragment}")
        for k in range(1, 10):
            k_str = str(k)
            print(
                f"  K-mer {k}: Environment Avg: {results[fragment][k_str]['environment_avg']}, Environment Var: {results[fragment][k_str]['environment_var']}, Taxonomy Avg: {results[fragment][k_str]['taxonomy_avg']}, Taxonomy Var: {results[fragment][k_str]['taxonomy_var']}")

    # Optionally, save the results to a file
    output_path = os.path.join(base_dir, f'statistical_results_{env}.json')
    with open(output_path, 'w') as file:
        json.dump(results, file, indent=4)

    print(f'Statistical results saved to {output_path}')

statistical_resutls('pH')
statistical_resutls('temperature')

# Load the previously saved statistical results
results_path = os.path.join(base_dir, 'statistical_results_pH.json')
with open(results_path, 'r') as file:
    results = json.load(file)

# Initialize a dictionary to store the maximum averages and their variances for each fragment size
max_averages = {
    fragment: {
        'max_env_avg': 0, 'max_tax_avg': 0,
        'var_env': 0, 'var_tax': 0,
        'kmer_env': 0, 'kmer_tax': 0
    } for fragment in results
}

# Iterate over each fragment size to find the maximum average accuracies
for fragment in results:
    for k in results[fragment]:
        current_env_avg = results[fragment][k]['environment_avg']
        current_tax_avg = results[fragment][k]['taxonomy_avg']

        # Update max average for environment if the current one is larger
        if current_env_avg > max_averages[fragment]['max_env_avg']:
            max_averages[fragment]['max_env_avg'] = current_env_avg
            max_averages[fragment]['var_env'] = results[fragment][k]['environment_var']
            max_averages[fragment]['kmer_env'] = k

        # Update max average for taxonomy if the current one is larger
        if current_tax_avg > max_averages[fragment]['max_tax_avg']:
            max_averages[fragment]['max_tax_avg'] = current_tax_avg
            max_averages[fragment]['var_tax'] = results[fragment][k]['taxonomy_var']
            max_averages[fragment]['kmer_tax'] = k

# Output the maximum averages and their variances for each fragment size
for fragment in max_averages:
    print(f"Fragment Size: {fragment}")
    print(
        f"  Max Environment Avg: {max_averages[fragment]['max_env_avg']} (Variance: {max_averages[fragment]['var_env']}, K-mer: {max_averages[fragment]['kmer_env']})")
    print(
        f"  Max Taxonomy Avg: {max_averages[fragment]['max_tax_avg']} (Variance: {max_averages[fragment]['var_tax']}, K-mer: {max_averages[fragment]['kmer_tax']})")

# save the max averages to a file
output_max_path = os.path.join(base_dir, 'max_averages_var_results_pH.json')
with open(output_max_path, 'w') as file:
    json.dump(max_averages, file, indent=4)

print(f'Maximum average results with variance saved to {output_max_path}')
