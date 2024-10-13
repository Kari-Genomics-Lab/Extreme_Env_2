import pandas as pd
import os
import matplotlib.pyplot as plt
import random
from scripts.src.utils import clean_sequence

plt.rcParams["font.family"] = "Times New Roman"
TAXA = 'Genus'
GT = "Domain"


def summary_fasta(filename, min_len):
    names, seqs, plasmids = [], [], []
    seq_id = ""
    lines = []

    with open(filename, "r") as file:
        for line in file:
            line = line.strip()

            if line.startswith('#'):
                continue

            if line.startswith('>'):
                if seq_id and lines:
                    process_sequence(lines, seq_id, min_len, names, seqs, plasmids)
                seq_id = line[1:-1]  # change this based on your label
                lines = []
            else:
                lines.append(line)

        # Handle the last sequence in the file
        if seq_id and lines:
            process_sequence(lines, seq_id, min_len, names, seqs, plasmids)

    return names, seqs, plasmids


def process_sequence(lines, seq_id, min_len, names, seqs, plasmids):
    sequence = "".join(lines)
    sequence = clean_sequence(sequence)
    if len(sequence) >= min_len:
        if 'plasmid' in seq_id.lower():
            plasmids.append(seq_id)
        else:
        # excluding plasmid in the fragment building
            names.append(seq_id)
            seqs.append(sequence)


def max_subfragments(L, l, max_num_part):
    for n in range(max_num_part, 1, -1):
        subfragment_length = l / n
        if subfragment_length * n <= L:
            return n
    return 1


def produce_fragment(names, seqs, min_len, is_whole_genome=False, max_num_part=10):
    if is_whole_genome:
        return ''.join(seq + 'N' for seq in seqs)

    random.shuffle(seqs)
    shuffled_whole_genome = ''.join(seq + 'N' for seq in seqs)[:-1]

    if len(shuffled_whole_genome) < min_len:
        print(len(shuffled_whole_genome))
        return shuffled_whole_genome

    clean_genome_length = len(shuffled_whole_genome.replace('N', ''))
    number_of_random_fragments = max_subfragments(clean_genome_length, min_len, max_num_part)

    if number_of_random_fragments < 2:
        start = random.randint(0, len(shuffled_whole_genome) - min_len)
        return shuffled_whole_genome[start:start + min_len]

    mini_fragments_length = (min_len + number_of_random_fragments - 1) // number_of_random_fragments  # ceil division
    potential_starting_points = list(
        range(0, len(shuffled_whole_genome) - mini_fragments_length + 1, mini_fragments_length))

    start_points = random.sample(potential_starting_points, number_of_random_fragments)
    fragments = [shuffled_whole_genome[start:start + mini_fragments_length] for start in start_points]
    fragment = 'N'.join(fragments)

    while len(fragment.replace('N', '')) < min_len:
        additional_needed = min_len - len(fragment.replace('N', ''))
        start = random.randint(0, len(shuffled_whole_genome) - additional_needed)
        fragment += 'N' + shuffled_whole_genome[start:start + additional_needed]

    if len(fragment.replace('N', '')) != min_len:
        print(len(fragment.replace('N', '')))

    return fragment

def run_fragment_builder(data_path, dataset_file, fragment_file, fragment_length, whole_genome, env):

    new_folder_path = os.path.join(data_path, fragment_file)
    os.makedirs(new_folder_path, exist_ok=True)
    env_folder = os.path.join(new_folder_path, env)
    os.makedirs(env_folder, exist_ok=True)

    dataset_path = os.path.join(data_path, dataset_file)
    dataset = pd.read_csv(dataset_path, delimiter='\t')

    # Filter dataset for missing assemblies
    dataset = filter_assemblies(dataset, data_path)


    # Filter dataset for non-null environmental values
    dataset_env = dataset[dataset[env].notnull()]

    process_environmental_data(dataset_env, data_path, env_folder, env, fragment_length, whole_genome)

def filter_assemblies(dataset, path):
    removed = []
    assemblies_path = os.path.join(path, 'Assemblies')
    for assembly in dataset['Assembly']:
        if not os.path.exists(os.path.join(assemblies_path, assembly)):
            removed.append(assembly)

    return dataset[~dataset['Assembly'].isin(removed)]


def process_environmental_data(dataset_env, path, env_folder, env_type, fragment_length, whole_genome):
    aux_dataset = []

    for index, row in dataset_env.iterrows():

        assembly_path = f'{path}/Assemblies/{row["Assembly"]}'
        seq_names, seqs = process_assembly(assembly_path, row, env_folder)
        if seq_names == None or seqs == None:
            print(f'Assembly {row["Assembly"]} not found')
        else:
            # creat the fragments for sequences that found
            assembly = row['Assembly']
            domain = row['Domain']

            fragment = produce_fragment(seq_names, seqs, fragment_length, is_whole_genome=whole_genome)

            acc = assembly
            aux_dataset.append([acc, assembly, domain, row[env_type], domain + '_' + env_type,
                                row[TAXA], row['Species'], row['tax_cluster_id'], len(fragment)])

            fasta = ('>%s\n%s\n' % (acc, fragment))
            with open(f'{env_folder}/Extremophiles_{env_type}.fas', 'a') as f:
                f.write(fasta)

    # Generate summary files
    summary_df = pd.DataFrame(aux_dataset, columns=[
        'sequence_id', 'Assembly', 'Domain', env_type, 'cluster_id', TAXA,
        'species', 'tax_cluster_id', 'len'])

    summary_path = os.path.join(env_folder, f'Extremophiles_{env_type}_Summary.tsv')
    summary_df.to_csv(summary_path, sep='\t')

    # Rename and export data for taxonomic and environmental summaries
    export_summary_data(summary_path, env_folder, env_type)


def process_assembly(assembly_path, row, env_folder):
    fasta_file = os.listdir(assembly_path)[0] if os.listdir(assembly_path) else None
    # use this to filter the sequence based on the length, the default is 0, so all sequences will be included
    _min = 0

    if fasta_file:
        fasta_path = os.path.join(assembly_path, fasta_file)
        seq_names, seqs, plasmid_names = summary_fasta(fasta_path, _min)

        if plasmid_names:
            with open(os.path.join(env_folder, f'plasmids.tsv'), 'a') as f:
                for plasmid in plasmid_names:
                    f.write(f'{row["Assembly"]} \t {plasmid} \n')

        if not seq_names:
            print('empty_fasta', row['Assembly'], 'min:', 0)

        return seq_names, seqs

    else:
        return None, None


def export_summary_data(summary_path, env_folder, env_type):
    df = pd.read_csv(summary_path, sep='\t', usecols=['sequence_id', 'Assembly', GT, TAXA])
    df.rename(columns={GT: 'cluster_id'}, inplace=True)
    df.to_csv(os.path.join(env_folder, f'Extremophiles_{env_type}_GT_Tax.tsv'), sep='\t')

    df = pd.read_csv(summary_path, sep='\t', usecols=['sequence_id', 'Assembly', env_type, TAXA])
    df.rename(columns={env_type: 'cluster_id'}, inplace=True)
    df.to_csv(os.path.join(env_folder, f'Extremophiles_{env_type}_GT_Env.tsv'), sep='\t')
