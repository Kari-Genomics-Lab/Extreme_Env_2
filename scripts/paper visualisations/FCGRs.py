import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import os
import scripts.src.CGR_utils as CGR_utils
import re
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt

ENVS = ["Temperature"]
k = 9
def draw_fcgrs(sequence):

    fcgr = CGR_utils.FCGR(k=k, bits=8)
    chaos1 = fcgr(sequence)/len(sequence)
    return fcgr.plot(chaos1)


def read_fasta(file_path):
    id_2_sequences = {}
    # Open the FASTA file
    with open(file_path, 'r') as file:
        # Iterate over each record
        for record in SeqIO.parse(file, 'fasta'):
            id_2_sequences[str(record.id)] = str(record.seq)

    return id_2_sequences


def clean_sequence(sequence):
    # Replace any character not A, C, G, T, or N with N
    return re.sub(r'[^ACGTN]', 'N', sequence)


def run():
    to_displays = { "Temperature":{"p1": ("GCA_000026225.1", "GCA_000026225.1")},
                    #"pH":{"p2": ("GCA_000145615.1", "GCA_000317795.1")}
                    }
    final_display = {}
    for env in ENVS:
    #     for fragment_length in FRAGMENT_LENGTHS:
        fragment_length = 100000
        data_folder = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'data'))
        result_folder = f"{data_folder}/fragments_{fragment_length}"
        fasta_file = os.path.join(result_folder, env, f'Extremophiles_{env}.fas')

        summary_path = f"{data_folder}/fragments_{fragment_length}/{env}/Extremophiles_{env}_Summary.tsv"
        summary_data = pd.read_csv(summary_path, sep='\t')

        id_2_sequences = read_fasta(fasta_file)

        print(len(id_2_sequences))
       # for id in id_2_sequences:
        display = to_displays[env]
        for key in display.keys():

                id1, id2 = display[key][0], display[key][1]
                print(id1, id2)
                # print(list(summary_data[summary_data['Assembly'] == id]))
                domain1 = str(list(summary_data[summary_data['Assembly'] == id1]['Domain'])[0])
                domain2 = str(list(summary_data[summary_data['Assembly'] == id2]['Domain'])[0])
                species1 = str(list(summary_data[summary_data['Assembly'] == id1]['species'])[0])
                species2 = str(list(summary_data[summary_data['Assembly'] == id2]['species'])[0])
                env_label1 = str(list(summary_data[summary_data['Assembly'] == id1][env])[0])
                env_label2 = str(list(summary_data[summary_data['Assembly'] == id2][env])[0])

                # print(f"Domain: {domain}, Env: {env_label}")
                sequence_id1 = clean_sequence(id_2_sequences[id1])
                sequence_id2 = clean_sequence(id_2_sequences[id2])
                img1 = draw_fcgrs(sequence_id1)
                img2 = draw_fcgrs(sequence_id2)
                final_display[key] = ((img1, img2), (id1, id2), (domain1, domain2), (env_label1, env_label2), (species1, species2))

    return final_display

to_displays = run()


# Create a figure
plt.figure(figsize=(8, 4))

# Define font properties
font_properties = {'fontname': 'Helvetica', 'fontsize': 18}

# Plot each image in a separate subplot
for i, key in enumerate(to_displays.keys()):
    for n in range(2):
        ax = plt.subplot(1, 2, 2 * i + n + 1)
        ax.imshow(to_displays[key][0][n], cmap="gray")

        # Title each subplot with letters and additional labels
        title_index = 2 * i + n
        label_prefix = 'a' if i == 0 else 'b'
        # plt.title(f"{label_prefix}-{n+1}) {to_displays[key][3][n]} - {to_displays[key][2][n]} - \n {to_displays[key][4][n]}", y=-0.15, **font_properties)
        print(f"{label_prefix}-{n+1}) {to_displays[key][3][n]} - {to_displays[key][2][n]} - \n {to_displays[key][4][n]}")
        plt.xticks([])  # Remove x ticks
        plt.yticks([])  # Remove y ticks

        # Coordinates of points to label (example points)
        points = [(-6,526),(-12,-10),(522,-10),(524,524)] # List of (x, y) tuples
        labels = ["A", "C", "G", "T"]  # Labels for each point

        # Annotating each point
        for (x, y), label in zip(points, labels):
            ax.text(x, y, label, color='black', ha='center', va='center', fontname= 'Helvetica', fontsize=14)

plt.tight_layout(rect=[0, 0.01, 1, 0.99])
result_folder = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'results'))
plt.savefig(f"{result_folder}/graphs/pairs_fcgr_poster_GCA_000026225.1.png", format='png', dpi=300, transparent=True)
plt.show()




# GCA_000026225.1
# GCA_000007185.1