from distance_analysis import read_fasta, calculate_dssim, descriptor_distance
import os
import pandas as pd
import lpips
import json
import scripts.src.CGR_utils as CGR_utils
import cv2
import numpy as np


# def fcgr_calculation(sequence):
#     k = 6
#     fcgr = CGR_utils.FCGR(k=k, bits=8)
#     chaos = fcgr(sequence)
#     img = np.array(fcgr.plot(chaos))
#     resized_image = cv2.resize(img, (64,64))
#     # Normalize the pixel values to the range [0, 1]
#     normalized_image = resized_image / (255.0)
#
#     return normalized_image

def fcgr_calculation(sequence):
    k = 6
    fcgr = CGR_utils.FCGR(k=k, bits=8)
    chaos = fcgr(sequence)
    img = np.array(fcgr.plot(chaos))

    return img




fragment_length = 100000
taxa = 'genus'
data_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'data'))
result_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'results'))

# loss_fn = lpips.LPIPS(net='alex')

distance = {}
id_2_sequences = {}
dist_metric = "MLDSP"


def get_value(row_label, column_label, df):
    try:
        value = df.loc[df['Unnamed: 0'] == row_label, column_label].values[0]
        return value
    except:
        return None

for env in ["Temperature", "pH"]:
    if dist_metric == "MLDSP":
        file_path = f"{result_path}/distances/dists100k6.xlsx"
        xls = pd.ExcelFile(file_path)
        print("Available sheets:", xls.sheet_names)
        if env == "Temperature":
            env1 = "temp"
        else:
            env1 = "pH"
        df = pd.read_excel(file_path, sheet_name=env1)
        summary_path = os.path.join(data_path, f'fragments_{fragment_length}/{env}/Extremophiles_{env}_Summary.tsv')
        summary_data = pd.read_csv(summary_path, sep='\t')
        for name1, group1 in summary_data.groupby(taxa):
            print("genus:", name1)
            if len(group1) > 1:
                for row in group1.iterrows():
                    id1 = row[1]['Assembly']
                    for row2 in group1.iterrows():
                        id2 = row2[1]['Assembly']
                        row_label, column_label = "'" + id1 + "'", "'" + id2 + "'"
                        value = get_value(row_label, column_label, df)
                        if name1 not in distance:
                            distance[name1] = {}
                        if id1 != id2:
                            if f"{id1}_{id2}" or f"{id2}_{id1}" not in distance[name1]:
                                if value != None:
                                    distance[name1][f"{id1}_{id2}"] = value

        json.dump(distance, open(f"{result_path}/distances/{taxa}_{env}_level_dist_MLDSP.json", "w"), indent=4)

    else:
        fasta_file = os.path.join(data_path, f'fragments_{fragment_length}', env, f'Extremophiles_{env}.fas')
        id_2_sequences = read_fasta(fasta_file, id_2_sequences)
        summary_path = os.path.join(data_path, f'fragments_{fragment_length}/{env}/Extremophiles_{env}_Summary.tsv')
        summary_data = pd.read_csv(summary_path, sep='\t')
        print(summary_data.groupby(taxa).count())
        for name1, group1 in summary_data.groupby(taxa):
            print("genus:", name1)
            if len(group1) > 1:
                for row in group1.iterrows():
                    id1 = row[1]['Assembly']
                    img1 = fcgr_calculation(id_2_sequences[id1])
                    for row2 in group1.iterrows():
                        id2 = row2[1]['Assembly']
                        img2 = fcgr_calculation(id_2_sequences[id2])
                        # dist = calculate_dssim(img1, img2)
                        dist = descriptor_distance(img1, img2, 2, [0, 8, 16])
                        if name1 not in distance:
                            distance[name1] = {}
                        if id1 != id2:
                            if f"{id1}_{id2}" or f"{id2}_{id1}" not in distance[name1]:
                                if dist != None:
                                    distance[name1][f"{id1}_{id2}"] = dist

        json.dump(distance, open(f"{result_path}/distances/{taxa}_{env}_level_dist_Desc.json", "w"), indent=4)




