import pandas as pd
import json
import matplotlib.pyplot as plt
from skimage.metrics import structural_similarity as ssim
import cv2
import lpips
import torch
import numpy as np
from PIL import Image
import scripts.src.CGR_utils as CGR_utils
import os
from Bio import SeqIO
import math
from skimage.metrics import normalized_root_mse


fragment_length = 100000
############################################################# candidats ids
# folder_path = '../results/non-parametric/candidates/fcgrs/all_6'
#
# # Get the list of all files in the folder
# file_list = os.listdir(folder_path)
#
# # Filter out directories, if any
# file_list = [f for f in file_list if os.path.isfile(os.path.join(folder_path, f))]
# candidates_ids = sorted([(i.split('*_*')[0].replace('.png','').strip(), i.split('*_*')[1].replace('.png','').strip()) for i in file_list])

############################################################# plot
def plot_dist(dist_list, metric):
    x_values = list(range(1, len(dist_list) + 1))
    # Create bar chart
    plt.bar(x_values, dist_list)
    # Create scatter plot
    # plt.scatter(x_values, gurjit_dist_list)
    plt.title(f'candidates bar chart plot - {metric}')
    plt.xlabel('Candidates')
    plt.ylabel('Distances')
    plt.savefig(f"../results/distances/k6_analysis/candidates_{metric}-dist.png")
    plt.grid(True)
    plt.show()

############################################################# MLDSP
def MLDP_run():
    file_path = "../results/distances/dists100k6.xlsx"
    def get_value(row_label, column_label, df_ph, df_temp):
        try:
            value = df_ph.loc[df_ph['Unnamed: 0'] == row_label, column_label].values[0]
            return value
        except:
            try:
                value = df_temp.loc[df_temp['Unnamed: 0'] == row_label, column_label].values[0]
                return value
            except:
                print("not found")
            return None

    xls = pd.ExcelFile(file_path)
    print("Available sheets:", xls.sheet_names)
      # Change this to the sheet you want to read
    df_ph = pd.read_excel(file_path, sheet_name='pH')
    df_temp = pd.read_excel(file_path, sheet_name='temp')

    mldsp_dist = {}
    for labels in candidates_ids:
        row_label, column_label = "'"+labels[0]+"'", "'"+labels[1]+"'"
        value = get_value(row_label, column_label, df_ph, df_temp)
        mldsp_dist[f"{row_label}_{column_label}"] = value
        # print(f"The value at row {row_label} and column {column_label} is: {value}")

    plot_dist(mldsp_dist.values(), "MLDSP")
    json.dump(mldsp_dist, open("../results/distances/k6_analysis/candidates_MLDSP_all.json", "w"), indent=4)

# MLDP_run()
############################################################# DSSIM



def fcgr_calculation(sequence):
    k = 6
    fcgr = CGR_utils.FCGR(k=k, bits=8)
    chaos = fcgr(sequence)
    img = np.array(fcgr.plot(chaos))
    resized_image = cv2.resize(img, (64,64))
    # Normalize the pixel values to the range [0, 1]
    normalized_image = resized_image / (255.0)

    return normalized_image


def calculate_dssim(image1, image2):
    # Ensure the images are the same size
    if image1.shape != image2.shape:
        raise ValueError("Images must have the same dimensions for DSSIM calculation")
    # Calculate SSIM
    ssim_index, _ = ssim(image1, image2, full=True, data_range=1.0)
    # Calculate DSSIM
    dssim_index = (1 - ssim_index)/2

    return dssim_index

def read_fasta(file_path, id_2_sequences):

    # Open the FASTA file
    with open(file_path, 'r') as file:
        # Iterate over each record
        for record in SeqIO.parse(file, 'fasta'):
            id_2_sequences[str(record.id)] = str(record.seq)

    return id_2_sequences


def creat_images():
    dssim_dist = {}
    id_2_sequences = {}
    for env in ["pH", "Temperature"]:
        result_folder = f"../data/fragments_{fragment_length}"
        fasta_file = os.path.join(result_folder, env, f'Extremophiles_{env}.fas')

        id_2_sequences = read_fasta(fasta_file, id_2_sequences)

    for id in candidates_ids:
        img1 = fcgr_calculation(id_2_sequences[id[0]])
        img2 = fcgr_calculation(id_2_sequences[id[1]])
        dist = calculate_dssim(img1, img2)
        dssim_dist[f"{id[0]}_{id[1]}"] = dist

    return dssim_dist


def DSSIM_run():
    dssim_dist = creat_images()
    plot_dist(dssim_dist.values(), "DSSIM")
    json.dump(dssim_dist, open("../results/distances/k6_analysis/candidates_DSSIM_all.json", "w"), indent=4)

# DSSIM_run()

############################################################# lpips


def load_image(path):
    image = Image.open(path).convert('L')  # Load image and convert to grayscale
    image = np.array(image)
    return image

def to_3channel(image):
    # Duplicate the single channel to create a 3-channel image
    return np.stack([image, image, image], axis=-1)

def normalize_image(image):
    # Convert image to range [-1, 1]
    return image / 127.5 - 1.0

def prepare_image_for_lpips(image):
    # Convert image to 3 channels and normalize
    image = to_3channel(image)
    image = normalize_image(image)
    # Convert to tensor and add batch dimension
    image = torch.tensor(image).permute(2, 0, 1).unsqueeze(0).float()
    return image

def fcgr_calculation(sequence):
    k = 6
    fcgr = CGR_utils.FCGR(k=k, bits=8)
    chaos = fcgr(sequence)
    img = np.array(fcgr.plot(chaos))

    return img

# Load and prepare images

def creat_images():
    loss_fn = lpips.LPIPS(net='alex')
    lpip_dist = {}
    id_2_sequences = {}
    for env in ["pH", "Temperature"]:
        result_folder = f"../data/fragments_{fragment_length}"
        fasta_file = os.path.join(result_folder, env, f'Extremophiles_{env}.fas')

        id_2_sequences = read_fasta(fasta_file, id_2_sequences)

    for id in candidates_ids:
        img1 = fcgr_calculation(id_2_sequences[id[0]])
        img2 = fcgr_calculation(id_2_sequences[id[1]])
        dist = loss_fn(img1, img2)
        lpip_dist[f"{id[0]}_{id[1]}"] = dist.item()


    return lpip_dist

def LPIPS_run():
    lpip_dist = creat_images()
    plot_dist(lpip_dist.values(), "LPIPS")
    json.dump(lpip_dist, open("../results/distances/k6_analysis/candidates_LPIPS_all.json", "w"), indent=4)

# LPIPS_run()

################################################################ Descriptor
def split_image(image, split_size=14):
    h, w = image.shape[0], image.shape[1]
    col_count = int(math.ceil(h / split_size))
    row_count = int(math.ceil(w / split_size))
    tiles_count = col_count * row_count
    tiles = np.zeros((tiles_count, split_size, split_size))
    for y in range(col_count):
        for x in range(row_count):
            ind = x + (y * row_count)
            tiles[ind:(ind + 1), :, :] = image[split_size * y: (y + 1) * split_size,
                                         split_size * x:(x + 1) * split_size]

    return tiles

def get_descriptor(patch, bin_bounds):
    descriptor = np.zeros(len(bin_bounds))

    for index, bin_point in enumerate(bin_bounds):
        if index < len(bin_bounds) - 1:
            low = bin_bounds[index]
            high = bin_bounds[index + 1]
        else:
            low = bin_bounds[index]
            high = np.inf

        descriptor[index] = ((low <= patch) & (patch < high)).sum()
    descriptor = descriptor / np.sum(descriptor)
    descriptor = list(descriptor)
    return descriptor


def descriptor_distance(img1, img2, m=2, bins_bound=None):
    p1 = split_image(img1, 2 ** m)
    p2 = split_image(img2, 2 ** m)

    sub_matrices = p1.shape[0]

    vec1 = []
    vec2 = []
    for i in range(sub_matrices):
        vec1 += get_descriptor(patch=p1[i], bin_bounds=bins_bound)
        vec2 += get_descriptor(patch=p2[i], bin_bounds=bins_bound)

    vec1 = np.asarray(vec1)
    vec2 = np.asarray(vec2)

    denom_1 = np.sqrt(np.mean((vec1 * vec1), dtype=np.float64))
    denom_2 = np.sqrt(np.mean((vec2 * vec2), dtype=np.float64))
    if denom_1 > denom_2:
        distance = normalized_root_mse(vec1, vec2, normalization='euclidean')
    else:
        distance = normalized_root_mse(vec2, vec1, normalization='euclidean')

    # distance = normalized_root_mse(vec1, vec2, normalization='euclidean')
    # distance = np.sqrt(np.sum((vec1 - vec2) ** 2))
    return distance

def creat_images():

    des_dist = {}
    id_2_sequences = {}
    for env in ["pH", "Temperature"]:
        result_folder = f"../data/fragments_{fragment_length}"
        fasta_file = os.path.join(result_folder, env, f'Extremophiles_{env}.fas')

        id_2_sequences = read_fasta(fasta_file, id_2_sequences)

    for id in candidates_ids:
        img1 = fcgr_calculation(id_2_sequences[id[0]])
        img2 = fcgr_calculation(id_2_sequences[id[1]])
        dist = descriptor_distance(img1, img2, 2, [0, 8, 16])
        des_dist[f"{id[0]}_{id[1]}"] = dist


    return des_dist

def Descriptor_run():
    des_dist = creat_images()
    plot_dist(des_dist.values(), "Descriptor")
    json.dump(des_dist, open("../results/distances/k6_analysis/candidates_Descriptor_all.json", "w"), indent=4)

# Descriptor_run()