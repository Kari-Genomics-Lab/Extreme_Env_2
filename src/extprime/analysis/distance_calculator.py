import json
from skimage.metrics import structural_similarity as ssim
from skimage.metrics import normalized_root_mse
import cv2
import lpips
import torch
import numpy as np
from PIL import Image
import extprime.utils.CGR_utils as CGR_utils
import os
from Bio import SeqIO
import math


class DistanceCalculator:
    def __init__(self, fragment_length=100000):
        self.fragment_length = fragment_length
        self.lpips_model = lpips.LPIPS(net='alex')

    def read_fasta(self, file_path):
        """Read FASTA file and return dictionary of sequences"""
        id_2_sequences = {}
        with open(file_path, 'r') as file:
            for record in SeqIO.parse(file, 'fasta'):
                id_2_sequences[str(record.id)] = str(record.seq)
        return id_2_sequences

    def fcgr_calculation(self, sequence, k=8):
        """Calculate FCGR for a sequence"""
        fcgr = CGR_utils.FCGR(k=k, bits=8)
        chaos = fcgr(sequence)
        img = np.array(fcgr.plot(chaos))
        resized_image = cv2.resize(img, (64, 64))
        normalized_image = resized_image / 255.0
        return normalized_image

    def calculate_dssim(self, image1, image2):
        """Calculate DSSIM between two images"""
        if image1.shape != image2.shape:
            raise ValueError("Images must have the same dimensions for DSSIM calculation")
        ssim_index, _ = ssim(image1, image2, full=True, data_range=1.0)
        dssim_index = (1 - ssim_index) / 2
        return dssim_index

    def prepare_image_for_lpips(self, image):
        """Prepare image for LPIPS calculation"""
        # Convert to 3 channels
        image_3ch = np.stack([image, image, image], axis=-1)
        # Normalize to [-1, 1]
        image_norm = image_3ch * 2.0 - 1.0
        # Convert to tensor
        image_tensor = torch.tensor(image_norm).permute(2, 0, 1).unsqueeze(0).float()
        return image_tensor

    def calculate_lpips(self, image1, image2):
        """Calculate LPIPS distance between two images"""
        img1_tensor = self.prepare_image_for_lpips(image1)
        img2_tensor = self.prepare_image_for_lpips(image2)
        distance = self.lpips_model(img1_tensor, img2_tensor)
        return distance.item()

    def split_image(self, image, split_size=4):
        """Split image into patches"""
        h, w = image.shape[0], image.shape[1]
        col_count = int(math.ceil(h / split_size))
        row_count = int(math.ceil(w / split_size))
        tiles_count = col_count * row_count
        tiles = np.zeros((tiles_count, split_size, split_size))

        for y in range(col_count):
            for x in range(row_count):
                ind = x + (y * row_count)
                y_start = split_size * y
                y_end = min((y + 1) * split_size, h)
                x_start = split_size * x
                x_end = min((x + 1) * split_size, w)

                patch = image[y_start:y_end, x_start:x_end]
                if patch.shape[0] < split_size or patch.shape[1] < split_size:
                    padded_patch = np.zeros((split_size, split_size))
                    padded_patch[:patch.shape[0], :patch.shape[1]] = patch
                    patch = padded_patch

                tiles[ind, :, :] = patch
        return tiles

    def get_descriptor(self, patch, bin_bounds=[0, 0.33, 0.66, 1.0]):
        """Get descriptor for a patch"""
        descriptor = np.zeros(len(bin_bounds) - 1)

        for i in range(len(bin_bounds) - 1):
            low = bin_bounds[i]
            high = bin_bounds[i + 1]
            descriptor[i] = ((low <= patch) & (patch < high)).sum()

        total = np.sum(descriptor)
        if total > 0:
            descriptor = descriptor / total
        return descriptor.tolist()

    def calculate_descriptor_distance(self, img1, img2, m=2):
        """Calculate descriptor distance between two images"""
        p1 = self.split_image(img1, 2 ** m)
        p2 = self.split_image(img2, 2 ** m)

        vec1 = []
        vec2 = []
        for i in range(p1.shape[0]):
            vec1 += self.get_descriptor(p1[i])
            vec2 += self.get_descriptor(p2[i])

        vec1 = np.asarray(vec1)
        vec2 = np.asarray(vec2)

        distance = normalized_root_mse(vec1, vec2, normalization='euclidean')
        return distance

    def load_data(self, env):
        """Load data for a specific environment"""
        result_folder = f"data/fragments_{self.fragment_length}"
        fasta_file = os.path.join(result_folder, env, f'Extremophiles_{env}.fas')

        if not os.path.exists(fasta_file):
            print(f"Warning: {fasta_file} not found")
            return {}

        return self.read_fasta(fasta_file)

    def normalize_distances(self, results):
        """Normalize all distances to range [0, 1] using min-max normalization"""
        print("\nNormalizing distances to [0, 1] range...")

        normalized_results = {}

        for k_val, k_results in results.items():
            normalized_results[k_val] = {}

            for method, distances in k_results.items():
                if not distances:  # Skip if no distances calculated
                    normalized_results[k_val][method] = {}
                    continue

                # Get all distance values
                distance_values = list(distances.values())

                # Find min and max values
                min_dist = min(distance_values)
                max_dist = max(distance_values)

                print(f"  {k_val} - {method}: min={min_dist:.6f}, max={max_dist:.6f}")

                # Normalize each distance value
                normalized_distances = {}

                if max_dist == min_dist:
                    # If all values are the same, set them all to 0
                    print(f"    Warning: All {method} distances are identical, setting to 0")
                    for pair_key in distances.keys():
                        normalized_distances[pair_key] = 0.0
                else:
                    # Apply min-max normalization: (x - min) / (max - min)
                    for pair_key, dist_value in distances.items():
                        normalized_value = (dist_value - min_dist) / (max_dist - min_dist)
                        normalized_distances[pair_key] = normalized_value

                normalized_results[k_val][method] = normalized_distances

                # Verify normalization
                if normalized_distances:
                    norm_values = list(normalized_distances.values())
                    norm_min = min(norm_values)
                    norm_max = max(norm_values)
                    print(f"    After normalization: min={norm_min:.6f}, max={norm_max:.6f}")

        return normalized_results

    def calculate_all_distances(self, k_values=[6, 8], environments=["pH", "Temperature"]):
        """Calculate all pairwise distances for all methods"""

        results = {}

        for k in k_values:
            print(f"\nProcessing k={k}")
            results[f"k_{k}"] = {
                "dssim": {},
                "lpips": {},
                "descriptor": {}
            }

            # Load all sequences from both environments
            all_sequences = {}
            for env in environments:
                print(f"Loading {env} data...")
                env_sequences = self.load_data(env)
                all_sequences.update(env_sequences)

            if not all_sequences:
                print("No sequences found!")
                continue

            sequence_ids = list(all_sequences.keys())
            total_pairs = len(sequence_ids) * (len(sequence_ids) - 1) // 2
            print(f"Total sequences: {len(sequence_ids)}")
            print(f"Total pairs to process: {total_pairs}")

            # Cache for images to avoid recalculation
            image_cache = {}

            processed_pairs = 0
            for i, id1 in enumerate(sequence_ids):
                for j, id2 in enumerate(sequence_ids[i + 1:], i + 1):

                    # Generate images if not cached
                    if id1 not in image_cache:
                        image_cache[id1] = self.fcgr_calculation(all_sequences[id1], k=k)
                    if id2 not in image_cache:
                        image_cache[id2] = self.fcgr_calculation(all_sequences[id2], k=k)

                    img1 = image_cache[id1]
                    img2 = image_cache[id2]

                    pair_key = f"{id1}_{id2}"

                    # Calculate DSSIM
                    dssim_dist = self.calculate_dssim(img1, img2)
                    results[f"k_{k}"]["dssim"][pair_key] = dssim_dist

                    # Calculate LPIPS
                    lpips_dist = self.calculate_lpips(img1, img2)
                    results[f"k_{k}"]["lpips"][pair_key] = lpips_dist

                    # Calculate Descriptor distance
                    desc_dist = self.calculate_descriptor_distance(img1, img2)
                    results[f"k_{k}"]["descriptor"][pair_key] = desc_dist

                    processed_pairs += 1
                    if processed_pairs % 100 == 0:
                        print(f"Processed {processed_pairs}/{total_pairs} pairs")

                print(f"Finished sequence {i + 1}/{len(sequence_ids)}: {id1}")

        return results

    def save_results(self, results, output_dir="results", suffix=""):
        """Save results to JSON files"""
        os.makedirs(output_dir, exist_ok=True)

        for k_val, k_results in results.items():
            for method, distances in k_results.items():
                filename = f"{k_val}_{method}_distances{suffix}.json"
                filepath = os.path.join(output_dir, filename)

                with open(filepath, 'w') as f:
                    json.dump(distances, f, indent=4)

                print(f"Saved {len(distances)} distances to {filepath}")

    def run_pipeline(self, k_values=[6, 8], environments=["pH", "Temperature"], output_dir="results"):
        """Run the complete pipeline"""
        print("Starting distance calculation pipeline...")
        print(f"K values: {k_values}")
        print(f"Environments: {environments}")

        # Calculate all distances
        results = self.calculate_all_distances(k_values, environments)

        # Save original (unnormalized) results
        print("\nSaving original distances...")
        self.save_results(results, output_dir, suffix="_original")

        # Normalize distances to [0, 1]
        normalized_results = self.normalize_distances(results)

        # Save normalized results
        print("\nSaving normalized distances...")
        self.save_results(normalized_results, output_dir, suffix="_normalized")

        print("\nPipeline completed!")
        print("Results saved:")
        print("  - Original distances: *_original.json")
        print("  - Normalized distances: *_normalized.json")

        return results, normalized_results


# Run the pipeline
if __name__ == "__main__":
    calculator = DistanceCalculator(fragment_length=100000)
    original_results, normalized_results = calculator.run_pipeline(
        k_values=[6, 7, 8, 9],
        environments=["pH", "Temperature"],
        output_dir="outputs"
    )