#!/usr/bin/env python3
import os
import json
from pathlib import Path
from typing import Dict, List
import numpy as np
import argparse

# Fixed defaults (simple & explicit)
EXPERIMENT_FOLDERS = [str(i) for i in range(10)]
DEFAULT_FRAGMENT_SIZES = ['10000', '50000', '100000', '250000', '500000', '1000000']
CLASSIFIER_KEY = "SVM"


def compute_stats(
    base_in: Path,
    env: str,
    max_k: int,
    fragments: List[str],
    whole_genome: bool = False,
) -> Dict[str, Dict[str, Dict[str, float]]]:
    data_storage: Dict[str, Dict[str, Dict[str, List[float]]]] = {
        fragment: {str(k): {'env_values': [], 'tax_values': []} for k in range(1, max_k + 1)}
        for fragment in fragments
    }

    for exp_folder in EXPERIMENT_FOLDERS:
        exp_folder_path = base_in / exp_folder

        for fragment in fragments:
            fragment_dir = (exp_folder_path / fragment) if whole_genome else (exp_folder_path / f"fragments_{fragment}")
            json_path = fragment_dir / f"Challenging_Supervised_Results_{env}.json"

            if not json_path.exists():
                print(f"[warn] missing: {json_path}")
                continue

            try:
                payload = json.loads(json_path.read_text())
            except json.JSONDecodeError as e:
                print(f"[warn] bad JSON: {json_path} → {e}")
                continue

            data_for_exp = payload.get(exp_folder, payload)

            for k in range(1, max_k + 1):
                k_str = str(k)
                entry = data_for_exp.get(k_str)
                if not entry:
                    continue
                vals = entry.get(CLASSIFIER_KEY)
                if (
                    isinstance(vals, list) and len(vals) >= 2
                    and isinstance(vals[0], (int, float)) and isinstance(vals[1], (int, float))
                ):
                    data_storage[fragment][k_str]['env_values'].append(float(vals[0]))
                    data_storage[fragment][k_str]['tax_values'].append(float(vals[1]))
                else:
                    print(f"[warn] {json_path} k={k_str} missing '{CLASSIFIER_KEY}' values")

    results: Dict[str, Dict[str, Dict[str, float]]] = {
        fragment: {str(k): {} for k in range(1, max_k + 1)} for fragment in fragments
    }

    for fragment in fragments:
        print(f"\nFragment: {fragment}")
        for k in range(1, max_k + 1):
            k_str = str(k)
            env_vals = data_storage[fragment][k_str]['env_values']
            tax_vals = data_storage[fragment][k_str]['tax_values']
            if env_vals and tax_vals:
                env_avg = float(np.mean(env_vals) * 100.0)
                env_var = float(np.var(env_vals) * 100.0)
                tax_avg = float(np.mean(tax_vals) * 100.0)
                tax_var = float(np.var(tax_vals) * 100.0)
                results[fragment][k_str] = {
                    "environment_avg": env_avg,
                    "environment_var": env_var,
                    "taxonomy_avg": tax_avg,
                    "taxonomy_var": tax_var,
                }
                print(
                    f"  k={k}: env avg={env_avg:.4f}, var={env_var:.4f} | "
                    f"tax avg={tax_avg:.4f}, var={tax_var:.4f}"
                )
            else:
                print(f"  k={k}: (no data)")
    return results


def compute_maxima(results: Dict[str, Dict[str, Dict[str, float]]]) -> Dict[str, Dict[str, object]]:
    maxima: Dict[str, Dict[str, object]] = {
        fragment: {'max_env_avg': 0.0, 'var_env': 0.0, 'kmer_env': None,
                   'max_tax_avg': 0.0, 'var_tax': 0.0, 'kmer_tax': None}
        for fragment in results
    }

    for fragment, kd in results.items():
        for k, vals in kd.items():
            if not vals:
                continue
            env_avg = vals.get('environment_avg', 0.0)
            tax_avg = vals.get('taxonomy_avg', 0.0)

            if env_avg > maxima[fragment]['max_env_avg']:
                maxima[fragment]['max_env_avg'] = env_avg
                maxima[fragment]['var_env'] = vals.get('environment_var', 0.0)
                maxima[fragment]['kmer_env'] = k

            if tax_avg > maxima[fragment]['max_tax_avg']:
                maxima[fragment]['max_tax_avg'] = tax_avg
                maxima[fragment]['var_tax'] = vals.get('taxonomy_var', 0.0)
                maxima[fragment]['kmer_tax'] = k

    print("\n=== Max Averages per Fragment ===")
    for fragment, md in maxima.items():
        print(f"\nFragment: {fragment}")
        print(f"  Taxonomy: max={md['max_tax_avg']:.4f}, var={md['var_tax']:.4f}, k={md['kmer_tax']}")
        print(f"  Environment: max={md['max_env_avg']:.4f}, var={md['var_env']:.4f}, k={md['kmer_env']}")
    return maxima


def run_pipeline(args):
    # args keys: env, max_k, whole_genome, data_root, output_root
    env = args['env']
    if env is None:
        raise SystemExit("[error] --env is required (e.g., pH or temperature)")

    max_k = int(args['max_k'] or 9)
    whole_genome = bool(args.get('whole_genome', False))

    base_in = Path(args['data_root']).expanduser().resolve()
    base_out = Path(args['output_root']).expanduser().resolve()
    base_out.mkdir(parents=True, exist_ok=True)

    if not base_in.exists():
        raise SystemExit(f"[error] input directory not found: {base_in}")

    fragments = ['whole_genome'] if whole_genome else DEFAULT_FRAGMENT_SIZES

    results = compute_stats(
        base_in=base_in,
        env=env,
        max_k=max_k,
        fragments=fragments,
        whole_genome=whole_genome,
    )

    stats_out = base_out / f"statistical_results_{env}.json"
    stats_out.write_text(json.dumps(results, indent=4))
    print(f"\n[ok] wrote {stats_out}")

    maxima = compute_maxima(results)
    max_out = base_out / f"max_averages_var_results_{env}.json"
    max_out.write_text(json.dumps(maxima, indent=4))
    print(f"[ok] wrote {max_out}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--env', action='store', type=str, help="Filename suffix, e.g., pH or temperature")
    parser.add_argument('--max_k', action='store', type=int, help="Max k-mer (default: 9)")
    parser.add_argument('--whole_genome', action='store_true', help="Look under <exp>/whole_genome instead of fragments_*")
    parser.add_argument('--data_root', default='data', help='input data root (directory containing experiment folders 0..9)')
    parser.add_argument('--output_root', default='outputs', help='output results root')
    args = vars(parser.parse_args())

    run_pipeline(args)


if __name__ == '__main__':
    main()
