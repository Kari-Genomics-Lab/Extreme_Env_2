#!/usr/bin/env python3
import argparse
import json
import os
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def load_json(file_paths: Dict[str, Path], exp_id: str) -> pd.DataFrame:
    all_data = pd.DataFrame()
    for file_label, file_path in file_paths.items():
        if not file_path.exists():
            print(f"[warn] missing: {file_path}")
            continue

        with open(file_path, "r") as f:
            data = json.load(f)

        # most dumps are nested under exp_id ('0'...'9'); fall back to root if not
        data_block = data.get(exp_id, data)

        results = []
        for k, models in data_block.items():
            try:
                k_int = int(k)
            except ValueError:
                continue

            for model, scores in models.items():
                # Expect scores like [env, taxa] in 0..1
                if (
                    isinstance(scores, list)
                    and len(scores) >= 2
                    and isinstance(scores[0], (int, float))
                    and isinstance(scores[1], (int, float))
                ):
                    results.append({
                        "file": file_label,
                        "k": k_int,
                        "model": model,
                        "env_score": float(f"{scores[0] * 100:.2f}"),
                        "taxa_score": float(f"{scores[1] * 100:.2f}"),
                    })

        if results:
            temp_df = pd.DataFrame(results)
            all_data = pd.concat([all_data, temp_df], ignore_index=True)

    return all_data


def find_best_k_for_each_fragment(df: pd.DataFrame, env: str, scenario: str, out_dir: Path) -> pd.DataFrame:
    best_ks = []
    for file_label, file_group in df.groupby("file"):
        for model, model_group in file_group.groupby("model"):
            max_acc_env = model_group["env_score"].max()
            best_k_env = list(model_group.loc[model_group["env_score"] == max_acc_env, "k"])
            max_acc_taxa = model_group["taxa_score"].max()
            best_k_taxa = list(model_group.loc[model_group["taxa_score"] == max_acc_taxa, "k"])

            best_ks.append({
                "length": file_label,
                "model": model,
                "best_k_env": best_k_env,
                "env_score_at_best_k": float(max_acc_env),
                "best_k_taxa": best_k_taxa,
                "taxa_score_at_best_k": float(max_acc_taxa),
            })

    best_ks_df = pd.DataFrame(best_ks)
    out_csv = out_dir / f"best_k_values_and_scores_{env}_{scenario}.csv"
    best_ks_df.to_csv(out_csv, index=False)
    print(f"[ok] wrote {out_csv}")
    return best_ks_df


def find_best_fragment_for_each_k(df: pd.DataFrame, env: str, scenario: str, out_dir: Path) -> pd.DataFrame:
    best_lengths = []
    for k_label, k_group in df.groupby("k"):
        for model, model_group in k_group.groupby("model"):
            max_acc_env = model_group["env_score"].max()
            best_len_env = list(model_group.loc[model_group["env_score"] == max_acc_env, "file"])
            max_acc_taxa = model_group["taxa_score"].max()
            best_len_taxa = list(model_group.loc[model_group["taxa_score"] == max_acc_taxa, "file"])

            best_lengths.append({
                "k": k_label,
                "model": model,
                "best_length_env": best_len_env,
                "env_score_at_best_length": float(max_acc_env),
                "best_length_taxa": best_len_taxa,
                "taxa_score_at_best_length": float(max_acc_taxa),
            })

    best_lengths_df = pd.DataFrame(best_lengths)
    out_csv = out_dir / f"best_lengths_values_and_scores_{env}_{scenario}.csv"
    best_lengths_df.to_csv(out_csv, index=False)
    print(f"[ok] wrote {out_csv}")
    return best_lengths_df


def create_plots(score_type: str, models: List[str], files: List[str], df: pd.DataFrame, env: str, scenario: str, out_dir: Path):
    # one subplot per model
    fig, axes = plt.subplots(nrows=1, ncols=len(models), figsize=(18, 6))

    if len(models) == 1:
        axes = [axes]

    markers = ["o", "^", "s", "p", "*", "h", "x"]
    colors = plt.cm.viridis(np.linspace(0, 1, len(files)))

    all_scores = df[score_type + "_score"].dropna().astype(float)
    if len(all_scores) == 0:
        print("[warn] no scores to plot")
        return

    min_score, max_score = all_scores.min(), all_scores.max()
    margin = (max_score - min_score) * 0.05 if max_score > min_score else 1.0
    y_min, y_max = min_score - margin, max_score + margin

    for ax, model in zip(axes, models):
        model_group = df[df["model"] == model]
        for file_idx, file_label in enumerate(files):
            file_group = model_group[model_group["file"] == file_label].sort_values("k")
            if file_group.empty:
                continue
            ax.plot(
                file_group["k"],
                file_group[score_type + "_score"].astype(float),
                marker=markers[file_idx % len(markers)],
                linestyle="-",
                label=f"{file_label}",
                color=colors[file_idx],
            )

        ax.set_title(f"{model} Model Classification", fontsize=14)
        ax.set_xlabel("k-mer Length", fontsize=14)
        ax.set_ylabel("Accuracy (%)", fontsize=14)
        ax.set_ylim(y_min, y_max)
        ax.grid(True, which="both", linestyle="--", linewidth=0.5)
        ax.legend(title="Fragment Length", fontsize=10)

    plt.subplots_adjust(left=0.06, right=0.97, top=0.90, bottom=0.15)
    out_pdf = out_dir / f"{score_type}_{env}_{scenario}.pdf"
    plt.savefig(out_pdf, format="pdf", dpi=300)
    plt.tight_layout()
    print(f"[ok] wrote {out_pdf}")
    # don't force plt.show() in headless/batch use
    plt.close(fig)


def build_file_map(data_root: Path, env: str, scenario: str) -> Dict[str, Path]:
    frag_sizes = ["10000", "50000", "100000", "250000", "500000", "1000000"]
    prefix = "" if scenario == "normal" else "Challenging_"
    name = f"{prefix}all_Results_{env}.json"
    return {fs: data_root / f"fragments_{fs}" / name for fs in frag_sizes}


def run(env: str, scenario: str, data_root: Path, output_root: Path, exp_id: str, max_k: int | None):
    # inputs
    file_paths = build_file_map(data_root, env, scenario)

    # outputs
    results_dir = output_root
    graphs_dir = output_root / "graph"
    results_dir.mkdir(parents=True, exist_ok=True)
    graphs_dir.mkdir(parents=True, exist_ok=True)

    # Load & optionally filter by max_k
    df = load_json(file_paths, exp_id=exp_id)
    if df.empty:
        print("[warn] no data loaded — check paths and filenames")
        return

    if max_k is not None:
        df = df[df["k"] <= int(max_k)]

    models = sorted(df["model"].unique().tolist())
    files = sorted(df["file"].unique().tolist())

    # CSVs
    find_best_k_for_each_fragment(df, env, scenario, results_dir)
    find_best_fragment_for_each_k(df, env, scenario, results_dir)

    # Plots (env only, like your script; uncomment taxa if desired)
    create_plots("env", models, files, df, env, scenario, graphs_dir)
    # create_plots("taxa", models, files, df, env, scenario, graphs_dir)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--env", action="store", type=str, required=True, help="Filename suffix, e.g., pH or Temperature")
    parser.add_argument("--scenario", action="store", type=str, choices=["normal", "challenging"], default="normal")
    parser.add_argument("--data_root", default="data", help="input data root (folder containing fragments_*/)")
    parser.add_argument("--output_root", default="outputs", help="output root for CSVs and graphs")
    parser.add_argument("--exp_id", default="0", help="experiment id key inside JSON (default: '0')")
    parser.add_argument("--max_k", action="store", type=int, help="cap k (optional)")
    args = parser.parse_args()

    data_root = Path(args.data_root).expanduser().resolve()
    output_root = Path(args.output_root).expanduser().resolve()

    if not data_root.exists():
        raise SystemExit(f"[error] input directory not found: {data_root}")

    run(
        env=args.env,
        scenario=args.scenario,
        data_root=data_root,
        output_root=output_root,
        exp_id=args.exp_id,
        max_k=args.max_k,
    )


if __name__ == "__main__":
    main()
