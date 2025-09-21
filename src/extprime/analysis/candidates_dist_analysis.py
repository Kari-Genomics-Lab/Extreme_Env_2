#!/usr/bin/env python3
import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional

import matplotlib.pyplot as plt

# Hard-coded (per your request)
METRICS = ["descriptor", "dssim", "lpips"]
THRESHOLDS = [0.190211, 0.501385, 0.17766]


def plot_dist(dist_list: List[float], dist_metric: str):
    x_values = list(range(1, len(dist_list) + 1))
    plt.bar(x_values, dist_list)
    plt.title(f"candidates bar chart plot / temp / {dist_metric}")
    plt.xlabel("Index")
    plt.ylabel("distances")
    plt.grid(True)
    plt.show()


def load_candidates_for_metric(results_root: Path, metric: str, env: Optional[str]) -> Dict[str, float]:
    """
    Load a single metric's JSON file.
    Priority: candidates_{metric}_{env}.json if --env is provided, else candidates_{metric}_all.json
    If the env-specific file is missing, fall back to *_all.json.
    """
    cand_dir = results_root / "distances"
    paths_to_try = []
    if env:
        paths_to_try.append(cand_dir / f"candidates_{metric}_{env}.json")
    paths_to_try.append(cand_dir / f"candidates_{metric}_all.json")

    for p in paths_to_try:
        if p.exists():
            with open(p, "r") as f:
                return json.load(f)
    print(f"[warn] no file found for metric '{metric}' (tried: {', '.join(str(p) for p in paths_to_try)})")
    return {}


def majority_candidates(
    results_root: Path,
    output_root: Path,
    majority_num: int,
    env: Optional[str],
) -> Dict[str, Dict[str, float]]:
    # load all metrics using hard-coded thresholds
    dist_data: Dict[str, Dict[str, float]] = {}
    limited_dist: Dict[str, set] = {}
    for metric, th in zip(METRICS, THRESHOLDS):
        data = load_candidates_for_metric(results_root, metric, env)
        dist_data[metric] = data
        kept = {cid for cid, val in data.items() if val <= th}
        limited_dist[metric] = kept
        print(f"[info] {metric}: kept {len(kept)} / {len(data)} with threshold ≤ {th}")

    # majority vote across metrics
    majority: Dict[str, Dict[str, float]] = {}
    all_ids = set().union(*[set(d.keys()) for d in dist_data.values()]) if dist_data else set()

    for cid in all_ids:
        count = sum(1 for m in METRICS if cid in limited_dist.get(m, set()))
        if count >= majority_num:
            majority[cid] = {m: dist_data[m][cid] for m in METRICS if cid in dist_data[m]}

    # write output
    out_dir = output_root / "candidates"
    out_dir.mkdir(parents=True, exist_ok=True)
    env_tag = env if env else "all"
    out_path = out_dir / f"majority_candidates_{env_tag}.json"
    with open(out_path, "w") as f:
        json.dump(majority, f, indent=4)
    print(f"[ok] wrote {out_path} (n={len(majority)})")

    return majority


def run_pipeline(args):
    env = args.get("env")
    results_root = Path(args["data_root"]).expanduser().resolve()   # expects distances/ under here
    output_root = Path(args["output_root"]).expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    majority_candidates(
        results_root=results_root,
        output_root=output_root,
        majority_num=int(args["majority"]),
        env=env,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--env", action="store", type=str, help="Optional suffix used in filenames (e.g., pH or Temperature)")
    parser.add_argument("--data_root", default="results", help="root containing distances/ (defaults to 'results')")
    parser.add_argument("--output_root", default="results", help="where to write candidates/ majority JSON (defaults to 'results')")
    parser.add_argument("--majority", type=int, default=3, help="minimum number of metrics to pass (default: 3)")
    # kept for compatibility with your prior script signature; not used here
    parser.add_argument("--fragment_length", type=int, default=100000, help="placeholder (unused)")
    args = vars(parser.parse_args())

    run_pipeline(args)


if __name__ == "__main__":
    main()
