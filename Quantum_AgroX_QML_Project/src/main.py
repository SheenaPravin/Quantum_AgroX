"""
Quantum_AgroX – main.py
End-to-end pipeline: prepare data, run classical baselines, run QML models.
"""

import os
import sys
import time
import yaml
import json
import argparse
from pathlib import Path
from datetime import datetime


def load_config(config_path: str = "config/project_config.yaml") -> dict:
    base_dir = Path(__file__).resolve().parent.parent
    with open(base_dir / config_path, "r") as f:
        return yaml.safe_load(f)


def run_pipeline(config, steps=None):
    base_dir = Path(__file__).resolve().parent.parent
    src_dir = base_dir / "src"
    all_steps = {
        "prepare": ("01_prepare_dataset", "Dataset Preparation"),
        "classical": ("02_classical_baselines", "Classical ML Baselines"),
        "qkernel_reg": ("03_quantum_kernel_regression", "Quantum Kernel Regression"),
        "qkernel_cls": ("04_quantum_kernel_classifier", "Quantum Kernel SVM Classifier"),
        "vqc": ("05_vqc_classifier", "Variational Quantum Classifier"),
    }
    if steps is None:
        steps = list(all_steps.keys())
    run_times = {}
    for step_key in steps:
        if step_key not in all_steps:
            print(f"  WARNING: Unknown step '{step_key}', skipping")
            continue
        module_name, description = all_steps[step_key]
        print(f"\n{'='*70}")
        print(f"STEP: {description} ({module_name}.py)")
        print(f"{'='*70}")
        start = time.time()
        try:
            module_path = src_dir / f"{module_name}.py"
            import importlib.util
            spec = importlib.util.spec_from_file_location(module_name, module_path)
            module = importlib.util.module_from_spec(spec)
            sys.modules[module_name] = module
            spec.loader.exec_module(module)
            elapsed = time.time() - start
            run_times[step_key] = {"status": "success", "time_seconds": round(elapsed, 2)}
            print(f"\n  {description} completed in {elapsed:.2f}s")
        except Exception as e:
            elapsed = time.time() - start
            run_times[step_key] = {"status": "failed", "error": str(e), "time_seconds": round(elapsed, 2)}
            print(f"\n  ERROR in {description}: {e}")
    return run_times


def print_summary(run_times: dict):
    print("\n" + "=" * 70)
    print("QUANTUM_AGROX PIPELINE SUMMARY")
    print("=" * 70)
    print(f"Run timestamp: {datetime.now().isoformat()}")
    print(f"\n{'Step':<25} {'Status':<12} {'Time (s)':<10}")
    print("-" * 50)
    for step, info in run_times.items():
        status = info["status"]
        t = info["time_seconds"]
        print(f"{step:<25} {status:<12} {t:<10.2f}")
    total_time = sum(info["time_seconds"] for info in run_times.values())
    n_success = sum(1 for info in run_times.values() if info["status"] == "success")
    print("-" * 50)
    print(f"{'TOTAL':<25} {n_success}/{len(run_times)} OK  {total_time:.2f}s")
    print("=" * 70)


def main():
    parser = argparse.ArgumentParser(description="Quantum_AgroX Pipeline Runner")
    parser.add_argument(
        "--steps",
        nargs="+",
        default=None,
        choices=["prepare", "classical", "qkernel_reg", "qkernel_cls", "vqc"],
        help="Pipeline steps to run (default: all)",
    )
    parser.add_argument(
        "--config",
        type=str,
        default="config/project_config.yaml",
        help="Path to config YAML",
    )
    args = parser.parse_args()
    print("=" * 70)
    print("Quantum_AgroX™ – Quantum ML Platform for Agrochemical Discovery")
    print("Pipeline Execution")
    print(f"Start: {datetime.now().isoformat()}")
    print("=" * 70)
    config = load_config(args.config)
    run_times = run_pipeline(config, steps=args.steps)
    print_summary(run_times)
    base_dir = Path(__file__).resolve().parent.parent
    results_dir = base_dir / "results" / "models"
    results_dir.mkdir(parents=True, exist_ok=True)
    with open(results_dir / "pipeline_run.json", "w") as f:
        json.dump({
            "timestamp": datetime.now().isoformat(),
            "steps": run_times,
        }, f, indent=2)
    print(f"\nPipeline log saved to {results_dir / 'pipeline_run.json'}")


if __name__ == "__main__":
    main()
