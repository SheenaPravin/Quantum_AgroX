"""
Quantum_AgroX – 01_prepare_dataset.py
Data ingestion, standardization, and QML-ready feature engineering.
"""

import os
import sys
import yaml
import numpy as np
import pandas as pd
from pathlib import Path


def load_config(config_path: str = "config/project_config.yaml") -> dict:
    base_dir = Path(__file__).resolve().parent.parent
    with open(base_dir / config_path, "r") as f:
        return yaml.safe_load(f)


def load_raw_data(config: dict) -> dict[str, pd.DataFrame]:
    base_dir = Path(__file__).resolve().parent.parent
    raw_dir = base_dir / config["data"]["raw_dir"]
    datasets = {}
    for key, fname in [
        ("mortality", config["data"]["mortality_dataset"]),
        ("docking", config["data"]["docking_candidates"]),
        ("validation", config["data"]["docking_validation"]),
        ("phytochemicals", config["data"]["phytochemicals"]),
        ("lc50", config["data"]["lc50_data"]),
    ]:
        fpath = raw_dir / fname
        if fpath.exists():
            datasets[key] = pd.read_csv(fpath)
            print(f"  Loaded {key}: {datasets[key].shape}")
        else:
            print(f"  WARNING: {fpath} not found")
    return datasets


def prepare_mortality_features(df: pd.DataFrame, config: dict) -> pd.DataFrame:
    feat_cols = config["features"]["numerical_cols"] + config["features"]["encoded_cols"]
    available = [c for c in feat_cols if c in df.columns]
    features = df[available].copy()
    target_col = config["features"]["target_mortality"]
    if target_col in df.columns:
        features[target_col] = df[target_col]
    high_thresh = config["classification_threshold"]["high_efficacy"]
    if target_col in features.columns:
        features["high_efficacy"] = (features[target_col] >= high_thresh).astype(int)
    meta_cols = ["Species", "Assay", "Extract", "Data_type"]
    for c in meta_cols:
        if c in df.columns:
            features[c] = df[c]
    features = features.dropna(subset=available)
    return features


def prepare_docking_features(df: pd.DataFrame) -> pd.DataFrame:
    features = df.copy()
    if "Binding_affinity_kcal_mol" in features.columns:
        features["binding_strength"] = features["Binding_affinity_kcal_mol"].abs()
    if "Affinity_strength" in features.columns:
        features["affinity_numeric"] = pd.to_numeric(
            features["Affinity_strength"], errors="coerce"
        )
    species_dummies = pd.get_dummies(features["Species"], prefix="species", dtype=int)
    extract_dummies = pd.get_dummies(features["Extract"], prefix="extract", dtype=int)
    features = pd.concat([features, species_dummies, extract_dummies], axis=1)
    return features


def normalize_features(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    df_norm = df.copy()
    for col in cols:
        if col in df_norm.columns:
            min_val = df_norm[col].min()
            max_val = df_norm[col].max()
            if max_val > min_val:
                df_norm[col] = (df_norm[col] - min_val) / (max_val - min_val)
            else:
                df_norm[col] = 0.0
    return df_norm


def main():
    print("=" * 70)
    print("Quantum_AgroX – Dataset Preparation")
    print("=" * 70)
    config = load_config()
    print("\n[1] Loading raw datasets...")
    datasets = load_raw_data(config)
    processed_dir = Path(__file__).resolve().parent.parent / config["data"]["processed_dir"]
    processed_dir.mkdir(parents=True, exist_ok=True)
    print("\n[2] Preparing mortality QML dataset...")
    if "mortality" in datasets:
        mortality_df = datasets["mortality"]
        num_cols = config["features"]["numerical_cols"]
        available_num = [c for c in num_cols if c in mortality_df.columns]
        norm_mortality = normalize_features(mortality_df, available_num)
        norm_mortality.to_csv(processed_dir / "mortality_normalized.csv", index=False)
        print(f"  Saved mortality_normalized.csv: {norm_mortality.shape}")
        features = prepare_mortality_features(mortality_df, config)
        features.to_csv(processed_dir / "mortality_features.csv", index=False)
        print(f"  Saved mortality_features.csv: {features.shape}")
        X_cols = config["features"]["numerical_cols"] + config["features"]["encoded_cols"]
        available_X = [c for c in X_cols if c in features.columns]
        X = features[available_X].values
        y_reg = features[config["features"]["target_mortality"]].values
        y_cls = features["high_efficacy"].values
        np.save(processed_dir / "X_mortality.npy", X)
        np.save(processed_dir / "y_mortality_regression.npy", y_reg)
        np.save(processed_dir / "y_mortality_classification.npy", y_cls)
        print(f"  Saved numpy arrays: X={X.shape}, y_reg={y_reg.shape}, y_cls={y_cls.shape}")
    print("\n[3] Preparing docking QML dataset...")
    if "docking" in datasets:
        docking_features = prepare_docking_features(datasets["docking"])
        docking_features.to_csv(processed_dir / "docking_features.csv", index=False)
        print(f"  Saved docking_features.csv: {docking_features.shape}")
    print("\n[4] Preparing docking validation dataset...")
    if "validation" in datasets:
        val_df = datasets["validation"]
        val_df.to_csv(processed_dir / "docking_validation_clean.csv", index=False)
        print(f"  Saved docking_validation_clean.csv: {val_df.shape}")
    print("\n[5] Dataset summary:")
    print(f"  Mortality records: {datasets.get('mortality', pd.DataFrame()).shape[0]}")
    print(f"  Docking candidates: {datasets.get('docking', pd.DataFrame()).shape[0]}")
    print(f"  Docking validation: {datasets.get('validation', pd.DataFrame()).shape[0]}")
    print(f"  Phytochemicals: {datasets.get('phytochemicals', pd.DataFrame()).shape[0]}")
    print(f"  LC50 entries: {datasets.get('lc50', pd.DataFrame()).shape[0]}")
    print("\n" + "=" * 70)
    print("Dataset preparation complete.")
    print("=" * 70)


if __name__ == "__main__":
    main()
