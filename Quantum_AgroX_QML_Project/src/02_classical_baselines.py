"""
Quantum_AgroX – 02_classical_baselines.py
Classical ML baselines: Random Forest, XGBoost, SVM for mortality prediction.
"""

import os
import sys
import yaml
import json
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.model_selection import cross_val_score, StratifiedKFold, KFold
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.svm import SVR, SVC
from sklearn.metrics import (
    mean_absolute_error, mean_squared_error, r2_score,
    accuracy_score, balanced_accuracy_score, f1_score, roc_auc_score,
)
from xgboost import XGBRegressor, XGBClassifier


def load_config(config_path: str = "config/project_config.yaml") -> dict:
    base_dir = Path(__file__).resolve().parent.parent
    with open(base_dir / config_path, "r") as f:
        return yaml.safe_load(f)


def load_processed_data(processed_dir: Path) -> dict:
    data = {}
    X_path = processed_dir / "X_mortality.npy"
    y_reg_path = processed_dir / "y_mortality_regression.npy"
    y_cls_path = processed_dir / "y_mortality_classification.npy"
    if X_path.exists() and y_reg_path.exists():
        data["X"] = np.load(X_path)
        data["y_reg"] = np.load(y_reg_path)
    if y_cls_path.exists():
        data["y_cls"] = np.load(y_cls_path)
    return data


def evaluate_regression(y_true, y_pred) -> dict:
    return {
        "MAE": round(mean_absolute_error(y_true, y_pred), 4),
        "RMSE": round(np.sqrt(mean_squared_error(y_true, y_pred)), 4),
        "R2": round(r2_score(y_true, y_pred), 4),
    }


def evaluate_classification(y_true, y_pred, y_prob=None) -> dict:
    metrics = {
        "accuracy": round(accuracy_score(y_true, y_pred), 4),
        "balanced_accuracy": round(balanced_accuracy_score(y_true, y_pred), 4),
        "f1": round(f1_score(y_true, y_pred, zero_division=0), 4),
    }
    if y_prob is not None and len(np.unique(y_true)) == 2:
        try:
            metrics["roc_auc"] = round(roc_auc_score(y_true, y_prob), 4)
        except ValueError:
            metrics["roc_auc"] = 0.0
    return metrics


def run_regression_baselines(X, y, config):
    print("\n--- REGRESSION BASELINES (Mortality Fraction Prediction) ---")
    models = {
        "RandomForest": RandomForestRegressor(
            n_estimators=config["models"]["random_forest"]["n_estimators"],
            max_depth=config["models"]["random_forest"]["max_depth"],
            random_state=config["models"]["random_forest"]["random_state"],
        ),
        "XGBoost": XGBRegressor(
            n_estimators=config["models"]["xgboost"]["n_estimators"],
            max_depth=config["models"]["xgboost"]["max_depth"],
            learning_rate=config["models"]["xgboost"]["learning_rate"],
            random_state=config["models"]["xgboost"]["random_state"],
            verbosity=0,
        ),
        "SVR": SVR(
            kernel=config["models"]["svm"]["kernel"],
            C=config["models"]["svm"]["C"],
            gamma=config["models"]["svm"]["gamma"],
        ),
    }
    results = {}
    n_splits = config["evaluation"]["cv_folds"]
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=config["evaluation"]["random_state"])
    for name, model in models.items():
        print(f"\n  Training {name}...")
        fold_metrics = {"MAE": [], "RMSE": [], "R2": []}
        for fold, (train_idx, test_idx) in enumerate(kf.split(X)):
            X_train, X_test = X[train_idx], X[test_idx]
            y_train, y_test = y[train_idx], y[test_idx]
            model.fit(X_train, y_train)
            y_pred = model.predict(X_test)
            fold_m = evaluate_regression(y_test, y_pred)
            for k in fold_metrics:
                fold_metrics[k].append(fold_m[k])
        avg_metrics = {k: round(np.mean(v), 4) for k, v in fold_metrics.items()}
        std_metrics = {k: round(np.std(v), 4) for k, v in fold_metrics.items()}
        results[name] = {"mean": avg_metrics, "std": std_metrics}
        print(f"    {name}: MAE={avg_metrics['MAE']}±{std_metrics['MAE']}, "
              f"RMSE={avg_metrics['RMSE']}±{std_metrics['RMSE']}, "
              f"R²={avg_metrics['R2']}±{std_metrics['R2']}")
    return results


def run_classification_baselines(X, y, config):
    print("\n--- CLASSIFICATION BASELINES (High vs Lower Efficacy) ---")
    models = {
        "RandomForest": RandomForestClassifier(
            n_estimators=config["models"]["random_forest"]["n_estimators"],
            max_depth=config["models"]["random_forest"]["max_depth"],
            random_state=config["models"]["random_forest"]["random_state"],
            class_weight="balanced",
        ),
        "XGBoost": XGBClassifier(
            n_estimators=config["models"]["xgboost"]["n_estimators"],
            max_depth=config["models"]["xgboost"]["max_depth"],
            learning_rate=config["models"]["xgboost"]["learning_rate"],
            random_state=config["models"]["xgboost"]["random_state"],
            verbosity=0,
            use_label_encoder=False,
            eval_metric="logloss",
        ),
        "SVC": SVC(
            kernel=config["models"]["svm"]["kernel"],
            C=config["models"]["svm"]["C"],
            gamma=config["models"]["svm"]["gamma"],
            probability=True,
            class_weight="balanced",
        ),
    }
    results = {}
    n_splits = config["evaluation"]["cv_folds"]
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=config["evaluation"]["random_state"])
    for name, model in models.items():
        print(f"\n  Training {name}...")
        fold_metrics = {"accuracy": [], "balanced_accuracy": [], "f1": [], "roc_auc": []}
        for fold, (train_idx, test_idx) in enumerate(skf.split(X, y)):
            X_train, X_test = X[train_idx], X[test_idx]
            y_train, y_test = y[train_idx], y[test_idx]
            model.fit(X_train, y_train)
            y_pred = model.predict(X_test)
            y_prob = model.predict_proba(X_test)[:, 1] if hasattr(model, "predict_proba") else None
            fold_m = evaluate_classification(y_test, y_pred, y_prob)
            for k in fold_metrics:
                fold_metrics[k].append(fold_m.get(k, 0.0))
        avg_metrics = {k: round(np.mean(v), 4) for k, v in fold_metrics.items()}
        std_metrics = {k: round(np.std(v), 4) for k, v in fold_metrics.items()}
        results[name] = {"mean": avg_metrics, "std": std_metrics}
        print(f"    {name}: Acc={avg_metrics['accuracy']}±{std_metrics['accuracy']}, "
              f"F1={avg_metrics['f1']}±{std_metrics['f1']}, "
              f"ROC-AUC={avg_metrics['roc_auc']}±{std_metrics['roc_auc']}")
    return results


def main():
    print("=" * 70)
    print("Quantum_AgroX – Classical ML Baselines")
    print("=" * 70)
    config = load_config()
    base_dir = Path(__file__).resolve().parent.parent
    processed_dir = base_dir / config["data"]["processed_dir"]
    data = load_processed_data(processed_dir)
    if not data:
        print("ERROR: Processed data not found. Run 01_prepare_dataset.py first.")
        sys.exit(1)
    X = data["X"]
    y_reg = data["y_reg"]
    y_cls = data["y_cls"]
    print(f"\nDataset: X={X.shape}, y_reg={y_reg.shape}, y_cls={y_cls.shape}")
    print(f"Class distribution: 0={np.sum(y_cls==0)}, 1={np.sum(y_cls==1)}")
    reg_results = run_regression_baselines(X, y_reg, config)
    cls_results = run_classification_baselines(X, y_cls, config)
    all_results = {"regression": reg_results, "classification": cls_results}
    results_dir = base_dir / "results" / "models"
    results_dir.mkdir(parents=True, exist_ok=True)
    with open(results_dir / "classical_baselines.json", "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nResults saved to {results_dir / 'classical_baselines.json'}")
    print("\n" + "=" * 70)
    print("Classical baselines complete.")
    print("=" * 70)


if __name__ == "__main__":
    main()
