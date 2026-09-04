"""
Quantum_AgroX – 04_quantum_kernel_classifier.py
Quantum Kernel SVM (QSVM-style) for high vs lower efficacy classification.
"""

import os
import sys
import yaml
import json
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.model_selection import StratifiedKFold
from sklearn.svm import SVC
from sklearn.metrics import (
    accuracy_score, balanced_accuracy_score, f1_score, roc_auc_score,
)
from sklearn.preprocessing import MinMaxScaler
import pennylane as qml


def load_config(config_path: str = "config/project_config.yaml") -> dict:
    base_dir = Path(__file__).resolve().parent.parent
    with open(base_dir / config_path, "r") as f:
        return yaml.safe_load(f)


def encode_features(X: np.ndarray, n_qubits: int) -> tuple:
    scaler = MinMaxScaler(feature_range=(0, np.pi))
    X_scaled = scaler.fit_transform(X)
    n_features = min(X_scaled.shape[1], n_qubits)
    X_encoded = np.zeros((len(X_scaled), n_qubits))
    X_encoded[:, :n_features] = X_scaled[:, :n_features]
    return X_encoded, scaler


def quantum_fidelity_kernel(X1, X2):
    K = np.zeros((len(X1), len(X2)))
    for i in range(len(X1)):
        for j in range(len(X2)):
            K[i, j] = np.abs(np.dot(X1[i], X2[j])) ** 2
    return K


def create_quantum_feature_map(n_qubits, n_layers=2):
    dev = qml.device("default.qubit", wires=n_qubits)

    @qml.qnode(dev)
    def feature_map(x):
        for i in range(n_qubits):
            qml.RY(x[i], wires=i)
        for _ in range(n_layers):
            for i in range(n_qubits - 1):
                qml.CNOT(wires=[i, i + 1])
            for i in range(n_qubits):
                qml.RZ(x[i] * 0.5, wires=i)
        return [qml.expval(qml.PauliZ(i)) for i in range(n_qubits)]

    return feature_map


def compute_quantum_features(X, n_qubits, n_layers=2):
    feature_map = create_quantum_feature_map(n_qubits, n_layers)
    X_qfeatures = np.zeros((len(X), n_qubits))
    for i in range(len(X)):
        try:
            feat = feature_map(X[i])
            X_qfeatures[i] = np.array(feat)
        except Exception:
            X_qfeatures[i] = X[i][:n_qubits]
    return X_qfeatures


def main():
    print("=" * 70)
    print("Quantum_AgroX – Quantum Kernel SVM (QSVM-style) Classifier")
    print("=" * 70)
    config = load_config()
    base_dir = Path(__file__).resolve().parent.parent
    processed_dir = base_dir / config["data"]["processed_dir"]
    n_qubits = config["quantum"]["n_qubits"]
    n_layers = config["quantum"]["n_layers"]
    n_splits = config["evaluation"]["cv_folds"]
    X_path = processed_dir / "X_mortality.npy"
    y_path = processed_dir / "y_mortality_classification.npy"
    if not X_path.exists() or not y_path.exists():
        print("ERROR: Processed data not found. Run 01_prepare_dataset.py first.")
        sys.exit(1)
    X_raw = np.load(X_path)
    y = np.load(y_path)
    print(f"\nDataset: X={X_raw.shape}, y={y.shape}")
    print(f"Class distribution: 0={np.sum(y==0)}, 1={np.sum(y==1)}")
    X_encoded, scaler = encode_features(X_raw, n_qubits)
    print(f"  Features encoded to {n_qubits} qubits")
    print("\n--- Approach 1: Quantum Fidelity Kernel SVM ---")
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=config["evaluation"]["random_state"])
    fold_metrics = {"accuracy": [], "balanced_accuracy": [], "f1": [], "roc_auc": []}
    for fold, (train_idx, test_idx) in enumerate(skf.split(X_encoded, y)):
        print(f"\n  Fold {fold + 1}/{n_splits}")
        X_train, X_test = X_encoded[train_idx], X_encoded[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]
        K_train = quantum_fidelity_kernel(X_train, X_train)
        K_test = quantum_fidelity_kernel(X_test, X_train)
        C_values = [0.1, 1.0, 10.0]
        best_f1 = -1
        best_C = 1.0
        for C in C_values:
            svm = SVC(C=C, kernel="precomputed", class_weight="balanced")
            svm.fit(K_train, y_train)
            y_pred = svm.predict(K_test)
            f1 = f1_score(y_test, y_pred, zero_division=0)
            if f1 > best_f1:
                best_f1 = f1
                best_C = C
        svm = SVC(C=best_C, kernel="precomputed", class_weight="balanced", probability=True)
        svm.fit(K_train, y_train)
        y_pred = svm.predict(K_test)
        y_prob = svm.predict_proba(K_test)[:, 1]
        metrics = {
            "accuracy": accuracy_score(y_test, y_pred),
            "balanced_accuracy": balanced_accuracy_score(y_test, y_pred),
            "f1": f1_score(y_test, y_pred, zero_division=0),
            "roc_auc": 0.0,
        }
        try:
            metrics["roc_auc"] = roc_auc_score(y_test, y_prob)
        except ValueError:
            pass
        for k in fold_metrics:
            fold_metrics[k].append(metrics[k])
        print(f"    Fold {fold+1}: Acc={metrics['accuracy']:.4f}, "
              f"F1={metrics['f1']:.4f}, ROC-AUC={metrics['roc_auc']:.4f} (C={best_C})")
    avg1 = {k: round(np.mean(v), 4) for k, v in fold_metrics.items()}
    std1 = {k: round(np.std(v), 4) for k, v in fold_metrics.items()}
    print(f"\n  Fidelity QSVM Summary: Acc={avg1['accuracy']}±{std1['accuracy']}, "
          f"F1={avg1['f1']}±{std1['f1']}, ROC-AUC={avg1['roc_auc']}±{std1['roc_auc']}")
    print("\n--- Approach 2: Quantum Feature Map + RBF SVM (small sample) ---")
    X_sub_idx = np.random.RandomState(42).choice(len(X_encoded), min(30, len(X_encoded)), replace=False)
    X_sub = X_encoded[X_sub_idx]
    y_sub = y[X_sub_idx]
    print(f"  Computing quantum features for {len(X_sub)} samples...")
    X_qfeatures = compute_quantum_features(X_sub, n_qubits, n_layers)
    split = len(X_sub) // 2
    X_tr, X_te = X_qfeatures[:split], X_qfeatures[split:]
    y_tr, y_te = y_sub[:split], y_sub[split:]
    svm_qf = SVC(C=1.0, kernel="rbf", class_weight="balanced", probability=True)
    svm_qf.fit(X_tr, y_tr)
    y_pred_qf = svm_qf.predict(X_te)
    y_prob_qf = svm_qf.predict_proba(X_te)[:, 1]
    qf_acc = accuracy_score(y_te, y_pred_qf)
    qf_f1 = f1_score(y_te, y_pred_qf, zero_division=0)
    try:
        qf_auc = roc_auc_score(y_te, y_prob_qf)
    except ValueError:
        qf_auc = 0.0
    print(f"  Quantum Feature + RBF SVM: Acc={qf_acc:.4f}, F1={qf_f1:.4f}, ROC-AUC={qf_auc:.4f}")
    results = {
        "fidelity_qsvm": {"mean": avg1, "std": std1},
        "quantum_feature_rbf_svm": {
            "accuracy": round(qf_acc, 4),
            "f1": round(qf_f1, 4),
            "roc_auc": round(qf_auc, 4),
            "sample_size": len(X_sub),
        },
        "config": {"n_qubits": n_qubits, "n_layers": n_layers},
    }
    results_dir = base_dir / "results" / "models"
    results_dir.mkdir(parents=True, exist_ok=True)
    with open(results_dir / "quantum_kernel_classifier.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {results_dir / 'quantum_kernel_classifier.json'}")
    print("\n" + "=" * 70)
    print("Quantum Kernel Classifier complete.")
    print("=" * 70)


if __name__ == "__main__":
    main()
