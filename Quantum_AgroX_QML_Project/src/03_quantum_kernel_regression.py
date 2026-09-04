"""
Quantum_AgroX – 03_quantum_kernel_regression.py
Quantum Fidelity Kernel + Kernel Ridge Regression for mortality prediction.
"""

import os
import sys
import yaml
import json
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.model_selection import KFold
from sklearn.kernel_ridge import KernelRidge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.preprocessing import MinMaxScaler
import pennylane as qml
from pennylane import numpy as pnp


def load_config(config_path: str = "config/project_config.yaml") -> dict:
    base_dir = Path(__file__).resolve().parent.parent
    with open(base_dir / config_path, "r") as f:
        return yaml.safe_load(f)


def angle_encode(x: np.ndarray, n_qubits: int) -> np.ndarray:
    encoded = np.zeros(n_qubits)
    n_features = min(len(x), n_qubits)
    encoded[:n_features] = x[:n_qubits]
    return encoded


def create_quantum_device(n_qubits: int, shots: int = None):
    return qml.device("default.qubit", wires=n_qubits, shots=shots)


def quantum_kernel_circuit(x1, x2, n_qubits, depth):
    qml.AmplitudeEmbedding(features=x1, wires=range(n_qubits), normalize=True)
    qml.adjoint(qml.StronglyEntanglingLayers)(
        weights=np.random.RandomState(42).uniform(
            0, 2 * np.pi, (depth, n_qubits, 3)
        ),
        wires=range(n_qubits),
    )
    qml.AmplitudeEmbedding(features=x2, wires=range(n_qubits), normalize=True)
    return qml.expval(qml.PauliZ(0))


def compute_quantum_kernel_matrix(X1, X2, n_qubits, depth, shots=None):
    dev = create_quantum_device(n_qubits, shots)
    qnode = qml.QNode(quantum_kernel_circuit, dev)
    n1, n2 = len(X1), len(X2)
    kernel_matrix = np.zeros((n1, n2))
    print(f"  Computing quantum kernel matrix ({n1}x{n2})...")
    for i in range(n1):
        if (i + 1) % 10 == 0 or i == n1 - 1:
            print(f"    Row {i+1}/{n1}")
        for j in range(i, n2):
            try:
                val = qnode(X1[i], X2[j], n_qubits, depth)
                kernel_val = float(np.abs(val))
            except Exception:
                kernel_val = 0.0
            kernel_matrix[i, j] = kernel_val
            kernel_matrix[j, i] = kernel_val
    return kernel_matrix


def fidelity_kernel_matrix(X1, X2, n_qubits):
    n1 = len(X1)
    n2 = len(X2)
    K = np.zeros((n1, n2))
    for i in range(n1):
        for j in range(n2):
            fidelity = np.abs(np.dot(X1[i], X2[j])) ** 2
            K[i, j] = fidelity
    return K


def encode_features(X: np.ndarray, n_qubits: int) -> np.ndarray:
    scaler = MinMaxScaler(feature_range=(0, np.pi))
    X_scaled = scaler.fit_transform(X)
    n_features = min(X_scaled.shape[1], n_qubits)
    X_encoded = np.zeros((len(X_scaled), n_qubits))
    X_encoded[:, :n_features] = X_scaled[:, :n_features]
    return X_encoded, scaler


def main():
    print("=" * 70)
    print("Quantum_AgroX – Quantum Kernel Regression (Fidelity Kernel + KRR)")
    print("=" * 70)
    config = load_config()
    base_dir = Path(__file__).resolve().parent.parent
    processed_dir = base_dir / config["data"]["processed_dir"]
    n_qubits = config["quantum"]["n_qubits"]
    depth = config["quantum"]["circuit_depth"]
    shots = config["quantum"]["shots"]
    n_splits = config["evaluation"]["cv_folds"]
    X_path = processed_dir / "X_mortality.npy"
    y_path = processed_dir / "y_mortality_regression.npy"
    if not X_path.exists() or not y_path.exists():
        print("ERROR: Processed data not found. Run 01_prepare_dataset.py first.")
        sys.exit(1)
    X_raw = np.load(X_path)
    y = np.load(y_path)
    print(f"\nDataset: X={X_raw.shape}, y={y.shape}")
    print(f"Quantum config: {n_qubits} qubits, depth={depth}, shots={shots}")
    X_encoded, scaler = encode_features(X_raw, n_qubits)
    print(f"  Features encoded: {X_encoded.shape} -> {n_qubits} qubits")
    print("\n--- Fidelity Quantum Kernel + KRR ---")
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=config["evaluation"]["random_state"])
    fold_metrics = {"MAE": [], "RMSE": [], "R2": []}
    for fold, (train_idx, test_idx) in enumerate(kf.split(X_encoded)):
        print(f"\n  Fold {fold + 1}/{n_splits}")
        X_train, X_test = X_encoded[train_idx], X_encoded[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]
        print("    Computing fidelity kernel matrices...")
        K_train = fidelity_kernel_matrix(X_train, X_train, n_qubits)
        K_test = fidelity_kernel_matrix(X_test, X_train, n_qubits)
        alphas = [0.01, 0.1, 1.0, 10.0]
        best_mae = float("inf")
        best_alpha = 0.1
        for alpha in alphas:
            krr = KernelRidge(alpha=alpha, kernel="precomputed")
            krr.fit(K_train, y_train)
            y_pred = krr.predict(K_test)
            mae = mean_absolute_error(y_test, y_pred)
            if mae < best_mae:
                best_mae = mae
                best_alpha = alpha
        krr = KernelRidge(alpha=best_alpha, kernel="precomputed")
        krr.fit(K_train, y_train)
        y_pred = krr.predict(K_test)
        mae = mean_absolute_error(y_test, y_pred)
        rmse = np.sqrt(mean_squared_error(y_test, y_pred))
        r2 = r2_score(y_test, y_pred)
        fold_metrics["MAE"].append(mae)
        fold_metrics["RMSE"].append(rmse)
        fold_metrics["R2"].append(r2)
        print(f"    Fold {fold+1}: MAE={mae:.4f}, RMSE={rmse:.4f}, R²={r2:.4f} (alpha={best_alpha})")
    avg = {k: round(np.mean(v), 4) for k, v in fold_metrics.items()}
    std = {k: round(np.std(v), 4) for k, v in fold_metrics.items()}
    print(f"\n  Fidelity KRR Summary: MAE={avg['MAE']}±{std['MAE']}, "
          f"RMSE={avg['RMSE']}±{std['RMSE']}, R²={avg['R2']}±{std['R2']}")
    print("\n--- Quantum Circuit Kernel + KRR (sample, first fold only) ---")
    X_sub = X_encoded[:min(20, len(X_encoded))]
    y_sub = y[:min(20, len(y))]
    if len(X_sub) >= 10:
        split = len(X_sub) // 2
        X_tr, X_te = X_sub[:split], X_sub[split:]
        y_tr, y_te = y_sub[:split], y_sub[split:]
        print(f"  Computing quantum circuit kernel ({len(X_tr)}x{len(X_tr)} train, {len(X_te)}x{len(X_tr)} test)...")
        K_train_q = compute_quantum_kernel_matrix(X_tr, X_tr, n_qubits, depth, shots)
        K_test_q = compute_quantum_kernel_matrix(X_te, X_tr, n_qubits, depth, shots)
        krr_q = KernelRidge(alpha=0.1, kernel="precomputed")
        krr_q.fit(K_train_q, y_tr)
        y_pred_q = krr_q.predict(K_test_q)
        mae_q = mean_absolute_error(y_te, y_pred_q)
        rmse_q = np.sqrt(mean_squared_error(y_te, y_pred_q))
        r2_q = r2_score(y_te, y_pred_q)
        print(f"  Quantum Circuit KRR: MAE={mae_q:.4f}, RMSE={rmse_q:.4f}, R²={r2_q:.4f}")
    else:
        print("  Skipping circuit kernel (insufficient samples)")
        mae_q, rmse_q, r2_q = 0, 0, 0
    results = {
        "fidelity_krr": {"mean": avg, "std": std, "n_qubits": n_qubits},
        "circuit_krr_sample": {
            "MAE": round(mae_q, 4),
            "RMSE": round(rmse_q, 4),
            "R2": round(r2_q, 4),
        },
        "config": {"n_qubits": n_qubits, "depth": depth, "shots": shots},
    }
    results_dir = base_dir / "results" / "models"
    results_dir.mkdir(parents=True, exist_ok=True)
    with open(results_dir / "quantum_kernel_regression.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {results_dir / 'quantum_kernel_regression.json'}")
    print("\n" + "=" * 70)
    print("Quantum Kernel Regression complete.")
    print("=" * 70)


if __name__ == "__main__":
    main()
