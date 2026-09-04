"""
Quantum_AgroX – 05_vqc_classifier.py
Variational Quantum Classifier (VQC) for acaricidal efficacy classification.
"""

import os
import sys
import yaml
import json
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score, roc_auc_score
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


def create_vqc(n_qubits, n_layers, n_classes=2):
    dev = qml.device("default.qubit", wires=n_qubits)

    weight_shapes = {"weights": (n_layers, n_qubits, 3)}

    @qml.qnode(dev)
    def vqc_circuit(inputs, weights):
        qml.AngleEmbedding(inputs, wires=range(n_qubits), rotation="Y")
        for layer in range(n_layers):
            for i in range(n_qubits):
                qml.RX(weights[layer, i, 0], wires=i)
                qml.RY(weights[layer, i, 1], wires=i)
                qml.RZ(weights[layer, i, 2], wires=i)
            for i in range(n_qubits - 1):
                qml.CNOT(wires=[i, i + 1])
            qml.CNOT(wires=[n_qubits - 1, 0])
        return qml.expval(qml.PauliZ(0))

    return vqc_circuit, weight_shapes


def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-np.clip(x, -500, 500)))


def train_vqc(X_train, y_train, n_qubits, n_layers, lr=0.1, epochs=30):
    dev = qml.device("default.qubit", wires=n_qubits)

    @qml.qnode(dev, diff_method="parameter-shift")
    def circuit(inputs, weights):
        qml.AngleEmbedding(inputs, wires=range(n_qubits), rotation="Y")
        for layer in range(n_layers):
            for i in range(n_qubits):
                qml.RX(weights[layer, i, 0], wires=i)
                qml.RY(weights[layer, i, 1], wires=i)
                qml.RZ(weights[layer, i, 2], wires=i)
            for i in range(n_qubits - 1):
                qml.CNOT(wires=[i, i + 1])
            qml.CNOT(wires=[n_qubits - 1, 0])
        return qml.expval(qml.PauliZ(0))

    weights = np.random.RandomState(42).uniform(-0.1, 0.1, (n_layers, n_qubits, 3), requires_grad=True)
    opt = qml.GradientDescentOptimizer(stepsize=lr)
    y_binary = 2 * y_train - 1.0

    def cost(w):
        total = 0.0
        for i in range(len(X_train)):
            out = circuit(X_train[i], w)
            total += (out - y_binary[i]) ** 2
        return total / len(X_train)

    for epoch in range(epochs):
        weights = opt.step(cost, weights)
        if (epoch + 1) % 10 == 0:
            loss = float(cost(weights))
            print(f"    Epoch {epoch+1}/{epochs}: loss={loss:.4f}")
    return weights, circuit


def predict_vqc(circuit, weights, X):
    preds = []
    for x in X:
        out = float(circuit(x, weights))
        prob_high = (1.0 + out) / 2.0
        preds.append(prob_high)
    return np.array(preds)


def main():
    print("=" * 70)
    print("Quantum_AgroX – Variational Quantum Classifier (VQC)")
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
    print(f"  Features encoded to {n_qubits} qubits (angle embedding)")
    print(f"  VQC: {n_layers} layers, {n_qubits} qubits")
    print("\n--- VQC Classification (2-fold CV, small sample due to simulation cost) ---")
    n_fold_vqc = min(2, n_splits)
    skf = StratifiedKFold(n_splits=n_fold_vqc, shuffle=True, random_state=config["evaluation"]["random_state"])
    fold_metrics = {"accuracy": [], "balanced_accuracy": [], "f1": [], "roc_auc": []}
    max_samples = min(40, len(X_encoded))
    X_sub_idx = np.random.RandomState(42).choice(len(X_encoded), max_samples, replace=False)
    X_sub = X_encoded[X_sub_idx]
    y_sub = y[X_sub_idx]
    print(f"  Using {len(X_sub)} samples for VQC training")
    for fold, (train_idx, test_idx) in enumerate(skf.split(X_sub, y_sub)):
        print(f"\n  Fold {fold + 1}/{n_fold_vqc}")
        X_train, X_test = X_sub[train_idx], X_sub[test_idx]
        y_train, y_test = y_sub[train_idx], y_sub[test_idx]
        weights, circuit = train_vqc(X_train, y_train, n_qubits, n_layers, lr=0.1, epochs=30)
        y_prob = predict_vqc(circuit, weights, X_test)
        y_pred = (y_prob >= 0.5).astype(int)
        metrics = {
            "accuracy": accuracy_score(y_test, y_pred),
            "balanced_accuracy": balanced_accuracy_score(y_test, y_pred),
            "f1": f1_score(y_test, y_pred, zero_division=0),
        }
        try:
            metrics["roc_auc"] = roc_auc_score(y_test, y_prob)
        except ValueError:
            metrics["roc_auc"] = 0.0
        for k in fold_metrics:
            fold_metrics[k].append(metrics[k])
        print(f"    Fold {fold+1}: Acc={metrics['accuracy']:.4f}, "
              f"F1={metrics['f1']:.4f}, ROC-AUC={metrics['roc_auc']:.4f}")
    avg = {k: round(np.mean(v), 4) for k, v in fold_metrics.items()}
    std = {k: round(np.std(v), 4) for k, v in fold_metrics.items()}
    print(f"\n  VQC Summary: Acc={avg['accuracy']}±{std['accuracy']}, "
          f"F1={avg['f1']}±{std['f1']}, ROC-AUC={avg['roc_auc']}±{std['roc_auc']}")
    print("\n--- VQC Regression Variant (mortality fraction) ---")
    y_reg_path = processed_dir / "y_mortality_regression.npy"
    y_reg = np.load(y_reg_path)
    y_reg_sub = y_reg[X_sub_idx]
    dev_reg = qml.device("default.qubit", wires=n_qubits)

    @qml.qnode(dev_reg, diff_method="parameter-shift")
    def vqc_regression(inputs, weights):
        qml.AngleEmbedding(inputs, wires=range(n_qubits), rotation="Y")
        for layer in range(n_layers):
            for i in range(n_qubits):
                qml.RX(weights[layer, i, 0], wires=i)
                qml.RY(weights[layer, i, 1], wires=i)
                qml.RZ(weights[layer, i, 2], wires=i)
            for i in range(n_qubits - 1):
                qml.CNOT(wires=[i, i + 1])
            qml.CNOT(wires=[n_qubits - 1, 0])
        return qml.expval(qml.PauliZ(0))

    split = len(X_sub) // 2
    X_tr_r, X_te_r = X_sub[:split], X_sub[split:]
    y_tr_r, y_te_r = y_reg_sub[:split], y_reg_sub[split:]
    weights_r = np.random.RandomState(42).uniform(-0.1, 0.1, (n_layers, n_qubits, 3), requires_grad=True)
    opt_r = qml.AdamOptimizer(stepsize=0.05)

    def cost_reg(w):
        total = 0.0
        for i in range(len(X_tr_r)):
            out = vqc_regression(X_tr_r[i], w)
            pred = (1.0 + out) / 2.0
            total += (pred - y_tr_r[i]) ** 2
        return total / len(X_tr_r)

    for epoch in range(20):
        weights_r = opt_r.step(cost_reg, weights_r)
        if (epoch + 1) % 5 == 0:
            loss = float(cost_reg(weights_r))
            print(f"    Reg Epoch {epoch+1}/20: MSE={loss:.4f}")
    y_pred_r = []
    for x in X_te_r:
        out = float(vqc_regression(x, weights_r))
        pred = (1.0 + out) / 2.0
        y_pred_r.append(pred)
    y_pred_r = np.array(y_pred_r)
    from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
    reg_mae = mean_absolute_error(y_te_r, y_pred_r)
    reg_rmse = np.sqrt(mean_squared_error(y_te_r, y_pred_r))
    reg_r2 = r2_score(y_te_r, y_pred_r)
    print(f"  VQC Regression: MAE={reg_mae:.4f}, RMSE={reg_rmse:.4f}, R²={reg_r2:.4f}")
    results = {
        "vqc_classifier": {"mean": avg, "std": std, "sample_size": len(X_sub)},
        "vqc_regressor": {
            "MAE": round(reg_mae, 4),
            "RMSE": round(reg_rmse, 4),
            "R2": round(reg_r2, 4),
        },
        "config": {
            "n_qubits": n_qubits,
            "n_layers": n_layers,
            "optimizer": "SGD (classifier) / Adam (regressor)",
        },
    }
    results_dir = base_dir / "results" / "models"
    results_dir.mkdir(parents=True, exist_ok=True)
    with open(results_dir / "vqc_classifier.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {results_dir / 'vqc_classifier.json'}")
    print("\n" + "=" * 70)
    print("VQC Classifier complete.")
    print("=" * 70)


if __name__ == "__main__":
    main()
