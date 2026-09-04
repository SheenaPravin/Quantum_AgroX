# Quantum_AgroX_QML

Quantum Machine Learning Platform for Botanical Acaricide & Agrochemical Discovery.

## Overview

Quantum_AgroX integrates classical machine learning with quantum machine learning (QML) models for efficacy prediction, target interaction analysis, and candidate ranking of botanical acaricides. This project uses the *Commiphora swynnertonii* benchmark dataset (GC-MS phytochemicals, tick bioassays, LC50/LC90, AChE modelling and docking).

## Project Structure

```
Quantum_AgroX_QML_Project/
├── data/
│   ├── raw/                  # Source CSV datasets
│   └── processed/            # QML-ready numpy arrays and normalized CSVs
├── src/
│   ├── 01_prepare_dataset.py
│   ├── 02_classical_baselines.py
│   ├── 03_quantum_kernel_regression.py
│   ├── 04_quantum_kernel_classifier.py
│   ├── 05_vqc_classifier.py
│   └── main.py
├── config/
│   └── project_config.yaml
├── results/
│   ├── models/               # Model results (JSON)
│   └── figures/              # Generated plots
└── requirements.txt
```

## Pipeline Steps

| Step | Script | Description |
|------|--------|-------------|
| 1 | `01_prepare_dataset.py` | Ingest raw CSVs, encode features, save QML-ready arrays |
| 2 | `02_classical_baselines.py` | RF, XGBoost, SVM baselines (regression + classification) |
| 3 | `03_quantum_kernel_regression.py` | Quantum fidelity kernel + Kernel Ridge Regression |
| 4 | `04_quantum_kernel_classifier.py` | Quantum kernel SVM (fidelity + feature map approaches) |
| 5 | `05_vqc_classifier.py` | Variational Quantum Classifier + regressor |

## Quick Start

```bash
# Install dependencies
pip install -r requirements.txt

# Run full pipeline
python src/main.py

# Run specific steps
python src/main.py --steps prepare classical

# Run individual scripts
python src/01_prepare_dataset.py
python src/02_classical_baselines.py
python src/03_quantum_kernel_regression.py
python src/04_quantum_kernel_classifier.py
python src/05_vqc_classifier.py
```

## Quantum Models

### 1. Quantum Fidelity Kernel + KRR
Uses the quantum fidelity kernel (|⟨ψ(x₁)|ψ(x₂)⟩|²) combined with Kernel Ridge Regression for mortality fraction prediction.

### 2. Quantum Kernel SVM (QSVM-style)
Two approaches:
- **Fidelity Kernel SVM**: Precomputed fidelity kernel with SVC
- **Quantum Feature Map + RBF SVM**: VQE-style feature map outputs fed to classical RBF SVM

### 3. Variational Quantum Classifier (VQC)
Parameterized quantum circuit with angle embedding, trainable RX/RY/RZ rotations, CNOT entanglement, and PauliZ measurement. Supports both classification and regression via expectation value mapping.

## Evaluation

- **Regression**: MAE, RMSE, R² (5-fold CV)
- **Classification**: Accuracy, Balanced Accuracy, F1, ROC-AUC (5-fold stratified CV)
- All QML results compared against classical baselines with identical data splits

## Configuration

All parameters are in `config/project_config.yaml`:
- Quantum: qubit count, circuit depth, layers, shots, feature map type
- Models: hyperparameters for RF, XGBoost, SVM, KRR
- Evaluation: CV folds, test size, metrics

## Important Safeguards

- Small benchmark dataset → for pipeline development, not quantum advantage claims
- Docking scores are not proof of AChE inhibition
- Aggregate means used; raw replicate-level data recommended for production
- Full provenance tracking from raw input to prediction

## Technology Stack

- Python, pandas, NumPy, scikit-learn, XGBoost
- PennyLane (quantum ML)
- RDKit (molecular descriptors)
- PyYAML (configuration)
