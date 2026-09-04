"""
Pathogenicity Classifier for Missense Variants.

Trains, serializes, and runs a Random Forest classifier on per-variant features
extracted from the ProteinScope pipeline stages. Uses ClinVar pathogenic/benign
labels as ground truth.

Supports:
  - Training from ClinVar-labeled variant features
  - Serialization/deserialization via joblib
  - Single-variant prediction with confidence score
  - Feature importance ranking for interpretability
  - Stratified train/test evaluation with full metrics
"""

import json
import logging
import os
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    accuracy_score,
    classification_report,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedKFold, cross_val_score, train_test_split
from sklearn.preprocessing import StandardScaler

logger = logging.getLogger("ProteinScope.Pathogenicity")


def _try_import_joblib():
    """Import joblib, falling back to sklearn's bundled version."""
    try:
        import joblib
        return joblib
    except ImportError:
        from sklearn.utils import _joblib
        return _joblib


def train_pathogenicity_classifier(
    feature_dicts: List[Dict[str, float]],
    labels: List[str],
    model_dir: str = "./data/models",
    n_estimators: int = 200,
    test_size: float = 0.2,
    random_state: int = 42,
) -> Dict[str, Any]:
    """
    Train a Random Forest classifier on labeled variant features.

    Args:
        feature_dicts: List of feature dictionaries (one per variant)
        labels: List of labels ('pathogenic' or 'benign')
        model_dir: Directory to save the trained model
        n_estimators: Number of trees in the Random Forest
        test_size: Fraction of data for the held-out test set
        random_state: Random seed for reproducibility

    Returns:
        Dict with training results including metrics, model path, and feature importances
    """
    logger.info(f"Training pathogenicity classifier on {len(feature_dicts)} samples")
    os.makedirs(model_dir, exist_ok=True)

    # Build feature matrix
    df = pd.DataFrame(feature_dicts)

    # Handle any NaN/inf values
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.fillna(0.0)

    feature_names = df.columns.tolist()
    X = df.values
    y = np.array([1 if label == "pathogenic" else 0 for label in labels])

    logger.info(
        f"Feature matrix shape: {X.shape}, "
        f"Pathogenic: {np.sum(y == 1)}, Benign: {np.sum(y == 0)}"
    )

    if len(set(y)) < 2:
        logger.error("Need at least 2 classes (pathogenic + benign) for training")
        return {"error": "Insufficient class diversity"}

    # Stratified train/test split
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=test_size, random_state=random_state, stratify=y
    )

    logger.info(
        f"Train set: {len(X_train)} samples "
        f"(path: {np.sum(y_train == 1)}, benign: {np.sum(y_train == 0)})"
    )
    logger.info(
        f"Test set:  {len(X_test)} samples "
        f"(path: {np.sum(y_test == 1)}, benign: {np.sum(y_test == 0)})"
    )

    # Train Random Forest with balanced class weights
    clf = RandomForestClassifier(
        n_estimators=n_estimators,
        class_weight="balanced",
        random_state=random_state,
        n_jobs=-1,
        max_depth=None,
        min_samples_split=5,
        min_samples_leaf=2,
    )
    clf.fit(X_train, y_train)

    # Evaluate on test set
    y_pred = clf.predict(X_test)
    y_prob = clf.predict_proba(X_test)[:, 1]

    accuracy = accuracy_score(y_test, y_pred)
    precision = precision_score(y_test, y_pred, zero_division=0)
    recall = recall_score(y_test, y_pred, zero_division=0)
    f1 = f1_score(y_test, y_pred, zero_division=0)

    try:
        auroc = roc_auc_score(y_test, y_prob)
    except ValueError:
        auroc = 0.0

    cm = confusion_matrix(y_test, y_pred)
    cm_dict = {
        "true_negatives": int(cm[0, 0]),
        "false_positives": int(cm[0, 1]),
        "false_negatives": int(cm[1, 0]),
        "true_positives": int(cm[1, 1]),
    }

    # Classification report
    cls_report = classification_report(
        y_test, y_pred, target_names=["benign", "pathogenic"], output_dict=True
    )

    # Feature importances
    importances = clf.feature_importances_
    importance_ranking = sorted(
        zip(feature_names, importances),
        key=lambda x: x[1],
        reverse=True,
    )
    top_features = {name: round(float(imp), 4) for name, imp in importance_ranking[:20]}

    # Cross-validation score (5-fold on full training data)
    try:
        cv_scores = cross_val_score(
            RandomForestClassifier(
                n_estimators=n_estimators,
                class_weight="balanced",
                random_state=random_state,
                n_jobs=-1,
            ),
            X,
            y,
            cv=StratifiedKFold(n_splits=5, shuffle=True, random_state=random_state),
            scoring="roc_auc",
        )
        cv_auroc_mean = round(float(np.mean(cv_scores)), 4)
        cv_auroc_std = round(float(np.std(cv_scores)), 4)
    except Exception as e:
        logger.warning(f"Cross-validation failed: {e}")
        cv_auroc_mean = 0.0
        cv_auroc_std = 0.0

    # Serialize model
    joblib = _try_import_joblib()
    model_path = os.path.join(model_dir, "tp53_pathogenicity_rf.joblib")
    model_bundle = {
        "classifier": clf,
        "feature_names": feature_names,
        "n_estimators": n_estimators,
        "random_state": random_state,
    }
    joblib.dump(model_bundle, model_path)
    logger.info(f"Saved trained model to {model_path}")

    # Also save a scaler (fit on training data) for optional use
    scaler = StandardScaler()
    scaler.fit(X_train)
    scaler_path = os.path.join(model_dir, "tp53_feature_scaler.joblib")
    joblib.dump(scaler, scaler_path)

    metrics = {
        "accuracy": round(accuracy, 4),
        "precision": round(precision, 4),
        "recall": round(recall, 4),
        "f1_score": round(f1, 4),
        "auroc": round(auroc, 4),
        "confusion_matrix": cm_dict,
        "classification_report": cls_report,
        "cv_5fold_auroc_mean": cv_auroc_mean,
        "cv_5fold_auroc_std": cv_auroc_std,
        "train_samples": len(X_train),
        "test_samples": len(X_test),
        "n_features": len(feature_names),
        "top_features": top_features,
        "model_path": model_path,
    }

    logger.info(
        f"Classifier trained — Accuracy: {accuracy:.3f}, AUROC: {auroc:.3f}, "
        f"F1: {f1:.3f}, 5-Fold CV AUROC: {cv_auroc_mean:.3f} ± {cv_auroc_std:.3f}"
    )

    return metrics


def load_trained_model(
    model_path: str = "./data/models/tp53_pathogenicity_rf.joblib",
) -> Optional[Dict[str, Any]]:
    """
    Load a previously trained and serialized pathogenicity classifier.

    Returns:
        Dict with 'classifier' (RandomForestClassifier) and 'feature_names' (list),
        or None if loading fails.
    """
    if not os.path.exists(model_path):
        logger.error(f"Trained model not found at {model_path}")
        return None

    try:
        joblib = _try_import_joblib()
        model_bundle = joblib.load(model_path)
        logger.info(f"Loaded trained model from {model_path}")
        return model_bundle
    except Exception as e:
        logger.error(f"Failed to load model from {model_path}: {e}")
        return None


def predict_variant(
    feature_dict: Dict[str, float],
    model_bundle: Optional[Dict[str, Any]] = None,
    model_path: str = "./data/models/tp53_pathogenicity_rf.joblib",
    confidence_threshold: float = 0.5,
) -> Dict[str, Any]:
    """
    Predict pathogenicity for a single variant given its feature vector.

    Args:
        feature_dict: Feature dictionary for the variant
        model_bundle: Pre-loaded model bundle (classifier + feature_names)
        model_path: Path to serialized model (used if model_bundle is None)
        confidence_threshold: Probability threshold for 'pathogenic' call

    Returns:
        Dict with prediction, confidence, and supporting details
    """
    if model_bundle is None:
        model_bundle = load_trained_model(model_path)

    if model_bundle is None:
        return {
            "prediction": "unknown",
            "confidence": 0.0,
            "error": "No trained model available",
        }

    clf = model_bundle["classifier"]
    feature_names = model_bundle["feature_names"]

    # Build feature vector in the correct order
    feature_vector = np.array(
        [feature_dict.get(name, 0.0) for name in feature_names]
    ).reshape(1, -1)

    # Handle NaN/inf
    feature_vector = np.nan_to_num(feature_vector, nan=0.0, posinf=0.0, neginf=0.0)

    # Predict
    prob = clf.predict_proba(feature_vector)[0]  # [prob_benign, prob_pathogenic]
    pathogenic_prob = float(prob[1])
    benign_prob = float(prob[0])

    prediction = "pathogenic" if pathogenic_prob >= confidence_threshold else "benign"
    confidence = pathogenic_prob if prediction == "pathogenic" else benign_prob

    # Feature importances for this prediction (top contributors)
    importances = clf.feature_importances_
    feature_contributions = sorted(
        zip(feature_names, importances, [feature_dict.get(n, 0.0) for n in feature_names]),
        key=lambda x: x[1],
        reverse=True,
    )
    top_contributors = [
        {
            "feature": name,
            "importance": round(float(imp), 4),
            "value": round(float(val), 4),
        }
        for name, imp, val in feature_contributions[:10]
    ]

    return {
        "prediction": prediction,
        "confidence": round(confidence, 4),
        "pathogenic_probability": round(pathogenic_prob, 4),
        "benign_probability": round(benign_prob, 4),
        "confidence_threshold": confidence_threshold,
        "top_feature_contributors": top_contributors,
    }


def generate_justification(
    feature_dict: Dict[str, float],
    prediction_result: Dict[str, Any],
    ref_aa: str,
    alt_aa: str,
    position: int,
) -> str:
    """
    Generate a plain-text structural/conservation justification for the
    pathogenicity prediction, built from the feature values.

    Returns a multi-line human-readable explanation.
    """
    lines = []

    # 1. Burial depth
    burial_class = feature_dict.get("burial_class", 0.0)
    burial_count = feature_dict.get("burial_neighbour_count", 0.0)
    if burial_class >= 2.0:
        lines.append(
            f"Residue {position} is BURIED ({int(burial_count)} neighbors within 8Å) "
            "— mutations at buried sites are more likely to disrupt protein folding"
        )
    elif burial_class >= 1.0:
        lines.append(
            f"Residue {position} is PARTIALLY BURIED ({int(burial_count)} neighbors within 8Å) "
            "— may affect local packing"
        )
    else:
        lines.append(
            f"Residue {position} is SURFACE-EXPOSED ({int(burial_count)} neighbors within 8Å) "
            "— surface mutations are often better tolerated"
        )

    # 2. Substitution properties
    grantham = feature_dict.get("grantham_distance", 0.0)
    if grantham > 100:
        lines.append(
            f"RADICAL substitution {ref_aa}→{alt_aa} (Grantham distance: {int(grantham)}) "
            "— large physicochemical difference"
        )
    elif grantham > 60:
        lines.append(
            f"NON-CONSERVATIVE substitution {ref_aa}→{alt_aa} (Grantham distance: {int(grantham)}) "
            "— moderate physicochemical change"
        )
    else:
        lines.append(
            f"Conservative substitution {ref_aa}→{alt_aa} (Grantham distance: {int(grantham)}) "
            "— relatively similar physicochemical properties"
        )

    # 3. Charge change
    charge_change = feature_dict.get("charge_change", 0.0)
    if abs(charge_change) > 0.5:
        if charge_change > 0:
            lines.append("CHARGE GAINED — introduces positive charge at this position")
        else:
            lines.append("CHARGE LOST — removes charge from this position")

    # 4. Conservation
    conservation = feature_dict.get("conservation_score", 0.0)
    diversity = feature_dict.get("homolog_diversity", 0.0)
    n_homologs = feature_dict.get("n_homologs_aligned", 0.0)
    if n_homologs > 0:
        if conservation > 0.8:
            lines.append(
                f"HIGH evolutionary conservation at position {position} "
                f"({conservation:.0%} of {int(n_homologs)} homologs have {ref_aa}) "
                "— strong constraint suggests functional importance"
            )
        elif conservation > 0.5:
            lines.append(
                f"MODERATE evolutionary conservation at position {position} "
                f"({conservation:.0%} of {int(n_homologs)} homologs have {ref_aa})"
            )
        else:
            lines.append(
                f"LOW evolutionary conservation at position {position} "
                f"({conservation:.0%} of {int(n_homologs)} homologs have {ref_aa}) "
                "— position is variable across homologs"
            )

        if diversity < 3:
            lines.append(
                f"LOW homolog tolerance — only {int(diversity)} distinct amino acids "
                "observed at this position across homologs"
            )
    else:
        lines.append("Conservation data unavailable (no BLAST homologs aligned at this position)")

    # 5. PTM proximity
    ptm_dist = feature_dict.get("min_ptm_distance", 999.0)
    is_ptm = feature_dict.get("is_ptm_site", 0.0)
    if is_ptm > 0:
        lines.append(
            f"Position {position} IS a known post-translational modification site "
            "— mutation may directly disrupt regulatory modification"
        )
    elif ptm_dist <= 5:
        lines.append(
            f"Position is {int(ptm_dist)} residues from a known PTM site "
            "— may affect nearby post-translational regulation"
        )

    # 6. Ramachandran classification
    rama_zone = feature_dict.get("ramachandran_zone", 0.0)
    if rama_zone == 0.0:
        lines.append(
            "Ramachandran zone: FAVORED — wild-type backbone conformation is well-defined"
        )
    elif rama_zone == 1.0:
        lines.append(
            "Ramachandran zone: ALLOWED — backbone conformation is acceptable but not optimal"
        )
    elif rama_zone >= 2.0:
        lines.append(
            "Ramachandran zone: OUTLIER — unusual backbone conformation at this position"
        )

    return "\n".join(f"• {line}" for line in lines)


def plot_confusion_matrix(
    cm_dict: Dict[str, int],
    output_path: str,
) -> None:
    """Generate and save a confusion matrix heatmap."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cm = np.array([
        [cm_dict["true_negatives"], cm_dict["false_positives"]],
        [cm_dict["false_negatives"], cm_dict["true_positives"]],
    ])

    fig, ax = plt.subplots(figsize=(7, 6), dpi=200)

    # Use a colormap
    im = ax.imshow(cm, interpolation="nearest", cmap="Blues")
    ax.set_title("Pathogenicity Classifier — Confusion Matrix", fontsize=14, fontweight="bold", pad=12)
    plt.colorbar(im, ax=ax, shrink=0.8)

    labels = ["Benign", "Pathogenic"]
    tick_marks = [0, 1]
    ax.set_xticks(tick_marks)
    ax.set_xticklabels(labels, fontsize=12)
    ax.set_yticks(tick_marks)
    ax.set_yticklabels(labels, fontsize=12)

    # Annotate cells
    thresh = cm.max() / 2.0
    for i in range(2):
        for j in range(2):
            ax.text(
                j, i, format(cm[i, j], "d"),
                ha="center", va="center",
                fontsize=16, fontweight="bold",
                color="white" if cm[i, j] > thresh else "black",
            )

    ax.set_ylabel("Actual Label", fontsize=12)
    ax.set_xlabel("Predicted Label", fontsize=12)
    plt.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)
    logger.info(f"Saved confusion matrix plot to {output_path}")


def plot_feature_importances(
    top_features: Dict[str, float],
    output_path: str,
) -> None:
    """Generate and save a horizontal bar chart of feature importances."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    names = list(top_features.keys())
    values = list(top_features.values())

    # Reverse for horizontal bar (top feature at top)
    names = names[::-1]
    values = values[::-1]

    fig, ax = plt.subplots(figsize=(10, 8), dpi=200)
    bars = ax.barh(range(len(names)), values, color="#2563eb", alpha=0.85, edgecolor="white")
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=10)
    ax.set_xlabel("Feature Importance (Gini)", fontsize=12)
    ax.set_title("Top Feature Importances — Pathogenicity Classifier", fontsize=14, fontweight="bold", pad=12)
    ax.grid(axis="x", alpha=0.3, linestyle="--")

    plt.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)
    logger.info(f"Saved feature importance plot to {output_path}")
