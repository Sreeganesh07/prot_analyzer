"""
ProteinScope Pathogenicity Benchmark Script.

Evaluates the missense variant pathogenicity classifier on ClinVar-labeled TP53
variants using stratified cross-validation and a held-out test set.

Usage:
    python benchmark.py
    python benchmark.py --accession P04637 --offline-blast
    python benchmark.py --output ./output --n-estimators 300
"""

import argparse
import json
import logging
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd

from proteinscope.protein_analyzer import ProteinAnalyzer
from proteinscope.stages.clinvar import load_clinvar_variants
from proteinscope.stages.variant_features import extract_feature_matrix
from proteinscope.stages.pathogenicity import (
    train_pathogenicity_classifier,
    plot_confusion_matrix,
    plot_feature_importances,
)

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("ProteinScope.Benchmark")


def run_benchmark(
    accession: str = "P04637",
    output_dir: str = "./output",
    config_path: str = "config/config.yaml",
    offline_blast: bool = False,
    offline: bool = False,
    n_estimators: int = 200,
    test_size: float = 0.2,
    random_state: int = 42,
) -> dict:
    """
    Full benchmark evaluation of the pathogenicity classifier.

    1. Runs the base pipeline for the accession (or uses cached output)
    2. Loads ClinVar labeled variants
    3. Extracts features for all labeled variants
    4. Trains a classifier with stratified train/test split
    5. Reports metrics: accuracy, AUROC, precision, recall, F1, confusion matrix
    6. Saves everything to output/<accession>/benchmark/
    """
    start_time = time.time()
    logger.info(f"Starting pathogenicity benchmark for {accession}")

    # Step 1: Run base pipeline (or use cached results)
    logger.info("Step 1: Running base pipeline stages...")
    analyzer = ProteinAnalyzer(
        accession=accession,
        output_dir=output_dir,
        config_path=config_path if os.path.exists(config_path) else None,
    )

    # Run base stages to populate internal state
    if not analyzer.validate_input():
        logger.error("Input validation failed")
        return {"error": "Input validation failed"}

    analyzer.fetch_from_uniprot()
    analyzer.analyze_sequence()
    analyzer.blast_search(hits=10, run_online=not offline_blast and not offline)
    analyzer.fetch_structures()
    analyzer.validate_structure()
    analyzer.comparative_analysis()

    sequence = analyzer.ingest_data.get("sequence", "")
    if not sequence:
        logger.error("No sequence available for feature extraction")
        return {"error": "No sequence available"}

    # Step 2: Load ClinVar variants
    logger.info("Step 2: Loading ClinVar variants...")
    clinvar_cfg = analyzer.config.get("clinvar", {})
    clinvar_df = load_clinvar_variants(
        gene_symbol=clinvar_cfg.get("gene_symbol", "TP53"),
        cache_dir=clinvar_cfg.get("cache_dir", "./data/clinvar"),
        offline=offline,
    )

    if clinvar_df.empty:
        logger.error("No ClinVar variants loaded — cannot benchmark")
        return {"error": "No ClinVar data"}

    logger.info(
        f"Loaded {len(clinvar_df)} ClinVar variants: "
        f"{(clinvar_df['label'] == 'pathogenic').sum()} pathogenic, "
        f"{(clinvar_df['label'] == 'benign').sum()} benign"
    )

    # Step 3: Extract features for all variants
    logger.info("Step 3: Extracting features for all ClinVar variants...")
    pdb_file = None
    best = analyzer.structure_data.get("best_structure")
    if best:
        pdb_file = best.get("file_path")

    blast_hits = analyzer.blast_data.get("hits", [])
    ptm_sites = analyzer.comparative_data.get("ptm_sites", [])
    if not ptm_sites:
        ptm_sites = analyzer.ingest_data.get("ptm_features", [])

    variants_list = clinvar_df.to_dict("records")

    feature_dicts, feature_names = extract_feature_matrix(
        variants=variants_list,
        query_sequence=sequence,
        pdb_file=pdb_file,
        blast_hits=blast_hits,
        ptm_sites=ptm_sites,
        include_ml_deltas=True,
    )

    # Match labels to features (some variants may have been skipped due to ref AA mismatch)
    labels = []
    feat_idx = 0
    for var in variants_list:
        pos = var.get("position", 0)
        ref = var.get("ref_aa", "")
        if sequence and 0 < pos <= len(sequence) and sequence[pos - 1] == ref:
            if feat_idx < len(feature_dicts):
                labels.append(var["label"])
                feat_idx += 1

    # Trim to consistent length
    min_len = min(len(feature_dicts), len(labels))
    feature_dicts = feature_dicts[:min_len]
    labels = labels[:min_len]

    logger.info(f"Feature extraction complete: {len(feature_dicts)} variants with {len(feature_names)} features")

    if len(feature_dicts) < 10:
        logger.error(f"Too few variants ({len(feature_dicts)}) for meaningful benchmark")
        return {"error": "Insufficient data for benchmark"}

    # Step 4: Train and evaluate
    logger.info("Step 4: Training and evaluating classifier...")
    model_dir = os.path.join("data", "models")

    metrics = train_pathogenicity_classifier(
        feature_dicts=feature_dicts,
        labels=labels,
        model_dir=model_dir,
        n_estimators=n_estimators,
        test_size=test_size,
        random_state=random_state,
    )

    # Step 5: Save benchmark results
    logger.info("Step 5: Saving benchmark results...")
    benchmark_dir = os.path.join(output_dir, accession, "benchmark")
    os.makedirs(benchmark_dir, exist_ok=True)

    # Save metrics JSON
    duration = round(time.time() - start_time, 2)
    metrics["benchmark_metadata"] = {
        "accession": accession,
        "total_variants": len(clinvar_df),
        "valid_variants_with_features": len(feature_dicts),
        "n_pathogenic": sum(1 for l in labels if l == "pathogenic"),
        "n_benign": sum(1 for l in labels if l == "benign"),
        "feature_count": len(feature_names),
        "n_estimators": n_estimators,
        "test_size": test_size,
        "random_state": random_state,
        "execution_time_seconds": duration,
    }

    metrics_path = os.path.join(benchmark_dir, "metrics.json")
    with open(metrics_path, "w", encoding="utf-8") as f:
        json.dump(metrics, f, indent=2, default=str)
    logger.info(f"Saved metrics to {metrics_path}")

    # Save plots
    if "confusion_matrix" in metrics:
        cm_path = os.path.join(benchmark_dir, "confusion_matrix.png")
        plot_confusion_matrix(metrics["confusion_matrix"], cm_path)

    if "top_features" in metrics:
        fi_path = os.path.join(benchmark_dir, "feature_importances.png")
        plot_feature_importances(metrics["top_features"], fi_path)

    # Save human-readable benchmark report
    report_path = os.path.join(benchmark_dir, "benchmark_report.txt")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("=" * 70 + "\n")
        f.write("     PROTEINSCOPE PATHOGENICITY CLASSIFIER BENCHMARK REPORT\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"Gene / Accession    : TP53 / {accession}\n")
        f.write(f"Total ClinVar Vars  : {len(clinvar_df)}\n")
        f.write(f"Valid (w/ features) : {len(feature_dicts)}\n")
        f.write(f"  Pathogenic        : {sum(1 for l in labels if l == 'pathogenic')}\n")
        f.write(f"  Benign            : {sum(1 for l in labels if l == 'benign')}\n")
        f.write(f"Feature Count       : {len(feature_names)}\n")
        f.write(f"RF Estimators       : {n_estimators}\n")
        f.write(f"Test Size           : {test_size:.0%}\n\n")

        f.write("── HELD-OUT TEST SET METRICS ──────────────────────────────────────\n")
        f.write(f"Accuracy            : {metrics.get('accuracy', 0):.4f}\n")
        f.write(f"AUROC               : {metrics.get('auroc', 0):.4f}\n")
        f.write(f"Precision           : {metrics.get('precision', 0):.4f}\n")
        f.write(f"Recall              : {metrics.get('recall', 0):.4f}\n")
        f.write(f"F1 Score            : {metrics.get('f1_score', 0):.4f}\n")
        f.write(f"5-Fold CV AUROC     : {metrics.get('cv_5fold_auroc_mean', 0):.4f} "
                f"± {metrics.get('cv_5fold_auroc_std', 0):.4f}\n\n")

        cm = metrics.get("confusion_matrix", {})
        f.write("── CONFUSION MATRIX ──────────────────────────────────────────────\n")
        f.write(f"                    Predicted Benign   Predicted Pathogenic\n")
        f.write(f"  Actual Benign     {cm.get('true_negatives', 0):>10}          {cm.get('false_positives', 0):>10}\n")
        f.write(f"  Actual Pathogenic {cm.get('false_negatives', 0):>10}          {cm.get('true_positives', 0):>10}\n\n")

        f.write("── TOP 15 FEATURE IMPORTANCES ────────────────────────────────────\n")
        top_feats = metrics.get("top_features", {})
        for i, (name, imp) in enumerate(top_feats.items(), 1):
            bar = "█" * int(imp * 200)
            f.write(f"  {i:2d}. {name:<35s} {imp:.4f}  {bar}\n")

        f.write(f"\n{'='*70}\n")
        f.write(f"Execution Time      : {duration} seconds\n")
        f.write(f"Results saved to    : {os.path.abspath(benchmark_dir)}\n")
        f.write(f"{'='*70}\n")

    logger.info(f"Saved benchmark report to {report_path}")

    # Print summary to console
    print(f"\n{'='*60}")
    print(f"  BENCHMARK RESULTS for TP53 ({accession})")
    print(f"{'='*60}")
    print(f"  Variants Used     : {len(feature_dicts)}")
    print(f"  Accuracy          : {metrics.get('accuracy', 0):.4f}")
    print(f"  AUROC             : {metrics.get('auroc', 0):.4f}")
    print(f"  F1 Score          : {metrics.get('f1_score', 0):.4f}")
    print(f"  5-Fold CV AUROC   : {metrics.get('cv_5fold_auroc_mean', 0):.4f} ± {metrics.get('cv_5fold_auroc_std', 0):.4f}")
    print(f"{'='*60}")
    print(f"  Full report: {report_path}")

    return metrics


def main():
    parser = argparse.ArgumentParser(
        description="ProteinScope Pathogenicity Classifier Benchmark",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--accession", type=str, default="P04637",
        help="UniProt accession code for the target protein"
    )
    parser.add_argument(
        "--output", type=str, default="./output",
        help="Base output directory"
    )
    parser.add_argument(
        "--config", type=str, default="config/config.yaml",
        help="Path to pipeline configuration file"
    )
    parser.add_argument(
        "--offline-blast", action="store_true",
        help="Skip online NCBI BLAST (use cached results)"
    )
    parser.add_argument(
        "--offline", action="store_true",
        help="Run fully offline (no API calls)"
    )
    parser.add_argument(
        "--n-estimators", type=int, default=200,
        help="Number of Random Forest estimators"
    )
    parser.add_argument(
        "--test-size", type=float, default=0.2,
        help="Fraction of data for held-out test set"
    )
    parser.add_argument(
        "--random-state", type=int, default=42,
        help="Random seed for reproducibility"
    )

    args = parser.parse_args()

    run_benchmark(
        accession=args.accession,
        output_dir=args.output,
        config_path=args.config,
        offline_blast=args.offline_blast,
        offline=args.offline,
        n_estimators=args.n_estimators,
        test_size=args.test_size,
        random_state=args.random_state,
    )


if __name__ == "__main__":
    main()
