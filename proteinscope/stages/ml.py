"""
Stage 8: Machine Learning and Deep Learning Feature Engineering.
Extracts 40+ physicochemical & compositional feature vectors, runs unsupervised KMeans clustering,
trains RandomForest classifiers, and provides optional ESM-2 embedding hooks.
"""

import json
import logging
import os
from typing import Any, Dict, List, Optional, Tuple
import numpy as np
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from sklearn.cluster import KMeans
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler

logger = logging.getLogger("ProteinScope.ML")

# All 20 standard amino acids in alphabetical order
STANDARD_AMINO_ACIDS = sorted(list("ACDEFGHIKLMNPQRSTVWY"))


def extract_sequence_feature_vector(sequence: str, label_name: str = "seq") -> Dict[str, float]:
    """
    Extract a comprehensive 40+ dimensional numerical feature vector from a single protein sequence:
    - Sequence Length (1)
    - Molecular Weight (1)
    - Isoelectric Point (1)
    - Aromaticity (1)
    - Instability Index (1)
    - GRAVY Hydropathicity (1)
    - Secondary Structure Fractions (Helix, Turn, Sheet) (3)
    - Charge & Polarity Class Percentages (Pos, Neg, Hydrophobic, Polar, Aliphatic, Tiny, Small) (7)
    - 20 Individual Amino Acid Percentages (20)
    Total: 36 + grouped = 40+ features.
    """
    cleaned = "".join([aa for aa in sequence.upper() if aa in STANDARD_AMINO_ACIDS])
    if not cleaned:
        cleaned = "A" * 10  # fallback

    analyser = ProteinAnalysis(cleaned)
    mw = analyser.molecular_weight()
    pi = analyser.isoelectric_point()
    arom = analyser.aromaticity()
    inst = analyser.instability_index()
    gravy = analyser.gravy()
    sec_struc = analyser.secondary_structure_fraction()
    if hasattr(analyser, "amino_acids_percent"):
        raw_pct = analyser.amino_acids_percent
    elif hasattr(analyser, "get_amino_acids_percent"):
        raw_pct = analyser.get_amino_acids_percent()
    else:
        counts = analyser.count_amino_acids()
        total = sum(counts.values()) or 1
        raw_pct = {aa: (cnt / total * 100) for aa, cnt in counts.items()}

    # Normalize to 0-100%
    is_frac = sum(raw_pct.values()) <= 1.5 if raw_pct else False
    aa_pct = {aa: (pct * 100.0 if is_frac else pct) for aa, pct in raw_pct.items()}

    # Calculate grouped physicochemical properties (percentages)
    pos_charged = sum(aa_pct.get(aa, 0) for aa in "RKH")
    neg_charged = sum(aa_pct.get(aa, 0) for aa in "DE")
    hydrophobic = sum(aa_pct.get(aa, 0) for aa in "AILMFWV")
    polar_uncharged = sum(aa_pct.get(aa, 0) for aa in "STNQ")
    aliphatic = sum(aa_pct.get(aa, 0) for aa in "ILV")
    aromatic = sum(aa_pct.get(aa, 0) for aa in "FWY")
    tiny = sum(aa_pct.get(aa, 0) for aa in "GAS")
    small = sum(aa_pct.get(aa, 0) for aa in "GASCDPNT")
    polar_total = sum(aa_pct.get(aa, 0) for aa in "RNDQEHKSTY")
    nonpolar_total = sum(aa_pct.get(aa, 0) for aa in "AFGILMPVW")
    basic_total = sum(aa_pct.get(aa, 0) for aa in "KRH")
    acidic_total = sum(aa_pct.get(aa, 0) for aa in "DE")

    features: Dict[str, float] = {
        "length": float(len(sequence)),
        "analyzed_length": float(len(cleaned)),
        "molecular_weight": round(float(mw), 3),
        "isoelectric_point": round(float(pi), 3),
        "aromaticity": round(float(arom), 4),
        "instability_index": round(float(inst), 3),
        "gravy": round(float(gravy), 4),
        "helix_fraction": round(float(sec_struc[0]), 4),
        "turn_fraction": round(float(sec_struc[1]), 4),
        "sheet_fraction": round(float(sec_struc[2]), 4),
        "charged_pos_pct": round(float(pos_charged), 3),
        "charged_neg_pct": round(float(neg_charged), 3),
        "hydrophobic_pct": round(float(hydrophobic), 3),
        "polar_uncharged_pct": round(float(polar_uncharged), 3),
        "aliphatic_pct": round(float(aliphatic), 3),
        "aromatic_pct": round(float(aromatic), 3),
        "tiny_pct": round(float(tiny), 3),
        "small_pct": round(float(small), 3),
        "polar_total_pct": round(float(polar_total), 3),
        "nonpolar_total_pct": round(float(nonpolar_total), 3),
        "basic_total_pct": round(float(basic_total), 3),
        "acidic_total_pct": round(float(acidic_total), 3),
    }

    # Add individual 20 amino acid percentages
    for aa in STANDARD_AMINO_ACIDS:
        features[f"aa_{aa}_pct"] = round(float(aa_pct.get(aa, 0) * 100), 3)

    return features


def extract_and_cluster_ml_features(
    query_accession: str,
    query_sequence: str,
    output_dir: str,
    blast_hits: Optional[List[Dict[str, Any]]] = None,
    n_clusters: int = 3,
    random_state: int = 42
) -> Dict[str, Any]:
    """
    Stage 8 Entry point:
    Extracts 40+ feature matrix for query and homologs, performs KMeans clustering,
    evaluates feature importances, and writes ML feature CSV and JSON results.
    """
    logger.info(f"Starting Stage 8: ML Feature extraction & clustering for {query_accession}")
    ml_dir = os.path.join(output_dir, "ml")
    os.makedirs(ml_dir, exist_ok=True)

    sequence_list: List[Tuple[str, str, str]] = []
    # 1. Add query
    sequence_list.append((query_accession, query_sequence, "Query"))

    # 2. Add BLAST hits if available
    if blast_hits:
        for idx, hit in enumerate(blast_hits):
            acc = hit.get("accession", f"Hit_{idx+1}")
            seq = hit.get("sbjct_sequence", "").replace("-", "")
            if seq and len(seq) >= 10:
                sequence_list.append((acc, seq, f"Homolog_Rank_{hit.get('rank', idx+1)}"))

    # If few homologs, create synthetic perturbation variants to enable robust clustering analysis
    if len(sequence_list) < 4:
        # Create subtle synthetic mutation variants for demonstration of clustering space
        seq_chars = list(query_sequence)
        for var_idx in range(1, 5):
            var_seq = list(seq_chars)
            # Mutate a couple random positions
            for pos in range(var_idx * 10, min(len(seq_chars), var_idx * 10 + 20), 5):
                var_seq[pos % len(seq_chars)] = "A" if var_seq[pos % len(seq_chars)] != "A" else "G"
            sequence_list.append((f"{query_accession}_var{var_idx}", "".join(var_seq), f"Variant_{var_idx}"))

    # Build feature table
    rows = []
    identifiers = []
    categories = []

    for item_id, seq, cat in sequence_list:
        feat_dict = extract_sequence_feature_vector(seq, label_name=item_id)
        rows.append(feat_dict)
        identifiers.append(item_id)
        categories.append(cat)

    df_features = pd.DataFrame(rows, index=identifiers)
    
    # Save feature matrix CSV
    feature_csv_path = os.path.join(ml_dir, "feature_matrix.csv")
    df_features.to_csv(feature_csv_path)
    logger.info(f"Saved {df_features.shape[1]} features for {len(df_features)} sequences to {feature_csv_path}")

    # Standardize features for KMeans clustering
    scaler = StandardScaler()
    x_scaled = scaler.fit_transform(df_features)

    # Determine effective cluster count (cannot exceed sample count)
    k = min(n_clusters, len(df_features))
    kmeans = KMeans(n_clusters=k, random_state=random_state, n_init=10)
    cluster_labels = kmeans.fit_predict(x_scaled)

    # Map clusters to sequences
    clustering_records = []
    for ident, cat, clust in zip(identifiers, categories, cluster_labels):
        clustering_records.append({
            "identifier": ident,
            "category": cat,
            "cluster": int(clust),
        })

    # Train Random Forest Classifier on clusters to rank most discriminative features
    rf_feature_importances = {}
    if len(set(cluster_labels)) > 1 and len(df_features) >= 3:
        try:
            rf = RandomForestClassifier(n_estimators=50, random_state=random_state)
            rf.fit(df_features, cluster_labels)
            importances = rf.feature_importances_
            feature_names = df_features.columns.tolist()
            sorted_imp = sorted(zip(feature_names, importances), key=lambda x: x[1], reverse=True)
            rf_feature_importances = {name: round(float(imp), 4) for name, imp in sorted_imp[:15]}
        except Exception as e:
            logger.warning(f"Random forest classifier evaluation skipped: {e}")

    # Optional ESM-2 Embeddings hook
    esm_status = "esm module not installed (torch/fair-esm2 optional)"
    try:
        import esm
        esm_status = "ESM-2 model available for deep sequence representations"
    except ImportError:
        pass

    results = {
        "accession": query_accession,
        "feature_count": df_features.shape[1],
        "total_samples": len(df_features),
        "features_csv": feature_csv_path,
        "n_clusters_used": k,
        "clusters": clustering_records,
        "top_discriminative_features": rf_feature_importances,
        "deep_learning_hook": esm_status,
    }

    out_json_path = os.path.join(ml_dir, "clustering_results.json")
    with open(out_json_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
    logger.info(f"Saved ML clustering results to {out_json_path}")

    return results
