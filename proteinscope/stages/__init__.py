"""
ProteinScope Analysis Stages (Stage 1 to Stage 9 + Pathogenicity).
"""

from .ingest import ingest_uniprot_data
from .sequence import analyze_protein_sequence
from .blast import run_blast_search
from .structure import retrieve_protein_structures
from .validate import validate_protein_structure
from .visualize import visualize_protein_structure
from .comparative import run_comparative_analysis
from .ml import extract_and_cluster_ml_features
from .ngs import process_ngs_data
from .clinvar import load_clinvar_variants
from .variant_features import extract_variant_features, extract_feature_matrix
from .pathogenicity import (
    train_pathogenicity_classifier,
    predict_variant,
    generate_justification,
)
from .stability import (
    analyze_protein_stability,
    calculate_net_charge,
    estimate_melting_temperature,
    evaluate_conformation_and_degradation,
)

__all__ = [
    "ingest_uniprot_data",
    "analyze_protein_sequence",
    "run_blast_search",
    "retrieve_protein_structures",
    "validate_protein_structure",
    "visualize_protein_structure",
    "analyze_protein_stability",
    "calculate_net_charge",
    "estimate_melting_temperature",
    "evaluate_conformation_and_degradation",
    "run_comparative_analysis",
    "extract_and_cluster_ml_features",
    "process_ngs_data",
    "load_clinvar_variants",
    "extract_variant_features",
    "extract_feature_matrix",
    "train_pathogenicity_classifier",
    "predict_variant",
    "generate_justification",
]

