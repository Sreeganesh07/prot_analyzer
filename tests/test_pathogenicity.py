"""
Unit tests for the new pathogenicity prediction modules:
  - ClinVar data loading (clinvar.py)
  - Variant feature extraction (variant_features.py)
  - Pathogenicity classifier (pathogenicity.py)
"""

import os
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from proteinscope.stages.clinvar import (
    parse_hgvs_missense,
    classify_clinical_significance,
)
from proteinscope.stages.variant_features import (
    get_grantham_distance,
    get_blosum62_score,
    get_hydrophobicity_delta,
    get_charge_change,
    get_size_delta,
    get_polarity_change,
    extract_ptm_features,
    extract_conservation_features,
    build_position_conservation,
    extract_variant_features,
)
from proteinscope.stages.pathogenicity import (
    train_pathogenicity_classifier,
    predict_variant,
    generate_justification,
)


class TestClinVar(unittest.TestCase):
    """Tests for ClinVar HGVS parsing and label classification."""

    def test_parse_hgvs_3letter(self):
        """p.Arg175His should parse to (R, 175, H)."""
        result = parse_hgvs_missense("p.Arg175His")
        self.assertIsNotNone(result)
        self.assertEqual(result, ("R", 175, "H"))

    def test_parse_hgvs_1letter(self):
        """p.R175H should parse to (R, 175, H)."""
        result = parse_hgvs_missense("p.R175H")
        self.assertIsNotNone(result)
        self.assertEqual(result, ("R", 175, "H"))

    def test_parse_hgvs_common_mutations(self):
        """Test a variety of well-known TP53 mutations."""
        cases = [
            ("p.Arg248Trp", ("R", 248, "W")),
            ("p.Arg273His", ("R", 273, "H")),
            ("p.Gly245Ser", ("G", 245, "S")),
            ("p.Pro72Arg", ("P", 72, "R")),
        ]
        for hgvs, expected in cases:
            with self.subTest(hgvs=hgvs):
                result = parse_hgvs_missense(hgvs)
                self.assertEqual(result, expected)

    def test_parse_hgvs_synonymous_returns_none(self):
        """Synonymous (silent) changes like p.R175R should return None."""
        result = parse_hgvs_missense("p.Arg175Arg")
        self.assertIsNone(result)

    def test_parse_hgvs_invalid_returns_none(self):
        """Invalid or empty HGVS should return None."""
        self.assertIsNone(parse_hgvs_missense(""))
        self.assertIsNone(parse_hgvs_missense(None))
        self.assertIsNone(parse_hgvs_missense("some random text"))
        self.assertIsNone(parse_hgvs_missense("p.Arg175del"))  # deletion, not missense

    def test_classify_pathogenic(self):
        self.assertEqual(classify_clinical_significance("Pathogenic"), "pathogenic")
        self.assertEqual(classify_clinical_significance("Likely pathogenic"), "pathogenic")
        self.assertEqual(classify_clinical_significance("Pathogenic/Likely pathogenic"), "pathogenic")

    def test_classify_benign(self):
        self.assertEqual(classify_clinical_significance("Benign"), "benign")
        self.assertEqual(classify_clinical_significance("Likely benign"), "benign")
        self.assertEqual(classify_clinical_significance("Benign/Likely benign"), "benign")

    def test_classify_vus_returns_none(self):
        self.assertIsNone(classify_clinical_significance("Uncertain significance"))
        self.assertIsNone(classify_clinical_significance("Conflicting classifications of pathogenicity"))
        self.assertIsNone(classify_clinical_significance(""))
        self.assertIsNone(classify_clinical_significance(None))


class TestVariantFeatures(unittest.TestCase):
    """Tests for physicochemical and structural feature extraction."""

    def test_grantham_distance_known_pair(self):
        """R→H is a small change (Grantham=29)."""
        dist = get_grantham_distance("R", "H")
        self.assertEqual(dist, 29.0)

    def test_grantham_distance_identical(self):
        """Same amino acid should have distance 0."""
        self.assertEqual(get_grantham_distance("R", "R"), 0.0)

    def test_grantham_distance_symmetric(self):
        """Grantham distance should be symmetric."""
        self.assertEqual(get_grantham_distance("R", "H"), get_grantham_distance("H", "R"))

    def test_grantham_radical_substitution(self):
        """Cys→Trp is a radical substitution (Grantham=215)."""
        dist = get_grantham_distance("C", "W")
        self.assertEqual(dist, 215.0)

    def test_blosum62_score_known(self):
        """R→H has BLOSUM62 score of 0."""
        self.assertEqual(get_blosum62_score("R", "H"), 0.0)

    def test_blosum62_identity_positive(self):
        """Same amino acid should have positive BLOSUM62 score."""
        self.assertGreater(get_blosum62_score("R", "R"), 0.0)

    def test_hydrophobicity_delta(self):
        """R is hydrophilic (-4.5), H is also hydrophilic (-3.2)."""
        delta = get_hydrophobicity_delta("R", "H")
        self.assertAlmostEqual(delta, -3.2 - (-4.5), places=1)

    def test_charge_change_r_to_h(self):
        """R is +1, H is +0.1, so charge loss."""
        change = get_charge_change("R", "H")
        self.assertLess(change, 0)

    def test_polarity_change_polar_to_nonpolar(self):
        """R (polar) → A (nonpolar) should flag polarity change."""
        self.assertEqual(get_polarity_change("R", "A"), 1)

    def test_polarity_change_same_class(self):
        """R (polar) → K (polar) should not flag polarity change."""
        self.assertEqual(get_polarity_change("R", "K"), 0)

    def test_ptm_features_at_ptm_site(self):
        ptms = [{"position": "175", "type": "phosphorylation"}]
        feats = extract_ptm_features(175, ptms)
        self.assertEqual(feats["is_ptm_site"], 1.0)
        self.assertEqual(feats["min_ptm_distance"], 0.0)

    def test_ptm_features_near_ptm(self):
        ptms = [{"position": "178", "type": "phosphorylation"}]
        feats = extract_ptm_features(175, ptms)
        self.assertEqual(feats["is_ptm_site"], 0.0)
        self.assertEqual(feats["min_ptm_distance"], 3.0)
        self.assertEqual(feats["n_ptms_within_5"], 1.0)

    def test_ptm_features_no_ptms(self):
        feats = extract_ptm_features(175, [])
        self.assertEqual(feats["min_ptm_distance"], 999.0)
        self.assertEqual(feats["is_ptm_site"], 0.0)

    def test_conservation_profile_basic(self):
        """Build conservation from simple aligned hits."""
        query = "MEEPQSDPSVEPPLSQE"  # 17 residues: Pos 6 is 'S'
        hits = [
            {"query_sequence": "MEEPQSDPSVEPPLSQE", "sbjct_sequence": "MEEPQSDPSVEPPLSQE"},
            {"query_sequence": "MEEPQSDPSVEPPLSQE", "sbjct_sequence": "MEEPQADPSVEPPLSQE"},
        ]
        profile = build_position_conservation(query, hits)
        # Position 6 (S in query): first hit has S, second has A -> conservation = 0.5
        self.assertEqual(profile[6]["conservation_score"], 0.5)
        # Position 1 (M in query): both hits have M -> conservation = 1.0
        self.assertEqual(profile[1]["conservation_score"], 1.0)

    def test_extract_variant_features_basic(self):
        """Extract features for a basic variant without PDB or BLAST data."""
        features = extract_variant_features(
            position=175,
            ref_aa="R",
            alt_aa="H",
            query_sequence="M" * 174 + "R" + "E" * 218,  # 393 residues with R at 175
        )
        self.assertIn("grantham_distance", features)
        self.assertIn("blosum62_score", features)
        self.assertIn("hydrophobicity_delta", features)
        self.assertIn("charge_change", features)
        self.assertIn("relative_position", features)
        self.assertGreater(len(features), 10)  # Should have many features


class TestPathogenicityClassifier(unittest.TestCase):
    """Tests for the classifier training and prediction pipeline."""

    def _make_synthetic_data(self, n=100):
        """Generate synthetic feature data for testing."""
        np.random.seed(42)
        feature_dicts = []
        labels = []
        for i in range(n):
            is_pathogenic = i < n // 2
            features = {
                "grantham_distance": np.random.uniform(100, 200) if is_pathogenic else np.random.uniform(10, 60),
                "blosum62_score": np.random.uniform(-4, -1) if is_pathogenic else np.random.uniform(-1, 4),
                "burial_class": np.random.choice([1.0, 2.0]) if is_pathogenic else np.random.choice([0.0, 1.0]),
                "conservation_score": np.random.uniform(0.7, 1.0) if is_pathogenic else np.random.uniform(0.0, 0.5),
                "polarity_change": float(np.random.choice([0, 1], p=[0.3, 0.7])) if is_pathogenic else float(np.random.choice([0, 1], p=[0.8, 0.2])),
                "charge_change": np.random.uniform(-1, 1),
                "relative_position": np.random.uniform(0.1, 0.9),
            }
            feature_dicts.append(features)
            labels.append("pathogenic" if is_pathogenic else "benign")
        return feature_dicts, labels

    def test_train_classifier(self):
        """Classifier trains and returns valid metrics on synthetic data."""
        feature_dicts, labels = self._make_synthetic_data(100)

        with tempfile.TemporaryDirectory() as tmpdir:
            metrics = train_pathogenicity_classifier(
                feature_dicts=feature_dicts,
                labels=labels,
                model_dir=tmpdir,
                n_estimators=20,
                test_size=0.3,
                random_state=42,
            )

        self.assertIn("accuracy", metrics)
        self.assertIn("auroc", metrics)
        self.assertIn("confusion_matrix", metrics)
        self.assertIn("top_features", metrics)
        self.assertGreater(metrics["accuracy"], 0.5)  # Should be better than random
        self.assertGreater(metrics["auroc"], 0.5)

    def test_predict_variant_with_model(self):
        """Prediction returns valid structure after training."""
        feature_dicts, labels = self._make_synthetic_data(80)

        with tempfile.TemporaryDirectory() as tmpdir:
            metrics = train_pathogenicity_classifier(
                feature_dicts=feature_dicts,
                labels=labels,
                model_dir=tmpdir,
                n_estimators=20,
            )

            model_path = metrics["model_path"]

            # Predict on a clearly pathogenic-looking variant
            pathogenic_features = {
                "grantham_distance": 180.0,
                "blosum62_score": -3.0,
                "burial_class": 2.0,
                "conservation_score": 0.95,
                "polarity_change": 1.0,
                "charge_change": -1.0,
                "relative_position": 0.5,
            }

            result = predict_variant(
                feature_dict=pathogenic_features,
                model_path=model_path,
            )

            self.assertIn("prediction", result)
            self.assertIn("confidence", result)
            self.assertIn(result["prediction"], ["pathogenic", "benign"])
            self.assertGreater(result["confidence"], 0.0)
            self.assertLessEqual(result["confidence"], 1.0)

    def test_generate_justification(self):
        """Justification text includes key feature descriptions."""
        features = {
            "burial_class": 2.0,
            "burial_neighbour_count": 22.0,
            "grantham_distance": 29.0,
            "charge_change": -0.9,
            "conservation_score": 0.9,
            "homolog_diversity": 2.0,
            "n_homologs_aligned": 10.0,
            "min_ptm_distance": 3.0,
            "is_ptm_site": 0.0,
            "ramachandran_zone": 0.0,
        }
        prediction = {"prediction": "pathogenic", "confidence": 0.94}

        text = generate_justification(features, prediction, "R", "H", 175)

        self.assertIn("BURIED", text)
        self.assertIn("175", text)
        self.assertIn("CHARGE", text)
        self.assertIn("conservation", text.lower())


if __name__ == "__main__":
    unittest.main()
