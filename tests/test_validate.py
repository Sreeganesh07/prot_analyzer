"""
Unit tests for Stage 5: Structure Validation & Stage 8: ML Feature Engineering.
"""

import os
import shutil
import tempfile
import unittest
import numpy as np

from proteinscope.stages.validate import (
    calc_dihedral_angle,
    classify_ramachandran_region,
)
from proteinscope.stages.ml import (
    extract_sequence_feature_vector,
    extract_and_cluster_ml_features,
)


class TestValidateAndML(unittest.TestCase):

    def setUp(self):
        self.test_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.test_dir, ignore_errors=True)

    def test_classify_ramachandran_region(self):
        # Alpha helix region
        self.assertEqual(classify_ramachandran_region(-60.0, -45.0, "ALA"), "favored")
        # Beta sheet region
        self.assertEqual(classify_ramachandran_region(-120.0, 135.0, "VAL"), "favored")
        # Outlier
        self.assertEqual(classify_ramachandran_region(100.0, -100.0, "LEU"), "outlier")
        # Glycine always favored
        self.assertEqual(classify_ramachandran_region(100.0, -100.0, "GLY"), "favored")

    def test_calc_dihedral_angle(self):
        # Planar cis conformation: (0,1,0) -> (0,0,0) -> (1,0,0) -> (1,1,0)
        p1 = np.array([0.0, 1.0, 0.0])
        p2 = np.array([0.0, 0.0, 0.0])
        p3 = np.array([1.0, 0.0, 0.0])
        p4 = np.array([1.0, 1.0, 0.0])
        angle_cis = calc_dihedral_angle(p1, p2, p3, p4)
        self.assertAlmostEqual(abs(angle_cis), 0.0, places=1)

        # Planar trans conformation: (0,1,0) -> (0,0,0) -> (1,0,0) -> (1,-1,0)
        p4_trans = np.array([1.0, -1.0, 0.0])
        angle_trans = calc_dihedral_angle(p1, p2, p3, p4_trans)
        self.assertAlmostEqual(abs(angle_trans), 180.0, places=1)

    def test_extract_sequence_feature_vector(self):
        seq = "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHS"
        feats = extract_sequence_feature_vector(seq)
        self.assertGreaterEqual(len(feats), 40)
        self.assertIn("molecular_weight", feats)
        self.assertIn("isoelectric_point", feats)
        self.assertIn("hydrophobic_pct", feats)
        self.assertIn("aa_A_pct", feats)
        self.assertIn("aa_Y_pct", feats)

    def test_extract_and_cluster_ml_features(self):
        seq = "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHS"
        res = extract_and_cluster_ml_features(
            query_accession="TEST01",
            query_sequence=seq,
            output_dir=self.test_dir,
            n_clusters=2
        )
        self.assertIn("feature_count", res)
        self.assertIn("clusters", res)
        self.assertTrue(os.path.exists(res["features_csv"]))


if __name__ == "__main__":
    unittest.main()
