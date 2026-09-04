"""
Unit tests for Stage 2: Sequence Analysis.
"""

import os
import shutil
import tempfile
import unittest

from proteinscope.stages.sequence import (
    analyze_protein_sequence,
    compute_physicochemical_properties,
    search_motifs,
)


class TestSequence(unittest.TestCase):

    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.sample_seq = "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFF"

    def tearDown(self):
        shutil.rmtree(self.test_dir, ignore_errors=True)

    def test_compute_physicochemical_properties(self):
        props = compute_physicochemical_properties(self.sample_seq)
        self.assertIn("molecular_weight", props)
        self.assertIn("isoelectric_point", props)
        self.assertIn("gravy", props)
        self.assertIn("secondary_structure_fraction", props)
        self.assertIn("amino_acid_percentages", props)
        
        self.assertGreater(props["molecular_weight"], 1000)
        self.assertGreaterEqual(props["isoelectric_point"], 1.0)
        self.assertLessEqual(props["isoelectric_point"], 14.0)
        
        # Check individual secondary structure fractions are valid probabilities in [0.0, 1.0]
        sec = props["secondary_structure_fraction"]
        self.assertGreaterEqual(sec["helix"], 0.0)
        self.assertLessEqual(sec["helix"], 1.0)
        self.assertGreaterEqual(sec["turn"], 0.0)
        self.assertLessEqual(sec["turn"], 1.0)
        self.assertGreaterEqual(sec["sheet"], 0.0)
        self.assertLessEqual(sec["sheet"], 1.0)

    def test_search_motifs(self):
        # Sequence containing N-glycosylation (NLT), PKC phospho (STQDLFLPFF -> S..E or T..E)
        test_seq = "MANLTTRTQSE"
        motifs = search_motifs(test_seq)
        
        self.assertIn("n_glycosylation", motifs)
        self.assertIn("signal_peptide", motifs)
        self.assertIn("pkc_phosphorylation", motifs)
        
        # NLT matches N[^P][ST]
        self.assertTrue(any(m["matched_sequence"] == "NLT" for m in motifs["n_glycosylation"]))

    def test_analyze_protein_sequence_file(self):
        # Save sample fasta
        seq_dir = os.path.join(self.test_dir, "sequences")
        os.makedirs(seq_dir, exist_ok=True)
        fasta_file = os.path.join(seq_dir, "TEST01.fasta")
        with open(fasta_file, "w") as f:
            f.write(f">sp|TEST01|Sample protein\n{self.sample_seq}\n")

        res = analyze_protein_sequence("TEST01", self.test_dir, fasta_path=fasta_file)
        self.assertEqual(res["accession"], "TEST01")
        self.assertEqual(res["sequence_length"], len(self.sample_seq))
        self.assertIn("physicochemical_properties", res)


if __name__ == "__main__":
    unittest.main()
