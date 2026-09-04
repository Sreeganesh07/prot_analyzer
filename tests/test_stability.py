"""
Unit tests for Stage 5b: Environmental Stability & Degradation (Alberts NBK26830).
"""

import json
import os
import shutil
import tempfile
import unittest

from proteinscope.stages.stability import (
    calculate_net_charge,
    generate_charge_titration_curve,
    calculate_salt_bridge_retention,
    estimate_melting_temperature,
    compute_two_state_unfolding,
    evaluate_conformation_and_degradation,
    analyze_protein_stability,
)


class TestStability(unittest.TestCase):

    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        # Typical protein sequence with mix of acidic, basic, polar, and hydrophobic residues
        self.sample_seq = "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFF"

    def tearDown(self):
        shutil.rmtree(self.test_dir, ignore_errors=True)

    def test_calculate_net_charge(self):
        # Extreme acid (pH 1.0) -> Asp and Glu protonated (uncharged), Lys/Arg/His/Nterm charged (+)
        q_acid = calculate_net_charge(self.sample_seq, ph=1.0)
        self.assertGreater(q_acid, 0.0)

        # Extreme base (pH 13.0) -> Lys and Arg deprotonated, Asp/Glu/Cys/Tyr/Cterm negative (-)
        q_base = calculate_net_charge(self.sample_seq, ph=13.0)
        self.assertLess(q_base, 0.0)

        # Physiological neutral (pH 7.4)
        q_phys = calculate_net_charge(self.sample_seq, ph=7.4)
        self.assertIsInstance(q_phys, float)

    def test_generate_charge_titration_curve(self):
        curve = generate_charge_titration_curve(self.sample_seq, step=1.0)
        self.assertGreaterEqual(len(curve), 14)
        self.assertEqual(curve[0]["ph"], 0.0)
        self.assertEqual(curve[-1]["ph"], 14.0)
        # Monotonically non-increasing net charge as pH increases
        for i in range(len(curve) - 1):
            self.assertGreaterEqual(curve[i]["net_charge"], curve[i+1]["net_charge"])

    def test_calculate_salt_bridge_retention(self):
        # Neutral pH should have optimal retention (> 90%)
        salt_neutral = calculate_salt_bridge_retention(7.4)
        self.assertGreaterEqual(salt_neutral, 90.0)

        # Acidic pH (pH 2.0) should dismantle salt bridges (< 20%)
        salt_acid = calculate_salt_bridge_retention(2.0)
        self.assertLess(salt_acid, 20.0)

        # Alkaline pH (pH 12.0) should dismantle salt bridges (< 20%)
        salt_base = calculate_salt_bridge_retention(12.0)
        self.assertLess(salt_base, 20.0)

    def test_estimate_melting_temperature(self):
        tm_without_ss = estimate_melting_temperature(self.sample_seq, ptm_sites=[])
        self.assertGreaterEqual(tm_without_ss, 45.0)
        self.assertLessEqual(tm_without_ss, 88.0)

        # Disulfide bonds should increase thermal stability
        fake_ptms = [{"type": "Disulfide bond", "position": "10-25"}, {"type": "Disulfide bond", "position": "30-50"}]
        tm_with_ss = estimate_melting_temperature(self.sample_seq, ptm_sites=fake_ptms)
        self.assertGreater(tm_with_ss, tm_without_ss)

    def test_compute_two_state_unfolding(self):
        tm = 60.0
        # Standard physiological temperature (37°C) -> majority folded
        thermo_norm = compute_two_state_unfolding(37.0, ph=7.4, tm_celsius=tm, net_charge=0.0)
        self.assertGreater(thermo_norm["fraction_folded"], 0.80)
        self.assertLess(thermo_norm["delta_g_folding_kcal_mol"], 0.0)

        # High temperature (85°C) -> thermally denatured / unfolded
        thermo_hot = compute_two_state_unfolding(85.0, ph=7.4, tm_celsius=tm, net_charge=0.0)
        self.assertLess(thermo_hot["fraction_folded"], 0.15)
        self.assertGreater(thermo_hot["delta_g_folding_kcal_mol"], 0.0)

    def test_evaluate_conformation_and_degradation(self):
        tm = 60.0

        # Standard physiological conditions -> Native
        res_phys = evaluate_conformation_and_degradation(37.0, 7.4, tm_celsius=tm, net_charge=0.0)
        self.assertEqual(res_phys["state"], "NATIVE")
        self.assertFalse(res_phys["is_degraded"])

        # Heat degradation (> 70°C) -> Degraded (aggregation)
        res_hot = evaluate_conformation_and_degradation(85.0, 7.4, tm_celsius=tm, net_charge=0.0)
        self.assertEqual(res_hot["state"], "DEGRADED")
        self.assertTrue(res_hot["is_degraded"])
        self.assertTrue(any("Thermal" in r for r in res_hot["degradation_reasons"]))

        # Acid degradation (pH 2.0) -> Degraded (hydrolysis)
        res_acid = evaluate_conformation_and_degradation(37.0, 2.0, tm_celsius=tm, net_charge=15.0)
        self.assertEqual(res_acid["state"], "DEGRADED")
        self.assertTrue(res_acid["is_degraded"])
        self.assertTrue(any("Acid" in r for r in res_acid["degradation_reasons"]))

        # Alkaline degradation (pH 12.5) -> Degraded
        res_alk = evaluate_conformation_and_degradation(37.0, 12.5, tm_celsius=tm, net_charge=-15.0)
        self.assertEqual(res_alk["state"], "DEGRADED")
        self.assertTrue(res_alk["is_degraded"])
        self.assertTrue(any("Alkaline" in r for r in res_alk["degradation_reasons"]))

    def test_custom_range_thresholds(self):
        tm = 60.0
        # Custom narrow stability range [20°C, 30°C]
        res = evaluate_conformation_and_degradation(
            temperature_celsius=35.0,
            ph=7.4,
            tm_celsius=tm,
            net_charge=0.0,
            custom_temp_range=(20.0, 30.0),
        )
        # 35°C is outside the custom range [20, 30] so it should be PERTURBED
        self.assertIn(res["state"], ["PERTURBED", "DENATURED"])

    def test_analyze_protein_stability_manifest(self):
        res = analyze_protein_stability(
            sequence=self.sample_seq,
            output_dir=self.test_dir,
            user_temperature=37.0,
            user_ph=7.4,
        )
        self.assertIn("estimated_melting_temperature_celsius", res)
        self.assertIn("standard_conditions", res)
        self.assertIn("physiological_37c_ph74", res["standard_conditions"])
        self.assertIn("invitro_25c_ph70", res["standard_conditions"])
        self.assertIn("active_evaluation", res)
        self.assertIn("titration_profile", res)

        # Check JSON manifest was saved
        manifest_path = os.path.join(self.test_dir, "stability", "stability_profile.json")
        self.assertTrue(os.path.exists(manifest_path))
        with open(manifest_path, "r", encoding="utf-8") as f:
            data = json.load(f)
            self.assertEqual(data["sequence_length"], len(self.sample_seq))


if __name__ == "__main__":
    unittest.main()
