"""
Stage 5b: Environmental Temperature & pH Stability and Denaturation / Degradation Modeling.
Based on Alberts et al., Molecular Biology of the Cell (4th ed.), Chapter 3:
"The Shape and Structure of Proteins" (NCBI Bookshelf NBK26830).

Provides biophysical modeling for:
- Standard physiological (37°C, pH 7.4) and in-vitro (25°C, pH 7.0) baseline stability.
- Net molecular charge titration via Henderson-Hasselbalch across pH 0-14.
- Sequence- and disulfide-dependent thermal melting temperature (Tm) prediction.
- Salt bridge electrostatic retention as a function of pH.
- Two-state thermodynamic folding stability: Delta G(T, pH) and fraction folded f_folded.
- Environmental range thresholds for structural perturbation, thermal/pH denaturation,
  and irreversible degradation (insoluble aggregation and hydrolytic breakdown).
"""

import json
import logging
import math
import os
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger("ProteinScope.Stability")

# Standard side-chain and terminal pKa values (Lehninger / Nelson & Cox scale)
PKA_VALUES: Dict[str, float] = {
    "C_term": 2.34,
    "D": 3.65,      # Aspartic acid
    "E": 4.25,      # Glutamic acid
    "H": 6.00,      # Histidine
    "C": 8.18,      # Cysteine
    "N_term": 9.69,
    "Y": 10.07,     # Tyrosine
    "K": 10.53,     # Lysine
    "R": 12.48,     # Arginine
}

# Standard Reference Conditions
STANDARD_CONDITIONS: Dict[str, Dict[str, Any]] = {
    "physiological": {
        "name": "Physiological Standard",
        "temperature_celsius": 37.0,
        "temperature_kelvin": 310.15,
        "ph": 7.4,
        "description": "Human physiological core temperature and blood/cytosolic pH",
    },
    "invitro": {
        "name": "Laboratory Standard (In-Vitro)",
        "temperature_celsius": 25.0,
        "temperature_kelvin": 298.15,
        "ph": 7.0,
        "description": "Standard room temperature and neutral experimental buffer",
    },
    "hypothermic": {
        "name": "Cold Storage / Hypothermic",
        "temperature_celsius": 4.0,
        "temperature_kelvin": 277.15,
        "ph": 7.4,
        "description": "Refrigerated / cold storage physiological buffer",
    },
}


def calculate_net_charge(sequence: str, ph: float, pka_table: Optional[Dict[str, float]] = None) -> float:
    """
    Calculate the net electrostatic charge of a protein sequence at a specific pH
    using the Henderson-Hasselbalch equation and side-chain pKa values.

    Positively charged groups: N-terminal, Lys, Arg, His
    Negatively charged groups: C-terminal, Asp, Glu, Cys, Tyr
    """
    if not sequence:
        return 0.0

    pka = pka_table or PKA_VALUES
    seq_upper = sequence.upper()

    # Counts of titratable residues
    counts = {
        "D": seq_upper.count("D"),
        "E": seq_upper.count("E"),
        "H": seq_upper.count("H"),
        "C": seq_upper.count("C"),
        "Y": seq_upper.count("Y"),
        "K": seq_upper.count("K"),
        "R": seq_upper.count("R"),
    }

    # N-terminal (+) and C-terminal (-)
    q_nterm = 1.0 / (1.0 + 10.0 ** (ph - pka["N_term"]))
    q_cterm = -1.0 / (1.0 + 10.0 ** (pka["C_term"] - ph))

    # Basic residues (+)
    q_k = counts["K"] / (1.0 + 10.0 ** (ph - pka["K"]))
    q_r = counts["R"] / (1.0 + 10.0 ** (ph - pka["R"]))
    q_h = counts["H"] / (1.0 + 10.0 ** (ph - pka["H"]))

    # Acidic residues (-)
    q_d = -counts["D"] / (1.0 + 10.0 ** (pka["D"] - ph))
    q_e = -counts["E"] / (1.0 + 10.0 ** (pka["E"] - ph))
    q_c = -counts["C"] / (1.0 + 10.0 ** (pka["C"] - ph))
    q_y = -counts["Y"] / (1.0 + 10.0 ** (pka["Y"] - ph))

    net_charge = q_nterm + q_cterm + q_k + q_r + q_h + q_d + q_e + q_c + q_y
    return round(net_charge, 2)


def generate_charge_titration_curve(sequence: str, step: float = 0.5) -> List[Dict[str, float]]:
    """
    Generate the full charge titration profile Q(pH) from pH 0.0 to 14.0.
    """
    profile: List[Dict[str, float]] = []
    current_ph = 0.0
    while current_ph <= 14.001:
        q = calculate_net_charge(sequence, current_ph)
        profile.append({"ph": round(current_ph, 1), "net_charge": q})
        current_ph += step
    return profile


def calculate_salt_bridge_retention(ph: float) -> float:
    """
    Compute fraction of intact electrostatic salt bridges (Asp-/Glu- with Lys+/Arg+/His+).
    At neutral pH (6.0 - 8.0), retention is near 100%. At extreme acid or base,
    neutralization of carboxylates or amines dismantles ionic bonds.
    """
    # Acidic deprotonation fraction (Asp pKa ~3.85, Glu ~4.25)
    f_acid_ionized = 1.0 / (1.0 + 10.0 ** (4.0 - ph))
    # Basic protonation fraction (Lys pKa ~10.5, Arg ~12.5)
    f_base_ionized = 1.0 / (1.0 + 10.0 ** (ph - 10.5))

    retention = f_acid_ionized * f_base_ionized
    return round(max(0.0, min(1.0, retention)) * 100.0, 1)


def estimate_melting_temperature(
    sequence: str,
    properties: Optional[Dict[str, Any]] = None,
    ptm_sites: Optional[List[Dict[str, Any]]] = None,
) -> float:
    """
    Estimate the sequence-specific melting temperature Tm (°C) using empirical biophysical correlations:
    - Base globular melting temperature: ~58.0°C
    - Aliphatic index (Ala, Val, Ile, Leu hydrophobic packing)
    - Proline rigidity in loops
    - Disulfide bonds (each disulfide covalent cross-link provides +3.0°C to +4.5°C stabilization)
    - ProtParam instability index penalty (instability index > 40 reduces Tm)
    """
    if not sequence:
        return 58.0

    seq_len = len(sequence)
    seq_u = sequence.upper()

    # Disulfide bond count from PTM annotations or Cys pairs
    disulfide_count = 0
    if ptm_sites:
        disulfide_count = sum(1 for p in ptm_sites if "disulfide" in str(p.get("type", "")).lower())
    if disulfide_count == 0:
        cys_count = seq_u.count("C")
        disulfide_count = cys_count // 2

    # Aliphatic index: X(Ala) + 2.9*X(Val) + 3.9*(X(Ile) + X(Leu))
    x_ala = (seq_u.count("A") / seq_len) * 100.0
    x_val = (seq_u.count("V") / seq_len) * 100.0
    x_ile = (seq_u.count("I") / seq_len) * 100.0
    x_leu = (seq_u.count("L") / seq_len) * 100.0
    aliphatic_index = x_ala + 2.9 * x_val + 3.9 * (x_ile + x_leu)

    # Proline fraction (backbone conformational entropy restriction)
    pro_pct = (seq_u.count("P") / seq_len) * 100.0

    # Instability index adjustment
    instab = 35.0
    if properties and "instability_index" in properties:
        try:
            instab = float(properties["instability_index"])
        except (ValueError, TypeError):
            pass

    # Baseline globular Tm
    tm_base = 56.5
    # Aliphatic hydrophobic effect
    aliph_contrib = (aliphatic_index - 75.0) * 0.12
    # Disulfide stabilization (+3.5°C per bond, capped at +18°C)
    ss_contrib = min(18.0, disulfide_count * 3.5)
    # Proline entropy restriction
    pro_contrib = (pro_pct - 5.0) * 0.35
    # Instability index penalty
    instab_penalty = (instab - 38.0) * 0.15

    tm = tm_base + aliph_contrib + ss_contrib + pro_contrib - instab_penalty
    # Bound to physically realistic protein denaturation ranges (42°C - 88°C)
    tm = max(45.0, min(88.0, tm))
    return round(tm, 1)


def compute_two_state_unfolding(
    temperature_celsius: float,
    ph: float,
    tm_celsius: float,
    net_charge: float,
    seq_len: int = 300,
) -> Dict[str, Any]:
    """
    Two-state Gibbs-Helmholtz thermodynamic unfolding model:
    Delta G_folding(T, pH) = Delta H_m * (1 - T/Tm) - Delta Cp * [(Tm - T) + T * ln(T/Tm)] - Delta G_charge(pH)
    Fraction folded: f_folded = 1 / (1 + exp(Delta G / RT)) (in folded convention, Delta G < 0)
    """
    t_kelvin = temperature_celsius + 273.15
    tm_kelvin = tm_celsius + 273.15
    r_const = 1.9872e-3  # kcal / (mol * K)

    # Enthalpy of unfolding at Tm (~1.2 kcal/mol per residue average)
    delta_h_m = 1.15 * max(50, min(600, seq_len))  # kcal / mol
    # Heat capacity change Delta Cp (~14 cal/mol/K per residue)
    delta_cp = 0.012 * max(50, min(600, seq_len))  # kcal / (mol * K)

    # Thermal Gibbs-Helmholtz term
    term1 = delta_h_m * (1.0 - (t_kelvin / tm_kelvin))
    term2 = delta_cp * ((tm_kelvin - t_kelvin) + t_kelvin * math.log(t_kelvin / tm_kelvin))
    delta_g_thermal = -(term1 - term2)  # delta_g_folding (negative when stably folded)

    # Electrostatic destabilization from extreme net charge repulsion
    # High net charge (|Q| >> 15) forces unfolding through mutual backbone repulsion
    charge_destabilization = 0.008 * (abs(net_charge) ** 1.8)
    delta_g_total = delta_g_thermal + charge_destabilization

    # Fraction folded f_folded
    # K_eq = [U] / [F] = exp(-Delta G_unfolding / RT) = exp(Delta G_folding / RT)
    exponent = max(-30.0, min(30.0, delta_g_total / (r_const * t_kelvin)))
    f_folded = 1.0 / (1.0 + math.exp(exponent))

    return {
        "temperature_celsius": temperature_celsius,
        "temperature_kelvin": round(t_kelvin, 2),
        "ph": ph,
        "delta_g_folding_kcal_mol": round(delta_g_total, 2),
        "fraction_folded": round(f_folded, 4),
        "fraction_folded_percent": round(f_folded * 100.0, 1),
    }


def evaluate_conformation_and_degradation(
    temperature_celsius: float,
    ph: float,
    tm_celsius: float,
    net_charge: float,
    custom_temp_range: Optional[Tuple[float, float]] = None,
    custom_ph_range: Optional[Tuple[float, float]] = None,
) -> Dict[str, Any]:
    """
    Classify the conformational and degradation state of the protein based on
    Alberts NBK26830 biophysical principles.

    States:
    1. NATIVE:
       Optimal folded state, active conformation, lowest free energy, intact noncovalent bonds.
    2. PERTURBED / MOLTEN GLOBULE:
       Partial loosening of loops, increased thermal fluctuations, pre-melting intermediate.
    3. DENATURED (REVERSIBLE / PARTIAL):
       Loss of tertiary and secondary structure, unfolded random coil, exposed hydrophobic core.
    4. DEGRADED (IRREVERSIBLE):
       Insoluble aggregation / amyloid-like precipitation (heat-induced degradation)
       or acid/base hydrolytic cleavage (extreme pH degradation).
    """
    # Thresholds: either custom or biophysically derived
    if custom_temp_range:
        temp_min_stable, temp_max_stable = custom_temp_range
    else:
        temp_min_stable = 15.0
        temp_max_stable = round(tm_celsius - 8.0, 1)

    if custom_ph_range:
        ph_min_stable, ph_max_stable = custom_ph_range
    else:
        ph_min_stable = 5.5
        ph_max_stable = 8.5

    # Critical irreversible degradation boundaries
    temp_degradation_threshold = round(max(70.0, tm_celsius + 10.0), 1)
    ph_acid_degradation_threshold = 3.0
    ph_alkaline_degradation_threshold = 11.5

    thermo = compute_two_state_unfolding(temperature_celsius, ph, tm_celsius, net_charge)
    f_folded = thermo["fraction_folded"]
    salt_bridges = calculate_salt_bridge_retention(ph)

    # Evaluate Degradation Triggers
    is_thermally_degraded = temperature_celsius >= temp_degradation_threshold
    is_cold_denatured = temperature_celsius < 4.0
    is_acid_degraded = ph <= ph_acid_degradation_threshold
    is_alkaline_degraded = ph >= ph_alkaline_degradation_threshold

    degradation_reasons = []
    if is_thermally_degraded:
        degradation_reasons.append(
            f"Thermal Aggregation: Temperature ({temperature_celsius:.1f}°C) exceeds degradation threshold "
            f"({temp_degradation_threshold:.1f}°C). Prolonged thermal excitation causes exposed hydrophobic residues "
            f"to coalesce irreversibly into insoluble aggregates (Alberts NBK26830)."
        )
    if is_acid_degraded:
        degradation_reasons.append(
            f"Acid Degradation & Hydrolysis: Extreme acidic environment (pH {ph:.1f} ≤ {ph_acid_degradation_threshold}) "
            f"protonates all carboxylate groups (Asp/Glu), abolishing all ionic salt bridges with intense positive "
            f"charge repulsion (+{net_charge:.1f} e), catalyzing peptide bond hydrolysis."
        )
    if is_alkaline_degraded:
        degradation_reasons.append(
            f"Alkaline Degradation: Extreme basic environment (pH {ph:.1f} ≥ {ph_alkaline_degradation_threshold}) "
            f"deprotonates amino sidechains (Lys/Arg/Tyr/Cys), generating massive negative charge repulsion "
            f"({net_charge:.1f} e) and promoting irreversible β-elimination of cystine disulfide linkages."
        )

    # Conformation Classification
    if degradation_reasons:
        state = "DEGRADED"
        badge_label = "Degraded / Irreversible Aggregation"
        badge_color = "#f43f5e"  # Rose / red
        status_code = "degraded"
    elif f_folded < 0.40 or is_cold_denatured or ph < (ph_min_stable - 1.5) or ph > (ph_max_stable + 1.5):
        state = "DENATURED"
        badge_label = "Denatured (Unfolded Polypeptide)"
        badge_color = "#f97316"  # Orange
        status_code = "denatured"
    elif f_folded < 0.85 or temperature_celsius > temp_max_stable or ph < ph_min_stable or ph > ph_max_stable:
        state = "PERTURBED"
        badge_label = "Perturbed / Molten Globule"
        badge_color = "#f59e0b"  # Amber
        status_code = "perturbed"
    else:
        state = "NATIVE"
        badge_label = "Native Folded State"
        badge_color = "#10b981"  # Emerald green
        status_code = "native"

    # Molecular Mechanism Description (Alberts NBK26830)
    if state == "NATIVE":
        mechanism = (
            f"At physiological/standard conditions ({temperature_celsius:.1f}°C, pH {ph:.1f}), the protein "
            f"resides in its thermodynamic free energy minimum (ΔG = {thermo['delta_g_folding_kcal_mol']:.1f} kcal/mol). "
            f"Noncovalent bonds—including backbone hydrogen bonds (C=O···H-N), buried hydrophobic van der Waals contacts, "
            f"and electrostatic salt bridges ({salt_bridges}% intact)—fully stabilize the functional 3D conformation."
        )
    elif state == "PERTURBED":
        mechanism = (
            f"Under moderate thermal or pH stress ({temperature_celsius:.1f}°C, pH {ph:.1f}), thermal kinetic energy "
            f"loosens flexible surface loops and begins to weaken noncovalent hydrogen bonds. The protein adopts a "
            f"'molten globule' intermediate: compact with secondary structure intact, but with fluctuating tertiary packing."
        )
    elif state == "DENATURED":
        mechanism = (
            f"The destabilizing environmental conditions ({temperature_celsius:.1f}°C, pH {ph:.1f}) have disrupted "
            f"the delicate balance of weak noncovalent bonds. Hydrogen bonds and hydrophobic interactions have collapsed, "
            f"converting the ordered 3D architecture into a flexible, denatured random coil polypeptide chain (Fraction folded: {thermo['fraction_folded_percent']}%)."
        )
    else:
        mechanism = " ".join(degradation_reasons)

    return {
        "status_code": status_code,
        "state": state,
        "badge_label": badge_label,
        "badge_color": badge_color,
        "temperature_celsius": temperature_celsius,
        "ph": ph,
        "net_charge": net_charge,
        "salt_bridge_retention_percent": salt_bridges,
        "fraction_folded_percent": thermo["fraction_folded_percent"],
        "delta_g_folding_kcal_mol": thermo["delta_g_folding_kcal_mol"],
        "is_degraded": len(degradation_reasons) > 0,
        "degradation_reasons": degradation_reasons,
        "molecular_mechanism": mechanism,
        "thresholds": {
            "tm_celsius": tm_celsius,
            "temp_min_stable": temp_min_stable,
            "temp_max_stable": temp_max_stable,
            "temp_degradation": temp_degradation_threshold,
            "ph_min_stable": ph_min_stable,
            "ph_max_stable": ph_max_stable,
            "ph_acid_degradation": ph_acid_degradation_threshold,
            "ph_alkaline_degradation": ph_alkaline_degradation_threshold,
        },
    }


def analyze_protein_stability(
    sequence: str,
    output_dir: str,
    properties: Optional[Dict[str, Any]] = None,
    ptm_sites: Optional[List[Dict[str, Any]]] = None,
    user_temperature: Optional[float] = None,
    user_ph: Optional[float] = None,
    custom_temp_range: Optional[Tuple[float, float]] = None,
    custom_ph_range: Optional[Tuple[float, float]] = None,
) -> Dict[str, Any]:
    """
    Stage 5b Entry Point:
    Compute full biophysical environmental stability analysis, baseline standards,
    isoelectric titration curve, melting temperature, and state evaluations.
    Saves environmental stability manifest to JSON.
    """
    logger.info("Executing Stage 5b: Environmental Stability and Degradation Analysis")

    stability_dir = os.path.join(output_dir, "stability")
    os.makedirs(stability_dir, exist_ok=True)

    # Sequence properties
    tm_est = estimate_melting_temperature(sequence, properties=properties, ptm_sites=ptm_sites)
    titration_curve = generate_charge_titration_curve(sequence, step=0.5)

    # Standard Physiological Baseline (37°C, pH 7.4)
    q_phys = calculate_net_charge(sequence, 7.4)
    phys_state = evaluate_conformation_and_degradation(
        temperature_celsius=37.0,
        ph=7.4,
        tm_celsius=tm_est,
        net_charge=q_phys,
        custom_temp_range=custom_temp_range,
        custom_ph_range=custom_ph_range,
    )

    # Standard In-Vitro Baseline (25°C, pH 7.0)
    q_invitro = calculate_net_charge(sequence, 7.0)
    invitro_state = evaluate_conformation_and_degradation(
        temperature_celsius=25.0,
        ph=7.0,
        tm_celsius=tm_est,
        net_charge=q_invitro,
        custom_temp_range=custom_temp_range,
        custom_ph_range=custom_ph_range,
    )

    # Active / User-Selected Evaluation
    eval_temp = user_temperature if user_temperature is not None else 37.0
    eval_ph = user_ph if user_ph is not None else 7.4
    q_eval = calculate_net_charge(sequence, eval_ph)
    active_eval = evaluate_conformation_and_degradation(
        temperature_celsius=eval_temp,
        ph=eval_ph,
        tm_celsius=tm_est,
        net_charge=q_eval,
        custom_temp_range=custom_temp_range,
        custom_ph_range=custom_ph_range,
    )

    result_manifest: Dict[str, Any] = {
        "sequence_length": len(sequence),
        "estimated_melting_temperature_celsius": tm_est,
        "standard_conditions": {
            "physiological_37c_ph74": phys_state,
            "invitro_25c_ph70": invitro_state,
        },
        "active_evaluation": active_eval,
        "titration_profile": titration_curve,
        "biophysical_reference": "Alberts et al., Molecular Biology of the Cell (4th ed.), Chapter 3: NBK26830",
    }

    manifest_path = os.path.join(stability_dir, "stability_profile.json")
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(result_manifest, f, indent=2)

    logger.info(f"Saved environmental stability manifest to {manifest_path}")
    return result_manifest
