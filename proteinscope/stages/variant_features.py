"""
Per-Variant Feature Extraction for Pathogenicity Prediction.

Takes a single missense variant (position, ref_aa, alt_aa) and the cached outputs
from existing pipeline stages to build a ~60-dimensional feature vector covering:
  - Structural burial depth and Ramachandran classification (Stage 5)
  - Conservation across BLAST homolog set (Stage 3/7)
  - PTM site proximity (Stage 7)
  - Amino acid substitution physicochemical properties
  - Wild-type vs mutant whole-sequence ML feature deltas (Stage 8)
"""

import logging
import os
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger("ProteinScope.VariantFeatures")

# ── Physicochemical Reference Data ──────────────────────────────────────────

# Grantham distance matrix (Grantham, 1974) — selected common pairs
# Full 20x20 matrix for all standard amino acid substitutions
GRANTHAM_MATRIX = {
    ("A", "R"): 112, ("A", "N"): 111, ("A", "D"): 126, ("A", "C"): 195,
    ("A", "Q"): 91,  ("A", "E"): 107, ("A", "G"): 60,  ("A", "H"): 86,
    ("A", "I"): 94,  ("A", "L"): 96,  ("A", "K"): 106, ("A", "M"): 84,
    ("A", "F"): 113, ("A", "P"): 27,  ("A", "S"): 99,  ("A", "T"): 58,
    ("A", "W"): 148, ("A", "Y"): 112, ("A", "V"): 64,
    ("R", "N"): 86,  ("R", "D"): 96,  ("R", "C"): 180, ("R", "Q"): 43,
    ("R", "E"): 54,  ("R", "G"): 125, ("R", "H"): 29,  ("R", "I"): 97,
    ("R", "L"): 102, ("R", "K"): 26,  ("R", "M"): 91,  ("R", "F"): 97,
    ("R", "P"): 103, ("R", "S"): 110, ("R", "T"): 71,  ("R", "W"): 101,
    ("R", "Y"): 77,  ("R", "V"): 96,
    ("N", "D"): 23,  ("N", "C"): 139, ("N", "Q"): 46,  ("N", "E"): 42,
    ("N", "G"): 80,  ("N", "H"): 68,  ("N", "I"): 149, ("N", "L"): 153,
    ("N", "K"): 94,  ("N", "M"): 142, ("N", "F"): 158, ("N", "P"): 91,
    ("N", "S"): 46,  ("N", "T"): 65,  ("N", "W"): 174, ("N", "Y"): 143,
    ("N", "V"): 133,
    ("D", "C"): 154, ("D", "Q"): 61,  ("D", "E"): 45,  ("D", "G"): 94,
    ("D", "H"): 81,  ("D", "I"): 168, ("D", "L"): 172, ("D", "K"): 101,
    ("D", "M"): 160, ("D", "F"): 177, ("D", "P"): 108, ("D", "S"): 65,
    ("D", "T"): 85,  ("D", "W"): 181, ("D", "Y"): 160, ("D", "V"): 152,
    ("C", "Q"): 154, ("C", "E"): 170, ("C", "G"): 159, ("C", "H"): 174,
    ("C", "I"): 198, ("C", "L"): 198, ("C", "K"): 202, ("C", "M"): 196,
    ("C", "F"): 205, ("C", "P"): 169, ("C", "S"): 112, ("C", "T"): 149,
    ("C", "W"): 215, ("C", "Y"): 194, ("C", "V"): 192,
    ("Q", "E"): 29,  ("Q", "G"): 87,  ("Q", "H"): 24,  ("Q", "I"): 109,
    ("Q", "L"): 113, ("Q", "K"): 53,  ("Q", "M"): 101, ("Q", "F"): 116,
    ("Q", "P"): 76,  ("Q", "S"): 68,  ("Q", "T"): 42,  ("Q", "W"): 130,
    ("Q", "Y"): 99,  ("Q", "V"): 96,
    ("E", "G"): 98,  ("E", "H"): 40,  ("E", "I"): 134, ("E", "L"): 138,
    ("E", "K"): 56,  ("E", "M"): 126, ("E", "F"): 140, ("E", "P"): 93,
    ("E", "S"): 80,  ("E", "T"): 65,  ("E", "W"): 152, ("E", "Y"): 122,
    ("E", "V"): 121,
    ("G", "H"): 98,  ("G", "I"): 135, ("G", "L"): 138, ("G", "K"): 127,
    ("G", "M"): 127, ("G", "F"): 153, ("G", "P"): 42,  ("G", "S"): 56,
    ("G", "T"): 59,  ("G", "W"): 184, ("G", "Y"): 147, ("G", "V"): 109,
    ("H", "I"): 94,  ("H", "L"): 99,  ("H", "K"): 32,  ("H", "M"): 87,
    ("H", "F"): 100, ("H", "P"): 77,  ("H", "S"): 89,  ("H", "T"): 47,
    ("H", "W"): 115, ("H", "Y"): 83,  ("H", "V"): 84,
    ("I", "L"): 5,   ("I", "K"): 102, ("I", "M"): 10,  ("I", "F"): 21,
    ("I", "P"): 95,  ("I", "S"): 142, ("I", "T"): 89,  ("I", "W"): 61,
    ("I", "Y"): 33,  ("I", "V"): 29,
    ("L", "K"): 107, ("L", "M"): 15,  ("L", "F"): 22,  ("L", "P"): 98,
    ("L", "S"): 145, ("L", "T"): 92,  ("L", "W"): 61,  ("L", "Y"): 36,
    ("L", "V"): 32,
    ("K", "M"): 95,  ("K", "F"): 102, ("K", "P"): 103, ("K", "S"): 121,
    ("K", "T"): 78,  ("K", "W"): 110, ("K", "Y"): 85,  ("K", "V"): 97,
    ("M", "F"): 28,  ("M", "P"): 87,  ("M", "S"): 135, ("M", "T"): 81,
    ("M", "W"): 67,  ("M", "Y"): 36,  ("M", "V"): 21,
    ("F", "P"): 114, ("F", "S"): 155, ("F", "T"): 103, ("F", "W"): 40,
    ("F", "Y"): 22,  ("F", "V"): 50,
    ("P", "S"): 74,  ("P", "T"): 38,  ("P", "W"): 147, ("P", "Y"): 110,
    ("P", "V"): 68,
    ("S", "T"): 58,  ("S", "W"): 177, ("S", "Y"): 144, ("S", "V"): 124,
    ("T", "W"): 128, ("T", "Y"): 92,  ("T", "V"): 69,
    ("W", "Y"): 37,  ("W", "V"): 88,
    ("Y", "V"): 55,
}

# Kyte-Doolittle hydrophobicity scale
KYTE_DOOLITTLE = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# Amino acid molecular weights (Da, monoisotopic residue mass)
AA_MASS = {
    "A": 71.04, "R": 156.10, "N": 114.04, "D": 115.03, "C": 103.01,
    "Q": 128.06, "E": 129.04, "G": 57.02, "H": 137.06, "I": 113.08,
    "L": 113.08, "K": 128.09, "M": 131.04, "F": 147.07, "P": 97.05,
    "S": 87.03, "T": 101.05, "W": 186.08, "Y": 163.06, "V": 99.07,
}

# Amino acid charge class at pH 7 (approx)
AA_CHARGE = {
    "R": 1, "K": 1, "H": 0.1,  # His is weakly positive
    "D": -1, "E": -1,
    "A": 0, "N": 0, "C": 0, "Q": 0, "G": 0,
    "I": 0, "L": 0, "M": 0, "F": 0, "P": 0,
    "S": 0, "T": 0, "W": 0, "Y": 0, "V": 0,
}

# Polarity classification
POLAR_AAS = set("RNDQEHKSTY")
NONPOLAR_AAS = set("AFGILMPVW")

# BLOSUM62 substitution matrix (subset for standard 20 AAs)
# Stored as dict of dicts for easy lookup
_BLOSUM62_RAW = {
    "A": {"A":4,"R":-1,"N":-2,"D":-2,"C":0,"Q":-1,"E":-1,"G":0,"H":-2,"I":-1,"L":-1,"K":-1,"M":-1,"F":-2,"P":-1,"S":1,"T":0,"W":-3,"Y":-2,"V":0},
    "R": {"A":-1,"R":5,"N":0,"D":-2,"C":-3,"Q":1,"E":0,"G":-2,"H":0,"I":-3,"L":-2,"K":2,"M":-1,"F":-3,"P":-2,"S":-1,"T":-1,"W":-3,"Y":-2,"V":-3},
    "N": {"A":-2,"R":0,"N":6,"D":1,"C":-3,"Q":0,"E":0,"G":0,"H":1,"I":-3,"L":-3,"K":0,"M":-2,"F":-3,"P":-2,"S":1,"T":0,"W":-4,"Y":-2,"V":-3},
    "D": {"A":-2,"R":-2,"N":1,"D":6,"C":-3,"Q":0,"E":2,"G":-1,"H":-1,"I":-3,"L":-4,"K":-1,"M":-3,"F":-3,"P":-1,"S":0,"T":-1,"W":-4,"Y":-3,"V":-3},
    "C": {"A":0,"R":-3,"N":-3,"D":-3,"C":9,"Q":-3,"E":-4,"G":-3,"H":-3,"I":-1,"L":-1,"K":-3,"M":-1,"F":-2,"P":-3,"S":-1,"T":-1,"W":-2,"Y":-2,"V":-1},
    "Q": {"A":-1,"R":1,"N":0,"D":0,"C":-3,"Q":5,"E":2,"G":-2,"H":0,"I":-3,"L":-2,"K":1,"M":0,"F":-3,"P":-1,"S":0,"T":-1,"W":-2,"Y":-1,"V":-2},
    "E": {"A":-1,"R":0,"N":0,"D":2,"C":-4,"Q":2,"E":5,"G":-2,"H":0,"I":-3,"L":-3,"K":1,"M":-2,"F":-3,"P":-1,"S":0,"T":-1,"W":-3,"Y":-2,"V":-2},
    "G": {"A":0,"R":-2,"N":0,"D":-1,"C":-3,"Q":-2,"E":-2,"G":6,"H":-2,"I":-4,"L":-4,"K":-2,"M":-3,"F":-3,"P":-2,"S":0,"T":-2,"W":-2,"Y":-3,"V":-3},
    "H": {"A":-2,"R":0,"N":1,"D":-1,"C":-3,"Q":0,"E":0,"G":-2,"H":8,"I":-3,"L":-3,"K":-1,"M":-2,"F":-1,"P":-2,"S":-1,"T":-2,"W":-2,"Y":2,"V":-3},
    "I": {"A":-1,"R":-3,"N":-3,"D":-3,"C":-1,"Q":-3,"E":-3,"G":-4,"H":-3,"I":4,"L":2,"K":-3,"M":1,"F":0,"P":-3,"S":-2,"T":-1,"W":-3,"Y":-1,"V":3},
    "L": {"A":-1,"R":-2,"N":-3,"D":-4,"C":-1,"Q":-2,"E":-3,"G":-4,"H":-3,"I":2,"L":4,"K":-2,"M":2,"F":0,"P":-3,"S":-2,"T":-1,"W":-2,"Y":-1,"V":1},
    "K": {"A":-1,"R":2,"N":0,"D":-1,"C":-3,"Q":1,"E":1,"G":-2,"H":-1,"I":-3,"L":-2,"K":5,"M":-1,"F":-3,"P":-1,"S":0,"T":-1,"W":-3,"Y":-2,"V":-2},
    "M": {"A":-1,"R":-1,"N":-2,"D":-3,"C":-1,"Q":0,"E":-2,"G":-3,"H":-2,"I":1,"L":2,"K":-1,"M":5,"F":0,"P":-2,"S":-1,"T":-1,"W":-1,"Y":-1,"V":1},
    "F": {"A":-2,"R":-3,"N":-3,"D":-3,"C":-2,"Q":-3,"E":-3,"G":-3,"H":-1,"I":0,"L":0,"K":-3,"M":0,"F":6,"P":-4,"S":-2,"T":-2,"W":1,"Y":3,"V":-1},
    "P": {"A":-1,"R":-2,"N":-2,"D":-1,"C":-3,"Q":-1,"E":-1,"G":-2,"H":-2,"I":-3,"L":-3,"K":-1,"M":-2,"F":-4,"P":7,"S":-1,"T":-1,"W":-4,"Y":-3,"V":-2},
    "S": {"A":1,"R":-1,"N":1,"D":0,"C":-1,"Q":0,"E":0,"G":0,"H":-1,"I":-2,"L":-2,"K":0,"M":-1,"F":-2,"P":-1,"S":4,"T":1,"W":-3,"Y":-2,"V":-2},
    "T": {"A":0,"R":-1,"N":0,"D":-1,"C":-1,"Q":-1,"E":-1,"G":-2,"H":-2,"I":-1,"L":-1,"K":-1,"M":-1,"F":-2,"P":-1,"S":1,"T":5,"W":-2,"Y":-2,"V":0},
    "W": {"A":-3,"R":-3,"N":-4,"D":-4,"C":-2,"Q":-2,"E":-3,"G":-2,"H":-2,"I":-3,"L":-2,"K":-3,"M":-1,"F":1,"P":-4,"S":-3,"T":-2,"W":11,"Y":2,"V":-3},
    "Y": {"A":-2,"R":-2,"N":-2,"D":-3,"C":-2,"Q":-1,"E":-2,"G":-3,"H":2,"I":-1,"L":-1,"K":-2,"M":-1,"F":3,"P":-3,"S":-2,"T":-2,"W":2,"Y":7,"V":-1},
    "V": {"A":0,"R":-3,"N":-3,"D":-3,"C":-1,"Q":-2,"E":-2,"G":-3,"H":-3,"I":3,"L":1,"K":-2,"M":1,"F":-1,"P":-2,"S":-2,"T":0,"W":-3,"Y":-1,"V":4},
}


def get_grantham_distance(aa1: str, aa2: str) -> float:
    """Return Grantham distance between two amino acids. 0 if identical."""
    if aa1 == aa2:
        return 0.0
    key = (aa1, aa2) if (aa1, aa2) in GRANTHAM_MATRIX else (aa2, aa1)
    return float(GRANTHAM_MATRIX.get(key, 100.0))  # 100 as default for unknown pairs


def get_blosum62_score(aa1: str, aa2: str) -> float:
    """Return BLOSUM62 substitution score for aa1 → aa2."""
    return float(_BLOSUM62_RAW.get(aa1, {}).get(aa2, -4))


def get_hydrophobicity_delta(ref: str, alt: str) -> float:
    """Kyte-Doolittle hydrophobicity change (alt - ref)."""
    return KYTE_DOOLITTLE.get(alt, 0.0) - KYTE_DOOLITTLE.get(ref, 0.0)


def get_charge_change(ref: str, alt: str) -> float:
    """Change in charge class: positive=gained charge, negative=lost charge."""
    return AA_CHARGE.get(alt, 0.0) - AA_CHARGE.get(ref, 0.0)


def get_size_delta(ref: str, alt: str) -> float:
    """Molecular weight difference (alt - ref) in Daltons."""
    return AA_MASS.get(alt, 100.0) - AA_MASS.get(ref, 100.0)


def get_polarity_change(ref: str, alt: str) -> int:
    """Binary: 1 if polarity class changed (polar↔nonpolar), 0 otherwise."""
    ref_polar = ref in POLAR_AAS
    alt_polar = alt in POLAR_AAS
    return 1 if ref_polar != alt_polar else 0


# ── Structural Feature Extraction ───────────────────────────────────────────

def extract_burial_features(
    position: int,
    validation_data: Dict[str, Any],
) -> Dict[str, float]:
    """
    Extract burial depth and Ramachandran zone for a specific residue position
    from Stage 5 validation output.

    Searches the per-residue burial_results and phi_psi_records for matching
    residue number.
    """
    features = {
        "burial_neighbour_count": 0.0,
        "burial_class": 0.0,  # 0=surface, 1=partial, 2=buried
        "ramachandran_zone": 0.0,  # 0=favored, 1=allowed, 2=outlier
        "phi_angle": 0.0,
        "psi_angle": 0.0,
    }

    if not validation_data:
        return features

    # Search burial results for matching residue position
    # burial_results have format: {"residue": "A:ARG_175", "neighbours_within_8A": 22, "environment": "buried"}
    # This data is not directly stored in the validation_data dict returned by validate.py,
    # but we need to reconstruct it. We'll re-parse the PDB if needed.
    burial_env = validation_data.get("burial_environment", {})

    # The validation stage doesn't store per-residue burial in its return dict,
    # so we need to use the raw data. We'll handle this by re-extracting from PDB
    # in the calling code and passing it here.

    return features


def extract_burial_from_pdb(
    position: int,
    pdb_file: str,
    burial_radius: float = 8.0,
) -> Dict[str, float]:
    """
    Extract burial depth and Ramachandran classification for a specific residue
    by directly parsing the PDB structure file.
    """
    import math
    from Bio.PDB import PDBParser, Polypeptide

    features = {
        "burial_neighbour_count": 0.0,
        "burial_class": 0.0,
        "ramachandran_zone": 0.0,
        "phi_angle": 0.0,
        "psi_angle": 0.0,
    }

    if not pdb_file or not os.path.exists(pdb_file):
        return features

    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", pdb_file)

        # Collect all CA coordinates for burial calculation
        all_ca_coords = []
        target_ca_coord = None

        for model in structure:
            for chain in model:
                for residue in chain:
                    if "CA" in residue:
                        resid = residue.get_id()[1]
                        ca_coord = residue["CA"].get_coord()
                        all_ca_coords.append(ca_coord)
                        if resid == position:
                            target_ca_coord = ca_coord

        # Burial: count CA neighbors within radius
        if target_ca_coord is not None and all_ca_coords:
            coords_array = np.array(all_ca_coords)
            dists = np.linalg.norm(coords_array - target_ca_coord, axis=1)
            neighbour_count = int(np.sum((dists < burial_radius) & (dists > 0.001)))
            features["burial_neighbour_count"] = float(neighbour_count)

            if neighbour_count > 15:
                features["burial_class"] = 2.0  # buried
            elif neighbour_count >= 8:
                features["burial_class"] = 1.0  # partial
            else:
                features["burial_class"] = 0.0  # surface

        # Ramachandran angles at the target position
        ppb = Polypeptide.PPBuilder()
        for pp in ppb.build_peptides(structure):
            phi_psi_list = pp.get_phi_psi_list()
            for res, (phi, psi) in zip(pp, phi_psi_list):
                if res.get_id()[1] == position:
                    if phi is not None:
                        features["phi_angle"] = round(math.degrees(phi), 2)
                    if psi is not None:
                        features["psi_angle"] = round(math.degrees(psi), 2)

                    if phi is not None and psi is not None:
                        phi_deg = math.degrees(phi)
                        psi_deg = math.degrees(psi)
                        resname = res.get_resname()

                        # Import the classifier from validate stage
                        from proteinscope.stages.validate import classify_ramachandran_region
                        zone = classify_ramachandran_region(phi_deg, psi_deg, resname)
                        zone_map = {"favored": 0.0, "allowed": 1.0, "outlier": 2.0}
                        features["ramachandran_zone"] = zone_map.get(zone, 1.0)
                    break

    except Exception as e:
        logger.warning(f"Error extracting burial features from PDB at position {position}: {e}")

    return features


# ── Conservation Feature Extraction ─────────────────────────────────────────

def build_position_conservation(
    query_sequence: str,
    blast_hits: List[Dict[str, Any]],
) -> Dict[int, Dict[str, Any]]:
    """
    Build a position-indexed conservation profile from BLAST pairwise alignments.

    For each query position, records:
      - residues observed across homologs at the aligned column
      - fraction of homologs with the same residue as the query (conservation_score)
      - number of distinct amino acids at this position (diversity)

    Returns dict keyed by 1-indexed query position.
    """
    if not query_sequence or not blast_hits:
        return {}

    n_homologs = 0
    # position → list of aligned residues from homologs
    position_residues: Dict[int, List[str]] = {}

    for hit in blast_hits:
        query_aln = hit.get("query_sequence", "")
        sbjct_aln = hit.get("sbjct_sequence", "")

        if not query_aln or not sbjct_aln:
            continue

        n_homologs += 1

        # Map aligned positions back to query sequence positions
        query_pos = 0  # 0-indexed position in the unaligned query
        for q_char, s_char in zip(query_aln, sbjct_aln):
            if q_char != "-":
                query_pos += 1  # Now 1-indexed
                if s_char != "-":
                    if query_pos not in position_residues:
                        position_residues[query_pos] = []
                    position_residues[query_pos].append(s_char.upper())

    # Build conservation profile
    conservation: Dict[int, Dict[str, Any]] = {}
    for pos in range(1, len(query_sequence) + 1):
        residues_at_pos = position_residues.get(pos, [])
        query_aa = query_sequence[pos - 1]

        if residues_at_pos:
            same_count = sum(1 for r in residues_at_pos if r == query_aa)
            conservation_score = same_count / len(residues_at_pos)
            diversity = len(set(residues_at_pos))
        else:
            conservation_score = 0.0  # No data → unknown
            diversity = 0

        conservation[pos] = {
            "query_aa": query_aa,
            "conservation_score": round(conservation_score, 4),
            "diversity": diversity,
            "n_homologs_aligned": len(residues_at_pos),
            "residues": residues_at_pos,
        }

    return conservation


def extract_conservation_features(
    position: int,
    ref_aa: str,
    conservation_profile: Dict[int, Dict[str, Any]],
) -> Dict[str, float]:
    """
    Extract conservation features for a specific position from the pre-built profile.
    """
    features = {
        "conservation_score": 0.0,
        "homolog_diversity": 0.0,
        "n_homologs_aligned": 0.0,
    }

    pos_data = conservation_profile.get(position)
    if pos_data:
        features["conservation_score"] = pos_data["conservation_score"]
        features["homolog_diversity"] = float(pos_data["diversity"])
        features["n_homologs_aligned"] = float(pos_data["n_homologs_aligned"])

    return features


# ── PTM Proximity Features ──────────────────────────────────────────────────

def extract_ptm_features(
    position: int,
    ptm_sites: List[Dict[str, Any]],
) -> Dict[str, float]:
    """
    Extract PTM proximity features for a specific residue position.

    Features:
      - min_ptm_distance: minimum sequence distance to any known PTM site
      - is_ptm_site: 1 if this exact position is a PTM site, 0 otherwise
      - n_ptms_within_5: count of PTM sites within ±5 residues
    """
    features = {
        "min_ptm_distance": 999.0,
        "is_ptm_site": 0.0,
        "n_ptms_within_5": 0.0,
    }

    if not ptm_sites:
        return features

    ptm_positions = []
    for ptm in ptm_sites:
        pos_str = str(ptm.get("position", ""))
        # Handle single position and range formats
        if pos_str.isdigit():
            ptm_positions.append(int(pos_str))
        elif "-" in pos_str:
            parts = pos_str.split("-")
            try:
                start, end = int(parts[0]), int(parts[1])
                ptm_positions.extend(range(start, end + 1))
            except (ValueError, IndexError):
                pass

    if not ptm_positions:
        return features

    distances = [abs(position - p) for p in ptm_positions]
    features["min_ptm_distance"] = float(min(distances))
    features["is_ptm_site"] = 1.0 if position in ptm_positions else 0.0
    features["n_ptms_within_5"] = float(sum(1 for d in distances if d <= 5))

    return features


# ── Whole-Sequence ML Feature Deltas ────────────────────────────────────────

def extract_ml_feature_delta(
    query_sequence: str,
    position: int,
    alt_aa: str,
) -> Dict[str, float]:
    """
    Compute the difference in the Stage 8 whole-sequence feature vector between
    the wild-type sequence and the mutant sequence (single amino acid substitution).

    Uses the existing extract_sequence_feature_vector() from ml.py.
    Returns a dict of 'delta_*' features.
    """
    from proteinscope.stages.ml import extract_sequence_feature_vector

    # Build mutant sequence
    if position < 1 or position > len(query_sequence):
        return {}

    mut_seq = query_sequence[:position - 1] + alt_aa + query_sequence[position:]

    # Extract feature vectors for both
    wt_features = extract_sequence_feature_vector(query_sequence, label_name="wt")
    mt_features = extract_sequence_feature_vector(mut_seq, label_name="mt")

    # Compute deltas
    deltas: Dict[str, float] = {}
    for key in wt_features:
        wt_val = wt_features[key]
        mt_val = mt_features.get(key, wt_val)
        deltas[f"delta_{key}"] = round(mt_val - wt_val, 6)

    return deltas


# ── Master Feature Extraction ───────────────────────────────────────────────

def extract_variant_features(
    position: int,
    ref_aa: str,
    alt_aa: str,
    query_sequence: str,
    pdb_file: Optional[str] = None,
    blast_hits: Optional[List[Dict[str, Any]]] = None,
    ptm_sites: Optional[List[Dict[str, Any]]] = None,
    conservation_profile: Optional[Dict[int, Dict[str, Any]]] = None,
    include_ml_deltas: bool = True,
) -> Dict[str, float]:
    """
    Master feature extraction for a single missense variant.

    Combines features from:
      1. Structural burial and Ramachandran (from PDB)
      2. Conservation (from BLAST homolog alignment)
      3. PTM proximity (from UniProt PTM annotations)
      4. Physicochemical substitution properties (Grantham, BLOSUM62, etc.)
      5. Wild-type vs mutant whole-sequence ML feature deltas

    Args:
        position: 1-indexed residue position in the protein sequence
        ref_aa: Reference (wild-type) amino acid single-letter code
        alt_aa: Alternate (mutant) amino acid single-letter code
        query_sequence: Full wild-type protein sequence
        pdb_file: Path to PDB structure file (for burial/Ramachandran)
        blast_hits: BLAST hit list from Stage 3 (for conservation)
        ptm_sites: PTM site list from Stage 7 / UniProt ingest (for proximity)
        conservation_profile: Pre-built conservation profile (if available)
        include_ml_deltas: Whether to include the ~40 ML feature deltas

    Returns:
        Dict of feature_name → float value (~60 features total)
    """
    features: Dict[str, float] = {}

    # 1. Physicochemical substitution features (always available, no external data needed)
    features["grantham_distance"] = get_grantham_distance(ref_aa, alt_aa)
    features["blosum62_score"] = get_blosum62_score(ref_aa, alt_aa)
    features["hydrophobicity_delta"] = get_hydrophobicity_delta(ref_aa, alt_aa)
    features["charge_change"] = get_charge_change(ref_aa, alt_aa)
    features["size_delta"] = get_size_delta(ref_aa, alt_aa)
    features["polarity_change"] = float(get_polarity_change(ref_aa, alt_aa))

    # Normalized position in the sequence (0.0 = N-term, 1.0 = C-term)
    seq_len = len(query_sequence) if query_sequence else 1
    features["relative_position"] = round(position / seq_len, 4)

    # 2. Structural features (from PDB)
    if pdb_file and os.path.exists(pdb_file):
        burial_feats = extract_burial_from_pdb(position, pdb_file)
        features.update(burial_feats)
    else:
        features.update({
            "burial_neighbour_count": 0.0,
            "burial_class": 0.0,
            "ramachandran_zone": 0.0,
            "phi_angle": 0.0,
            "psi_angle": 0.0,
        })

    # 3. Conservation features
    if conservation_profile:
        cons_feats = extract_conservation_features(position, ref_aa, conservation_profile)
    elif blast_hits:
        profile = build_position_conservation(query_sequence, blast_hits)
        cons_feats = extract_conservation_features(position, ref_aa, profile)
    else:
        cons_feats = {
            "conservation_score": 0.0,
            "homolog_diversity": 0.0,
            "n_homologs_aligned": 0.0,
        }
    features.update(cons_feats)

    # 4. PTM proximity features
    ptm_feats = extract_ptm_features(position, ptm_sites or [])
    features.update(ptm_feats)

    # 5. ML feature deltas (wild-type vs mutant whole-sequence)
    if include_ml_deltas and query_sequence:
        try:
            ml_deltas = extract_ml_feature_delta(query_sequence, position, alt_aa)
            features.update(ml_deltas)
        except Exception as e:
            logger.warning(f"Failed to compute ML feature deltas: {e}")

    return features


def extract_feature_matrix(
    variants: List[Dict[str, Any]],
    query_sequence: str,
    pdb_file: Optional[str] = None,
    blast_hits: Optional[List[Dict[str, Any]]] = None,
    ptm_sites: Optional[List[Dict[str, Any]]] = None,
    include_ml_deltas: bool = True,
) -> Tuple[List[Dict[str, float]], List[str]]:
    """
    Extract feature vectors for a list of variants.

    Args:
        variants: List of dicts with keys: position, ref_aa, alt_aa
        query_sequence: Full wild-type protein sequence
        pdb_file: Path to PDB structure file
        blast_hits: BLAST hit list from Stage 3
        ptm_sites: PTM site list from Stage 7
        include_ml_deltas: Whether to include ML feature deltas

    Returns:
        Tuple of (feature_dicts_list, feature_names_list)
    """
    # Pre-build conservation profile once (expensive for many homologs)
    conservation_profile = None
    if blast_hits:
        logger.info("Building conservation profile from BLAST homologs...")
        conservation_profile = build_position_conservation(query_sequence, blast_hits)

    feature_dicts = []
    n_total = len(variants)

    for i, var in enumerate(variants):
        pos = var.get("position", 0)
        ref = var.get("ref_aa", "")
        alt = var.get("alt_aa", "")

        if not pos or not ref or not alt:
            continue

        # Validate that the reference AA matches the sequence
        if query_sequence and 0 < pos <= len(query_sequence):
            actual_aa = query_sequence[pos - 1]
            if actual_aa != ref:
                logger.debug(
                    f"Ref AA mismatch at position {pos}: expected {ref}, "
                    f"found {actual_aa} in sequence. Skipping variant {ref}{pos}{alt}."
                )
                continue

        features = extract_variant_features(
            position=pos,
            ref_aa=ref,
            alt_aa=alt,
            query_sequence=query_sequence,
            pdb_file=pdb_file,
            blast_hits=blast_hits,
            ptm_sites=ptm_sites,
            conservation_profile=conservation_profile,
            include_ml_deltas=include_ml_deltas,
        )
        feature_dicts.append(features)

        if (i + 1) % 100 == 0:
            logger.info(f"  Extracted features for {i + 1}/{n_total} variants...")

    feature_names = list(feature_dicts[0].keys()) if feature_dicts else []
    logger.info(
        f"Extracted {len(feature_names)} features for {len(feature_dicts)}/{n_total} variants"
    )

    return feature_dicts, feature_names
