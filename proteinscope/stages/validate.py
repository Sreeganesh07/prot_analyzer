"""
Stage 5: 3D Structure Validation.
Calculates Ramachandran dihedral angles, performs CA-CB bond length geometry checks,
and evaluates residue burial environments (Verify3D-style).
"""

import json
import logging
import math
import os
from typing import Any, Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")  # Non-interactive backend
import matplotlib.pyplot as plt

import numpy as np

from Bio.PDB import PDBParser, Polypeptide

logger = logging.getLogger("ProteinScope.Validate")


def calc_dihedral_angle(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray, p4: np.ndarray) -> float:
    """
    Calculate dihedral angle between 4 points in 3D space according to standard IUPAC conventions.
    Returns angle in degrees within [-180, 180].
    """
    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3

    # Normal vectors to the two planes
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)

    n1_norm = np.linalg.norm(n1)
    n2_norm = np.linalg.norm(n2)
    b2_norm = np.linalg.norm(b2)

    if n1_norm == 0 or n2_norm == 0 or b2_norm == 0:
        return 0.0

    n1 /= n1_norm
    n2 /= n2_norm
    u_b2 = b2 / b2_norm

    m1 = np.cross(n1, u_b2)
    x = np.dot(n1, n2)
    y = np.dot(m1, n2)
    angle = np.degrees(np.arctan2(y, x))
    return float(angle)


def classify_ramachandran_region(phi: float, psi: float, resname: str = "") -> str:
    """
    Classify (phi, psi) into favored, allowed, or outlier based on standard conformational zones.
    - Alpha-helix: phi in [-160, -40], psi in [-70, 50]
    - Beta-sheet: phi in [-180, -45], psi in [60, 180] or [-180, -150]
    - Left-handed helix: phi in [30, 90], psi in [10, 80]
    """
    if resname == "GLY":
        # Glycine has high conformational freedom
        return "favored"

    # Favored zones
    is_alpha = (-160 <= phi <= -40) and (-70 <= psi <= 50)
    is_beta = (-180 <= phi <= -45) and ((60 <= psi <= 180) or (-180 <= psi <= -150))
    is_left_alpha = (30 <= phi <= 90) and (10 <= psi <= 80)

    if is_alpha or is_beta or is_left_alpha:
        return "favored"

    # Allowed buffer zones (+/- 20 degrees)
    is_alpha_allowed = (-180 <= phi <= -20) and (-90 <= psi <= 70)
    is_beta_allowed = (-180 <= phi <= -30) and ((40 <= psi <= 180) or (-180 <= psi <= -130))
    is_left_allowed = (10 <= phi <= 110) and (-10 <= psi <= 100)

    if is_alpha_allowed or is_beta_allowed or is_left_allowed:
        return "allowed"

    return "outlier"


def plot_ramachandran_figure(
    phi_psi_data: List[Dict[str, Any]],
    output_png: str
) -> None:
    """Generate high-quality Ramachandran plot PNG using matplotlib."""
    if plt is None:
        logger.warning("Matplotlib not available; skipping Ramachandran figure generation.")
        return

    fig, ax = plt.subplots(figsize=(8, 8), dpi=300)

    # Background favored/allowed regions
    # Alpha region
    ax.fill_between([-160, -40], -70, 50, color="#d0e1fd", alpha=0.6, label="Favored (Alpha-Helix)")
    # Beta region
    ax.fill_between([-180, -45], 60, 180, color="#fef3c7", alpha=0.6, label="Favored (Beta-Sheet)")
    ax.fill_between([-180, -45], -180, -150, color="#fef3c7", alpha=0.6)
    # Left-handed alpha
    ax.fill_between([30, 90], 10, 80, color="#dcfce7", alpha=0.6, label="Favored (Left-Helix)")

    fav_phi = [d["phi"] for d in phi_psi_data if d["zone"] == "favored"]
    fav_psi = [d["psi"] for d in phi_psi_data if d["zone"] == "favored"]
    all_phi = [d["phi"] for d in phi_psi_data if d["zone"] == "allowed"]
    all_psi = [d["psi"] for d in phi_psi_data if d["zone"] == "allowed"]
    out_phi = [d["phi"] for d in phi_psi_data if d["zone"] == "outlier"]
    out_psi = [d["psi"] for d in phi_psi_data if d["zone"] == "outlier"]

    if fav_phi:
        ax.scatter(fav_phi, fav_psi, c="#2563eb", s=18, alpha=0.7, label=f"Favored ({len(fav_phi)})")
    if all_phi:
        ax.scatter(all_phi, all_psi, c="#16a34a", s=20, marker="^", alpha=0.8, label=f"Allowed ({len(all_phi)})")
    if out_phi:
        ax.scatter(out_phi, out_psi, c="#dc2626", s=32, marker="x", alpha=0.9, label=f"Outlier ({len(out_phi)})")

    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.set_xticks(np.arange(-180, 181, 60))
    ax.set_yticks(np.arange(-180, 181, 60))
    ax.axhline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.axvline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)

    ax.set_title(r"Ramachandran Dihedral Angle Distribution ($\phi$ vs $\psi$)", fontsize=14, fontweight="bold", pad=12)
    ax.set_xlabel(r"$\phi$ Dihedral Angle (degrees)", fontsize=12)
    ax.set_ylabel(r"$\psi$ Dihedral Angle (degrees)", fontsize=12)
    ax.legend(loc="upper right", frameon=True, facecolor="white", framealpha=0.9)
    ax.grid(True, linestyle=":", alpha=0.4)

    plt.tight_layout()
    fig.savefig(output_png)
    plt.close(fig)
    logger.info(f"Saved Ramachandran plot to {output_png}")


def _parse_pdb_atoms_pure_python(pdb_file: str) -> Tuple[List[Tuple[str, int, np.ndarray]], Dict[Tuple[str, int], Dict[str, np.ndarray]]]:
    """Extract CA coords and atom dictionary from PDB file in pure Python."""
    ca_coords = []
    residues = {}
    with open(pdb_file, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("ATOM  ") or line.startswith("HETATM"):
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain_id = line[21:22].strip() or "A"
                try:
                    res_seq = int(line[22:26].strip())
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    coord = np.array([x, y, z], dtype=float)
                    
                    key = (chain_id, res_seq)
                    if key not in residues:
                        residues[key] = {"resname": res_name}
                    residues[key][atom_name] = coord

                    if atom_name == "CA":
                        ca_coords.append((f"{chain_id}:{res_name}_{res_seq}", res_seq, coord))
                except (ValueError, IndexError):
                    continue
    return ca_coords, residues


def validate_protein_structure(
    pdb_file: str,
    output_dir: str,
    ca_cb_standard: float = 1.52,
    ca_cb_tolerance: float = 0.05,
    burial_radius: float = 8.0
) -> Dict[str, Any]:
    """
    Stage 5 Entry point:
    Validates PDB structure geometry, Ramachandran distribution, and residue burial.
    """
    logger.info(f"Starting Stage 5: Validating structure {pdb_file}")
    val_dir = os.path.join(output_dir, "validation")
    fig_dir = os.path.join(output_dir, "figures")
    os.makedirs(val_dir, exist_ok=True)
    os.makedirs(fig_dir, exist_ok=True)

    if not os.path.exists(pdb_file):
        logger.error(f"PDB file not found: {pdb_file}")
        return {"error": "PDB file not found"}

    phi_psi_records: List[Dict[str, Any]] = []
    ca_coords: List[Tuple[str, int, np.ndarray]] = []
    cb_violations: List[Dict[str, Any]] = []
    total_ca_cb_checked = 0

    parsed_biopython = False
    if PDBParser is not None and Polypeptide is not None:
        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("protein", pdb_file)

            # 1. Traverse polypeptides for Ramachandran and geometry checks
            polypeptide_builder = Polypeptide.PPBuilder()
            for pp in polypeptide_builder.build_peptides(structure):
                phi_psi_list = pp.get_phi_psi_list()
                for i, (res, (phi, psi)) in enumerate(zip(pp, phi_psi_list)):
                    resname = res.get_resname()
                    resid = res.get_id()[1]
                    chain = res.get_parent().get_id()

                    # Record CA coords for burial check
                    if "CA" in res:
                        ca_coords.append((f"{chain}:{resname}_{resid}", resid, res["CA"].get_coord()))

                    # Check CA-CB bond length (for non-GLY)
                    if resname != "GLY" and "CA" in res and "CB" in res:
                        total_ca_cb_checked += 1
                        ca_pos = res["CA"].get_coord()
                        cb_pos = res["CB"].get_coord()
                        dist = float(np.linalg.norm(ca_pos - cb_pos))
                        dev = abs(dist - ca_cb_standard)
                        if dev > ca_cb_tolerance:
                            cb_violations.append({
                                "chain": chain,
                                "resname": resname,
                                "resid": resid,
                                "bond_length": round(dist, 3),
                                "deviation": round(dev, 3),
                                "standard": ca_cb_standard,
                            })

                    # Calculate Ramachandran if angles exist
                    if phi is not None and psi is not None:
                        phi_deg = math.degrees(phi)
                        psi_deg = math.degrees(psi)
                        zone = classify_ramachandran_region(phi_deg, psi_deg, resname)
                        phi_psi_records.append({
                            "chain": chain,
                            "resid": resid,
                            "resname": resname,
                            "phi": round(phi_deg, 2),
                            "psi": round(psi_deg, 2),
                            "zone": zone,
                        })
            parsed_biopython = True
        except Exception as e:
            logger.warning(f"Biopython PDB parser encountered issue, falling back to pure Python parser: {e}")

    if not parsed_biopython:
        # Fallback pure Python PDB parser
        ca_coords, residues = _parse_pdb_atoms_pure_python(pdb_file)
        for (chain, resid), atoms in residues.items():
            resname = atoms.get("resname", "UNK")
            if resname != "GLY" and "CA" in atoms and "CB" in atoms:
                total_ca_cb_checked += 1
                ca_pos = atoms["CA"]
                cb_pos = atoms["CB"]
                dist = float(np.linalg.norm(ca_pos - cb_pos))
                dev = abs(dist - ca_cb_standard)
                if dev > ca_cb_tolerance:
                    cb_violations.append({
                        "chain": chain,
                        "resname": resname,
                        "resid": resid,
                        "bond_length": round(dist, 3),
                        "deviation": round(dev, 3),
                        "standard": ca_cb_standard,
                    })

    # 2. Ramachandran Statistics
    total_angles = len(phi_psi_records)
    favored_count = sum(1 for d in phi_psi_records if d["zone"] == "favored")
    allowed_count = sum(1 for d in phi_psi_records if d["zone"] == "allowed")
    outlier_count = sum(1 for d in phi_psi_records if d["zone"] == "outlier")

    ramachandran_stats = {
        "total_evaluated_residues": total_angles,
        "favored_count": favored_count,
        "favored_percent": round((favored_count / total_angles * 100), 2) if total_angles > 0 else 0.0,
        "allowed_count": allowed_count,
        "allowed_percent": round((allowed_count / total_angles * 100), 2) if total_angles > 0 else 0.0,
        "outlier_count": outlier_count,
        "outlier_percent": round((outlier_count / total_angles * 100), 2) if total_angles > 0 else 0.0,
    }

    # Plot Ramachandran diagram
    ram_png_path = os.path.join(fig_dir, "ramachandran.png")
    plot_ramachandran_figure(phi_psi_records, ram_png_path)

    # 3. Verify3D-Style Burial Check (CA Neighbour counting within 8 Angstroms)
    burial_results = []
    buried_count = 0
    partial_count = 0
    surface_count = 0

    if ca_coords:
        coords_array = np.array([c[2] for c in ca_coords])
        for i, (label, resid, coord) in enumerate(ca_coords):
            dists = np.linalg.norm(coords_array - coord, axis=1)
            # Count neighbours within burial_radius (excluding self distance = 0)
            neighbour_count = int(np.sum((dists < burial_radius) & (dists > 0.001)))
            
            if neighbour_count > 15:
                classification = "buried"
                buried_count += 1
            elif neighbour_count >= 8:
                classification = "partial"
                partial_count += 1
            else:
                classification = "surface"
                surface_count += 1

            burial_results.append({
                "residue": label,
                "neighbours_within_8A": neighbour_count,
                "environment": classification,
            })

    total_ca = len(ca_coords)
    burial_stats = {
        "total_ca_evaluated": total_ca,
        "buried_count": buried_count,
        "buried_percent": round((buried_count / total_ca * 100), 2) if total_ca > 0 else 0.0,
        "partial_count": partial_count,
        "partial_percent": round((partial_count / total_ca * 100), 2) if total_ca > 0 else 0.0,
        "surface_count": surface_count,
        "surface_percent": round((surface_count / total_ca * 100), 2) if total_ca > 0 else 0.0,
    }

    validation_data = {
        "pdb_file": os.path.basename(pdb_file),
        "ramachandran_statistics": ramachandran_stats,
        "geometry_ca_cb": {
            "total_bonds_checked": total_ca_cb_checked,
            "violations_count": len(cb_violations),
            "pass_rate_percent": round(((total_ca_cb_checked - len(cb_violations)) / total_ca_cb_checked * 100), 2) if total_ca_cb_checked > 0 else 100.0,
            "violations": cb_violations[:20],  # Sample of top 20
        },
        "burial_environment": burial_stats,
        "ramachandran_plot_path": ram_png_path,
    }

    out_json = os.path.join(val_dir, "structure_validation.json")
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(validation_data, f, indent=2)
    logger.info(f"Saved structure validation summary to {out_json}")

    return validation_data
