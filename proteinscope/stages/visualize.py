"""
Stage 6: Structure Visualization and PyMOL Scripting.
Generates PyMOL (.pml) commands, session (.pse), calculates RMSD superimposition,
and exports publication-quality structure figures (cartoon, surface, domain views).
"""

import json
import logging
import os
import shutil
import subprocess
from typing import Any, Dict, List, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import numpy as np

from Bio.PDB import PDBParser, Superimposer

logger = logging.getLogger("ProteinScope.Visualize")


def calculate_rmsd_superimposition(fixed_pdb: str, moving_pdb: str) -> Optional[float]:
    """
    Superimpose two PDB structures on matching CA atoms and calculate RMSD.
    """
    if PDBParser is None or Superimposer is None:
        logger.warning("Biopython Bio.PDB Superimposer unavailable; skipping RMSD superimposition.")
        return None

    try:
        parser = PDBParser(QUIET=True)
        struct_fix = parser.get_structure("fixed", fixed_pdb)
        struct_mov = parser.get_structure("moving", moving_pdb)

        ca_fix = [atom for atom in struct_fix.get_atoms() if atom.get_name() == "CA"]
        ca_mov = [atom for atom in struct_mov.get_atoms() if atom.get_name() == "CA"]

        min_len = min(len(ca_fix), len(ca_mov))
        if min_len < 3:
            return None

        sup = Superimposer()
        sup.set_atoms(ca_fix[:min_len], ca_mov[:min_len])
        rmsd = float(sup.rms)
        logger.info(f"Superimposed {len(ca_fix[:min_len])} CA atoms. RMSD: {rmsd:.3f} A")
        return round(rmsd, 3)
    except Exception as e:
        logger.warning(f"Superimposition failed between {fixed_pdb} and {moving_pdb}: {e}")
        return None


def generate_pymol_script(
    pdb_path: str,
    output_dir: str,
    domains: Optional[List[Dict[str, Any]]] = None
) -> str:
    """Generate a complete, executable PyMOL script (.pml) for rendering and session saving."""
    pdb_filename = os.path.basename(pdb_path)
    abs_pdb_path = os.path.abspath(pdb_path).replace("\\", "/")
    fig_dir = os.path.abspath(os.path.join(output_dir, "figures")).replace("\\", "/")
    struct_dir = os.path.abspath(os.path.join(output_dir, "structures")).replace("\\", "/")

    script_lines = [
        "# ProteinScope Automated PyMOL Visualization Script",
        "reinitialize",
        f"load {abs_pdb_path}, protein",
        "bg_color white",
        "hide everything, all",
        "",
        "# Color Scheme: Secondary Structure",
        "show cartoon, protein",
        "color cyan, ss h        # Alpha helices",
        "color yellow, ss s      # Beta sheets",
        "color white, ss l+''    # Loops / Coils",
        "set cartoon_fancy_helices, 1",
        "set ray_shadows, 0",
        "zoom protein",
        "orient protein",
        f"png {fig_dir}/structure_full.png, width=1600, height=1200, dpi=300, ray=1",
        "",
        "# Surface View",
        "show surface, protein",
        "set transparency, 0.2",
        f"png {fig_dir}/structure_surface.png, width=1600, height=1200, dpi=300, ray=1",
        "hide surface, protein",
        "",
        "# Domain View",
    ]

    palette = ["salmon", "forest", "marine", "orange", "magenta", "teal", "warmpink", "slate"]
    if domains:
        for idx, dom in enumerate(domains[:8]):
            color_name = palette[idx % len(palette)]
            locs = dom.get("locations", [])
            for loc in locs:
                start = loc.get("start")
                end = loc.get("end")
                if start and end:
                    dom_id = f"dom_{idx}_{start}_{end}"
                    script_lines.append(f"select {dom_id}, protein and resi {start}-{end}")
                    script_lines.append(f"color {color_name}, {dom_id}")
    else:
        script_lines.append("color spectrum, protein")

    script_lines.extend([
        f"png {fig_dir}/structure_domains.png, width=1600, height=1200, dpi=300, ray=1",
        "",
        f"save {struct_dir}/session.pse",
        "quit"
    ])

    pml_content = "\n".join(script_lines)
    pml_path = os.path.join(output_dir, "structures", "render_script.pml")
    with open(pml_path, "w", encoding="utf-8") as f:
        f.write(pml_content)
    logger.info(f"Generated PyMOL script at {pml_path}")
    return pml_path


def render_matplotlib_3d_fallbacks(
    pdb_path: str,
    output_dir: str,
    domains: Optional[List[Dict[str, Any]]] = None
) -> None:
    """
    Render 3D cartoon, surface/contact, and domain visualizations using matplotlib 3D
    to guarantee high-resolution figures when headless or PyMOL is not available.
    """
    fig_dir = os.path.join(output_dir, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("model", pdb_path)
    except Exception as e:
        logger.error(f"Error parsing PDB for fallback rendering: {e}")
        return

    ca_atoms = []
    res_indices = []
    res_names = []

    for model in structure:
        for chain in model:
            for residue in chain:
                if "CA" in residue:
                    ca_atoms.append(residue["CA"].get_coord())
                    res_indices.append(residue.get_id()[1])
                    res_names.append(residue.get_resname())

    if len(ca_atoms) < 2:
        return

    coords = np.array(ca_atoms)

    # 1. Full Structure Backbone / Secondary Structure View
    fig = plt.figure(figsize=(9, 7), dpi=300)
    ax = fig.add_subplot(111, projection="3d")
    # Color by sequence progression / secondary structure gradient
    colors = plt.cm.viridis(np.linspace(0, 1, len(coords)))
    ax.plot(coords[:, 0], coords[:, 1], coords[:, 2], color="#3b82f6", alpha=0.8, linewidth=2.0, label="C-alpha Backbone")
    ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], c=colors, s=15, alpha=0.9)
    ax.set_title(f"3D Structure Cartoon View ({os.path.basename(pdb_path)})", fontsize=13, fontweight="bold", pad=10)
    ax.set_xlabel("X (Å)")
    ax.set_ylabel("Y (Å)")
    ax.set_zlabel("Z (Å)")
    ax.grid(True, linestyle=":", alpha=0.3)
    plt.tight_layout()
    fig.savefig(os.path.join(fig_dir, "structure_full.png"))
    plt.close(fig)

    # 2. Surface / Density Envelope View
    fig = plt.figure(figsize=(9, 7), dpi=300)
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], c="#0284c7", s=120, alpha=0.4, edgecolors="none")
    ax.plot(coords[:, 0], coords[:, 1], coords[:, 2], color="#0f172a", alpha=0.7, linewidth=1.2)
    ax.set_title(f"3D Structure Surface Envelope View", fontsize=13, fontweight="bold", pad=10)
    ax.set_xlabel("X (Å)")
    ax.set_ylabel("Y (Å)")
    ax.set_zlabel("Z (Å)")
    ax.grid(True, linestyle=":", alpha=0.3)
    plt.tight_layout()
    fig.savefig(os.path.join(fig_dir, "structure_surface.png"))
    plt.close(fig)

    # 3. Domain Annotated View
    fig = plt.figure(figsize=(9, 7), dpi=300)
    ax = fig.add_subplot(111, projection="3d")
    domain_colors = np.array(["#94a3b8"] * len(coords), dtype=object)  # default gray
    
    palette = ["#ef4444", "#10b981", "#3b82f6", "#f59e0b", "#8b5cf6", "#ec4899"]
    if domains:
        for idx, dom in enumerate(domains[:6]):
            col = palette[idx % len(palette)]
            for loc in dom.get("locations", []):
                start = loc.get("start", 0)
                end = loc.get("end", 0)
                for i, r_idx in enumerate(res_indices):
                    if start <= r_idx <= end:
                        domain_colors[i] = col

    ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], c=domain_colors, s=35, alpha=0.9)
    ax.plot(coords[:, 0], coords[:, 1], coords[:, 2], color="#64748b", alpha=0.5, linewidth=1.5)
    ax.set_title(f"3D Structure Functional Domain Annotation", fontsize=13, fontweight="bold", pad=10)
    ax.set_xlabel("X (Å)")
    ax.set_ylabel("Y (Å)")
    ax.set_zlabel("Z (Å)")
    ax.grid(True, linestyle=":", alpha=0.3)
    plt.tight_layout()
    fig.savefig(os.path.join(fig_dir, "structure_domains.png"))
    plt.close(fig)

    logger.info("Rendered 3D structure figures via matplotlib fallback engine")


def visualize_protein_structure(
    structures_manifest: Dict[str, Any],
    output_dir: str,
    domains: Optional[List[Dict[str, Any]]] = None
) -> Dict[str, Any]:
    """
    Stage 6 Entry point:
    Finds best structure, attempts PyMOL rendering or falls back to Matplotlib 3D,
    calculates superimposition RMSD if multiple structures exist, and saves session.
    """
    logger.info("Starting Stage 6: Visualization and PyMOL Scripting")
    fig_dir = os.path.join(output_dir, "figures")
    struct_dir = os.path.join(output_dir, "structures")
    os.makedirs(fig_dir, exist_ok=True)
    os.makedirs(struct_dir, exist_ok=True)

    best = structures_manifest.get("best_structure")
    if not best or not os.path.exists(best.get("file_path", "")):
        logger.warning("No structure available to visualize")
        return {"error": "No structure available"}

    best_pdb = best["file_path"]

    # 1. Superimpose if both experimental and AlphaFold exist
    rmsd = None
    exp_info = structures_manifest.get("experimental")
    af_info = structures_manifest.get("alphafold")

    if exp_info and af_info:
        exp_file = exp_info.get("file_path")
        af_file = af_info.get("file_path")
        if exp_file and af_file and os.path.exists(exp_file) and os.path.exists(af_file):
            rmsd = calculate_rmsd_superimposition(exp_file, af_file)

    # 2. Generate PyMOL PML script
    pml_script_path = generate_pymol_script(best_pdb, output_dir, domains=domains)

    # 3. Try executing PyMOL if installed
    pymol_executed = False
    pymol_bin = shutil.which("pymol")
    if pymol_bin:
        try:
            cmd = [pymol_bin, "-c", "-q", pml_script_path]
            subprocess.run(cmd, timeout=30, check=True)
            pymol_executed = True
            logger.info("Successfully executed PyMOL batch rendering")
        except Exception as e:
            logger.warning(f"PyMOL execution failed or timed out: {e}")

    # 4. If PyMOL was not run or figures missing, render matplotlib fallbacks
    full_png = os.path.join(fig_dir, "structure_full.png")
    if not os.path.exists(full_png) or not pymol_executed:
        render_matplotlib_3d_fallbacks(best_pdb, output_dir, domains=domains)

    # Create dummy / placeholder session.pse if not generated by pymol
    pse_path = os.path.join(struct_dir, "session.pse")
    if not os.path.exists(pse_path):
        with open(pse_path, "w", encoding="utf-8") as f:
            f.write("# ProteinScope PyMOL Session Placeholder\n# Load render_script.pml in PyMOL GUI to view interactive 3D session.\n")

    return {
        "best_structure_used": best["source"],
        "pdb_file": os.path.basename(best_pdb),
        "rmsd_superimposition_angstroms": rmsd,
        "pymol_script_path": pml_script_path,
        "pymol_session_path": pse_path,
        "figures": {
            "full_view": os.path.join(fig_dir, "structure_full.png"),
            "surface_view": os.path.join(fig_dir, "structure_surface.png"),
            "domains_view": os.path.join(fig_dir, "structure_domains.png"),
        }
    }


def generate_interactive_html_report(report_data: Dict[str, Any], output_dir: str) -> str:
    """
    Generate a modern, responsive, standalone HTML report with an interactive 3Dmol.js protein viewer.
    Embeds downloaded local PDB coordinates directly so rendering never fails or requires external IDs.
    """
    report_dir = os.path.join(output_dir, "report")
    os.makedirs(report_dir, exist_ok=True)
    html_path = os.path.join(report_dir, "protein_viewer.html")

    accession = report_data.get("accession", "UNKNOWN")
    protein_name = report_data.get("protein_name", "Unknown Protein")
    organism = report_data.get("organism", "Unknown Organism")
    seq_len = report_data.get("sequence_length", 0)
    duration = report_data.get("run_metadata", {}).get("execution_duration_seconds", 0.0)

    props = report_data.get("properties", {})
    mw = props.get("molecular_weight", 0.0)
    pi = props.get("isoelectric_point", 0.0)
    instab = props.get("instability_index", 0.0)
    is_stable = props.get("is_stable", True)
    gravy = props.get("gravy", 0.0)
    aroma = props.get("aromaticity", 0.0)
    sec_struct = props.get("secondary_structure_fraction", {})
    helix_pct = round(sec_struct.get("helix", 0.0) * 100, 1)
    sheet_pct = round(sec_struct.get("sheet", 0.0) * 100, 1)
    turn_pct = round(sec_struct.get("turn", 0.0) * 100, 1)

    grouped = props.get("grouped_composition") or {}
    pos_pct = grouped.get("positively_charged_percent", 0.0)
    neg_pct = grouped.get("negatively_charged_percent", 0.0)
    hydro_pct = grouped.get("hydrophobic_percent", 0.0)

    val = report_data.get("validation") or {}
    ram = val.get("ramachandran") or {}
    fav_pct = ram.get("favored_percent", 0.0)
    allow_pct = ram.get("allowed_percent", 0.0)
    out_pct = ram.get("outlier_percent", 0.0)
    ca_cb_pass = val.get("ca_cb_bond_pass_rate", 100.0)
    burial = val.get("burial") or {}
    surf_pct = burial.get("surface_percent", 0.0)
    part_pct = burial.get("partial_percent", 0.0)
    bur_pct = burial.get("buried_percent", 0.0)

    ptm_sites = report_data.get("ptm_sites", [])
    glyco_sites = [p for p in ptm_sites if "glycosylation" in str(p.get("type", "")).lower()]
    disulf_bonds = [p for p in ptm_sites if "disulfide" in str(p.get("type", "")).lower()]

    struct_source = report_data.get("structure_source", "None")
    struct_file = report_data.get("structure_best_file", "")
    rmsd = report_data.get("structure_rmsd_angstroms", "N/A")

    # Read all available local PDB files to embed directly in HTML
    struct_dir = os.path.join(output_dir, "structures")
    struct_models: Dict[str, str] = {}
    struct_labels: Dict[str, str] = {}

    # Check manifest if available
    manifest_p = os.path.join(struct_dir, "structures_manifest.json")
    manifest_data = {}
    if os.path.exists(manifest_p):
        try:
            with open(manifest_p, "r", encoding="utf-8") as mf:
                manifest_data = json.load(mf)
        except Exception:
            pass

    # 1. Experimental PDB
    exp_info = manifest_data.get("experimental")
    if exp_info and exp_info.get("file_path") and os.path.exists(exp_info["file_path"]):
        try:
            with open(exp_info["file_path"], "r", encoding="utf-8", errors="ignore") as f:
                struct_models["experimental"] = f.read()
                struct_labels["experimental"] = f"Experimental ({exp_info.get('pdb_id', 'PDB')})"
        except Exception:
            pass

    # 2. AlphaFold PDB
    af_info = manifest_data.get("alphafold")
    if af_info and af_info.get("file_path") and os.path.exists(af_info["file_path"]):
        try:
            with open(af_info["file_path"], "r", encoding="utf-8", errors="ignore") as f:
                struct_models["alphafold"] = f.read()
                struct_labels["alphafold"] = "AlphaFold Model"
        except Exception:
            pass

    # 3. Homology SWISS-MODEL PDB
    hom_info = manifest_data.get("homology")
    if hom_info and hom_info.get("file_path") and os.path.exists(hom_info["file_path"]):
        try:
            with open(hom_info["file_path"], "r", encoding="utf-8", errors="ignore") as f:
                struct_models["homology"] = f.read()
                tmpl = hom_info.get("template", "Homology")
                struct_labels["homology"] = f"SWISS-MODEL ({tmpl.split('.')[0].upper()})"
        except Exception:
            pass

    # Determine default initial model key
    default_model_key = "experimental"
    if "experimental" not in struct_models:
        if "alphafold" in struct_models:
            default_model_key = "alphafold"
        elif "homology" in struct_models:
            default_model_key = "homology"
        elif struct_models:
            default_model_key = list(struct_models.keys())[0]
        else:
            default_model_key = "none"

    ml = report_data.get("ml_features", {})
    top_features = ml.get("top_features", {})

    # Pathogenicity section
    patho_data = report_data.get("pathogenicity", {})
    mutation_name = report_data.get("mutation", "")
    patho_html = ""
    if patho_data:
        verdict = patho_data.get("prediction", "unknown").upper()
        conf = patho_data.get("confidence", 0.0)
        conf_pct = round(conf * 100, 1)
        path_prob = round(patho_data.get("pathogenic_probability", 0.0) * 100, 1)
        justification = patho_data.get("justification", "")
        is_path = verdict == "PATHOGENIC"
        badge_bg = "rgba(244,63,94,0.15)" if is_path else "rgba(16,185,129,0.15)"
        badge_border = "var(--accent-rose)" if is_path else "var(--accent-emerald)"
        badge_color = "#fb7185" if is_path else "#34d399"

        just_items = ""
        for line in justification.split("\n"):
            line = line.strip()
            if line.startswith("•"):
                line = line[1:].strip()
            if line:
                just_items += f'<li style="margin-bottom:8px; line-height:1.5; font-size:13px; color:var(--text-secondary);">{line}</li>'

        patho_html = f'''
        <div class="card" style="grid-column: 1 / -1; border-left: 4px solid {badge_border}; background: linear-gradient(180deg, var(--bg-card) 0%, rgba(26,31,53,0.8) 100%);">
            <div style="display:flex; justify-content:space-between; align-items:flex-start; flex-wrap:wrap; gap:12px; margin-bottom:16px;">
                <div style="display:flex; align-items:center; gap:12px;">
                    <div class="card-icon" style="background:{badge_bg}; color:{badge_color}; font-size:20px;">🛡️</div>
                    <div>
                        <div class="card-title" style="font-size:18px;">Pathogenicity Prediction: <span style="color:{badge_color};">{verdict}</span></div>
                        <div class="card-subtitle">Mutation Target: <strong style="color:var(--text-primary); font-family:\'JetBrains Mono\', monospace;">{mutation_name}</strong></div>
                    </div>
                </div>
                <div style="display:flex; gap:12px; align-items:center;">
                    <div style="text-align:right;">
                        <div style="font-size:11px; color:var(--text-muted); text-transform:uppercase;">Confidence Score</div>
                        <div style="font-size:20px; font-weight:800; font-family:\'JetBrains Mono\', monospace; color:{badge_color};">{conf_pct}%</div>
                    </div>
                    <div style="text-align:right; border-left:1px solid var(--border); padding-left:12px;">
                        <div style="font-size:11px; color:var(--text-muted); text-transform:uppercase;">Pathogenic Probability</div>
                        <div style="font-size:20px; font-weight:800; font-family:\'JetBrains Mono\', monospace; color:var(--text-primary);">{path_prob}%</div>
                    </div>
                </div>
            </div>
            <div style="background:rgba(10,14,26,0.6); border:1px solid var(--border); border-radius:var(--radius-sm); padding:16px;">
                <div style="font-size:12px; font-weight:700; color:var(--accent-blue); text-transform:uppercase; letter-spacing:1px; margin-bottom:10px;">Structural & Evolutionary Justification</div>
                <ul style="padding-left:20px; margin:0;">
                    {just_items}
                </ul>
            </div>
        </div>'''

    glyco_tags = "".join([f'<span class="ptm-tag glyco">{p.get("residue", "Res")}{p.get("position", "")}</span>' for p in glyco_sites[:24]])
    if not glyco_tags:
        glyco_tags = '<span style="color:var(--text-muted); font-size:12px;">No glycosylation sites detected</span>'

    disulf_tags = "".join([f'<span class="ptm-tag disulfide">{p.get("position", "")}</span>' for p in disulf_bonds[:16]])
    if not disulf_tags:
        disulf_tags = '<span style="color:var(--text-muted); font-size:12px;">No disulfide bonds detected</span>'

    ml_rows = ""
    for feat, imp in list(top_features.items())[:5]:
        val_pct = round(imp * 100, 1)
        ml_rows += f'''
        <div class="progress-row">
            <div class="progress-label"><span class="name">{feat}</span><span class="val">{val_pct}%</span></div>
            <div class="progress-bar"><div class="progress-fill purple" data-width="{min(100, val_pct * 4)}"></div></div>
        </div>'''

    # Build source tabs HTML
    tabs_html = ""
    for key, label in struct_labels.items():
        is_active = "active" if key == default_model_key else ""
        tabs_html += f'<button class="source-tab {is_active}" onclick="loadSource(\'{key}\')">{label}</button> '

    # Sequence and Titratable Residues for Real-time Henderson-Hasselbalch Simulation
    sequence = report_data.get("sequence", "")
    if not sequence:
        fasta_path = os.path.join(output_dir, "sequences", f"{accession}.fasta")
        if os.path.exists(fasta_path):
            try:
                with open(fasta_path, "r", encoding="utf-8") as f:
                    sequence = "".join(f.read().split("\n")[1:]).strip()
            except Exception:
                pass

    seq_u = sequence.upper() if sequence else ""
    titratable_counts = {
        "D": seq_u.count("D"),
        "E": seq_u.count("E"),
        "H": seq_u.count("H"),
        "C": seq_u.count("C"),
        "Y": seq_u.count("Y"),
        "K": seq_u.count("K"),
        "R": seq_u.count("R"),
        "len": len(seq_u) or seq_len or 300,
    }
    counts_json = json.dumps(titratable_counts)

    # Environmental Stability Data (Alberts NBK26830)
    stab_data = report_data.get("environmental_stability", {})
    if not stab_data and sequence:
        try:
            from proteinscope.stages.stability import analyze_protein_stability
            stab_data = analyze_protein_stability(
                sequence=sequence,
                output_dir=output_dir,
                properties=props,
                ptm_sites=ptm_sites,
            )
        except Exception:
            stab_data = {}

    tm_est = stab_data.get("estimated_melting_temperature_celsius", 62.0)
    active_eval = stab_data.get("active_evaluation", {})
    init_temp = active_eval.get("temperature_celsius", 37.0)
    init_ph = active_eval.get("ph", 7.4)
    init_state = active_eval.get("state", "NATIVE")
    init_badge = active_eval.get("badge_label", "Native Folded State")
    init_color = active_eval.get("badge_color", "#10b981")
    init_mechanism = active_eval.get("molecular_mechanism", "Optimal folded conformation in standard conditions.")
    init_folded_pct = active_eval.get("fraction_folded_percent", 98.8)
    init_q = active_eval.get("net_charge", 0.0)
    init_salt = active_eval.get("salt_bridge_retention_percent", 99.0)
    init_dg = active_eval.get("delta_g_folding_kcal_mol", -12.5)

    thresholds = active_eval.get("thresholds", {
        "tm_celsius": tm_est,
        "temp_min_stable": 15.0,
        "temp_max_stable": round(tm_est - 8.0, 1),
        "temp_degradation": round(max(70.0, tm_est + 10.0), 1),
        "ph_min_stable": 5.5,
        "ph_max_stable": 8.5,
        "ph_acid_degradation": 3.0,
        "ph_alkaline_degradation": 11.5,
    })
    thresholds_json = json.dumps(thresholds)

    models_json = json.dumps(struct_models)
    labels_json = json.dumps(struct_labels)

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{accession} // {protein_name} — ProteinScope Workstation</title>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800&family=JetBrains+Mono:wght@400;500;600;700&display=swap" rel="stylesheet">
    <style>
        :root {{
            --bg-canvas: #090b10;
            --bg-panel: #10141d;
            --bg-subpanel: #0a0d14;
            --border-technical: #1d2433;
            --border-highlight: #0284c7;
            --text-title: #f8fafc;
            --text-body: #94a3b8;
            --text-muted: #526077;
            --accent-cyan: #06b6d4;
            --accent-emerald: #10b981;
            --accent-amber: #f59e0b;
            --accent-rose: #f43f5e;
            --accent-blue: #3b82f6;
            --radius: 8px;
            --radius-sm: 4px;
        }}
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
            background: var(--bg-canvas);
            background-image: radial-gradient(#182030 1px, transparent 1px);
            background-size: 24px 24px;
            color: var(--text-body);
            min-height: 100vh;
            padding: 24px 20px;
        }}
        .container {{ max-width: 1420px; margin: 0 auto; }}

        /* ── WORKSTATION TOP BAR ── */
        .workstation-topbar {{
            display: flex; justify-content: space-between; align-items: center;
            border-bottom: 1px solid var(--border-technical); padding-bottom: 16px; margin-bottom: 24px;
            flex-wrap: wrap; gap: 12px;
        }}
        .brand-block {{ display: flex; align-items: center; gap: 12px; }}
        .brand-logo {{
            width: 32px; height: 32px; border-radius: var(--radius-sm);
            background: #0284c7; color: #ffffff; display: flex;
            align-items: center; justify-content: center; font-weight: 900;
            font-family: 'JetBrains Mono', monospace; font-size: 14px;
            letter-spacing: -1px;
        }}
        .brand-meta {{ display: flex; flex-direction: column; }}
        .brand-sys {{ font-size: 11px; font-weight: 700; color: var(--accent-cyan); letter-spacing: 1.5px; text-transform: uppercase; }}
        .brand-sub {{ font-size: 11px; color: var(--text-muted); font-family: 'JetBrains Mono', monospace; }}

        .workstation-telemetry-pill {{
            display: flex; gap: 16px; font-family: 'JetBrains Mono', monospace; font-size: 11px;
            background: var(--bg-panel); border: 1px solid var(--border-technical);
            padding: 8px 16px; border-radius: var(--radius-sm);
        }}
        .w-item {{ display: flex; gap: 6px; }}
        .w-lbl {{ color: var(--text-muted); }}
        .w-val {{ color: var(--text-title); font-weight: 600; }}

        /* ── HERO METADATA CARD ── */
        .hero-banner {{
            background: var(--bg-panel); border: 1px solid var(--border-technical);
            border-radius: var(--radius); padding: 24px; margin-bottom: 20px;
            position: relative;
        }}
        .hero-title-row {{ display: flex; justify-content: space-between; align-items: flex-start; flex-wrap: wrap; gap: 16px; }}
        .hero-accession {{
            display: inline-block; font-family: 'JetBrains Mono', monospace;
            background: rgba(6,182,212,0.1); color: var(--accent-cyan);
            border: 1px solid rgba(6,182,212,0.25); padding: 3px 10px;
            border-radius: var(--radius-sm); font-size: 12px; font-weight: 700;
            letter-spacing: 1px; margin-bottom: 8px;
        }}
        .hero-title {{
            font-size: 30px; font-weight: 800; color: var(--text-title);
            letter-spacing: -0.5px; margin-bottom: 4px;
        }}
        .hero-organism {{ font-size: 14px; color: var(--text-muted); font-style: italic; }}

        .hero-chips {{ display: flex; gap: 12px; flex-wrap: wrap; margin-top: 20px; }}
        .h-chip {{
            background: var(--bg-subpanel); border: 1px solid var(--border-technical);
            border-radius: var(--radius-sm); padding: 8px 14px; display: flex; flex-direction: column; gap: 2px;
        }}
        .h-chip .lbl {{ font-size: 10px; color: var(--text-muted); text-transform: uppercase; font-weight: 700; }}
        .h-chip .val {{ font-size: 13px; font-weight: 700; color: var(--text-title); font-family: 'JetBrains Mono', monospace; }}

        /* ── 3D MOLECULAR VIEWER DECK ── */
        .viewer-card {{
            background: var(--bg-panel); border: 1px solid var(--border-technical);
            border-radius: var(--radius); overflow: hidden; margin-bottom: 20px;
            transition: border-color 0.25s ease;
        }}
        .viewer-card.degraded-glow {{
            border-color: var(--accent-rose) !important;
        }}
        .viewer-toolbar {{
            display: flex; align-items: center; justify-content: space-between;
            padding: 12px 18px; border-bottom: 1px solid var(--border-technical);
            background: #0d1017; flex-wrap: wrap; gap: 10px;
        }}
        .toolbar-title {{
            font-size: 12px; font-weight: 700; color: var(--text-title);
            letter-spacing: 1px; text-transform: uppercase; font-family: 'JetBrains Mono', monospace;
        }}
        .btn-group {{ display: flex; gap: 4px; flex-wrap: wrap; }}
        .btn {{
            padding: 6px 12px; border-radius: var(--radius-sm); border: 1px solid var(--border-technical);
            background: var(--bg-panel); color: var(--text-body); font-size: 11px;
            font-weight: 600; cursor: pointer; transition: all 0.15s; font-family: 'JetBrains Mono', monospace;
        }}
        .btn:hover {{ border-color: var(--accent-cyan); color: var(--text-title); }}
        .btn.active {{ background: #0284c7; color: #ffffff; border-color: #0284c7; }}
        .source-tab {{
            padding: 5px 12px; border-radius: var(--radius-sm); border: 1px solid var(--border-technical);
            background: transparent; color: var(--text-muted); font-size: 11px; font-weight: 600; cursor: pointer;
            font-family: 'JetBrains Mono', monospace;
        }}
        .source-tab.active {{ background: #1e293b; color: var(--accent-cyan); border-color: var(--accent-cyan); }}

        .viewer-container {{ position: relative; height: 500px; background: #07090e; }}
        #viewer-3d {{ width: 100%; height: 100%; }}
        .hud-reticle {{
            position: absolute; color: rgba(255,255,255,0.15); font-family: 'JetBrains Mono', monospace;
            font-size: 14px; pointer-events: none; user-select: none;
        }}
        .hud-tl {{ top: 12px; left: 12px; }}
        .hud-tr {{ top: 12px; right: 12px; }}
        .hud-bl {{ bottom: 12px; left: 12px; }}
        .hud-br {{ bottom: 12px; right: 12px; }}

        .viewer-info-overlay {{
            position: absolute; bottom: 14px; left: 14px;
            background: rgba(10,13,20,0.92); backdrop-filter: blur(8px);
            border: 1px solid var(--border-technical); border-radius: var(--radius-sm);
            padding: 8px 14px; font-size: 11px; color: var(--text-muted); font-family: 'JetBrains Mono', monospace;
            display: flex; gap: 14px; align-items: center;
        }}
        .viewer-info-overlay strong {{ color: var(--text-title); }}

        .loading-overlay {{
            position: absolute; inset: 0; display: flex; flex-direction: column;
            align-items: center; justify-content: center; background: rgba(9,11,16,0.92); z-index: 10;
        }}
        .loading-overlay.hidden {{ display: none; }}
        .spinner {{
            width: 36px; height: 36px; border: 2px solid var(--border-technical);
            border-top-color: var(--accent-cyan); border-radius: 50%; animation: spin 0.8s linear infinite; margin-bottom: 12px;
        }}
        @keyframes spin {{ to {{ transform: rotate(360deg); }} }}

        /* ── ENVIRONMENTAL CALIBRATION & BIOPHYSICAL CHAMBER ── */
        .env-card {{
            background: var(--bg-panel); border: 1px solid var(--border-technical);
            border-radius: var(--radius); padding: 22px; margin-bottom: 24px;
        }}
        .env-header {{
            display: flex; justify-content: space-between; align-items: center;
            flex-wrap: wrap; gap: 14px; margin-bottom: 18px; border-bottom: 1px solid var(--border-technical);
            padding-bottom: 14px;
        }}
        .env-tagline {{ display: flex; flex-direction: column; gap: 2px; }}
        .env-title {{
            font-size: 15px; font-weight: 800; color: var(--text-title);
            letter-spacing: 0.5px; text-transform: uppercase; font-family: 'JetBrains Mono', monospace;
        }}
        .env-subtitle {{ font-size: 12px; color: var(--text-muted); }}

        .env-status-lamp {{
            display: flex; align-items: center; gap: 8px; font-family: 'JetBrains Mono', monospace;
            font-size: 12px; font-weight: 700; padding: 6px 14px; border-radius: var(--radius-sm);
            background: var(--bg-subpanel); border: 1px solid var(--border-technical);
        }}
        .lamp-dot {{ width: 8px; height: 8px; border-radius: 50%; }}

        .env-presets-strip {{
            display: flex; gap: 6px; flex-wrap: wrap; margin-bottom: 20px;
        }}
        .preset-chip {{
            padding: 6px 12px; border-radius: var(--radius-sm); border: 1px solid var(--border-technical);
            background: var(--bg-subpanel); color: var(--text-body); font-size: 11px;
            font-weight: 600; cursor: pointer; transition: all 0.15s; font-family: 'JetBrains Mono', monospace;
        }}
        .preset-chip:hover {{ border-color: var(--accent-cyan); color: var(--text-title); }}
        .preset-chip.active {{
            background: #0369a1; color: #ffffff; border-color: #0369a1;
        }}

        .instrument-dials-grid {{
            display: grid; grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
            gap: 16px; margin-bottom: 20px;
        }}
        .dial-box {{
            background: var(--bg-subpanel); border: 1px solid var(--border-technical);
            border-radius: var(--radius-sm); padding: 16px;
        }}
        .dial-top {{ display: flex; justify-content: space-between; align-items: baseline; margin-bottom: 12px; }}
        .dial-label {{ font-size: 11px; font-weight: 700; color: var(--text-muted); text-transform: uppercase; font-family: 'JetBrains Mono', monospace; }}
        .dial-readout {{ font-size: 22px; font-weight: 800; font-family: 'JetBrains Mono', monospace; color: var(--text-title); }}

        .slider-bracket {{ margin-bottom: 6px; }}
        .env-slider {{
            -webkit-appearance: none; width: 100%; height: 6px;
            border-radius: 3px; outline: none; cursor: pointer;
        }}
        .env-slider::-webkit-slider-thumb {{
            -webkit-appearance: none; appearance: none; width: 18px; height: 18px;
            border-radius: 2px; background: #ffffff; border: 2px solid #0284c7;
            cursor: pointer; transition: transform 0.1s;
        }}
        .env-slider::-webkit-slider-thumb:hover {{ transform: scale(1.15); }}

        #temp-slider {{
            background: linear-gradient(90deg, #38bdf8 0%, #10b981 35%, #f59e0b 60%, #f43f5e 100%);
        }}
        #ph-slider {{
            background: linear-gradient(90deg, #f43f5e 0%, #f59e0b 25%, #10b981 50%, #0ea5e9 75%, #a855f7 100%);
        }}

        .graduations {{
            display: flex; justify-content: space-between; font-size: 10px;
            color: var(--text-muted); font-family: 'JetBrains Mono', monospace; margin-top: 4px;
        }}

        /* Range Tuning Drawer */
        .drawer-toggle {{
            display: inline-flex; align-items: center; gap: 8px; font-size: 11px;
            font-weight: 700; color: var(--accent-cyan); font-family: 'JetBrains Mono', monospace;
            cursor: pointer; margin-bottom: 14px; text-transform: uppercase;
        }}
        .drawer-panel {{
            background: #080a0f; border: 1px dashed var(--border-technical);
            border-radius: var(--radius-sm); padding: 14px; margin-bottom: 20px;
        }}
        .drawer-grid {{
            display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
            gap: 10px; margin-bottom: 12px;
        }}
        .d-field {{ display: flex; flex-direction: column; gap: 3px; }}
        .d-field label {{ font-size: 10px; color: var(--text-muted); font-family: 'JetBrains Mono', monospace; text-transform: uppercase; }}
        .d-field input {{
            background: var(--bg-panel); border: 1px solid var(--border-technical);
            border-radius: var(--radius-sm); padding: 5px 8px; color: var(--text-title);
            font-family: 'JetBrains Mono', monospace; font-size: 12px;
        }}

        /* Telemetry Instrument Matrix */
        .telemetry-rack {{
            display: grid; grid-template-columns: repeat(auto-fit, minmax(170px, 1fr));
            gap: 10px; margin-bottom: 16px;
        }}
        .rack-module {{
            background: var(--bg-subpanel); border: 1px solid var(--border-technical);
            border-radius: var(--radius-sm); padding: 10px 14px;
        }}
        .rack-label {{ font-size: 10px; color: var(--text-muted); text-transform: uppercase; font-family: 'JetBrains Mono', monospace; margin-bottom: 2px; }}
        .rack-value {{ font-size: 16px; font-weight: 700; font-family: 'JetBrains Mono', monospace; color: var(--text-title); }}

        /* Degradation Warning Console */
        .degradation-console {{
            background: rgba(244,63,94,0.08); border: 1px solid var(--accent-rose);
            border-left: 4px solid var(--accent-rose); border-radius: var(--radius-sm);
            padding: 12px 16px; margin-bottom: 16px; font-family: 'JetBrains Mono', monospace;
        }}
        .degradation-console.hidden {{ display: none; }}
        .deg-tag {{ font-size: 11px; font-weight: 800; color: var(--accent-rose); margin-bottom: 2px; }}
        .deg-text {{ font-size: 12px; color: #fecdd3; line-height: 1.4; }}

        /* Scientific Dossier Box */
        .dossier-box {{
            background: var(--bg-subpanel); border: 1px solid var(--border-technical);
            border-radius: var(--radius-sm); padding: 12px 16px; font-size: 12px;
            line-height: 1.6; color: var(--text-body);
        }}
        .dossier-box strong {{ color: var(--accent-cyan); }}

        /* ── SECTION DASHBOARD ── */
        .dashboard-grid {{
            display: grid; grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
            gap: 16px; margin-bottom: 30px;
        }}
        .panel-card {{
            background: var(--bg-panel); border: 1px solid var(--border-technical);
            border-radius: var(--radius); padding: 18px;
        }}
        .panel-title {{
            font-size: 11px; font-weight: 700; color: var(--accent-cyan);
            letter-spacing: 1px; text-transform: uppercase; font-family: 'JetBrains Mono', monospace;
            margin-bottom: 12px; padding-bottom: 8px; border-bottom: 1px solid var(--border-technical);
        }}

        .data-table {{ width: 100%; border-collapse: collapse; }}
        .data-table td {{ padding: 6px 0; font-size: 12px; border-bottom: 1px solid rgba(29,36,51,0.6); }}
        .data-table td:first-child {{ color: var(--text-muted); }}
        .data-table td:last-child {{ color: var(--text-title); font-family: 'JetBrains Mono', monospace; font-weight: 600; text-align: right; }}

        .metric-bar-group {{ margin-bottom: 10px; }}
        .mb-head {{ display: flex; justify-content: space-between; font-size: 11px; margin-bottom: 3px; font-family: 'JetBrains Mono', monospace; }}
        .mb-track {{ height: 5px; background: #1a2233; border-radius: 2px; overflow: hidden; }}
        .mb-fill {{ height: 100%; width: 0; transition: width 0.8s ease; }}
        .mb-fill.cyan {{ background: #06b6d4; }}
        .mb-fill.blue {{ background: #3b82f6; }}
        .mb-fill.purple {{ background: #8b5cf6; }}

        .tag-cloud {{ display: flex; flex-wrap: wrap; gap: 4px; }}
        .tag-item {{
            padding: 2px 7px; border-radius: 3px; font-size: 11px; font-family: 'JetBrains Mono', monospace;
            background: var(--bg-subpanel); border: 1px solid var(--border-technical); color: var(--text-title);
        }}
        .tag-item.glyco {{ border-color: rgba(6,182,212,0.4); color: var(--accent-cyan); }}
        .tag-item.disulf {{ border-color: rgba(245,158,11,0.4); color: var(--accent-amber); }}

        .footer {{
            text-align: center; padding: 24px; color: var(--text-muted);
            font-size: 11px; font-family: 'JetBrains Mono', monospace;
            border-top: 1px solid var(--border-technical); margin-top: 40px;
        }}
        .footer a {{ color: var(--accent-cyan); text-decoration: none; }}
    </style>
</head>
<body>
    <div class="container">
        <!-- WORKSTATION TOP BAR -->
        <header class="workstation-topbar">
            <div class="brand-block">
                <div class="brand-logo">PS</div>
                <div class="brand-meta">
                    <span class="brand-sys">ProteinScope // Structural Biophysics Workstation</span>
                    <span class="brand-sub">Biophysical Model · Alberts et al. (NBK26830)</span>
                </div>
            </div>
            <div class="workstation-telemetry-pill">
                <div class="w-item"><span class="w-lbl">ACC:</span><span class="w-val">{accession}</span></div>
                <div class="w-item"><span class="w-lbl">SOURCE:</span><span class="w-val">{struct_source.split(' ')[0]}</span></div>
                <div class="w-item"><span class="w-lbl">RMSD:</span><span class="w-val">{rmsd if isinstance(rmsd, str) else f"{rmsd} Å"}</span></div>
                <div class="w-item"><span class="w-lbl">EXEC:</span><span class="w-val">{duration}s</span></div>
            </div>
        </header>

        <!-- HERO HEADER -->
        <div class="hero-banner">
            <div class="hero-title-row">
                <div>
                    <span class="hero-accession">{accession}</span>
                    <h1 class="hero-title">{protein_name}</h1>
                    <div class="hero-organism">{organism}</div>
                </div>
            </div>
            <div class="hero-chips">
                <div class="h-chip"><span class="lbl">Sequence Length</span><span class="val">{seq_len:,} aa</span></div>
                <div class="h-chip"><span class="lbl">Molecular Mass</span><span class="val">{mw:,.2f} Da</span></div>
                <div class="h-chip"><span class="lbl">Isoelectric Point</span><span class="val">pI {pi}</span></div>
                <div class="h-chip"><span class="lbl">Instability Index</span><span class="val">{instab} ({'Stable' if is_stable else 'Unstable'})</span></div>
                <div class="h-chip"><span class="lbl">Estimated Tm</span><span class="val" style="color:var(--accent-amber);">{tm_est:.1f}°C</span></div>
            </div>
        </div>

        <!-- 3D MOLECULAR GRAPHICS WORKSTATION -->
        <div class="viewer-card" id="main-viewer-card">
            <div class="viewer-toolbar">
                <div class="toolbar-title">3D Coordinate Engine // {struct_source}</div>
                <div style="display:flex; gap:8px; align-items:center; flex-wrap:wrap;">
                    <div class="btn-group" id="source-tabs-container">
                        {tabs_html}
                    </div>
                    <div class="btn-group">
                        <button class="btn active" onclick="setStyle('cartoon')">Cartoon</button>
                        <button class="btn" onclick="setStyle('stick')">Stick</button>
                        <button class="btn" onclick="setStyle('sphere')">Sphere</button>
                        <button class="btn" onclick="setStyle('surface')">Surface</button>
                    </div>
                    <div class="btn-group">
                        <button class="btn active" onclick="colorBy('spectrum')">Rainbow</button>
                        <button class="btn" onclick="colorBy('chain')">Chain</button>
                        <button class="btn" onclick="colorBy('ss')">Sec. Struct</button>
                    </div>
                </div>
            </div>
            <div class="viewer-container">
                <div id="viewer-3d"></div>
                <div class="hud-reticle hud-tl">[+]</div>
                <div class="hud-reticle hud-tr">[+]</div>
                <div class="hud-reticle hud-bl">[+]</div>
                <div class="hud-reticle hud-br">[+]</div>

                <div class="loading-overlay" id="loading">
                    <div class="spinner"></div>
                    <div style="color:var(--text-muted); font-size:12px; font-family:'JetBrains Mono',monospace;">STREAMING PDB GEOMETRY...</div>
                </div>
                <div class="viewer-info-overlay" id="viewer-info">
                    <div>MODEL: <strong id="active-struct-name">{struct_source}</strong></div>
                    <div>STATE: <span id="viewer-state-badge" style="color:{init_color}; font-weight:700;">{init_badge}</span></div>
                    <button class="btn" style="padding:2px 8px; font-size:10px; display:none;" id="renature-btn" onclick="applyPreset('phys')">RESTORE STANDARD STATE</button>
                </div>
            </div>
        </div>

        <!-- ENVIRONMENTAL STABILITY & CONFORMATION CALIBRATION CHAMBER -->
        <div class="env-card">
            <div class="env-header">
                <div class="env-tagline">
                    <div class="env-title">01 // Environmental Stability & Conformation Dynamics Chamber</div>
                    <div class="env-subtitle">Thermodynamics, net charge titration, and degradation boundaries (Alberts NBK26830)</div>
                </div>
                <div class="env-status-lamp">
                    <span class="lamp-dot" id="lamp-dot" style="background:{init_color};"></span>
                    <span id="conformation-badge" style="color:{init_color};">{init_badge}</span>
                </div>
            </div>

            <!-- Standard Reference Calibration Presets -->
            <div class="env-presets-strip">
                <button class="preset-chip active" id="preset-phys" onclick="applyPreset('phys')">37.0°C / pH 7.4 Physiological</button>
                <button class="preset-chip" id="preset-invitro" onclick="applyPreset('invitro')">25.0°C / pH 7.0 In-Vitro Lab</button>
                <button class="preset-chip" id="preset-cold" onclick="applyPreset('cold')">4.0°C Cryogenic Buffer</button>
                <button class="preset-chip" id="preset-stress" onclick="applyPreset('stress')">58.0°C Thermal Pre-Melt</button>
                <button class="preset-chip" id="preset-heat-deg" onclick="applyPreset('heat_deg')">85.0°C Thermal Aggregation</button>
                <button class="preset-chip" id="preset-acid-deg" onclick="applyPreset('acid_deg')">37.0°C / pH 2.0 Acid Hydrolysis</button>
                <button class="preset-chip" id="preset-alk-deg" onclick="applyPreset('alk_deg')">37.0°C / pH 12.5 Base Cleavage</button>
            </div>

            <!-- Dual Calibrated Environmental Sliders -->
            <div class="instrument-dials-grid">
                <!-- Temperature Module -->
                <div class="dial-box">
                    <div class="dial-top">
                        <span class="dial-label">Temperature [T]</span>
                        <div style="text-align:right;">
                            <span class="dial-readout" id="temp-display" style="color:var(--accent-amber);">{init_temp:.1f} °C</span>
                            <div style="font-size:11px; color:var(--text-muted); font-family:'JetBrains Mono',monospace;" id="temp-kelvin">{(init_temp + 273.15):.2f} K</div>
                        </div>
                    </div>
                    <div class="slider-bracket">
                        <input type="range" class="env-slider" id="temp-slider" min="0" max="100" step="0.5" value="{init_temp}" oninput="onEnvInput()">
                    </div>
                    <div class="graduations">
                        <span>0°C (Ice)</span>
                        <span>25°C Lab</span>
                        <span>37°C Core</span>
                        <span>Tm ~{tm_est:.0f}°C</span>
                        <span>100°C Boil</span>
                    </div>
                </div>

                <!-- pH Level Module -->
                <div class="dial-box">
                    <div class="dial-top">
                        <span class="dial-label">Proton Activity [pH]</span>
                        <div style="text-align:right;">
                            <span class="dial-readout" id="ph-display" style="color:var(--accent-cyan);">pH {init_ph:.2f}</span>
                            <div style="font-size:11px; color:var(--text-muted); font-family:'JetBrains Mono',monospace;" id="ph-status-pill">Physiological Neutral</div>
                        </div>
                    </div>
                    <div class="slider-bracket">
                        <input type="range" class="env-slider" id="ph-slider" min="1.0" max="14.0" step="0.1" value="{init_ph}" oninput="onEnvInput()">
                    </div>
                    <div class="graduations">
                        <span>pH 1.0 (Gastric)</span>
                        <span>pH 4.0</span>
                        <span>pH 7.4 (Blood)</span>
                        <span>pH 10.0</span>
                        <span>pH 14.0 (Alkaline)</span>
                    </div>
                </div>
            </div>

            <!-- Range Boundary Tuning Drawer -->
            <div class="drawer-toggle" onclick="toggleTuningPanel()">
                <span id="tuning-arrow">▸</span>
                <span>Configure Degradation & Stability Threshold Cutoffs</span>
            </div>
            <div class="drawer-panel" id="tuning-panel" style="display:none;">
                <div style="font-size:11px; color:var(--text-muted); margin-bottom:10px; font-family:'JetBrains Mono',monospace;">
                    // CUSTOM THRESHOLD CALIBRATION (MOVING PAST CUTOFF TRIGGERS STRUCTURAL DEGRADATION):
                </div>
                <div class="drawer-grid">
                    <div class="d-field">
                        <label>Min Stable Temp (°C)</label>
                        <input type="number" id="cfg-temp-min" value="{thresholds['temp_min_stable']}" step="1" onchange="onBoundaryChange()">
                    </div>
                    <div class="d-field">
                        <label>Max Stable Temp (°C)</label>
                        <input type="number" id="cfg-temp-max" value="{thresholds['temp_max_stable']}" step="1" onchange="onBoundaryChange()">
                    </div>
                    <div class="d-field">
                        <label>Heat Degrade Temp (°C)</label>
                        <input type="number" id="cfg-temp-deg" value="{thresholds['temp_degradation']}" step="1" onchange="onBoundaryChange()">
                    </div>
                    <div class="d-field">
                        <label>Min Stable pH</label>
                        <input type="number" id="cfg-ph-min" value="{thresholds['ph_min_stable']}" step="0.1" onchange="onBoundaryChange()">
                    </div>
                    <div class="d-field">
                        <label>Max Stable pH</label>
                        <input type="number" id="cfg-ph-max" value="{thresholds['ph_max_stable']}" step="0.1" onchange="onBoundaryChange()">
                    </div>
                    <div class="d-field">
                        <label>Acid Degrade pH</label>
                        <input type="number" id="cfg-ph-acid-deg" value="{thresholds['ph_acid_degradation']}" step="0.1" onchange="onBoundaryChange()">
                    </div>
                    <div class="d-field">
                        <label>Base Degrade pH</label>
                        <input type="number" id="cfg-ph-alk-deg" value="{thresholds['ph_alkaline_degradation']}" step="0.1" onchange="onBoundaryChange()">
                    </div>
                </div>
                <button class="btn" onclick="resetBoundaries()">Reset Thresholds to Sequence Defaults</button>
            </div>

            <!-- Degradation Fail-Safe Banner -->
            <div class="degradation-console hidden" id="degradation-banner">
                <div class="deg-tag">[CRITICAL ALERT] PROTEIN DEGRADATION DETECTED</div>
                <div class="deg-text" id="degradation-desc">Structure collapsed into irreversible insoluble aggregates or underwent acid/base cleavage.</div>
            </div>

            <!-- Real-time Instrument Telemetry Rack -->
            <div class="telemetry-rack">
                <div class="rack-module">
                    <div class="rack-label">Fraction Folded [f]</div>
                    <div class="rack-value" id="t-folded" style="color:var(--accent-emerald);">{init_folded_pct:.1f}%</div>
                </div>
                <div class="rack-module">
                    <div class="rack-label">Net Charge [Q(pH)]</div>
                    <div class="rack-value" id="t-charge">{init_q:+.1f} e</div>
                </div>
                <div class="rack-module">
                    <div class="rack-label">Salt Bridges [Intact]</div>
                    <div class="rack-value" id="t-salt" style="color:var(--accent-cyan);">{init_salt:.1f}%</div>
                </div>
                <div class="rack-module">
                    <div class="rack-label">Free Energy [ΔG_fold]</div>
                    <div class="rack-value" id="t-dg">{init_dg:.1f} kcal/mol</div>
                </div>
                <div class="rack-module">
                    <div class="rack-label">Melting Temp [Tm]</div>
                    <div class="rack-value" id="t-tm" style="color:var(--accent-amber);">{tm_est:.1f}°C</div>
                </div>
            </div>

            <!-- Molecular Mechanism Dossier -->
            <div class="dossier-box" id="mechanism-text">
                {init_mechanism}
            </div>
        </div>

        <!-- SCIENTIFIC ANALYSIS PANELS -->
        <div class="dashboard-grid">
            {patho_html}

            <div class="panel-card">
                <div class="panel-title">02 // Physicochemical Characterization</div>
                <table class="data-table">
                    <tr><td>Molecular Mass</td><td>{mw:,.2f} Da</td></tr>
                    <tr><td>Isoelectric Point (pI)</td><td>{pi}</td></tr>
                    <tr><td>Instability Index</td><td>{instab}</td></tr>
                    <tr><td>GRAVY Hydropathicity</td><td>{gravy}</td></tr>
                    <tr><td>Aromaticity Index</td><td>{aroma}</td></tr>
                    <tr><td>Cationic (+) Residues</td><td>{pos_pct}%</td></tr>
                    <tr><td>Anionic (−) Residues</td><td>{neg_pct}%</td></tr>
                    <tr><td>Hydrophobic Core</td><td>{hydro_pct}%</td></tr>
                </table>
            </div>

            <div class="panel-card">
                <div class="panel-title">03 // Secondary Structure Decomposition</div>
                <div class="metric-bar-group">
                    <div class="mb-head"><span>α-Helix</span><span style="color:var(--accent-cyan);">{helix_pct}%</span></div>
                    <div class="mb-track"><div class="mb-fill cyan" data-width="{helix_pct}"></div></div>
                </div>
                <div class="metric-bar-group">
                    <div class="mb-head"><span>β-Sheet</span><span style="color:var(--accent-blue);">{sheet_pct}%</span></div>
                    <div class="mb-track"><div class="mb-fill blue" data-width="{sheet_pct}"></div></div>
                </div>
                <div class="metric-bar-group">
                    <div class="mb-head"><span>Turn / Coil</span><span style="color:var(--accent-amber);">{turn_pct}%</span></div>
                    <div class="mb-track"><div class="mb-fill purple" data-width="{turn_pct}"></div></div>
                </div>
            </div>

            <div class="panel-card">
                <div class="panel-title">04 // Stereochemical Ramachandran Validation</div>
                <table class="data-table" style="margin-bottom:12px;">
                    <tr><td>Core Favored Dihedrals</td><td style="color:var(--accent-emerald);">{fav_pct}%</td></tr>
                    <tr><td>Allowed Dihedrals</td><td>{allow_pct}%</td></tr>
                    <tr><td>Conformational Outliers</td><td>{out_pct}%</td></tr>
                    <tr><td>CA–CB Geometry Pass Rate</td><td style="color:var(--accent-emerald);">{ca_cb_pass}%</td></tr>
                    <tr><td>Coordinate Superimposition</td><td>{rmsd if isinstance(rmsd, str) else f"{rmsd} Å"}</td></tr>
                </table>
            </div>

            <div class="panel-card">
                <div class="panel-title">05 // Post-Translational Modification Matrix</div>
                <div style="font-size:10px; color:var(--text-muted); text-transform:uppercase; font-family:'JetBrains Mono',monospace; margin-bottom:4px;">Glycosylation Residues</div>
                <div class="tag-cloud" style="margin-bottom:12px;">{glyco_tags}</div>
                <div style="font-size:10px; color:var(--text-muted); text-transform:uppercase; font-family:'JetBrains Mono',monospace; margin-bottom:4px;">Disulfide Linkages / Modifications</div>
                <div class="tag-cloud">{disulf_tags}</div>
            </div>

            <div class="panel-card">
                <div class="panel-title">06 // Machine Learning Vector Importances</div>
                {ml_rows if ml_rows else '<div style="color:var(--text-muted); font-size:11px; font-family:\'JetBrains Mono\',monospace;">NO VARIANT DESCRIPTORS RECORDED</div>'}
            </div>
        </div>

        <footer class="footer">
            PROTEINSCOPE BIOINFORMATICS PIPELINE · UNIPROT <a href="https://www.uniprot.org/uniprotkb/{accession}" target="_blank">{accession} ↗</a> · MOLECULAR BIOLOGY OF THE CELL (ALBERTS NBK26830)
        </footer>
    </div>

    <script>
        let viewer;
        let currentStyle = 'cartoon';
        let currentColor = 'spectrum';
        const structModels = {models_json};
        const structLabels = {labels_json};
        let activeKey = "{default_model_key}";

        const titratable = {counts_json};
        const defaultThresholds = {thresholds_json};
        let userThresholds = Object.assign({{}}, defaultThresholds);
        let currentConformation = "{init_state}";

        function initViewer() {{
            const el = document.getElementById('viewer-3d');
            viewer = $3Dmol.createViewer(el, {{ backgroundColor: '#07090e', antialias: true }});
            loadSource(activeKey);
        }}

        function loadSource(sourceKey) {{
            activeKey = sourceKey;
            const loading = document.getElementById('loading');
            const infoLabel = document.getElementById('active-struct-name');
            loading.classList.remove('hidden');
            viewer.clear();

            document.querySelectorAll('.source-tab').forEach(t => t.classList.remove('active'));
            if (event && event.target && event.target.classList.contains('source-tab')) {{
                event.target.classList.add('active');
            }}

            if (infoLabel && structLabels[sourceKey]) {{
                infoLabel.innerText = structLabels[sourceKey];
            }}

            if (structModels[sourceKey] && structModels[sourceKey].length > 50) {{
                viewer.addModel(structModels[sourceKey], "pdb");
                applyStyle();
                viewer.zoomTo();
                viewer.zoom(0.85);
                loading.classList.add('hidden');
            }} else {{
                loading.classList.add('hidden');
            }}
        }}

        function setStyle(style) {{
            currentStyle = style;
            document.querySelectorAll('.btn-group:nth-child(2) .btn').forEach(b => b.classList.remove('active'));
            if (event && event.target) event.target.classList.add('active');
            applyStyle();
        }}

        function colorBy(color) {{
            currentColor = color;
            document.querySelectorAll('.btn-group:nth-child(3) .btn').forEach(b => b.classList.remove('active'));
            if (event && event.target) event.target.classList.add('active');
            applyStyle();
        }}

        function applyStyle() {{
            if (!viewer) return;
            viewer.setStyle({{}}, {{}});
            let cs = {{ color: 'spectrum' }};
            if (currentColor === 'chain') cs = {{ colorscheme: 'chainHetatm' }};
            if (currentColor === 'ss') cs = {{ colorscheme: 'ssJmol' }};

            if (currentConformation === 'DEGRADED') {{
                viewer.setStyle({{}}, {{
                    sphere: {{ color: '#f43f5e', scale: 0.28, opacity: 0.6 }},
                    stick: {{ color: '#881337', radius: 0.08, opacity: 0.4 }}
                }});
            }} else if (currentConformation === 'DENATURED') {{
                viewer.setStyle({{}}, {{
                    cartoon: {{ style: 'trace', color: '#ea580c', opacity: 0.5 }},
                    stick: {{ color: '#c2410c', radius: 0.12, opacity: 0.55 }}
                }});
            }} else if (currentConformation === 'PERTURBED') {{
                if (currentStyle === 'cartoon') {{
                    viewer.setStyle({{}}, {{ cartoon: {{ ...cs, opacity: 0.85 }} }});
                }} else if (currentStyle === 'stick') {{
                    viewer.setStyle({{}}, {{ stick: {{ ...cs, radius: 0.18 }} }});
                }} else if (currentStyle === 'sphere') {{
                    viewer.setStyle({{}}, {{ sphere: {{ ...cs, scale: 0.3 }} }});
                }} else if (currentStyle === 'surface') {{
                    viewer.setStyle({{}}, {{ cartoon: {{ ...cs, opacity: 0.4 }} }});
                    viewer.addSurface($3Dmol.SurfaceType.VDW, {{ opacity: 0.6, color: '#f59e0b' }});
                }}
            }} else {{
                if (currentStyle === 'cartoon') viewer.setStyle({{}}, {{ cartoon: cs }});
                else if (currentStyle === 'stick') viewer.setStyle({{}}, {{ stick: cs }});
                else if (currentStyle === 'sphere') viewer.setStyle({{}}, {{ sphere: {{ ...cs, scale: 0.3 }} }});
                else if (currentStyle === 'surface') {{
                    viewer.setStyle({{}}, {{ cartoon: {{ ...cs, opacity: 0.4 }} }});
                    viewer.addSurface($3Dmol.SurfaceType.VDW, {{ opacity: 0.7, ...cs }});
                }}
            }}
            viewer.render();
        }}

        function computeNetCharge(ph) {{
            const qNterm = 1.0 / (1.0 + Math.pow(10, ph - 9.69));
            const qCterm = -1.0 / (1.0 + Math.pow(10, 2.34 - ph));
            const qK = (titratable.K || 0) / (1.0 + Math.pow(10, ph - 10.53));
            const qR = (titratable.R || 0) / (1.0 + Math.pow(10, ph - 12.48));
            const qH = (titratable.H || 0) / (1.0 + Math.pow(10, ph - 6.00));
            const qD = -(titratable.D || 0) / (1.0 + Math.pow(10, 3.65 - ph));
            const qE = -(titratable.E || 0) / (1.0 + Math.pow(10, 4.25 - ph));
            const qC = -(titratable.C || 0) / (1.0 + Math.pow(10, 8.18 - ph));
            const qY = -(titratable.Y || 0) / (1.0 + Math.pow(10, 10.07 - ph));
            return qNterm + qCterm + qK + qR + qH + qD + qE + qC + qY;
        }}

        function computeSaltBridges(ph) {{
            const fAcid = 1.0 / (1.0 + Math.pow(10, 4.0 - ph));
            const fBase = 1.0 / (1.0 + Math.pow(10, ph - 10.5));
            return Math.max(0, Math.min(100, fAcid * fBase * 100));
        }}

        function computeThermodynamics(T_celsius, ph, netCharge) {{
            const T_kelvin = T_celsius + 273.15;
            const Tm_celsius = userThresholds.tm_celsius || 62.0;
            const Tm_kelvin = Tm_celsius + 273.15;
            const R_gas = 1.9872e-3;
            const N_res = titratable.len || 300;

            const deltaHm = 1.15 * Math.min(500, Math.max(50, N_res));
            const deltaCp = 0.012 * Math.min(500, Math.max(50, N_res));

            const term1 = deltaHm * (1.0 - (T_kelvin / Tm_kelvin));
            const term2 = deltaCp * ((Tm_kelvin - T_kelvin) + T_kelvin * Math.log(T_kelvin / Tm_kelvin));
            const deltaG_thermal = -(term1 - term2);
            const chargeDestab = 0.008 * Math.pow(Math.abs(netCharge), 1.8);
            const deltaG_total = deltaG_thermal + chargeDestab;

            const expVal = Math.max(-30, Math.min(30, deltaG_total / (R_gas * T_kelvin)));
            const f_folded = 1.0 / (1.0 + Math.exp(expVal));
            return {{ deltaG: deltaG_total, fFolded: f_folded }};
        }}

        function onEnvInput() {{
            const temp = parseFloat(document.getElementById('temp-slider').value);
            const ph = parseFloat(document.getElementById('ph-slider').value);
            updateBiophysics(temp, ph);
        }}

        function updateBiophysics(temp, ph) {{
            document.getElementById('temp-display').innerText = temp.toFixed(1) + ' °C';
            document.getElementById('temp-kelvin').innerText = (temp + 273.15).toFixed(2) + ' K';
            document.getElementById('ph-display').innerText = 'pH ' + ph.toFixed(2);

            const netQ = computeNetCharge(ph);
            const saltPct = computeSaltBridges(ph);
            const thermo = computeThermodynamics(temp, ph, netQ);
            const fFoldedPct = thermo.fFolded * 100;

            const isHeatDegraded = temp >= userThresholds.temp_degradation;
            const isAcidDegraded = ph <= userThresholds.ph_acid_degradation;
            const isAlkDegraded = ph >= userThresholds.ph_alkaline_degradation;
            const isColdDenatured = temp < 4.0;
            const isDegraded = isHeatDegraded || isAcidDegraded || isAlkDegraded;

            let state = 'NATIVE';
            let badgeLabel = 'Native Folded State';
            let badgeColor = '#10b981';
            let mechanism = '';

            if (isDegraded) {{
                state = 'DEGRADED';
                badgeLabel = 'Degraded / Irreversible Aggregation';
                badgeColor = '#f43f5e';
                let degMsg = [];
                if (isHeatDegraded) degMsg.push(`Thermal Aggregation: Kinetic energy (${{temp.toFixed(1)}}°C ≥ ${{userThresholds.temp_degradation}}°C) has completely dismantled backbone hydrogen bonding. Hydrophobic residues form insoluble aggregates (Alberts NBK26830).`);
                if (isAcidDegraded) degMsg.push(`Acid Hydrolysis: pH ${{ph.toFixed(1)}} ≤ ${{userThresholds.ph_acid_degradation}} protonates Asp/Glu, breaking salt bridges with strong positive repulsion (+${{netQ.toFixed(1)}} e), catalyzing peptide cleavage.`);
                if (isAlkDegraded) degMsg.push(`Alkaline Degradation: pH ${{ph.toFixed(1)}} ≥ ${{userThresholds.ph_alkaline_degradation}} deprotonates basic residues and hydrolyzes disulfides via β-elimination, destroying covalent integrity.`);
                mechanism = degMsg.join(' ');
                document.getElementById('degradation-desc').innerText = mechanism;
                document.getElementById('degradation-banner').classList.remove('hidden');
                document.getElementById('main-viewer-card').classList.add('degraded-glow');
                document.getElementById('renature-btn').style.display = 'inline-block';
            }} else if (thermo.fFolded < 0.40 || isColdDenatured || ph < (userThresholds.ph_min_stable - 1.5) || ph > (userThresholds.ph_max_stable + 1.5)) {{
                state = 'DENATURED';
                badgeLabel = 'Denatured (Unfolded Chain)';
                badgeColor = '#ea580c';
                mechanism = `Severe stress (${{temp.toFixed(1)}}°C, pH ${{ph.toFixed(1)}}) has disrupted tertiary architecture into a flexible random coil (Fraction folded: ${{fFoldedPct.toFixed(1)}}%, Salt bridges: ${{saltPct.toFixed(1)}}%).`;
                document.getElementById('degradation-banner').classList.add('hidden');
                document.getElementById('main-viewer-card').classList.remove('degraded-glow');
                document.getElementById('renature-btn').style.display = 'inline-block';
            }} else if (thermo.fFolded < 0.85 || temp > userThresholds.temp_max_stable || ph < userThresholds.ph_min_stable || ph > userThresholds.ph_max_stable) {{
                state = 'PERTURBED';
                badgeLabel = 'Perturbed / Molten Globule';
                badgeColor = '#f59e0b';
                mechanism = `Moderate environmental stress (${{temp.toFixed(1)}}°C, pH ${{ph.toFixed(1)}}) weakens weak hydrogen bonds and loop interactions. The protein adopts a fluctuating 'molten globule' intermediate with preserved secondary structure but loosened packing.`;
                document.getElementById('degradation-banner').classList.add('hidden');
                document.getElementById('main-viewer-card').classList.remove('degraded-glow');
                document.getElementById('renature-btn').style.display = 'none';
            }} else {{
                state = 'NATIVE';
                badgeLabel = 'Native Folded State';
                badgeColor = '#10b981';
                mechanism = `Under standard conditions (${{temp.toFixed(1)}}°C, pH ${{ph.toFixed(1)}}), the protein resides in its lowest free energy conformation (ΔG = ${{thermo.deltaG.toFixed(1)}} kcal/mol). Backbone hydrogen bonds, buried hydrophobic core, and salt bridges (${{saltPct.toFixed(1)}}% intact) maintain native function (Alberts NBK26830).`;
                document.getElementById('degradation-banner').classList.add('hidden');
                document.getElementById('main-viewer-card').classList.remove('degraded-glow');
                document.getElementById('renature-btn').style.display = 'none';
            }}

            currentConformation = state;

            const confBadge = document.getElementById('conformation-badge');
            confBadge.innerText = badgeLabel;
            confBadge.style.color = badgeColor;
            document.getElementById('lamp-dot').style.background = badgeColor;

            const viewerBadge = document.getElementById('viewer-state-badge');
            if (viewerBadge) {{
                viewerBadge.innerText = badgeLabel;
                viewerBadge.style.color = badgeColor;
            }}

            document.getElementById('t-folded').innerText = fFoldedPct.toFixed(1) + '%';
            document.getElementById('t-folded').style.color = (fFoldedPct > 80) ? '#10b981' : (fFoldedPct > 40 ? '#f59e0b' : '#f43f5e');
            document.getElementById('t-charge').innerText = (netQ > 0 ? '+' : '') + netQ.toFixed(1) + ' e';
            document.getElementById('t-salt').innerText = saltPct.toFixed(1) + '%';
            document.getElementById('t-dg').innerText = thermo.deltaG.toFixed(1) + ' kcal/mol';
            document.getElementById('mechanism-text').innerHTML = mechanism;

            const phPill = document.getElementById('ph-status-pill');
            if (ph <= userThresholds.ph_acid_degradation) {{
                phPill.innerText = 'Acid Hydrolysis';
                phPill.style.color = '#f43f5e';
            }} else if (ph >= userThresholds.ph_alkaline_degradation) {{
                phPill.innerText = 'Alkaline Cleavage';
                phPill.style.color = '#f43f5e';
            }} else if (ph < userThresholds.ph_min_stable || ph > userThresholds.ph_max_stable) {{
                phPill.innerText = 'Perturbed Buffer';
                phPill.style.color = '#f59e0b';
            }} else {{
                phPill.innerText = 'Physiological Neutral';
                phPill.style.color = '#10b981';
            }}

            applyStyle();
        }}

        function applyPreset(presetName) {{
            document.querySelectorAll('.preset-chip').forEach(b => b.classList.remove('active'));
            let t = 37.0, p = 7.4;
            if (presetName === 'phys') {{
                t = 37.0; p = 7.4;
                document.getElementById('preset-phys').classList.add('active');
            }} else if (presetName === 'invitro') {{
                t = 25.0; p = 7.0;
                document.getElementById('preset-invitro').classList.add('active');
            }} else if (presetName === 'cold') {{
                t = 4.0; p = 7.4;
                document.getElementById('preset-cold').classList.add('active');
            }} else if (presetName === 'stress') {{
                t = 58.0; p = 7.4;
                document.getElementById('preset-stress').classList.add('active');
            }} else if (presetName === 'heat_deg') {{
                t = 85.0; p = 7.4;
                document.getElementById('preset-heat-deg').classList.add('active');
            }} else if (presetName === 'acid_deg') {{
                t = 37.0; p = 2.0;
                document.getElementById('preset-acid-deg').classList.add('active');
            }} else if (presetName === 'alk_deg') {{
                t = 37.0; p = 12.5;
                document.getElementById('preset-alk-deg').classList.add('active');
            }}

            document.getElementById('temp-slider').value = t;
            document.getElementById('ph-slider').value = p;
            updateBiophysics(t, p);
        }}

        function toggleTuningPanel() {{
            const panel = document.getElementById('tuning-panel');
            const arrow = document.getElementById('tuning-arrow');
            if (panel.style.display === 'none') {{
                panel.style.display = 'block';
                arrow.innerText = '▾';
            }} else {{
                panel.style.display = 'none';
                arrow.innerText = '▸';
            }}
        }}

        function onBoundaryChange() {{
            userThresholds.temp_min_stable = parseFloat(document.getElementById('cfg-temp-min').value) || 15.0;
            userThresholds.temp_max_stable = parseFloat(document.getElementById('cfg-temp-max').value) || 45.0;
            userThresholds.temp_degradation = parseFloat(document.getElementById('cfg-temp-deg').value) || 70.0;
            userThresholds.ph_min_stable = parseFloat(document.getElementById('cfg-ph-min').value) || 5.5;
            userThresholds.ph_max_stable = parseFloat(document.getElementById('cfg-ph-max').value) || 8.5;
            userThresholds.ph_acid_degradation = parseFloat(document.getElementById('cfg-ph-acid-deg').value) || 3.0;
            userThresholds.ph_alkaline_degradation = parseFloat(document.getElementById('cfg-ph-alk-deg').value) || 11.5;
            onEnvInput();
        }}

        function resetBoundaries() {{
            userThresholds = Object.assign({{}}, defaultThresholds);
            document.getElementById('cfg-temp-min').value = userThresholds.temp_min_stable;
            document.getElementById('cfg-temp-max').value = userThresholds.temp_max_stable;
            document.getElementById('cfg-temp-deg').value = userThresholds.temp_degradation;
            document.getElementById('cfg-ph-min').value = userThresholds.ph_min_stable;
            document.getElementById('cfg-ph-max').value = userThresholds.ph_max_stable;
            document.getElementById('cfg-ph-acid-deg').value = userThresholds.ph_acid_degradation;
            document.getElementById('cfg-ph-alk-deg').value = userThresholds.ph_alkaline_degradation;
            onEnvInput();
        }}

        window.addEventListener('DOMContentLoaded', function() {{
            document.querySelectorAll('.mb-fill').forEach(bar => {{
                bar.style.width = (bar.dataset.width || 0) + '%';
            }});
            initViewer();
            onEnvInput();
        }});
    </script>
</body>
</html>
"""
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(html_content)
    logger.info(f"Generated standalone interactive 3D HTML report at {html_path}")
    return html_path



