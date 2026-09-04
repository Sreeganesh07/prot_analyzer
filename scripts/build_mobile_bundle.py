"""
Bundle the 9 core proteins with their PDB structures and biophysical manifests
into a self-contained JavaScript module for the mobile application.
"""

import json
import os

OUTPUT_DIR = os.path.abspath("output")
MOBILE_JS_DIR = os.path.abspath(os.path.join("mobile_app", "js"))
os.makedirs(MOBILE_JS_DIR, exist_ok=True)

CORE_PROTEINS = [
    ("P04637", "Tumor Antigen p53", "Homo sapiens"),
    ("P0DTC2", "Spike Glycoprotein", "SARS-CoV-2"),
    ("P01308", "Insulin", "Homo sapiens"),
    ("P68871", "Hemoglobin Subunit Beta", "Homo sapiens"),
    ("P00533", "EGFR (Epidermal Growth Factor Receptor)", "Homo sapiens"),
    ("P00734", "Prothrombin (Thrombin)", "Homo sapiens"),
    ("P42212", "Green Fluorescent Protein (GFP)", "Aequorea victoria"),
    ("P11021", "BiP Endoplasmic Chaperone", "Homo sapiens"),
    ("P02769", "Serum Albumin", "Homo sapiens"),
]

def build_bundle():
    bundle = {}
    for acc, fallback_name, fallback_org in CORE_PROTEINS:
        rep_p = os.path.join(OUTPUT_DIR, acc, "report", "report.json")
        mf_p = os.path.join(OUTPUT_DIR, acc, "structures", "structures_manifest.json")

        if not os.path.exists(rep_p):
            print(f"Warning: {rep_p} not found")
            continue

        with open(rep_p, "r", encoding="utf-8") as f:
            data = json.load(f)

        # Load best PDB
        pdb_content = ""
        if os.path.exists(mf_p):
            with open(mf_p, "r", encoding="utf-8") as mf:
                mfd = json.load(mf)
                pdb_file = mfd.get("best_structure", {}).get("file_path", "")
                if pdb_file and os.path.exists(pdb_file):
                    with open(pdb_file, "r", encoding="utf-8", errors="ignore") as pf:
                        # Extract ATOM/HETATM/TER lines to keep bundle lightweight
                        lines = [l for l in pf.readlines() if l.startswith(("ATOM", "HETATM", "TER", "HELIX", "SHEET"))]
                        pdb_content = "".join(lines)

        stab = data.get("environmental_stability", {})
        active = stab.get("active_evaluation", {})
        physico = data.get("properties", {})

        # Count titratable amino acids
        seq = data.get("sequence", "")
        seq_u = seq.upper()
        titratable_counts = {
            "D": seq_u.count("D"),
            "E": seq_u.count("E"),
            "H": seq_u.count("H"),
            "C": seq_u.count("C"),
            "Y": seq_u.count("Y"),
            "K": seq_u.count("K"),
            "R": seq_u.count("R"),
            "len": len(seq_u),
        }

        bundle[acc] = {
            "accession": acc,
            "name": data.get("protein_name", fallback_name),
            "organism": data.get("organism", fallback_org),
            "sequence": seq,
            "sequence_length": len(seq),
            "titratable_counts": titratable_counts,
            "molecular_weight": physico.get("molecular_weight", 0),
            "isoelectric_point": physico.get("isoelectric_point", 7.0),
            "instability_index": physico.get("instability_index", 40.0),
            "estimated_tm": stab.get("estimated_melting_temperature_celsius", 65.0),
            "active_evaluation": active,
            "thresholds": active.get("thresholds", {
                "tm_celsius": stab.get("estimated_melting_temperature_celsius", 65.0),
                "temp_min_stable": 15.0,
                "temp_max_stable": 55.0,
                "temp_degradation": 75.0,
                "ph_min_stable": 5.5,
                "ph_max_stable": 8.5,
                "ph_acid_degradation": 3.0,
                "ph_alkaline_degradation": 11.5,
            }),
            "titration_profile": stab.get("titration_profile", []),
            "pdb_content": pdb_content,
        }
        print(f"Bundled {acc}: {bundle[acc]['name']} (PDB size: {len(pdb_content):,} chars)")

    out_js = os.path.join(MOBILE_JS_DIR, "cached_data.js")
    with open(out_js, "w", encoding="utf-8") as f:
        f.write("// ProteinScope Mobile Cached Offline Protein Database (9 Core Targets)\n")
        f.write("window.PROTEINSCOPE_CACHED_DB = ")
        json.dump(bundle, f, indent=2)
        f.write(";\n")
    print(f"\nSuccessfully generated {out_js} with {len(bundle)} proteins.")

if __name__ == "__main__":
    build_bundle()
