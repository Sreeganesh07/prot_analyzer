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

CORE_ANNOTATIONS = {
    "P04637": {
        "gene_name": "TP53",
        "gene_synonyms": ["P53", "BCC7", "LFS1"],
        "subcellular_location": ["Nucleus", "Cytoplasm", "Mitochondrion"],
        "function_summary": "Master cellular tumor suppressor. Induces cell cycle arrest, DNA repair, senescence, or apoptosis in response to diverse cellular stresses. Acts as a tetrameric sequence-specific transcription factor regulating hundreds of stress-response genes.",
        "disease_associations": "Li-Fraumeni syndrome (LFS), Adrenocortical carcinoma, Choroid plexus papilloma, and over 50% of all human somatic cancers.",
        "pdb_cross_references": ["1TUP", "1OLG", "2OCJ", "3KMD"],
    },
    "P0DTC2": {
        "gene_name": "S",
        "gene_synonyms": ["Spike", "Surface Glycoprotein"],
        "subcellular_location": ["Virion membrane", "Host cell surface", "Endoplasmic reticulum-Golgi intermediate compartment"],
        "function_summary": "Trimeric class I viral fusion glycoprotein. Binds host ACE2 receptor via its Receptor-Binding Domain (RBD) and undergoes proteolytic cleavage at S1/S2 to drive viral-host membrane fusion and cellular entry.",
        "disease_associations": "Coronavirus disease 2019 (COVID-19), severe acute respiratory distress syndrome.",
        "pdb_cross_references": ["6VXX", "6VYB", "7KRR", "7C2L"],
    },
    "P01308": {
        "gene_name": "INS",
        "gene_synonyms": ["Insulin", "IRDN", "IDDM2"],
        "subcellular_location": ["Secreted extracellular space"],
        "function_summary": "Essential anabolic peptide hormone regulating carbohydrate, lipid, and protein metabolism. Promotes cellular glucose uptake via GLUT4 translocation, stimulates glycogenesis and lipogenesis, and inhibits hepatic gluconeogenesis.",
        "disease_associations": "Diabetes mellitus type 1 (IDDM), Permanent neonatal diabetes mellitus (PNDM), Hyperproinsulinemia.",
        "pdb_cross_references": ["4INS", "1TRZ", "2KQP", "3I40"],
    },
    "P68871": {
        "gene_name": "HBB",
        "gene_synonyms": ["Beta-Globin", "CD113t-C"],
        "subcellular_location": ["Erythrocyte cytoplasm"],
        "function_summary": "Forms the beta-subunit of heterotetrameric adult hemoglobin (alpha-2, beta-2). Binds four heme iron cofactors to cooperatively bind and transport oxygen from lungs to peripheral tissues, modulated by allosteric effectors like 2,3-BPG, pH (Bohr effect), and CO2.",
        "disease_associations": "Sickle cell anemia (Glu6Val mutation causing polymerization), Beta-thalassemia, Erythrocytosis.",
        "pdb_cross_references": ["4HHB", "1HHO", "2HHB", "1A3N"],
    },
    "P00533": {
        "gene_name": "EGFR",
        "gene_synonyms": ["ERBB", "ERBB1", "HER1"],
        "subcellular_location": ["Cell membrane", "Endosome", "Nucleus"],
        "function_summary": "Transmembrane receptor tyrosine kinase. Ligand binding (EGF, TGF-alpha) induces receptor homodimerization and autophosphorylation, activating Ras-Raf-MEK-ERK and PI3K-Akt cascades governing cell proliferation, survival, and differentiation.",
        "disease_associations": "Non-small cell lung cancer (NSCLC), Glioblastoma multiforme, Colorectal carcinoma, Squamous cell carcinoma.",
        "pdb_cross_references": ["1IVO", "2GS6", "3VJO", "4I22"],
    },
    "P00734": {
        "gene_name": "F2",
        "gene_synonyms": ["Prothrombin", "Coagulation Factor II", "Thrombin"],
        "subcellular_location": ["Secreted blood plasma"],
        "function_summary": "Central vitamin K-dependent serine protease of the blood coagulation cascade. Cleaved by Factor Xa into active thrombin, which proteolytically converts soluble fibrinogen into insoluble fibrin mesh and activates platelets via PAR receptors.",
        "disease_associations": "Prothrombin G20210A thrombophilia, Deep vein thrombosis, Dysprothrombinemia.",
        "pdb_cross_references": ["1HAP", "1PPB", "2HGT", "3KCG"],
    },
    "P42212": {
        "gene_name": "GFP",
        "gene_synonyms": ["Green Fluorescent Protein"],
        "subcellular_location": ["Cytoplasm"],
        "function_summary": "11-stranded beta-barrel fluorophore from Aequorea victoria. The internal Ser65-Tyr66-Gly67 tripeptide undergoes spontaneous post-translational cyclization and oxidation to yield an intensely fluorescent chromophore (excitation 395/475 nm, emission 509 nm).",
        "disease_associations": "Widely utilized biological reporter and biosensor in molecular medicine and protein dynamics research.",
        "pdb_cross_references": ["1EMA", "1GFL", "2B3P", "1QY3"],
    },
    "P11021": {
        "gene_name": "HSPA5",
        "gene_synonyms": ["BiP", "GRP78", "MIF2"],
        "subcellular_location": ["Endoplasmic reticulum lumen"],
        "function_summary": "Major Hsp70-family molecular chaperone in the ER lumen. Binds hydrophobic patches of nascent polypeptides to facilitate folding, prevent misfolded aggregation, and acts as the master sensor controlling IRE1, PERK, and ATF6 in the Unfolded Protein Response (UPR).",
        "disease_associations": "ER stress neurodegeneration, Tumor chemotherapeutic resistance adaptation, Neurodevelopmental disorders.",
        "pdb_cross_references": ["5E84", "3LDN", "6ASY", "3QFD"],
    },
    "P02769": {
        "gene_name": "ALB",
        "gene_synonyms": ["Albumin", "Serum Albumin"],
        "subcellular_location": ["Secreted blood plasma"],
        "function_summary": "Most abundant globular protein in human blood plasma (~50 g/L). Provides 80% of colloidal intravascular osmotic oncotic pressure and serves as the primary transport vehicle for non-esterified fatty acids, bilirubin, calcium, steroid hormones, and therapeutic drugs.",
        "disease_associations": "Analbuminemia, Familial dysalbuminemic hyperthyroxinemia, Hypoalbuminemia in liver failure.",
        "pdb_cross_references": ["1AO6", "1E78", "4F5S", "1UOR"],
    },
}

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

        # Extinction coefficient Pace et al.
        w_count = seq_u.count("W")
        y_count = seq_u.count("Y")
        c_count = seq_u.count("C")
        extinction_coef = (w_count * 5500) + (y_count * 1490) + (c_count * 125)

        annot = CORE_ANNOTATIONS.get(acc, {})

        bundle[acc] = {
            "accession": acc,
            "name": data.get("protein_name", fallback_name),
            "organism": data.get("organism", fallback_org),
            "gene_name": annot.get("gene_name", "UNKNOWN"),
            "gene_synonyms": annot.get("gene_synonyms", []),
            "subcellular_location": annot.get("subcellular_location", ["Cytoplasm"]),
            "function_summary": annot.get("function_summary", "Biological macromolecule investigated under physiological and thermal stress."),
            "disease_associations": annot.get("disease_associations", "No clinical pathology registered."),
            "pdb_cross_references": annot.get("pdb_cross_references", []),
            "extinction_coefficient": extinction_coef,
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
