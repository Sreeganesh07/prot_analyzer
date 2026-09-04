# ProteinScope Snakemake Workflow (Stage 10 Automation)
# Configuration loaded from config/config.yaml

configfile: "config/config.yaml"

ACCESSION = config.get("accession", "P0DTC2")
OUTDIR = config.get("pipeline", {}).get("default_output_dir", "./output")
ACC_OUT = f"{OUTDIR}/{ACCESSION}"

rule all:
    input:
        f"{ACC_OUT}/report/report.json",
        f"{ACC_OUT}/report/summary.txt"

rule validate_input:
    output:
        log=f"{ACC_OUT}/run.log"
    params:
        acc=ACCESSION,
        out=OUTDIR
    shell:
        """
        python -c "from proteinscope import ProteinAnalyzer; p = ProteinAnalyzer('{params.acc}', '{params.out}'); p.validate_input()"
        """

rule fetch_data:
    input:
        f"{ACC_OUT}/run.log"
    output:
        fasta=f"{ACC_OUT}/sequences/{ACCESSION}.fasta",
        raw_xml=f"{ACC_OUT}/raw_data/{ACCESSION}.xml",
        raw_json=f"{ACC_OUT}/raw_data/{ACCESSION}.json"
    params:
        acc=ACCESSION,
        out=OUTDIR
    shell:
        """
        python -c "from proteinscope import ProteinAnalyzer; p = ProteinAnalyzer('{params.acc}', '{params.out}'); p.fetch_from_uniprot()"
        """

rule sequence_analysis:
    input:
        fasta=f"{ACC_OUT}/sequences/{ACCESSION}.fasta"
    output:
        json=f"{ACC_OUT}/analysis/sequence_properties.json"
    params:
        acc=ACCESSION,
        out=OUTDIR
    shell:
        """
        python -c "from proteinscope import ProteinAnalyzer; p = ProteinAnalyzer('{params.acc}', '{params.out}'); p.analyze_sequence()"
        """

rule blast_search:
    input:
        fasta=f"{ACC_OUT}/sequences/{ACCESSION}.fasta"
    output:
        json=f"{ACC_OUT}/analysis/blast_results.json"
    params:
        acc=ACCESSION,
        out=OUTDIR
    shell:
        """
        python -c "from proteinscope import ProteinAnalyzer; p = ProteinAnalyzer('{params.acc}', '{params.out}'); p.blast_search(hits=10, run_online=False)"
        """

rule fetch_structures:
    input:
        f"{ACC_OUT}/sequences/{ACCESSION}.fasta"
    output:
        manifest=f"{ACC_OUT}/structures/structures_manifest.json"
    params:
        acc=ACCESSION,
        out=OUTDIR
    shell:
        """
        python -c "from proteinscope import ProteinAnalyzer; p = ProteinAnalyzer('{params.acc}', '{params.out}'); p.fetch_structures()"
        """

rule validate_structure:
    input:
        manifest=f"{ACC_OUT}/structures/structures_manifest.json"
    output:
        val_json=f"{ACC_OUT}/validation/structure_validation.json",
        ram_png=f"{ACC_OUT}/figures/ramachandran.png"
    params:
        acc=ACCESSION,
        out=OUTDIR
    shell:
        """
        python -c "from proteinscope import ProteinAnalyzer; p = ProteinAnalyzer('{params.acc}', '{params.out}'); p.fetch_structures(); p.validate_structure()"
        """

rule visualize:
    input:
        manifest=f"{ACC_OUT}/structures/structures_manifest.json",
        seq_json=f"{ACC_OUT}/analysis/sequence_properties.json"
    output:
        full_png=f"{ACC_OUT}/figures/structure_full.png",
        surf_png=f"{ACC_OUT}/figures/structure_surface.png",
        dom_png=f"{ACC_OUT}/figures/structure_domains.png",
        script=f"{ACC_OUT}/structures/render_script.pml"
    params:
        acc=ACCESSION,
        out=OUTDIR
    shell:
        """
        python -c "from proteinscope import ProteinAnalyzer; p = ProteinAnalyzer('{params.acc}', '{params.out}'); p.fetch_structures(); p.analyze_sequence(); p.visualize()"
        """

rule comparative:
    input:
        fasta=f"{ACC_OUT}/sequences/{ACCESSION}.fasta",
        blast_json=f"{ACC_OUT}/analysis/blast_results.json"
    output:
        comp_json=f"{ACC_OUT}/analysis/comparative.json",
        dist_csv=f"{ACC_OUT}/analysis/distance_matrix.csv"
    params:
        acc=ACCESSION,
        out=OUTDIR
    shell:
        """
        python -c "from proteinscope import ProteinAnalyzer; p = ProteinAnalyzer('{params.acc}', '{params.out}'); p.fetch_from_uniprot(); p.comparative_analysis()"
        """

rule ml_features:
    input:
        fasta=f"{ACC_OUT}/sequences/{ACCESSION}.fasta",
        blast_json=f"{ACC_OUT}/analysis/blast_results.json"
    output:
        features_csv=f"{ACC_OUT}/ml/feature_matrix.csv",
        clustering_json=f"{ACC_OUT}/ml/clustering_results.json"
    params:
        acc=ACCESSION,
        out=OUTDIR
    shell:
        """
        python -c "from proteinscope import ProteinAnalyzer; p = ProteinAnalyzer('{params.acc}', '{params.out}'); p.fetch_from_uniprot(); p.extract_ml_features()"
        """

rule ngs:
    output:
        ngs_json=f"{ACC_OUT}/ngs/ngs_summary.json"
    params:
        acc=ACCESSION,
        out=OUTDIR
    shell:
        """
        python -c "from proteinscope import ProteinAnalyzer; p = ProteinAnalyzer('{params.acc}', '{params.out}'); p.process_ngs()"
        """

rule compile_report:
    input:
        seq_json=f"{ACC_OUT}/analysis/sequence_properties.json",
        blast_json=f"{ACC_OUT}/analysis/blast_results.json",
        val_json=f"{ACC_OUT}/validation/structure_validation.json",
        comp_json=f"{ACC_OUT}/analysis/comparative.json",
        ml_json=f"{ACC_OUT}/ml/clustering_results.json",
        ngs_json=f"{ACC_OUT}/ngs/ngs_summary.json"
    output:
        rep_json=f"{ACC_OUT}/report/report.json",
        rep_txt=f"{ACC_OUT}/report/summary.txt"
    params:
        acc=ACCESSION,
        out=OUTDIR
    shell:
        """
        python -c "from proteinscope import ProteinAnalyzer; p = ProteinAnalyzer('{params.acc}', '{params.out}'); p.save_report()"
        """

# ── Pathogenicity Prediction Rules ──────────────────────────────────────

MUTATION = config.get("mutation", "R175H")

rule clinvar_data:
    output:
        tsv="data/clinvar/clinvar_tp53_variants.tsv"
    shell:
        """
        python -c "from proteinscope.stages.clinvar import load_clinvar_variants; load_clinvar_variants()"
        """

rule train_classifier:
    input:
        tsv="data/clinvar/clinvar_tp53_variants.tsv",
        fasta=f"{ACC_OUT}/sequences/{ACCESSION}.fasta",
        blast_json=f"{ACC_OUT}/analysis/blast_results.json",
        val_json=f"{ACC_OUT}/validation/structure_validation.json"
    output:
        model="data/models/tp53_pathogenicity_rf.joblib"
    params:
        acc=ACCESSION,
        out=OUTDIR
    shell:
        """
        python -c "
from proteinscope import ProteinAnalyzer
p = ProteinAnalyzer('{params.acc}', '{params.out}', mutation='R175H')
p.fetch_from_uniprot()
p.blast_search(run_online=False)
p.fetch_structures()
p.validate_structure()
p.comparative_analysis()
p.predict_pathogenicity(train=True)
"
        """

rule predict_variant:
    input:
        model="data/models/tp53_pathogenicity_rf.joblib",
        fasta=f"{ACC_OUT}/sequences/{ACCESSION}.fasta"
    output:
        rep_json=f"{ACC_OUT}/report/pathogenicity_report.json"
    params:
        acc=ACCESSION,
        out=OUTDIR,
        mut=MUTATION
    shell:
        """
        python main.py {params.acc} --mutation {params.mut} --offline-blast --output {params.out}
        """

rule benchmark:
    input:
        tsv="data/clinvar/clinvar_tp53_variants.tsv",
        fasta=f"{ACC_OUT}/sequences/{ACCESSION}.fasta"
    output:
        metrics=f"{ACC_OUT}/benchmark/metrics.json",
        report=f"{ACC_OUT}/benchmark/benchmark_report.txt"
    params:
        acc=ACCESSION,
        out=OUTDIR
    shell:
        """
        python benchmark.py --accession {params.acc} --output {params.out} --offline-blast
        """
