"""
ProteinScope Main Engine: ProteinAnalyzer Class.
Coordinates the complete 11-stage automated bioinformatics pipeline,
with optional missense variant pathogenicity prediction.
"""

import json
import logging
import os
import re
import shutil
import sys
import time
from typing import Any, Dict, List, Optional, Tuple
import yaml

from proteinscope.stages.ingest import ingest_uniprot_data, validate_accession
from proteinscope.stages.sequence import analyze_protein_sequence
from proteinscope.stages.blast import run_blast_search
from proteinscope.stages.structure import retrieve_protein_structures
from proteinscope.stages.validate import validate_protein_structure
from proteinscope.stages.stability import analyze_protein_stability
from proteinscope.stages.visualize import visualize_protein_structure, generate_interactive_html_report
from proteinscope.stages.comparative import run_comparative_analysis
from proteinscope.stages.ml import extract_and_cluster_ml_features
from proteinscope.stages.ngs import process_ngs_data
from proteinscope.stages.clinvar import load_clinvar_variants
from proteinscope.stages.variant_features import (
    extract_variant_features,
    extract_feature_matrix,
    build_position_conservation,
)
from proteinscope.stages.pathogenicity import (
    train_pathogenicity_classifier,
    load_trained_model,
    predict_variant,
    generate_justification,
    plot_confusion_matrix,
    plot_feature_importances,
)


class ProteinAnalyzer:
    """
    Object-oriented protein bioinformatics analyzer for ProteinScope.
    Executes an automated multi-stage pipeline from a single UniProt accession code,
    with optional missense variant pathogenicity prediction.
    """

    @staticmethod
    def parse_mutation(mutation_str: str) -> Optional[Tuple[str, int, str]]:
        """
        Parse a mutation string like 'R175H' into (ref_aa, position, alt_aa).
        Returns None if the string doesn't match the expected format.
        """
        if not mutation_str:
            return None
        mutation_str = mutation_str.strip().upper()
        m = re.match(r"^([A-Z])(\d+)([A-Z])$", mutation_str)
        if m:
            return (m.group(1), int(m.group(2)), m.group(3))
        return None

    def __init__(
        self,
        accession: str,
        output_dir: str = "./output",
        config_path: Optional[str] = None,
        mutation: Optional[str] = None,
        temperature: Optional[float] = None,
        ph: Optional[float] = None,
        temp_range: Optional[Tuple[float, float]] = None,
        ph_range: Optional[Tuple[float, float]] = None,
    ) -> None:
        """
        Initialize the ProteinAnalyzer pipeline instance.

        Args:
            accession: UniProt accession code (e.g. 'P04637')
            output_dir: Base output directory
            config_path: Path to YAML configuration file
            mutation: Optional missense mutation string (e.g. 'R175H')
            temperature: Optional environmental temperature in °C (default 37.0)
            ph: Optional environmental pH (default 7.4)
            temp_range: Optional tuple of (min_temp, max_temp) stability cutoff
            ph_range: Optional tuple of (min_ph, max_ph) stability cutoff
        """
        self.accession = accession.strip().upper()
        self.base_output_dir = output_dir
        self.output_dir = os.path.join(output_dir, self.accession)
        self.config = self._load_config(config_path)
        self.start_time: Optional[float] = None
        self.end_time: Optional[float] = None

        # Mutation state
        self.mutation_str = mutation
        self.mutation_parsed = self.parse_mutation(mutation) if mutation else None

        # Environmental stability parameters (Alberts NBK26830)
        self.temperature = temperature
        self.ph = ph
        self.temp_range = temp_range
        self.ph_range = ph_range

        # Internal state store
        self.ingest_data: Dict[str, Any] = {}
        self.sequence_data: Dict[str, Any] = {}
        self.blast_data: Dict[str, Any] = {}
        self.structure_data: Dict[str, Any] = {}
        self.validation_data: Dict[str, Any] = {}
        self.stability_data: Dict[str, Any] = {}
        self.visualization_data: Dict[str, Any] = {}
        self.comparative_data: Dict[str, Any] = {}
        self.ml_data: Dict[str, Any] = {}
        self.ngs_data: Dict[str, Any] = {}
        self.clinvar_data: Optional[Any] = None  # pandas DataFrame
        self.pathogenicity_data: Dict[str, Any] = {}
        self.variant_features: Dict[str, float] = {}
        self.final_report: Dict[str, Any] = {}

        # Setup logging
        self._setup_logging()

    def _load_config(self, config_path: Optional[str]) -> Dict[str, Any]:
        """Load YAML configuration or use robust defaults."""
        default_config = {
            "pipeline": {"timeout_seconds": 15, "log_level": "INFO"},
            "blast": {"max_hits": 10, "matrix": "BLOSUM62", "query_max_length": 500},
            "validation": {"ca_cb_bond_standard": 1.52, "ca_cb_bond_tolerance": 0.05, "burial_contact_radius": 8.0},
            "comparative": {"top_homologs": 5, "string_species_taxid": 9606},
            "ml": {"n_clusters": 3, "random_state": 42},
            "ngs": {"input_dir": "./input/ngs"},
        }

        if config_path and os.path.exists(config_path):
            try:
                with open(config_path, "r", encoding="utf-8") as f:
                    user_cfg = yaml.safe_load(f)
                    if isinstance(user_cfg, dict):
                        default_config.update(user_cfg)
            except Exception as e:
                print(f"[Warning] Failed to load config from {config_path}: {e}")
        return default_config

    def _setup_logging(self) -> None:
        """Configure dual console and file logging."""
        os.makedirs(self.output_dir, exist_ok=True)
        log_file = os.path.join(self.output_dir, "run.log")

        self.logger = logging.getLogger(f"ProteinScope.{self.accession}")
        self.logger.setLevel(logging.INFO)

        # Clear existing handlers to prevent duplicate lines
        if self.logger.hasHandlers():
            self.logger.handlers.clear()

        formatter = logging.Formatter(
            "[%(asctime)s] [%(levelname)s] [%(name)s] %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S"
        )

        # File handler
        fh = logging.FileHandler(log_file, mode="a", encoding="utf-8")
        fh.setLevel(logging.INFO)
        fh.setFormatter(formatter)
        self.logger.addHandler(fh)

        # Console handler
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(logging.INFO)
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)

    # ── STAGE 0: INPUT VALIDATION ─────────────────────────────────────────────
    def validate_input(self) -> bool:
        """
        Stage 0: Validate UniProt accession code format and create directory hierarchy.
        """
        self.start_time = time.time()
        self.logger.info("=" * 60)
        self.logger.info(f"Starting ProteinScope Pipeline for accession: {self.accession}")
        self.logger.info("=" * 60)

        is_valid = validate_accession(self.accession)
        if not is_valid:
            self.logger.error(
                f"Invalid UniProt accession format: '{self.accession}'. "
                "Expected pattern: [O,P,Q][0-9][A-Z0-9]{3}[0-9] or [A-Z][0-9][A-Z0-9]{3}[0-9]"
            )
            return False

        # Create output directories
        subdirs = ["sequences", "structures", "analysis", "validation", "figures", "ml", "ngs", "report", "raw_data"]
        for sub in subdirs:
            os.makedirs(os.path.join(self.output_dir, sub), exist_ok=True)

        self.logger.info(f"Stage 0 passed: Accession format validated and output directories created at {self.output_dir}")
        return True

    # ── STAGE 1: DATA INGESTION ───────────────────────────────────────────────
    def fetch_from_uniprot(self) -> Dict[str, Any]:
        """
        Stage 1: Fetch FASTA, XML, and JSON metadata from UniProt REST API.
        """
        timeout = self.config.get("pipeline", {}).get("timeout_seconds", 15)
        self.ingest_data = ingest_uniprot_data(self.accession, self.output_dir, timeout=timeout)
        return self.ingest_data

    # ── STAGE 2: SEQUENCE ANALYSIS ────────────────────────────────────────────
    def analyze_sequence(self) -> Dict[str, Any]:
        """
        Stage 2: Physicochemical properties, motif search, and InterPro domain annotation.
        """
        timeout = self.config.get("pipeline", {}).get("timeout_seconds", 15)
        fasta_path = self.ingest_data.get("fasta_path")
        self.sequence_data = analyze_protein_sequence(
            self.accession,
            self.output_dir,
            fasta_path=fasta_path,
            timeout=timeout
        )
        return self.sequence_data

    # ── STAGE 3: BLAST SEARCH ─────────────────────────────────────────────────
    def blast_search(self, hits: int = 10, run_online: bool = True) -> Dict[str, Any]:
        """
        Stage 3: Run BLASTp similarity search against NCBI nr database.
        """
        seq = self.ingest_data.get("sequence", "")
        if not seq and self.sequence_data:
            # Check fasta
            fasta_p = os.path.join(self.output_dir, "sequences", f"{self.accession}.fasta")
            if os.path.exists(fasta_p):
                with open(fasta_p, "r", encoding="utf-8") as f:
                    seq = "".join(f.read().split("\n")[1:]).strip()

        blast_cfg = self.config.get("blast", {})
        max_hits = blast_cfg.get("max_hits", hits)
        matrix = blast_cfg.get("matrix", "BLOSUM62")
        q_len = blast_cfg.get("query_max_length", 500)

        self.blast_data = run_blast_search(
            sequence=seq,
            output_dir=self.output_dir,
            max_hits=max_hits,
            run_online=run_online,
            matrix_name=matrix,
            query_max_len=q_len,
        )
        return self.blast_data

    # ── STAGE 4: STRUCTURE RETRIEVAL ──────────────────────────────────────────
    def fetch_structures(self) -> Dict[str, Any]:
        """
        Stage 4: Retrieve structures across experimental (PDB), AlphaFold, and SWISS-MODEL paths.
        """
        timeout = self.config.get("pipeline", {}).get("timeout_seconds", 15)
        pdb_refs = self.ingest_data.get("pdb_cross_references", [])
        self.structure_data = retrieve_protein_structures(
            self.accession,
            self.output_dir,
            pdb_cross_refs=pdb_refs,
            timeout=timeout
        )
        return self.structure_data

    # ── STAGE 5: STRUCTURE VALIDATION ─────────────────────────────────────────
    def validate_structure(self, pdb_file: Optional[str] = None) -> Dict[str, Any]:
        """
        Stage 5: Ramachandran dihedral distribution, bond geometry checks, and burial analysis.
        """
        target_pdb = pdb_file
        if not target_pdb:
            best = self.structure_data.get("best_structure")
            if best:
                target_pdb = best.get("file_path")

        if not target_pdb or not os.path.exists(target_pdb):
            self.logger.warning("No structure available for validation stage.")
            return {"error": "No structure file provided"}

        val_cfg = self.config.get("validation", {})
        standard = val_cfg.get("ca_cb_bond_standard", 1.52)
        tol = val_cfg.get("ca_cb_bond_tolerance", 0.05)
        radius = val_cfg.get("burial_contact_radius", 8.0)

        self.validation_data = validate_protein_structure(
            pdb_file=target_pdb,
            output_dir=self.output_dir,
            ca_cb_standard=standard,
            ca_cb_tolerance=tol,
            burial_radius=radius,
        )
        return self.validation_data

    # ── STAGE 5b: ENVIRONMENTAL STABILITY & DEGRADATION (Alberts NBK26830) ────
    def analyze_stability(self) -> Dict[str, Any]:
        """
        Stage 5b: Biophysical temperature and pH stability profile,
        thermodynamic folding stability, net charge titration, and
        denaturation/degradation threshold modeling based on Alberts NBK26830.
        """
        seq = self.ingest_data.get("sequence", "")
        if not seq and self.sequence_data:
            seq = self.sequence_data.get("sequence", "")

        physico = self.sequence_data.get("physicochemical_properties", {})
        ptms = self.ingest_data.get("ptm_features", [])

        self.stability_data = analyze_protein_stability(
            sequence=seq,
            output_dir=self.output_dir,
            properties=physico,
            ptm_sites=ptms,
            user_temperature=self.temperature,
            user_ph=self.ph,
            custom_temp_range=self.temp_range,
            custom_ph_range=self.ph_range,
        )
        return self.stability_data

    # ── STAGE 6: VISUALIZATION ────────────────────────────────────────────────
    def visualize(self) -> Dict[str, Any]:
        """
        Stage 6: PyMOL scripting, superimposition RMSD, and cartoon/surface/domain exports.
        """
        domains = self.sequence_data.get("domains", [])
        self.visualization_data = visualize_protein_structure(
            structures_manifest=self.structure_data,
            output_dir=self.output_dir,
            domains=domains
        )
        return self.visualization_data

    # ── STAGE 7: COMPARATIVE ANALYSIS ─────────────────────────────────────────
    def comparative_analysis(self) -> Dict[str, Any]:
        """
        Stage 7: Homolog pairwise alignment distance matrix, PTM mapping, and STRING PPI network.
        """
        seq = self.ingest_data.get("sequence", "")
        blast_hits = self.blast_data.get("hits", [])
        ptms = self.ingest_data.get("ptm_features", [])
        taxid = self.ingest_data.get("taxid", 9606)
        timeout = self.config.get("pipeline", {}).get("timeout_seconds", 15)

        self.comparative_data = run_comparative_analysis(
            accession=self.accession,
            sequence=seq,
            output_dir=self.output_dir,
            blast_hits=blast_hits,
            ptm_features=ptms,
            taxid=taxid,
            timeout=timeout,
        )
        return self.comparative_data

    # ── STAGE 8: ML FEATURE ENGINEERING ───────────────────────────────────────
    def extract_ml_features(self) -> Dict[str, Any]:
        """
        Stage 8: Extract 40+ physicochemical feature vectors, KMeans clustering, and feature importance.
        """
        seq = self.ingest_data.get("sequence", "")
        blast_hits = self.blast_data.get("hits", [])
        ml_cfg = self.config.get("ml", {})
        n_clusters = ml_cfg.get("n_clusters", 3)
        rand_state = ml_cfg.get("random_state", 42)

        self.ml_data = extract_and_cluster_ml_features(
            query_accession=self.accession,
            query_sequence=seq,
            output_dir=self.output_dir,
            blast_hits=blast_hits,
            n_clusters=n_clusters,
            random_state=rand_state,
        )
        return self.ml_data

    # ── STAGE 9: NGS INTEGRATION ──────────────────────────────────────────────
    def process_ngs(self, ngs_dir: Optional[str] = None) -> Dict[str, Any]:
        """
        Stage 9: Parse FASTQ, BAM, and VCF mutation mapping if present.
        """
        target_ngs = ngs_dir or self.config.get("ngs", {}).get("input_dir", "./input/ngs")
        self.ngs_data = process_ngs_data(output_dir=self.output_dir, ngs_input_dir=target_ngs)
        return self.ngs_data

    # ── CLINVAR DATA LOADING ───────────────────────────────────────────────────
    def load_clinvar(self, offline: bool = False) -> Any:
        """
        Load ClinVar missense variant labels for the target gene.
        """
        clinvar_cfg = self.config.get("clinvar", {})
        gene_symbol = clinvar_cfg.get("gene_symbol", "TP53")
        cache_dir = clinvar_cfg.get("cache_dir", "./data/clinvar")
        clinvar_url = clinvar_cfg.get(
            "ftp_url",
            "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
        )
        min_stars = clinvar_cfg.get("min_review_stars", 1)

        self.logger.info(f"Loading ClinVar variants for gene: {gene_symbol}")
        self.clinvar_data = load_clinvar_variants(
            gene_symbol=gene_symbol,
            cache_dir=cache_dir,
            clinvar_url=clinvar_url,
            min_review_stars=min_stars,
            offline=offline,
        )
        return self.clinvar_data

    # ── VARIANT FEATURE EXTRACTION ────────────────────────────────────────────
    def extract_variant_feature_vector(self) -> Dict[str, float]:
        """
        Extract the feature vector for the current mutation using cached
        pipeline stage outputs.
        """
        if not self.mutation_parsed:
            self.logger.warning("No mutation specified for feature extraction")
            return {}

        ref_aa, position, alt_aa = self.mutation_parsed
        self.logger.info(f"Extracting features for variant {ref_aa}{position}{alt_aa}")

        # Get best PDB file path
        pdb_file = None
        best = self.structure_data.get("best_structure")
        if best:
            pdb_file = best.get("file_path")

        # Get sequence
        sequence = self.ingest_data.get("sequence", "")

        # Get BLAST hits
        blast_hits = self.blast_data.get("hits", [])

        # Get PTM sites
        ptm_sites = self.comparative_data.get("ptm_sites", [])
        if not ptm_sites:
            ptm_sites = self.ingest_data.get("ptm_features", [])

        self.variant_features = extract_variant_features(
            position=position,
            ref_aa=ref_aa,
            alt_aa=alt_aa,
            query_sequence=sequence,
            pdb_file=pdb_file,
            blast_hits=blast_hits,
            ptm_sites=ptm_sites,
            include_ml_deltas=True,
        )

        return self.variant_features

    # ── PATHOGENICITY PREDICTION ──────────────────────────────────────────────
    def predict_pathogenicity(self, train: bool = False, offline: bool = False) -> Dict[str, Any]:
        """
        Run the pathogenicity prediction sub-pipeline.

        If train=True, trains a new classifier on ClinVar data first.
        Otherwise, loads a pre-trained model.
        """
        if not self.mutation_parsed:
            self.logger.warning("No mutation specified for pathogenicity prediction")
            return {}

        ref_aa, position, alt_aa = self.mutation_parsed
        self.logger.info(f"Running pathogenicity prediction for {ref_aa}{position}{alt_aa}")

        patho_cfg = self.config.get("pathogenicity", {})
        model_dir = os.path.dirname(
            patho_cfg.get("model_path", "./data/models/tp53_pathogenicity_rf.joblib")
        )
        model_path = patho_cfg.get("model_path", "./data/models/tp53_pathogenicity_rf.joblib")
        confidence_threshold = patho_cfg.get("confidence_threshold", 0.5)

        # Train if requested or if no model exists
        if train or not os.path.exists(model_path):
            self.logger.info("Training pathogenicity classifier...")

            # Load ClinVar data if not already loaded
            if self.clinvar_data is None or (hasattr(self.clinvar_data, 'empty') and self.clinvar_data.empty):
                self.load_clinvar(offline=offline)

            if self.clinvar_data is None or len(self.clinvar_data) == 0:
                self.logger.error("No ClinVar data available for training")
                return {"error": "No ClinVar training data"}

            # Extract features for all ClinVar variants
            sequence = self.ingest_data.get("sequence", "")
            pdb_file = None
            best = self.structure_data.get("best_structure")
            if best:
                pdb_file = best.get("file_path")

            blast_hits = self.blast_data.get("hits", [])
            ptm_sites = self.comparative_data.get("ptm_sites", [])
            if not ptm_sites:
                ptm_sites = self.ingest_data.get("ptm_features", [])

            variants_list = self.clinvar_data.to_dict("records")

            self.logger.info(f"Extracting features for {len(variants_list)} ClinVar variants...")
            feature_dicts, feature_names = extract_feature_matrix(
                variants=variants_list,
                query_sequence=sequence,
                pdb_file=pdb_file,
                blast_hits=blast_hits,
                ptm_sites=ptm_sites,
                include_ml_deltas=True,
            )

            # Match labels to successfully extracted features
            labels = []
            valid_variants = []
            feat_idx = 0
            for var in variants_list:
                pos = var.get("position", 0)
                ref = var.get("ref_aa", "")
                if sequence and 0 < pos <= len(sequence) and sequence[pos - 1] == ref:
                    if feat_idx < len(feature_dicts):
                        labels.append(var["label"])
                        valid_variants.append(var)
                        feat_idx += 1

            if len(feature_dicts) != len(labels):
                # Trim to minimum length
                min_len = min(len(feature_dicts), len(labels))
                feature_dicts = feature_dicts[:min_len]
                labels = labels[:min_len]

            n_estimators = patho_cfg.get("n_estimators", 200)
            test_size = patho_cfg.get("test_size", 0.2)
            random_state = patho_cfg.get("random_state", 42)

            training_metrics = train_pathogenicity_classifier(
                feature_dicts=feature_dicts,
                labels=labels,
                model_dir=model_dir,
                n_estimators=n_estimators,
                test_size=test_size,
                random_state=random_state,
            )

            # Save training metrics
            benchmark_dir = os.path.join(self.output_dir, "benchmark")
            os.makedirs(benchmark_dir, exist_ok=True)
            metrics_path = os.path.join(benchmark_dir, "metrics.json")
            with open(metrics_path, "w", encoding="utf-8") as f:
                json.dump(training_metrics, f, indent=2, default=str)
            self.logger.info(f"Saved training metrics to {metrics_path}")

            # Generate plots
            if "confusion_matrix" in training_metrics:
                cm_path = os.path.join(benchmark_dir, "confusion_matrix.png")
                plot_confusion_matrix(training_metrics["confusion_matrix"], cm_path)

            if "top_features" in training_metrics:
                fi_path = os.path.join(benchmark_dir, "feature_importances.png")
                plot_feature_importances(training_metrics["top_features"], fi_path)

            self.pathogenicity_data["training_metrics"] = training_metrics

        # Extract features for the target variant
        if not self.variant_features:
            self.extract_variant_feature_vector()

        # Predict
        prediction = predict_variant(
            feature_dict=self.variant_features,
            model_path=model_path,
            confidence_threshold=confidence_threshold,
        )

        # Generate justification
        justification = generate_justification(
            feature_dict=self.variant_features,
            prediction_result=prediction,
            ref_aa=ref_aa,
            alt_aa=alt_aa,
            position=position,
        )
        prediction["justification"] = justification

        self.pathogenicity_data["prediction"] = prediction
        self.pathogenicity_data["mutation"] = f"{ref_aa}{position}{alt_aa}"
        self.pathogenicity_data["variant_features"] = self.variant_features

        self.logger.info(
            f"Pathogenicity prediction for {ref_aa}{position}{alt_aa}: "
            f"{prediction['prediction'].upper()} "
            f"(confidence: {prediction['confidence']:.2%})"
        )

        return self.pathogenicity_data

    # ── STAGE 11: FINAL REPORT GENERATION ─────────────────────────────────────
    def save_report(self) -> Dict[str, Any]:
        """
        Stage 11: Compile consolidated report.json and human-readable summary.txt.
        When a mutation is specified, leads with pathogenicity verdict and justification.
        Copies all generated figures into report/figures/.
        """
        self.logger.info("Starting Stage 11: Compiling final multi-stage report")
        report_dir = os.path.join(self.output_dir, "report")
        report_fig_dir = os.path.join(report_dir, "figures")
        os.makedirs(report_dir, exist_ok=True)
        os.makedirs(report_fig_dir, exist_ok=True)

        # Copy generated figures into report/figures/
        for source_dir_name in ["figures", "benchmark"]:
            source_dir = os.path.join(self.output_dir, source_dir_name)
            if os.path.exists(source_dir):
                for fig_file in os.listdir(source_dir):
                    if fig_file.endswith(".png") or fig_file.endswith(".pdf"):
                        shutil.copy2(
                            os.path.join(source_dir, fig_file),
                            os.path.join(report_fig_dir, fig_file)
                        )

        self.end_time = time.time()
        duration = round(self.end_time - (self.start_time or self.end_time), 2)

        # 1. Structured JSON Report
        best_struct = self.structure_data.get("best_structure") or {}
        physico = self.sequence_data.get("physicochemical_properties") or {}
        ram_stats = self.validation_data.get("ramachandran_statistics") or {}

        self.final_report = {
            "accession": self.accession,
            "organism": self.ingest_data.get("organism", "Unknown"),
            "protein_name": self.ingest_data.get("protein_name", "Unknown"),
            "sequence_length": self.ingest_data.get("sequence_length", 0),
        }

        # Add pathogenicity data if a mutation was analyzed
        if self.mutation_parsed and self.pathogenicity_data:
            self.final_report["mutation"] = self.pathogenicity_data.get("mutation")
            self.final_report["pathogenicity"] = self.pathogenicity_data.get("prediction", {})

        self.final_report.update({
            "properties": physico,
            "motifs_found": {k: len(v) for k, v in self.sequence_data.get("motifs", {}).items()},
            "domains": self.sequence_data.get("domains", []),
            "blast_top5": self.blast_data.get("hits", [])[:5],
            "structure_source": best_struct.get("source") if best_struct else "None",
            "structure_best_file": os.path.basename(best_struct.get("file_path", "")) if best_struct else None,
            "structure_rmsd_angstroms": self.visualization_data.get("rmsd_superimposition_angstroms"),
            "validation": {
                "ramachandran": ram_stats,
                "ca_cb_bond_pass_rate": self.validation_data.get("geometry_ca_cb", {}).get("pass_rate_percent"),
                "burial": self.validation_data.get("burial_environment"),
            },
            "ptm_sites": self.comparative_data.get("ptm_sites", []),
            "ppi_partners": self.comparative_data.get("ppi_interactions", []),
            "environmental_stability": self.stability_data,
            "ml_features": {
                "feature_count": self.ml_data.get("feature_count"),
                "clusters": self.ml_data.get("clusters"),
                "top_features": self.ml_data.get("top_discriminative_features"),
            },
            "ngs_summary": self.ngs_data,
            "run_metadata": {
                "pipeline_version": "2.0.0",
                "execution_duration_seconds": duration,
                "output_directory": os.path.abspath(self.output_dir),
            }
        })

        json_path = os.path.join(report_dir, "report.json")
        with open(json_path, "w", encoding="utf-8") as f:
            json.dump(self.final_report, f, indent=2, default=str)

        # 2. Text Summary Report
        summary_txt_path = os.path.join(report_dir, "summary.txt")
        with open(summary_txt_path, "w", encoding="utf-8") as f:
            # ── Header ──
            if self.mutation_parsed:
                ref_aa, pos, alt_aa = self.mutation_parsed
                # Amino acid full names for readability
                aa_names = {
                    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
                    "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
                    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
                    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
                }
                ref_name = aa_names.get(ref_aa, ref_aa)
                alt_name = aa_names.get(alt_aa, alt_aa)

                f.write("=" * 70 + "\n")
                f.write(f"       PROTEINSCOPE PATHOGENICITY PREDICTION REPORT\n")
                f.write("=" * 70 + "\n\n")
                f.write(f"Accession       : {self.accession}\n")
                f.write(f"Protein Name    : {self.ingest_data.get('protein_name', 'Unknown')}\n")
                f.write(f"Organism        : {self.ingest_data.get('organism', 'Unknown')}\n")
                f.write(f"Mutation        : {ref_aa}{pos}{alt_aa} ({ref_name} \u2192 {alt_name} at position {pos})\n")
                f.write(f"Sequence Length : {self.ingest_data.get('sequence_length', 0)} amino acids\n")
                f.write(f"Execution Time  : {duration} seconds\n\n")

                # ── Pathogenicity Verdict ──
                pred = self.pathogenicity_data.get("prediction", {})
                verdict = pred.get("prediction", "unknown").upper()
                confidence = pred.get("confidence", 0.0)
                path_prob = pred.get("pathogenic_probability", 0.0)

                f.write("\u2500\u2500 PATHOGENICITY VERDICT \u2500" * 3 + "\u2500" * 22 + "\n")
                f.write(f"Prediction      : {verdict}\n")
                f.write(f"Confidence      : {confidence:.4f} ({confidence:.1%})\n")
                f.write(f"Path. Prob.     : {path_prob:.4f}\n")
                f.write(f"Threshold       : {pred.get('confidence_threshold', 0.5)}\n\n")

                # ── Structural Justification ──
                justification = pred.get("justification", "")
                if justification:
                    f.write("\u2500\u2500 STRUCTURAL JUSTIFICATION \u2500" * 2 + "\u2500" * 22 + "\n")
                    f.write(justification + "\n\n")

            else:
                f.write("=" * 70 + "\n")
                f.write(f"              PROTEINSCOPE AUTOMATED ANALYSIS REPORT\n")
                f.write("=" * 70 + "\n\n")
                f.write(f"Accession       : {self.accession}\n")
                f.write(f"Protein Name    : {self.ingest_data.get('protein_name', 'Unknown')}\n")
                f.write(f"Organism        : {self.ingest_data.get('organism', 'Unknown')}\n")
                f.write(f"Sequence Length : {self.ingest_data.get('sequence_length', 0)} amino acids\n")
                f.write(f"Execution Time  : {duration} seconds\n\n")

            # ── Supporting Pipeline Data (always included) ──
            f.write("\u2500\u2500 SUPPORTING PIPELINE DATA \u2500" * 2 + "\u2500" * 22 + "\n\n")

            f.write("  1. PHYSICOCHEMICAL PROPERTIES\n")
            f.write(f"  Molecular Weight      : {physico.get('molecular_weight', 'N/A')} Da\n")
            f.write(f"  Isoelectric Point (pI): {physico.get('isoelectric_point', 'N/A')}\n")
            f.write(f"  Instability Index     : {physico.get('instability_index', 'N/A')} ({'Stable' if physico.get('is_stable') else 'Unstable'})\n")
            f.write(f"  GRAVY Hydropathicity  : {physico.get('gravy', 'N/A')}\n")
            f.write(f"  Aromaticity           : {physico.get('aromaticity', 'N/A')}\n")
            sec = physico.get("secondary_structure_fraction", {})
            f.write(f"  Secondary Structure   : Helix: {sec.get('helix', 0)*100:.1f}%, Sheet: {sec.get('sheet', 0)*100:.1f}%, Turn: {sec.get('turn', 0)*100:.1f}%\n\n")

            f.write("  2. STRUCTURE & VALIDATION\n")
            f.write(f"  Best Structure Source : {best_struct.get('source', 'None')}\n")
            f.write(f"  Ramachandran Favored  : {ram_stats.get('favored_percent', 'N/A')}%\n")
            f.write(f"  Ramachandran Allowed  : {ram_stats.get('allowed_percent', 'N/A')}%\n")
            f.write(f"  Ramachandran Outliers : {ram_stats.get('outlier_percent', 'N/A')}%\n")
            f.write(f"  CA-CB Bond Pass Rate  : {self.validation_data.get('geometry_ca_cb', {}).get('pass_rate_percent', 'N/A')}%\n\n")

            f.write("  2b. ENVIRONMENTAL STABILITY & DEGRADATION (Alberts NBK26830)\n")
            active_stab = self.stability_data.get("active_evaluation", {})
            f.write(f"  Condition Evaluated   : {active_stab.get('temperature_celsius', 37.0)}°C, pH {active_stab.get('ph', 7.4)}\n")
            f.write(f"  Conformation State    : {active_stab.get('state', 'NATIVE')} ({active_stab.get('badge_label', 'Folded')})\n")
            f.write(f"  Estimated Tm (Melting): {self.stability_data.get('estimated_melting_temperature_celsius', 'N/A')}°C\n")
            f.write(f"  Fraction Folded       : {active_stab.get('fraction_folded_percent', 'N/A')}%\n")
            f.write(f"  Net Charge Q(pH)      : {active_stab.get('net_charge', 'N/A')} e\n")
            f.write(f"  Salt Bridge Retention : {active_stab.get('salt_bridge_retention_percent', 'N/A')}%\n")
            f.write(f"  Folding Free Energy dG: {active_stab.get('delta_g_folding_kcal_mol', 'N/A')} kcal/mol\n")
            if active_stab.get("is_degraded"):
                f.write("  Degradation Alert     : CRITICAL - Protein degraded / aggregated\n")
            f.write("\n")

            f.write("  3. HOMOLOGY & EVOLUTIONARY (BLAST + PPI)\n")
            f.write(f"  BLAST Homologs Found  : {len(self.blast_data.get('hits', []))}\n")
            f.write(f"  STRING PPI Interactors: {len(self.comparative_data.get('ppi_interactions', []))}\n")
            f.write(f"  Identified PTM Sites  : {len(self.comparative_data.get('ptm_sites', []))}\n\n")

            f.write("  4. MACHINE LEARNING & NGS\n")
            f.write(f"  ML Features Extracted : {self.ml_data.get('feature_count', 0)} numeric descriptors\n")
            clusters = self.ml_data.get('clusters', [{}])
            cluster_val = clusters[0].get('cluster', 'N/A') if clusters else 'N/A'
            f.write(f"  KMeans Cluster Group  : Cluster {cluster_val}\n")
            f.write(f"  NGS Integration Status: {'Detected' if self.ngs_data.get('ngs_data_detected') else 'No NGS datasets supplied (skipped)'}\n\n")

            f.write("=" * 70 + "\n")
            f.write(f"Outputs generated at: {os.path.abspath(self.output_dir)}\n")
            f.write("=" * 70 + "\n")

        self.logger.info(f"Final reports saved: {json_path} and {summary_txt_path}")

        # 3. Interactive 3D HTML Report
        try:
            html_path = generate_interactive_html_report(self.final_report, self.output_dir)
            self.logger.info(f"Interactive 3D HTML Report generated at: {html_path}")
        except Exception as e:
            self.logger.warning(f"Failed to generate interactive HTML report: {e}")

        return self.final_report

    # ── RUN ALL WORKFLOW ──────────────────────────────────────────────────────
    def run_all(
        self,
        run_blast_online: bool = True,
        train_classifier: bool = False,
        offline: bool = False,
    ) -> Dict[str, Any]:
        """
        Execute all stages of the ProteinScope pipeline sequentially.
        When a mutation is specified, includes pathogenicity prediction.

        Args:
            run_blast_online: Whether to submit BLAST queries to NCBI
            train_classifier: Whether to train the classifier (vs using pre-trained)
            offline: If True, skip all online API calls
        """
        self.logger.info("Executing Complete ProteinScope Pipeline Workflow")

        if self.mutation_parsed:
            ref, pos, alt = self.mutation_parsed
            self.logger.info(f"Mutation mode: analyzing {ref}{pos}{alt}")

        # Stage 0: Input Validation
        if not self.validate_input():
            return {"error": "Input validation failed"}

        # Stage 1: Data Ingestion
        self.fetch_from_uniprot()

        # Stage 2: Sequence Analysis
        self.analyze_sequence()

        # Stage 3: BLAST Search
        self.blast_search(hits=10, run_online=run_blast_online and not offline)

        # Stage 4: Structure Retrieval
        self.fetch_structures()

        # Stage 5: Structure Validation
        self.validate_structure()

        # Stage 5b: Environmental Stability & Degradation (Alberts NBK26830)
        self.analyze_stability()

        # Stage 6: Visualization
        self.visualize()

        # Stage 7: Comparative Analysis
        self.comparative_analysis()

        # Stage 8: ML Feature Engineering
        self.extract_ml_features()

        # Stage 9: NGS Integration
        self.process_ngs()

        # Pathogenicity prediction (when mutation is specified)
        if self.mutation_parsed:
            self.predict_pathogenicity(
                train=train_classifier,
                offline=offline,
            )

        # Stage 11: Save Reports
        final = self.save_report()

        # Print manifest to console
        self._print_manifest()
        return final

    def _print_manifest(self) -> None:
        """Print file manifest and run summary to standard output."""
        print("\n" + "=" * 70)
        print("          PROTEINSCOPE PIPELINE EXECUTION COMPLETE")
        print("=" * 70)
        print(f"Target Accession : {self.accession}")
        print(f"Output Directory : {os.path.abspath(self.output_dir)}")
        print("\nGenerated File Manifest:")
        
        for root, dirs, files in os.walk(self.output_dir):
            rel_root = os.path.relpath(root, self.output_dir)
            if rel_root == ".":
                prefix = ""
            else:
                prefix = f"  [{rel_root}]"
                print(prefix)
            for f in sorted(files):
                f_path = os.path.join(root, f)
                sz = os.path.getsize(f_path)
                print(f"    • {f} ({sz:,} bytes)")
        print("=" * 70 + "\n")
