"""
ProteinScope: Automated 11-Stage Protein Bioinformatics Pipeline CLI Entry Point.

Usage:
    # Standard pipeline (no mutation)
    python main.py P0DTC2
    python main.py P0DTC2 --config config/config.yaml --output ./output

    # Pathogenicity prediction for a specific mutation
    python main.py P04637 --mutation R175H
    python main.py P04637 --mutation R175H --train

    # Offline mode (use cached data only)
    python main.py P04637 --mutation R175H --offline-blast --offline
"""

import argparse
import os
import sys

# Ensure current directory is in sys.path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from proteinscope.protein_analyzer import ProteinAnalyzer


def main() -> None:
    parser = argparse.ArgumentParser(
        description="ProteinScope: Full-Stack Protein Bioinformatics Analysis Pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "accession",
        type=str,
        help="UniProt or NCBI Protein accession code (e.g., P04637, NP_000537 for TP53)",
    )
    parser.add_argument(
        "--config",
        type=str,
        default="config/config.yaml",
        help="Path to pipeline YAML configuration file",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="./output",
        help="Base directory for pipeline outputs",
    )
    parser.add_argument(
        "--offline-blast",
        action="store_true",
        help="Skip online NCBI QBLAST if internet is restricted or offline",
    )
    parser.add_argument(
        "--mutation",
        type=str,
        default=None,
        help=(
            "Missense mutation to analyze, in format REF_POS_ALT "
            "(e.g., R175H for Arg->His at position 175)"
        ),
    )
    parser.add_argument(
        "--train",
        action="store_true",
        help=(
            "Train (or retrain) the pathogenicity classifier on ClinVar data "
            "before running prediction. Required on first run."
        ),
    )
    parser.add_argument(
        "--offline",
        action="store_true",
        help="Run in fully offline mode — use only cached/local data (no API calls)",
    )
    parser.add_argument(
        "--temperature",
        type=float,
        default=37.0,
        help="Environmental temperature in °C (default: 37.0 standard physiological)",
    )
    parser.add_argument(
        "--ph",
        type=float,
        default=7.4,
        help="Environmental pH level (default: 7.4 standard physiological)",
    )
    parser.add_argument(
        "--temp-range",
        type=str,
        default=None,
        help="Permissible stability temperature range as MIN:MAX in °C (e.g. 15.0:45.0)",
    )
    parser.add_argument(
        "--ph-range",
        type=str,
        default=None,
        help="Permissible stability pH range as MIN:MAX (e.g. 5.5:8.5)",
    )

    args = parser.parse_args()

    accession = args.accession.strip().upper()
    mutation = args.mutation.strip().upper() if args.mutation else None

    # Parse custom ranges if provided
    temp_range = None
    if args.temp_range:
        try:
            parts = [float(p) for p in args.temp_range.split(":")]
            if len(parts) == 2:
                temp_range = (parts[0], parts[1])
        except ValueError:
            pass

    ph_range = None
    if args.ph_range:
        try:
            parts = [float(p) for p in args.ph_range.split(":")]
            if len(parts) == 2:
                ph_range = (parts[0], parts[1])
        except ValueError:
            pass

    if mutation:
        print(f"\n[ProteinScope] Initializing pathogenicity analysis for {accession} mutation {mutation}")
    else:
        print(f"\n[ProteinScope] Initializing analysis for accession: {accession}")
    print(f"[ProteinScope] Target environment: {args.temperature}°C, pH {args.ph}")

    analyzer = ProteinAnalyzer(
        accession=accession,
        output_dir=args.output,
        config_path=args.config if os.path.exists(args.config) else None,
        mutation=mutation,
        temperature=args.temperature,
        ph=args.ph,
        temp_range=temp_range,
        ph_range=ph_range,
    )

    # Determine if classifier training is needed
    train_classifier = args.train
    if mutation and not train_classifier:
        # Auto-train if no pre-trained model exists
        patho_cfg = analyzer.config.get("pathogenicity", {})
        model_path = patho_cfg.get("model_path", "./data/models/tp53_pathogenicity_rf.joblib")
        if not os.path.exists(model_path):
            print("[ProteinScope] No pre-trained model found — will train classifier on ClinVar data")
            train_classifier = True

    # Run complete workflow
    results = analyzer.run_all(
        run_blast_online=(not args.offline_blast and not args.offline),
        train_classifier=train_classifier,
        offline=args.offline,
    )

    # Print key result to console
    if mutation and "pathogenicity" in results:
        pred = results.get("pathogenicity", {})
        print(f"\n{'='*60}")
        print(f"  PATHOGENICITY RESULT: {pred.get('prediction', 'unknown').upper()}")
        print(f"  Confidence: {pred.get('confidence', 0):.1%}")
        print(f"{'='*60}")

    # Print Environmental Stability Result
    if "environmental_stability" in results:
        stab_res = results["environmental_stability"]
        active_stab = stab_res.get("active_evaluation", {})
        print(f"\n{'='*60}")
        print(f"  BIOPHYSICAL STABILITY (Alberts NBK26830)")
        print(f"  Condition        : {active_stab.get('temperature_celsius', 37.0)} C, pH {active_stab.get('ph', 7.4)}")
        print(f"  State            : {active_stab.get('badge_label', 'Native')}")
        print(f"  Estimated Tm     : {stab_res.get('estimated_melting_temperature_celsius', 'N/A')} C")
        print(f"  Fraction Folded  : {active_stab.get('fraction_folded_percent', 'N/A')}%")
        print(f"  Net Charge Q(pH) : {active_stab.get('net_charge', 'N/A')} e")
        print(f"  Salt Bridges     : {active_stab.get('salt_bridge_retention_percent', 'N/A')}% intact")
        if active_stab.get("is_degraded"):
            print("  ALERT            : [CRITICAL] PROTEIN DEGRADED / IRREVERSIBLY INACTIVATED")
        print(f"{'='*60}")

    summary_path = os.path.join(args.output, accession, "report", "summary.txt")
    print(f"\nDone. Full report at {summary_path}")


if __name__ == "__main__":
    main()
