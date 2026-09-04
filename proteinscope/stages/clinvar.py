"""
ClinVar Data Loading for Missense Variant Pathogenicity Labels.

Downloads and parses NCBI ClinVar variant_summary.txt.gz, filters to a target gene
(default: TP53), extracts missense variants with pathogenic/benign labels, and returns
a pandas DataFrame ready for classifier training.

Supports:
  - Online download from NCBI FTP (cached locally after first fetch)
  - Offline fallback from a manually-placed TSV file
  - HGVS protein notation parsing (e.g. p.Arg175His → R175H)
"""

import gzip
import logging
import os
import re
import shutil
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
import requests

logger = logging.getLogger("ProteinScope.ClinVar")

# Standard 3-letter to 1-letter amino acid mapping
AA_3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}

# Regex for HGVS protein missense notation: p.Arg175His or p.R175H
HGVS_MISSENSE_3LETTER = re.compile(
    r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})"
)
HGVS_MISSENSE_1LETTER = re.compile(
    r"p\.([A-Z])(\d+)([A-Z])"
)

# ClinVar FTP URL for the tab-delimited variant summary
DEFAULT_CLINVAR_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
)


def parse_hgvs_missense(hgvs: str) -> Optional[Tuple[str, int, str]]:
    """
    Parse HGVS protein-level missense notation into (ref_aa, position, alt_aa).

    Supports both 3-letter (p.Arg175His) and 1-letter (p.R175H) formats.
    Returns None if the notation is not a recognizable missense change.
    """
    if not hgvs or not isinstance(hgvs, str):
        return None

    # Try 3-letter format first (more common in ClinVar)
    m = HGVS_MISSENSE_3LETTER.search(hgvs)
    if m:
        ref_3 = m.group(1)
        pos = int(m.group(2))
        alt_3 = m.group(3)
        ref_1 = AA_3TO1.get(ref_3)
        alt_1 = AA_3TO1.get(alt_3)
        if ref_1 and alt_1 and ref_1 != alt_1:
            return (ref_1, pos, alt_1)

    # Try 1-letter format
    m = HGVS_MISSENSE_1LETTER.search(hgvs)
    if m:
        ref_1 = m.group(1)
        pos = int(m.group(2))
        alt_1 = m.group(3)
        if ref_1 != alt_1:
            return (ref_1, pos, alt_1)

    return None


def classify_clinical_significance(clin_sig: str) -> Optional[str]:
    """
    Map ClinVar clinical significance string to a binary label.

    Returns:
      'pathogenic' for Pathogenic, Likely pathogenic, Pathogenic/Likely pathogenic
      'benign'     for Benign, Likely benign, Benign/Likely benign
      None         for VUS, Conflicting, Not provided, or anything else (dropped)
    """
    if not clin_sig or not isinstance(clin_sig, str):
        return None

    sig_lower = clin_sig.strip().lower()

    pathogenic_terms = [
        "pathogenic",
        "likely pathogenic",
        "pathogenic/likely pathogenic",
    ]
    benign_terms = [
        "benign",
        "likely benign",
        "benign/likely benign",
    ]

    for term in pathogenic_terms:
        if sig_lower == term:
            return "pathogenic"

    for term in benign_terms:
        if sig_lower == term:
            return "benign"

    return None


def download_clinvar_summary(
    url: str = DEFAULT_CLINVAR_URL,
    cache_dir: str = "./data/clinvar",
    timeout: int = 120,
) -> Optional[str]:
    """
    Download variant_summary.txt.gz from NCBI ClinVar FTP and cache locally.
    Returns path to the decompressed TSV file, or None on failure.
    """
    os.makedirs(cache_dir, exist_ok=True)

    gz_path = os.path.join(cache_dir, "variant_summary.txt.gz")
    tsv_path = os.path.join(cache_dir, "variant_summary.txt")

    # Use cached decompressed file if it exists and has content
    if os.path.exists(tsv_path) and os.path.getsize(tsv_path) > 1000:
        logger.info(f"Using cached ClinVar summary at {tsv_path}")
        return tsv_path

    # Use cached gzip file if it exists
    if os.path.exists(gz_path) and os.path.getsize(gz_path) > 1000:
        logger.info(f"Found cached gzip at {gz_path}, decompressing...")
        try:
            with gzip.open(gz_path, "rt", encoding="utf-8", errors="replace") as f_in:
                with open(tsv_path, "w", encoding="utf-8") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            logger.info(f"Decompressed ClinVar summary to {tsv_path}")
            return tsv_path
        except Exception as e:
            logger.error(f"Failed to decompress cached gzip: {e}")

    # Download from NCBI FTP
    logger.info(f"Downloading ClinVar variant summary from {url} ...")
    try:
        resp = requests.get(url, timeout=timeout, stream=True)
        resp.raise_for_status()

        with open(gz_path, "wb") as f:
            for chunk in resp.iter_content(chunk_size=1024 * 1024):
                f.write(chunk)

        logger.info(f"Downloaded ClinVar summary to {gz_path} ({os.path.getsize(gz_path):,} bytes)")

        # Decompress
        with gzip.open(gz_path, "rt", encoding="utf-8", errors="replace") as f_in:
            with open(tsv_path, "w", encoding="utf-8") as f_out:
                shutil.copyfileobj(f_in, f_out)

        logger.info(f"Decompressed ClinVar summary to {tsv_path}")
        return tsv_path

    except Exception as e:
        logger.error(f"Failed to download ClinVar summary: {e}")
        return None


def load_clinvar_variants(
    gene_symbol: str = "TP53",
    cache_dir: str = "./data/clinvar",
    clinvar_url: str = DEFAULT_CLINVAR_URL,
    min_review_stars: int = 1,
    offline: bool = False,
) -> pd.DataFrame:
    """
    Load and filter ClinVar missense variants for a specific gene.

    Args:
        gene_symbol: HGNC gene symbol to filter (default: TP53)
        cache_dir: Directory for cached ClinVar data files
        clinvar_url: URL for the ClinVar variant summary download
        min_review_stars: Minimum ClinVar review status star count (0-4)
        offline: If True, only use locally-cached files (no download)

    Returns:
        DataFrame with columns:
            position (int), ref_aa (str), alt_aa (str), mutation (str),
            hgvs (str), label (str: 'pathogenic'|'benign'),
            clinvar_id (int), clinical_significance (str),
            review_stars (int)
    """
    logger.info(f"Loading ClinVar variants for gene: {gene_symbol}")
    os.makedirs(cache_dir, exist_ok=True)

    # Check for pre-processed gene-specific TSV first (fastest path)
    gene_tsv = os.path.join(cache_dir, f"clinvar_{gene_symbol.lower()}_variants.tsv")
    if os.path.exists(gene_tsv) and os.path.getsize(gene_tsv) > 100:
        logger.info(f"Loading pre-processed gene variants from {gene_tsv}")
        try:
            df = pd.read_csv(gene_tsv, sep="\t")
            if len(df) > 0 and "label" in df.columns:
                logger.info(f"Loaded {len(df)} pre-processed variants for {gene_symbol}")
                return df
        except Exception as e:
            logger.warning(f"Failed to read pre-processed TSV: {e}")

    # Download or load full ClinVar summary
    tsv_path = None
    if not offline:
        tsv_path = download_clinvar_summary(
            url=clinvar_url, cache_dir=cache_dir
        )

    if tsv_path is None:
        # Fallback: look for manually-placed full summary
        fallback = os.path.join(cache_dir, "variant_summary.txt")
        if os.path.exists(fallback):
            tsv_path = fallback
            logger.info(f"Using manually-placed ClinVar summary at {fallback}")
        else:
            logger.error(
                f"No ClinVar data available. Place variant_summary.txt in {cache_dir} "
                "or run with internet access."
            )
            return pd.DataFrame()

    # Parse the full ClinVar TSV, filtering to target gene
    logger.info(f"Parsing ClinVar summary for gene {gene_symbol} ...")
    variants: List[Dict[str, Any]] = []

    # Map review status strings to star counts
    review_star_map = {
        "practice guideline": 4,
        "reviewed by expert panel": 3,
        "criteria provided, multiple submitters, no conflicts": 2,
        "criteria provided, conflicting classifications": 1,
        "criteria provided, conflicting interpretations": 1,
        "criteria provided, single submitter": 1,
        "no assertion for the individual variant": 0,
        "no assertion criteria provided": 0,
        "no classification provided": 0,
        "no classification for the individual variant": 0,
    }

    try:
        # Read in chunks to handle the large file efficiently
        chunk_iter = pd.read_csv(
            tsv_path,
            sep="\t",
            dtype=str,
            chunksize=50000,
            on_bad_lines="skip",
            low_memory=False,
        )

        for chunk in chunk_iter:
            # Normalize column names (ClinVar TSV can have '#' prefix)
            chunk.columns = [c.lstrip("#").strip() for c in chunk.columns]

            # Filter to target gene
            gene_col = None
            for candidate in ["GeneSymbol", "Gene symbol", "Gene(s)"]:
                if candidate in chunk.columns:
                    gene_col = candidate
                    break

            if gene_col is None:
                logger.warning(f"Could not find gene symbol column. Available: {list(chunk.columns)[:10]}")
                continue

            gene_mask = chunk[gene_col].str.upper().fillna("") == gene_symbol.upper()
            gene_chunk = chunk[gene_mask]

            if len(gene_chunk) == 0:
                continue

            # Filter to single nucleotide variants (missense)
            type_col = None
            for candidate in ["Type", "VariantType"]:
                if candidate in chunk.columns:
                    type_col = candidate
                    break

            if type_col:
                snv_mask = gene_chunk[type_col].str.lower().fillna("").str.contains("single nucleotide")
                gene_chunk = gene_chunk[snv_mask]

            # Extract HGVS protein notation
            hgvs_col = None
            for candidate in ["Name", "name", "HGVS(p)", "ProteinChange"]:
                if candidate in chunk.columns:
                    hgvs_col = candidate
                    break

            # Clinical significance column
            sig_col = None
            for candidate in [
                "ClinicalSignificance",
                "Clinical significance (Last reviewed)",
                "ClinSigSimple",
                "Clinical significance",
                "ClassType",
            ]:
                if candidate in chunk.columns:
                    sig_col = candidate
                    break

            # Review status column
            review_col = None
            for candidate in ["ReviewStatus", "Review status"]:
                if candidate in chunk.columns:
                    review_col = candidate
                    break

            # Variant ID column
            id_col = None
            for candidate in ["VariationID", "#VariationID", "AlleleID"]:
                if candidate in chunk.columns:
                    id_col = candidate
                    break

            for _, row in gene_chunk.iterrows():
                # Parse clinical significance
                clin_sig_raw = str(row.get(sig_col, "")) if sig_col else ""
                label = classify_clinical_significance(clin_sig_raw)
                if label is None:
                    continue  # Skip VUS, conflicting, etc.

                # Parse HGVS protein change
                hgvs_raw = str(row.get(hgvs_col, "")) if hgvs_col else ""
                # Also try the Name column which often contains the full HGVS
                name_raw = str(row.get("Name", "")) if "Name" in row.index else ""

                parsed = parse_hgvs_missense(hgvs_raw)
                if parsed is None:
                    parsed = parse_hgvs_missense(name_raw)
                if parsed is None:
                    continue  # Not a parseable missense variant

                ref_aa, position, alt_aa = parsed

                # Review status → star count
                review_raw = str(row.get(review_col, "")).strip().lower() if review_col else ""
                stars = review_star_map.get(review_raw, 0)

                if stars < min_review_stars:
                    continue

                # Variant ID
                var_id = 0
                if id_col:
                    try:
                        var_id = int(row.get(id_col, 0))
                    except (ValueError, TypeError):
                        var_id = 0

                variants.append({
                    "position": position,
                    "ref_aa": ref_aa,
                    "alt_aa": alt_aa,
                    "mutation": f"{ref_aa}{position}{alt_aa}",
                    "hgvs": hgvs_raw if hgvs_raw != "nan" else name_raw,
                    "label": label,
                    "clinvar_id": var_id,
                    "clinical_significance": clin_sig_raw,
                    "review_stars": stars,
                })

    except Exception as e:
        logger.error(f"Error parsing ClinVar summary: {e}")
        return pd.DataFrame()

    if not variants:
        logger.warning(f"No labeled missense variants found for {gene_symbol} in ClinVar data")
        return pd.DataFrame()

    df = pd.DataFrame(variants)

    # Deduplicate: keep one entry per unique mutation, preferring higher review stars
    df = df.sort_values("review_stars", ascending=False)
    df = df.drop_duplicates(subset=["position", "ref_aa", "alt_aa"], keep="first")
    df = df.sort_values("position").reset_index(drop=True)

    logger.info(
        f"Loaded {len(df)} unique labeled missense variants for {gene_symbol}: "
        f"{(df['label'] == 'pathogenic').sum()} pathogenic, "
        f"{(df['label'] == 'benign').sum()} benign"
    )

    # Cache the processed gene-specific TSV for fast reload
    df.to_csv(gene_tsv, sep="\t", index=False)
    logger.info(f"Cached processed variants to {gene_tsv}")

    return df
