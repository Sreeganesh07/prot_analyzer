"""
Stage 7: Comparative, Evolutionary, PTM, and Protein-Protein Interaction (STRING DB) Analysis.
Aligns top homologs, constructs distance matrix, maps PTMs, and retrieves PPI interactors.
"""

import json
import logging
import os
from typing import Any, Dict, List, Optional, Tuple
import numpy as np
import pandas as pd
import requests

from Bio import Align
from Bio.Seq import Seq

logger = logging.getLogger("ProteinScope.Comparative")


def fetch_string_ppi_network(accession: str, species_taxid: int = 9606, timeout: int = 15) -> List[Dict[str, Any]]:
    """
    Query STRING DB REST API for known and predicted protein-protein interactions.
    Endpoint: https://string-db.org/api/json/network?identifiers={accession}&species={taxid}
    """
    url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": accession,
        "species": species_taxid,
        "caller_identity": "proteinscope_pipeline",
    }
    interactions = []

    try:
        resp = requests.get(url, params=params, timeout=timeout)
        if resp.status_code == 200:
            data = resp.json()
            # Sort by combined score descending
            if isinstance(data, list):
                data.sort(key=lambda x: x.get("score", 0.0), reverse=True)
                for item in data[:10]:
                    partner_a = item.get("preferredName_A", item.get("stringId_A"))
                    partner_b = item.get("preferredName_B", item.get("stringId_B"))
                    # Identify the interacting partner that is not query if possible
                    partner = partner_b if accession in str(partner_a) else partner_a
                    interactions.append({
                        "partner_a": partner_a,
                        "partner_b": partner_b,
                        "preferred_partner_name": partner,
                        "combined_score": round(float(item.get("score", 0.0)), 3),
                        "escore_experiments": round(float(item.get("escore", 0.0)), 3),
                        "dscore_database": round(float(item.get("dscore", 0.0)), 3),
                        "tscore_textmining": round(float(item.get("tscore", 0.0)), 3),
                    })
        else:
            logger.warning(f"STRING DB returned status code {resp.status_code} for {accession}")
    except Exception as e:
        logger.error(f"Error querying STRING DB for {accession}: {e}")

    return interactions


def _simple_pairwise_identity(seq1: str, seq2: str) -> Tuple[int, int, float]:
    """Fallback simple pairwise identity counter when Bio.Align is unavailable."""
    min_len = min(len(seq1), len(seq2))
    max_len = max(len(seq1), len(seq2))
    if max_len == 0:
        return 0, 0, 0.0
    matches = sum(1 for a, b in zip(seq1[:min_len], seq2[:min_len]) if a == b)
    identity_pct = (matches / max_len) * 100.0
    return matches, max_len, identity_pct


def perform_pairwise_homolog_alignments(
    query_accession: str,
    query_sequence: str,
    blast_hits: List[Dict[str, Any]]
) -> Tuple[List[Dict[str, Any]], pd.DataFrame]:
    """
    Perform pairwise alignment of query sequence against top BLAST hit sequences.
    Builds NxN distance matrix based on (1.0 - identity_fraction).
    """
    aligner = None
    if Align is not None:
        try:
            aligner = Align.PairwiseAligner()
            aligner.mode = "global"
            aligner.open_gap_score = -10.0
            aligner.extend_gap_score = -0.5
        except Exception as e:
            logger.warning(f"Could not initialize Bio.Align.PairwiseAligner: {e}")
            aligner = None

    sequences = [(query_accession, query_sequence)]
    alignments_summary = []

    for hit in blast_hits[:5]:
        hit_acc = hit.get("accession", f"Hit_{hit.get('rank')}")
        # Use sbjct_sequence from BLAST HSP or fetch full sequence
        hit_seq = hit.get("sbjct_sequence", "").replace("-", "")
        if hit_seq:
            sequences.append((hit_acc, hit_seq))

    n = len(sequences)
    labels = [s[0] for s in sequences]
    distance_mat = np.zeros((n, n), dtype=float)

    # 1. Pairwise alignments against query
    for i in range(1, n):
        target_label, target_seq = sequences[i]
        try:
            if aligner is not None and Seq is not None:
                alns = aligner.align(Seq(query_sequence), Seq(target_seq))
                if alns:
                    top_aln = alns[0]
                    identities = top_aln.counts().identities
                    aln_len = top_aln.length
                    identity_pct = (identities / aln_len * 100) if aln_len > 0 else 0.0
                    score = float(top_aln.score)
                    dist = round(1.0 - (identities / aln_len), 4) if aln_len > 0 else 1.0
                else:
                    identities, aln_len, identity_pct = _simple_pairwise_identity(query_sequence, target_seq)
                    score = float(identities)
                    dist = round(1.0 - (identity_pct / 100.0), 4)
            else:
                identities, aln_len, identity_pct = _simple_pairwise_identity(query_sequence, target_seq)
                score = float(identities)
                dist = round(1.0 - (identity_pct / 100.0), 4)

            alignments_summary.append({
                "query": query_accession,
                "target": target_label,
                "score": round(score, 2),
                "alignment_length": aln_len,
                "identities": identities,
                "identity_percent": round(identity_pct, 2),
                "evolutionary_distance": dist,
            })
        except Exception as e:
            logger.error(f"Pairwise alignment failed for {target_label}: {e}")

    # 2. Build full NxN distance matrix
    for i in range(n):
        for j in range(n):
            if i == j:
                distance_mat[i, j] = 0.0
            elif i < j:
                try:
                    if aligner is not None and Seq is not None:
                        s1 = Seq(sequences[i][1])
                        s2 = Seq(sequences[j][1])
                        alns = aligner.align(s1, s2)
                        if alns:
                            top_aln = alns[0]
                            identities = top_aln.counts().identities
                            aln_len = top_aln.length
                            dist = round(1.0 - (identities / aln_len), 4) if aln_len > 0 else 1.0
                        else:
                            _, _, id_pct = _simple_pairwise_identity(sequences[i][1], sequences[j][1])
                            dist = round(1.0 - (id_pct / 100.0), 4)
                    else:
                        _, _, id_pct = _simple_pairwise_identity(sequences[i][1], sequences[j][1])
                        dist = round(1.0 - (id_pct / 100.0), 4)

                    distance_mat[i, j] = dist
                    distance_mat[j, i] = dist
                except Exception:
                    distance_mat[i, j] = 1.0
                    distance_mat[j, i] = 1.0

    df_dist = pd.DataFrame(distance_mat, index=labels, columns=labels)
    return alignments_summary, df_dist


def run_comparative_analysis(
    accession: str,
    sequence: str,
    output_dir: str,
    blast_hits: List[Dict[str, Any]],
    ptm_features: Optional[List[Dict[str, Any]]] = None,
    taxid: Optional[Any] = 9606,
    timeout: int = 15
) -> Dict[str, Any]:
    """
    Stage 7 Entry point:
    Executes evolutionary distance matrix, PTM site mapping, and STRING DB network query.
    """
    logger.info(f"Starting Stage 7: Comparative and evolutionary analysis for {accession}")
    analysis_dir = os.path.join(output_dir, "analysis")
    os.makedirs(analysis_dir, exist_ok=True)

    # 1. Pairwise Homolog Alignments and Distance Matrix
    aln_summary, df_dist = perform_pairwise_homolog_alignments(accession, sequence, blast_hits)
    
    # Save distance matrix CSV
    dist_csv_path = os.path.join(analysis_dir, "distance_matrix.csv")
    df_dist.to_csv(dist_csv_path)
    logger.info(f"Saved phylogenetic distance matrix to {dist_csv_path}")

    # 2. PTM Analysis and Site Mapping
    ptm_sites = []
    if ptm_features:
        for ptm in ptm_features:
            ptm_type = ptm.get("type", "Unknown")
            pos = ptm.get("position", "Unknown")
            desc = ptm.get("description", "")
            
            # Map amino acid at position if single numeric position
            aa_at_pos = None
            if pos and pos.isdigit():
                idx = int(pos) - 1
                if 0 <= idx < len(sequence):
                    aa_at_pos = sequence[idx]

            ptm_sites.append({
                "type": ptm_type,
                "position": pos,
                "residue": aa_at_pos,
                "description": desc,
            })

    # 3. PPI Network via STRING DB
    species = int(taxid) if taxid and str(taxid).isdigit() else 9606
    ppi_interactions = fetch_string_ppi_network(accession, species_taxid=species, timeout=timeout)

    comparative_results = {
        "accession": accession,
        "pairwise_alignments_top5": aln_summary,
        "distance_matrix_csv": dist_csv_path,
        "ptm_sites": ptm_sites,
        "ptm_count": len(ptm_sites),
        "ppi_interactions": ppi_interactions,
        "ppi_partner_count": len(ppi_interactions),
    }

    comp_json_path = os.path.join(analysis_dir, "comparative.json")
    with open(comp_json_path, "w", encoding="utf-8") as f:
        json.dump(comparative_results, f, indent=2)
    logger.info(f"Saved comparative analysis to {comp_json_path}")

    return comparative_results
