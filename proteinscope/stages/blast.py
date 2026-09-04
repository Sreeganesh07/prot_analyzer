"""
Stage 3: Similarity Search using NCBI BLASTp and PSI-BLAST.
Performs homology searching, parses BLAST XML, extracts top hits, and interprets e-values and identity zones.
"""

import io
import json
import logging
import os
from typing import Any, Dict, List, Optional

from Bio.Blast import NCBIWWW, NCBIXML

logger = logging.getLogger("ProteinScope.BLAST")


def interpret_blast_hit(evalue: float, identity_pct: float) -> Dict[str, str]:
    """
    Interpret biological significance based on E-value threshold and Sequence Identity zone.
    Zones:
      - Safe homology zone: Identity > 30%
      - Twilight zone: 20% <= Identity <= 30%
      - Midnight zone: Identity < 20%
    Confidence:
      - High confidence homology: E-value < 1e-10
      - Moderate homology: 1e-10 <= E-value <= 1e-3
      - Low / Insignificant: E-value > 1e-3
    """
    if evalue < 1e-10:
        conf = "High confidence homology (statistically significant)"
    elif evalue <= 1e-3:
        conf = "Moderate confidence homology (likely related superfamily)"
    else:
        conf = "Weak / Low confidence (possible spurious alignment)"

    if identity_pct >= 30.0:
        zone = "Safe homology zone (common structural fold guaranteed)"
    elif identity_pct >= 20.0:
        zone = "Twilight zone (structural similarity possible, requires structural verification)"
    else:
        zone = "Midnight zone (low sequence similarity)"

    return {
        "confidence_level": conf,
        "identity_zone": zone,
    }


def parse_blast_xml_content(xml_content: str, query_len: int, max_hits: int = 10) -> List[Dict[str, Any]]:
    """Parse BLAST XML string into structured list of top homolog hits."""
    hits = []
    if not xml_content or not xml_content.strip():
        return hits

    # 1. Try standard Biopython NCBIXML parser if available
    if NCBIXML is not None:
        try:
            blast_records = NCBIXML.parse(io.StringIO(xml_content))
            for record in blast_records:
                for alignment in record.alignments:
                    if len(hits) >= max_hits:
                        break
                
                # Best HSP (High-scoring Segment Pair) for this alignment
                if alignment.hsps:
                    best_hsp = alignment.hsps[0]
                    hsp_len = best_hsp.align_length
                    identities = best_hsp.identities
                    positives = best_hsp.positives
                    identity_pct = round((identities / hsp_len * 100), 2) if hsp_len > 0 else 0.0
                    positives_pct = round((positives / hsp_len * 100), 2) if hsp_len > 0 else 0.0
                    coverage_pct = round((hsp_len / query_len * 100), 2) if query_len > 0 else 0.0

                    interp = interpret_blast_hit(best_hsp.expect, identity_pct)

                    # Extract accession from hit title / id
                    hit_id = alignment.hit_id
                    hit_def = alignment.hit_def
                    acc = alignment.accession if hasattr(alignment, "accession") and alignment.accession else hit_id

                    hits.append({
                        "rank": len(hits) + 1,
                        "hit_id": hit_id,
                        "accession": acc,
                        "title": hit_def,
                        "length": alignment.length,
                        "score": round(best_hsp.score, 2),
                        "bit_score": round(best_hsp.bits, 2),
                        "evalue": best_hsp.expect,
                        "evalue_formatted": f"{best_hsp.expect:.2e}",
                        "identities": identities,
                        "positives": positives,
                        "alignment_length": hsp_len,
                        "identity_percent": identity_pct,
                        "positives_percent": positives_pct,
                        "query_coverage_percent": coverage_pct,
                        "query_sequence": best_hsp.query,
                        "match_string": best_hsp.match,
                        "sbjct_sequence": best_hsp.sbjct,
                        "interpretation": interp,
                    })
        except Exception as e:
            logger.warning(f"Biopython NCBIXML parser encountered issue, trying direct XML parser: {e}")

    # 2. Fallback to ElementTree if NCBIXML returned 0 hits
    if not hits:
        try:
            import xml.etree.ElementTree as ET
            root = ET.fromstring(xml_content)
            for hit_elem in root.findall(".//Hit"):
                if len(hits) >= max_hits:
                    break
                hit_id = hit_elem.findtext("Hit_id", "")
                hit_def = hit_elem.findtext("Hit_def", "")
                hit_acc = hit_elem.findtext("Hit_accession", hit_id)
                hit_len = int(hit_elem.findtext("Hit_len", "0"))

                hsp = hit_elem.find(".//Hsp")
                if hsp is not None:
                    bit_score = float(hsp.findtext("Hsp_bit-score", "0"))
                    score = float(hsp.findtext("Hsp_score", "0"))
                    evalue = float(hsp.findtext("Hsp_evalue", "10.0"))
                    identities = int(hsp.findtext("Hsp_identity", "0"))
                    positives = int(hsp.findtext("Hsp_positive", "0"))
                    align_len = int(hsp.findtext("Hsp_align-len", "0"))
                    qseq = hsp.findtext("Hsp_qseq", "")
                    hseq = hsp.findtext("Hsp_hseq", "")
                    midline = hsp.findtext("Hsp_midline", "")

                    identity_pct = round((identities / align_len * 100), 2) if align_len > 0 else 0.0
                    positives_pct = round((positives / align_len * 100), 2) if align_len > 0 else 0.0
                    coverage_pct = round((align_len / query_len * 100), 2) if query_len > 0 else 0.0
                    interp = interpret_blast_hit(evalue, identity_pct)

                    hits.append({
                        "rank": len(hits) + 1,
                        "hit_id": hit_id,
                        "accession": hit_acc,
                        "title": hit_def,
                        "length": hit_len,
                        "score": round(score, 2),
                        "bit_score": round(bit_score, 2),
                        "evalue": evalue,
                        "evalue_formatted": f"{evalue:.2e}",
                        "identities": identities,
                        "positives": positives,
                        "alignment_length": align_len,
                        "identity_percent": identity_pct,
                        "positives_percent": positives_pct,
                        "query_coverage_percent": coverage_pct,
                        "query_sequence": qseq,
                        "match_string": midline,
                        "sbjct_sequence": hseq,
                        "interpretation": interp,
                    })
        except Exception as e:
            logger.error(f"Error in fallback BLAST XML ElementTree parsing: {e}")

    return hits


def run_blast_search(
    sequence: str,
    output_dir: str,
    max_hits: int = 10,
    run_online: bool = True,
    matrix_name: str = "BLOSUM62",
    query_max_len: int = 500
) -> Dict[str, Any]:
    """
    Stage 3 Entry point:
    Submits query sequence (first query_max_len residues) to NCBI QBLAST.
    Saves raw XML and parsed JSON results.
    """
    logger.info(f"Starting Stage 3: BLAST similarity search (max_hits={max_hits}, online={run_online})")
    analysis_dir = os.path.join(output_dir, "analysis")
    os.makedirs(analysis_dir, exist_ok=True)

    xml_path = os.path.join(analysis_dir, "blast_results.xml")
    json_path = os.path.join(analysis_dir, "blast_results.json")

    query_seq = sequence[:query_max_len].strip().upper()
    query_len = len(query_seq)

    xml_content = ""

    # If previous XML results already exist in output directory, reuse them
    if os.path.exists(xml_path) and os.path.getsize(xml_path) > 0:
        logger.info(f"Found existing BLAST XML at {xml_path}, loading...")
        with open(xml_path, "r", encoding="utf-8") as f:
            xml_content = f.read()

    if not xml_content and run_online:
        if NCBIWWW is not None:
            logger.info(f"Submitting BLASTp query (length {query_len}) to NCBI QBLAST...")
            try:
                # Query NCBI QBLAST service
                result_handle = NCBIWWW.qblast(
                    program="blastp",
                    database="nr",
                    sequence=query_seq,
                    hitlist_size=max_hits,
                    matrix_name=matrix_name,
                    expect=10.0,
                )
                xml_content = result_handle.read()
                result_handle.close()

                with open(xml_path, "w", encoding="utf-8") as f:
                    f.write(xml_content)
                logger.info(f"Saved raw BLAST XML to {xml_path}")
            except Exception as e:
                logger.error(f"NCBI BLAST search failed or timed out: {e}")
        else:
            logger.warning("Biopython Bio.Blast (NCBIWWW) is unavailable; skipping online QBLAST.")

    # Parse hits
    hits = parse_blast_xml_content(xml_content, query_len=query_len, max_hits=max_hits)

    results_data = {
        "program": "blastp",
        "database": "nr",
        "matrix": matrix_name,
        "query_length_submitted": query_len,
        "total_hits_found": len(hits),
        "xml_path": xml_path if os.path.exists(xml_path) else None,
        "hits": hits,
    }

    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(results_data, f, indent=2)
    logger.info(f"Saved BLAST results JSON to {json_path}")

    return results_data
