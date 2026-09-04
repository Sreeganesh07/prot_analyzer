"""
Stage 2: Sequence Analysis using Biopython ProtParam and InterPro REST API.
Computes physicochemical properties, amino acid percentages, motif occurrences, and domain annotations.
"""

import json
import logging
import os
import re
from typing import Any, Dict, List, Optional
import requests
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

logger = logging.getLogger("ProteinScope.Sequence")


def search_motifs(sequence: str) -> Dict[str, List[Dict[str, Any]]]:
    """
    Search sequence for key biological motifs using regular expressions:
    - N-glycosylation: N[^P][ST]
    - Signal peptide: M[^P]{15,30}[LIVMF]{4,6}
    - PKC phosphorylation: [ST]..E
    """
    motifs = {
        "n_glycosylation": r"N[^P][ST]",
        "signal_peptide": r"M[^P]{15,30}[LIVMF]{4,6}",
        "pkc_phosphorylation": r"[ST]..E",
    }
    results: Dict[str, List[Dict[str, Any]]] = {}

    for name, pattern in motifs.items():
        matches = []
        # Use finditer with overlapping consideration if needed
        for m in re.finditer(f"(?=({pattern}))", sequence):
            start = m.start()
            matched_str = m.group(1)
            end = start + len(matched_str)
            matches.append({
                "start": start + 1,  # 1-indexed for bio conventions
                "end": end,
                "pattern": pattern,
                "matched_sequence": matched_str,
            })
        results[name] = matches

    return results


def fetch_interpro_domains(accession: str, timeout: int = 15) -> List[Dict[str, Any]]:
    """
    Query InterPro REST API for functional domains and signatures.
    Endpoint: https://www.ebi.ac.uk/interpro/api/protein/UniProt/{accession}
    """
    url = f"https://www.ebi.ac.uk/interpro/api/protein/UniProt/{accession}"
    domains = []
    headers = {"Accept": "application/json"}

    try:
        response = requests.get(url, headers=headers, timeout=timeout)
        if response.status_code == 200:
            data = response.json()
            results = data.get("results", [])
            for res in results:
                metadata = res.get("metadata", {})
                entry_id = metadata.get("accession", "")
                name = metadata.get("name", "")
                entry_type = metadata.get("type", "")
                source_db = metadata.get("source_database", "")
                
                # Extract locations
                locations = []
                for loc in res.get("locations", []):
                    for fragment in loc.get("fragments", []):
                        locations.append({
                            "start": fragment.get("start"),
                            "end": fragment.get("end"),
                            "score": loc.get("score")
                        })

                domains.append({
                    "id": entry_id,
                    "name": name,
                    "type": entry_type,
                    "source_database": source_db,
                    "locations": locations
                })
        else:
            logger.warning(f"InterPro API returned status {response.status_code} for {accession}")
    except Exception as e:
        logger.error(f"Error fetching InterPro domains for {accession}: {e}")

    return domains


def compute_physicochemical_properties(sequence: str) -> Dict[str, Any]:
    """
    Compute molecular weight, isoelectric point, aromaticity, instability index,
    GRAVY, secondary structure fraction, and amino acid composition using ProtParam.
    """
    # Clean sequence of non-standard amino acids for ProtParam compatibility
    standard_aa = set("ACDEFGHIKLMNPQRSTVWY")
    cleaned_seq = "".join([aa for aa in sequence.upper() if aa in standard_aa])

    if not cleaned_seq:
        return {}

    analyser = ProteinAnalysis(cleaned_seq)
    
    # Secondary structure fraction returns tuple of (helix, turn, sheet)
    sec_struc = analyser.secondary_structure_fraction()
    # Support both property and method access across Biopython versions
    if hasattr(analyser, "amino_acids_percent"):
        aa_percent = analyser.amino_acids_percent
    elif hasattr(analyser, "get_amino_acids_percent"):
        aa_percent = analyser.get_amino_acids_percent()
    else:
        counts = analyser.count_amino_acids()
        total = sum(counts.values()) or 1
        aa_percent = {aa: (cnt / total * 100) for aa, cnt in counts.items()}

    # Format amino acid percentages to standard percentage scale (0-100%)
    aa_percent_formatted = {}
    is_fraction = sum(aa_percent.values()) <= 1.5 if aa_percent else False
    for aa, pct in aa_percent.items():
        val = pct * 100.0 if is_fraction else pct
        aa_percent_formatted[aa] = round(val, 3)

    # Calculate grouped properties
    pos_charged = aa_percent_formatted.get("R", 0) + aa_percent_formatted.get("K", 0) + aa_percent_formatted.get("H", 0)
    neg_charged = aa_percent_formatted.get("D", 0) + aa_percent_formatted.get("E", 0)
    hydrophobic = sum(aa_percent_formatted.get(aa, 0) for aa in "AILMFWV")
    aromatic = sum(aa_percent_formatted.get(aa, 0) for aa in "FWY")

    return {
        "sequence_length": len(sequence),
        "analyzed_length": len(cleaned_seq),
        "molecular_weight": round(analyser.molecular_weight(), 2),
        "isoelectric_point": round(analyser.isoelectric_point(), 2),
        "aromaticity": round(analyser.aromaticity(), 4),
        "instability_index": round(analyser.instability_index(), 2),
        "is_stable": analyser.instability_index() <= 40.0,
        "gravy": round(analyser.gravy(), 3),
        "secondary_structure_fraction": {
            "helix": round(sec_struc[0], 4),
            "turn": round(sec_struc[1], 4),
            "sheet": round(sec_struc[2], 4),
        },
        "amino_acid_percentages": aa_percent_formatted,
        "grouped_composition": {
            "positively_charged_percent": round(pos_charged, 2),
            "negatively_charged_percent": round(neg_charged, 2),
            "hydrophobic_percent": round(hydrophobic, 2),
            "aromatic_percent": round(aromatic, 2),
        }
    }


def analyze_protein_sequence(
    accession: str,
    output_dir: str,
    fasta_path: Optional[str] = None,
    timeout: int = 15
) -> Dict[str, Any]:
    """
    Stage 2 Entry point:
    Loads FASTA, computes physicochemical properties, identifies motifs, fetches InterPro domains,
    and saves sequence_properties.json.
    """
    logger.info(f"Starting Stage 2: Sequence analysis for {accession}")
    analysis_dir = os.path.join(output_dir, "analysis")
    os.makedirs(analysis_dir, exist_ok=True)

    sequence = ""
    header = ""

    if fasta_path and os.path.exists(fasta_path):
        try:
            record = SeqIO.read(fasta_path, "fasta")
            sequence = str(record.seq).upper()
            header = record.description
        except Exception as e:
            logger.error(f"Error reading FASTA from {fasta_path}: {e}")

    if not sequence:
        # Try finding in default sequences folder
        default_fasta = os.path.join(output_dir, "sequences", f"{accession}.fasta")
        if os.path.exists(default_fasta):
            try:
                record = SeqIO.read(default_fasta, "fasta")
                sequence = str(record.seq).upper()
                header = record.description
            except Exception as e:
                logger.error(f"Error reading FASTA from {default_fasta}: {e}")

    if not sequence:
        logger.warning(f"No sequence found to analyze for accession {accession}")
        return {"accession": accession, "error": "Sequence not found"}

    # 1. Physicochemical analysis
    properties = compute_physicochemical_properties(sequence)

    # 2. Motif search
    motifs = search_motifs(sequence)

    # 3. InterPro domain annotations
    domains = fetch_interpro_domains(accession, timeout=timeout)

    result = {
        "accession": accession,
        "header": header,
        "sequence_length": len(sequence),
        "physicochemical_properties": properties,
        "motifs": motifs,
        "domains": domains,
        "domain_count": len(domains),
    }

    # Save to analysis/sequence_properties.json
    out_file = os.path.join(analysis_dir, "sequence_properties.json")
    with open(out_file, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)
    logger.info(f"Saved sequence properties to {out_file}")

    return result
