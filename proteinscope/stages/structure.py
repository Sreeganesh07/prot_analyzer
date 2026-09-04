"""
Stage 4: 3D Structure Retrieval through 3 parallel pathways.
Path A: Experimental structures from RCSB PDB
Path B: Predicted structures from AlphaFold DB
Path C: Homology models from SWISS-MODEL Repository
"""

import json
import logging
import os
from typing import Any, Dict, List, Optional
import requests

from Bio.PDB import PDBList

logger = logging.getLogger("ProteinScope.Structure")


def retrieve_rcsb_experimental(
    accession: str,
    output_dir: str,
    pdb_cross_refs: Optional[List[Dict[str, Any]]] = None,
    timeout: int = 15
) -> Optional[Dict[str, Any]]:
    """
    Path A: Retrieve best-resolution experimental structure from RCSB PDB.
    """
    exp_dir = os.path.join(output_dir, "structures", "experimental")
    os.makedirs(exp_dir, exist_ok=True)

    pdb_candidates = []

    # Check cross-references from UniProt XML if available
    if pdb_cross_refs:
        for ref in pdb_cross_refs:
            pdb_id = ref.get("pdb_id")
            res_str = ref.get("resolution", "")
            try:
                # Parse numeric resolution (e.g. '2.80 A' -> 2.80)
                res_val = float(res_str.replace("A", "").strip()) if res_str else 999.0
            except ValueError:
                res_val = 999.0
            if pdb_id:
                pdb_candidates.append({"pdb_id": pdb_id.lower(), "resolution": res_val, "method": ref.get("method", "")})

    # If no cross refs, query RCSB search API
    if not pdb_candidates:
        search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
        query_payload = {
            "query": {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                    "operator": "exact_match",
                    "value": accession
                }
            },
            "return_type": "entry",
            "request_options": {
                "paginate": {"start": 0, "rows": 10}
            }
        }
        try:
            resp = requests.post(search_url, json=query_payload, timeout=timeout)
            if resp.status_code == 200:
                result_json = resp.json()
                for hit in result_json.get("result_set", []):
                    pdb_id = hit.get("identifier", "").lower()
                    if pdb_id:
                        pdb_candidates.append({"pdb_id": pdb_id, "resolution": 999.0, "method": "PDB Search"})
        except Exception as e:
            logger.warning(f"RCSB PDB search query failed for {accession}: {e}")

    if not pdb_candidates:
        logger.info(f"No experimental PDB structures found for {accession}")
        return None

    # Sort by resolution (best/lowest resolution first)
    pdb_candidates.sort(key=lambda x: x["resolution"])
    best_candidate = pdb_candidates[0]
    best_pdb_id = best_candidate["pdb_id"].upper()

    logger.info(f"Selected best experimental PDB candidate: {best_pdb_id} (res: {best_candidate['resolution']} A)")

    # Download PDB file
    pdb_filename = f"{best_pdb_id.lower()}.pdb"
    saved_path = os.path.join(exp_dir, pdb_filename)

    # 1. Try direct HTTPS download from files.rcsb.org
    try:
        download_url = f"https://files.rcsb.org/download/{best_pdb_id}.pdb"
        r = requests.get(download_url, timeout=timeout)
        if r.status_code == 200 and len(r.content) > 100:
            with open(saved_path, "wb") as f:
                f.write(r.content)
            logger.info(f"Downloaded experimental PDB structure to {saved_path}")
            return {
                "source": "Experimental (RCSB PDB)",
                "pdb_id": best_pdb_id,
                "file_path": saved_path,
                "resolution": best_candidate["resolution"],
                "method": best_candidate["method"],
                "file_size": len(r.content),
            }
    except Exception as e:
        logger.warning(f"Direct download for PDB {best_pdb_id} failed: {e}")

    # 2. Fallback to Bio.PDB.PDBList if available
    if PDBList is not None:
        try:
            pdbl = PDBList(verbose=False)
            fetched_file = pdbl.retrieve_pdb_file(best_pdb_id, pdir=exp_dir, file_format="pdb")
            if fetched_file and os.path.exists(fetched_file):
                return {
                    "source": "Experimental (RCSB PDB)",
                    "pdb_id": best_pdb_id,
                    "file_path": fetched_file,
                    "resolution": best_candidate["resolution"],
                    "method": best_candidate["method"],
                    "file_size": os.path.getsize(fetched_file),
                }
        except Exception as e:
            logger.error(f"Biopython PDBList download failed for {best_pdb_id}: {e}")

    return None


def retrieve_alphafold_model(
    accession: str,
    output_dir: str,
    timeout: int = 15
) -> Optional[Dict[str, Any]]:
    """
    Path B: Retrieve predicted AlphaFold model.
    URL: https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v4.pdb
    """
    af_dir = os.path.join(output_dir, "structures", "alphafold")
    os.makedirs(af_dir, exist_ok=True)
    af_path = os.path.join(af_dir, f"AF-{accession}.pdb")

    # 1. Try direct v4 URL
    v4_url = f"https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v4.pdb"
    try:
        r = requests.get(v4_url, timeout=timeout)
        if r.status_code == 200 and len(r.content) > 100:
            with open(af_path, "wb") as f:
                f.write(r.content)
            logger.info(f"Downloaded AlphaFold model to {af_path}")
            return {
                "source": "AlphaFold DB",
                "model_id": f"AF-{accession}-F1",
                "file_path": af_path,
                "version": "v4",
                "file_size": len(r.content),
            }
    except Exception as e:
        logger.warning(f"Direct AlphaFold v4 fetch failed for {accession}: {e}")

    # 2. Query AlphaFold API
    try:
        api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{accession}"
        resp = requests.get(api_url, timeout=timeout)
        if resp.status_code == 200:
            data = resp.json()
            if data and isinstance(data, list) and len(data) > 0:
                pdb_url = data[0].get("pdbUrl")
                if pdb_url:
                    r2 = requests.get(pdb_url, timeout=timeout)
                    if r2.status_code == 200:
                        with open(af_path, "wb") as f:
                            f.write(r2.content)
                        logger.info(f"Downloaded AlphaFold model via API to {af_path}")
                        return {
                            "source": "AlphaFold DB",
                            "model_id": data[0].get("entryId", f"AF-{accession}"),
                            "file_path": af_path,
                            "plddt_url": data[0].get("plddtReportUrl"),
                            "file_size": len(r2.content),
                        }
    except Exception as e:
        logger.warning(f"AlphaFold API fetch failed for {accession}: {e}")

    return None


def retrieve_swissmodel_homology(
    accession: str,
    output_dir: str,
    timeout: int = 15
) -> Optional[Dict[str, Any]]:
    """
    Path C: Retrieve homology model from SWISS-MODEL repository.
    URL: https://swissmodel.expasy.org/repository/uniprot/{accession}.json
    """
    homology_dir = os.path.join(output_dir, "structures", "homology")
    os.makedirs(homology_dir, exist_ok=True)
    homology_path = os.path.join(homology_dir, f"{accession}_model.pdb")

    repo_url = f"https://swissmodel.expasy.org/repository/uniprot/{accession}.json"
    headers = {"Accept": "application/json"}

    try:
        resp = requests.get(repo_url, headers=headers, timeout=timeout)
        if resp.status_code == 200:
            data = resp.json()
            structures = data.get("result", {}).get("structures", [])
            if structures:
                # Sort models by GMQE score (higher is better)
                structures.sort(key=lambda s: s.get("gmqe", 0.0), reverse=True)
                best_struct = structures[0]
                coordinates_url = best_struct.get("coordinates")
                gmqe = best_struct.get("gmqe")
                qmean = best_struct.get("qmean", {}).get("qmean4_norm")
                template = best_struct.get("template")

                if coordinates_url:
                    coord_resp = requests.get(coordinates_url, timeout=timeout)
                    if coord_resp.status_code == 200:
                        with open(homology_path, "wb") as f:
                            f.write(coord_resp.content)
                        logger.info(f"Downloaded SWISS-MODEL homology structure to {homology_path}")
                        return {
                            "source": "SWISS-MODEL Homology",
                            "template": template,
                            "gmqe_score": gmqe,
                            "qmean_score": qmean,
                            "file_path": homology_path,
                            "file_size": len(coord_resp.content),
                        }
    except Exception as e:
        logger.warning(f"SWISS-MODEL repository query failed for {accession}: {e}")

    return None


def retrieve_protein_structures(
    accession: str,
    output_dir: str,
    pdb_cross_refs: Optional[List[Dict[str, Any]]] = None,
    timeout: int = 15
) -> Dict[str, Any]:
    """
    Stage 4 Entry point:
    Executes Path A, Path B, and Path C structure retrievals, evaluates results,
    and designates the best available structure (Experimental > AlphaFold > Homology).
    """
    logger.info(f"Starting Stage 4: Structure retrieval for accession {accession}")
    struct_summary: Dict[str, Any] = {
        "accession": accession,
        "experimental": None,
        "alphafold": None,
        "homology": None,
        "best_structure": None,
    }

    # Path A: Experimental RCSB PDB
    exp_res = retrieve_rcsb_experimental(accession, output_dir, pdb_cross_refs=pdb_cross_refs, timeout=timeout)
    struct_summary["experimental"] = exp_res

    # Path B: AlphaFold DB
    af_res = retrieve_alphafold_model(accession, output_dir, timeout=timeout)
    struct_summary["alphafold"] = af_res

    # Path C: SWISS-MODEL
    sm_res = retrieve_swissmodel_homology(accession, output_dir, timeout=timeout)
    struct_summary["homology"] = sm_res

    # Priority hierarchy: Experimental > AlphaFold > Homology
    if exp_res and os.path.exists(exp_res["file_path"]):
        struct_summary["best_structure"] = exp_res
    elif af_res and os.path.exists(af_res["file_path"]):
        struct_summary["best_structure"] = af_res
    elif sm_res and os.path.exists(sm_res["file_path"]):
        struct_summary["best_structure"] = sm_res
    else:
        logger.warning(f"No 3D structures could be retrieved for {accession}")

    # Save summary
    summary_path = os.path.join(output_dir, "structures", "structures_manifest.json")
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(struct_summary, f, indent=2)
    logger.info(f"Saved structure manifest to {summary_path}")

    return struct_summary
