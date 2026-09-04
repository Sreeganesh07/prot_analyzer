"""
Stage 1: Data Ingestion from UniProt REST API.
Retrieves FASTA, XML, and JSON metadata for a given accession.
"""

import json
import logging
import os
import re
import xml.etree.ElementTree as ET
from typing import Any, Dict, Optional
import requests

logger = logging.getLogger("ProteinScope.Ingest")


def validate_accession(accession: str) -> bool:
    """
    Validate UniProt or NCBI accession format.
    Standard UniProt: [O,P,Q][0-9][A-Z0-9]{3}[0-9] or 6/10-char format
    NCBI RefSeq: [A-Z]{2}_[0-9]{5,} (e.g., NP_000537, XP_012345)
    NCBI GenBank: [A-Z]{3}[0-9]{5} (e.g., AAA12345)
    """
    if not accession or not isinstance(accession, str):
        return False
    accession = accession.strip().upper()
    uniprot_regex = re.compile(
        r"^([O,P,Q][0-9][A-Z0-9]{3}[0-9]|[A-N,R-Z][0-9][A-Z0-9]{3}[0-9]([A-Z][A-Z0-9]{2}[0-9])?|[A-Z][0-9][A-Z0-9]{3}[0-9]([A-Z0-9]{4})?)$"
    )
    ncbi_refseq_regex = re.compile(r"^[A-Z]{2}_[0-9]{5,}(\.[0-9]+)?$")
    ncbi_genbank_regex = re.compile(r"^[A-Z]{3}[0-9]{5}(\.[0-9]+)?$")

    return bool(
        uniprot_regex.match(accession)
        or ncbi_refseq_regex.match(accession)
        or ncbi_genbank_regex.match(accession)
    )


def fetch_ncbi_protein(query: str, timeout: int = 15) -> Optional[Dict[str, Any]]:
    """
    Query the NCBI Protein Entrez database via E-utilities API.
    Retrieves FASTA sequence and descriptive metadata for any NCBI accession,
    RefSeq identifier (e.g., NP_000537), or protein name.
    """
    clean_query = query.strip()
    try:
        search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        params = {
            "db": "protein",
            "term": clean_query,
            "retmode": "json",
            "retmax": 1,
        }
        resp = requests.get(search_url, params=params, timeout=timeout)
        if resp.status_code != 200:
            logger.warning(f"NCBI esearch returned HTTP {resp.status_code} for query '{clean_query}'")
            return None

        search_data = resp.json()
        id_list = search_data.get("esearchresult", {}).get("idlist", [])
        if not id_list:
            logger.warning(f"No NCBI protein entries found for query '{clean_query}'")
            return None

        ncbi_id = id_list[0]

        fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        fetch_params = {
            "db": "protein",
            "id": ncbi_id,
            "rettype": "fasta",
            "retmode": "text",
        }
        fasta_resp = requests.get(fetch_url, params=fetch_params, timeout=timeout)
        if fasta_resp.status_code != 200 or not fasta_resp.text.startswith(">"):
            logger.warning(f"NCBI efetch FASTA failed for ID {ncbi_id}")
            return None

        fasta_text = fasta_resp.text.strip()
        lines = fasta_text.split("\n")
        header = lines[0] if lines else ""
        sequence = "".join(lines[1:]).replace(" ", "").upper()

        summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        sum_params = {"db": "protein", "id": ncbi_id, "retmode": "json"}
        sum_resp = requests.get(summary_url, params=sum_params, timeout=timeout)

        protein_name = clean_query
        organism = "Unknown"
        acc = clean_query
        taxid = None

        if sum_resp.status_code == 200:
            try:
                sum_json = sum_resp.json()
                res_obj = sum_json.get("result", {}).get(str(ncbi_id), {})
                title = res_obj.get("title", "")
                if title:
                    protein_name = title.split("[")[0].strip()
                    if "[" in title and "]" in title:
                        organism = title.split("[")[-1].split("]")[0].strip()
                acc = res_obj.get("caption", clean_query)
                taxid = str(res_obj.get("taxid")) if res_obj.get("taxid") else None
            except Exception as e:
                logger.debug(f"Error parsing NCBI summary JSON: {e}")

        return {
            "accession": acc,
            "ncbi_id": ncbi_id,
            "header": header,
            "sequence": sequence,
            "protein_name": protein_name,
            "organism": organism,
            "taxid": taxid,
            "fasta_text": fasta_text,
        }
    except Exception as e:
        logger.error(f"Error during NCBI protein retrieval for '{clean_query}': {e}")
        return None


def get_fasta(accession: str, timeout: int = 15) -> Optional[str]:
    """Fetch raw FASTA sequence from UniProt REST API."""
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    try:
        response = requests.get(url, timeout=timeout)
        if response.status_code == 200 and response.text.startswith(">"):
            return response.text
        logger.warning(f"UniProt FASTA returned status code {response.status_code} for {accession}")
    except Exception as e:
        logger.error(f"Error fetching FASTA for {accession}: {e}")
    return None


def get_xml(accession: str, timeout: int = 15) -> Optional[str]:
    """Fetch raw UniProt XML content."""
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.xml"
    try:
        response = requests.get(url, timeout=timeout)
        if response.status_code == 200:
            return response.text
        logger.warning(f"UniProt XML returned status code {response.status_code} for {accession}")
    except Exception as e:
        logger.error(f"Error fetching XML for {accession}: {e}")
    return None


def get_json(accession: str, timeout: int = 15) -> Optional[Dict[str, Any]]:
    """Fetch raw UniProt JSON metadata."""
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    try:
        response = requests.get(url, timeout=timeout)
        if response.status_code == 200:
            return response.json()
        logger.warning(f"UniProt JSON returned status code {response.status_code} for {accession}")
    except Exception as e:
        logger.error(f"Error fetching JSON for {accession}: {e}")
    return None


def get_metadata(accession: str, timeout: int = 15) -> Dict[str, Any]:
    """
    Fetch and consolidate UniProt metadata (XML and JSON) for an accession.
    Returns parsed dictionary containing organism, function, PTMs, isoforms, and full metadata.
    """
    xml_content = get_xml(accession, timeout=timeout)
    json_content = get_json(accession, timeout=timeout)
    parsed = parse_uniprot_xml(xml_content) if xml_content else {}
    parsed["raw_json"] = json_content
    return parsed


def parse_uniprot_xml(xml_content: str) -> Dict[str, Any]:
    """Parse key biological features, organism, function, PTMs, and isoforms from XML."""
    parsed: Dict[str, Any] = {
        "organism": "Unknown",
        "taxid": None,
        "gene_names": [],
        "protein_name": "Unknown",
        "function_comments": [],
        "ptm_features": [],
        "isoforms": [],
        "pdb_cross_references": [],
    }

    if not xml_content:
        return parsed

    try:
        root = ET.fromstring(xml_content)
        # XML namespace in UniProt
        ns = {"up": "http://uniprot.org/uniprot"}
        entry = root.find("up:entry", ns)
        if entry is None:
            entry = root

        # Organism
        org_elem = entry.find("up:organism/up:name[@type='scientific']", ns)
        if org_elem is not None and org_elem.text:
            parsed["organism"] = org_elem.text
        
        tax_elem = entry.find("up:organism/up:dbReference[@type='NCBI Taxonomy']", ns)
        if tax_elem is not None:
            parsed["taxid"] = tax_elem.attrib.get("id")

        # Protein Recommended Name
        rec_name = entry.find("up:protein/up:recommendedName/up:fullName", ns)
        if rec_name is not None and rec_name.text:
            parsed["protein_name"] = rec_name.text
        else:
            sub_name = entry.find("up:protein/up:submittedName/up:fullName", ns)
            if sub_name is not None and sub_name.text:
                parsed["protein_name"] = sub_name.text

        # Gene Names
        for gene in entry.findall("up:gene/up:name", ns):
            if gene.text and gene.text not in parsed["gene_names"]:
                parsed["gene_names"].append(gene.text)

        # Comments: Function
        for comment in entry.findall("up:comment[@type='function']/up:text", ns):
            if comment.text:
                parsed["function_comments"].append(comment.text.strip())

        # Comments: Isoforms
        for iso in entry.findall("up:comment[@type='alternative products']/up:isoform/up:name", ns):
            if iso.text:
                parsed["isoforms"].append(iso.text.strip())

        # Features: PTMs, Modified residues, Glycosylation, Disulfide bonds
        ptm_types = {"modified residue", "glycosylation site", "disulfide bond", "cross-link", "phosphorylation"}
        for feat in entry.findall("up:feature", ns):
            ft_type = feat.attrib.get("type", "").lower()
            if any(t in ft_type for t in ptm_types) or ft_type in ptm_types:
                desc = feat.attrib.get("description", "")
                loc = feat.find("up:location/up:position", ns)
                pos = loc.attrib.get("position") if loc is not None else None
                if not pos:
                    begin = feat.find("up:location/up:begin", ns)
                    end = feat.find("up:location/up:end", ns)
                    if begin is not None and end is not None:
                        pos = f"{begin.attrib.get('position')}-{end.attrib.get('position')}"
                
                parsed["ptm_features"].append({
                    "type": ft_type,
                    "description": desc,
                    "position": pos or "Unknown"
                })

        # PDB Cross References
        for db_ref in entry.findall("up:dbReference[@type='PDB']", ns):
            pdb_id = db_ref.attrib.get("id")
            method = None
            resolution = None
            for prop in db_ref.findall("up:property", ns):
                if prop.attrib.get("type") == "method":
                    method = prop.attrib.get("value")
                elif prop.attrib.get("type") == "resolution":
                    resolution = prop.attrib.get("value")
            if pdb_id:
                parsed["pdb_cross_references"].append({
                    "pdb_id": pdb_id,
                    "method": method,
                    "resolution": resolution
                })

    except Exception as e:
        logger.error(f"Error parsing UniProt XML: {e}")

    return parsed


def ingest_uniprot_data(accession: str, output_dir: str, timeout: int = 15) -> Dict[str, Any]:
    """
    Stage 1 Entry point:
    Fetches FASTA, XML, and JSON from UniProt, saves raw files, and returns parsed metadata.
    """
    accession = accession.strip().upper()
    logger.info(f"Starting Stage 1: Ingesting data for UniProt accession {accession}")

    seq_dir = os.path.join(output_dir, "sequences")
    raw_dir = os.path.join(output_dir, "raw_data")
    os.makedirs(seq_dir, exist_ok=True)
    os.makedirs(raw_dir, exist_ok=True)

    fasta_content = get_fasta(accession, timeout=timeout)
    xml_content = get_xml(accession, timeout=timeout)
    json_content = get_json(accession, timeout=timeout)
    parsed_xml = parse_uniprot_xml(xml_content) if xml_content else {}

    # If not found directly in UniProt, search and fetch from NCBI Protein database
    if not fasta_content:
        logger.info(f"Accession {accession} not directly found in UniProt; searching NCBI Protein Entrez database...")
        ncbi_data = fetch_ncbi_protein(accession, timeout=timeout)
        if ncbi_data and ncbi_data.get("fasta_text"):
            fasta_content = ncbi_data["fasta_text"]
            logger.info(f"Retrieved sequence for {accession} from NCBI ({ncbi_data.get('protein_name')})")
            
            # Check UniProt KB for potential structural/PTM cross-references
            try:
                up_search = requests.get(
                    f"https://rest.uniprot.org/uniprotkb/search?query={accession}&format=json&size=1",
                    timeout=timeout
                )
                if up_search.status_code == 200:
                    results = up_search.json().get("results", [])
                    if results:
                        alt_acc = results[0].get("primaryAccession")
                        if alt_acc:
                            logger.info(f"Found corresponding UniProt cross-reference: {alt_acc}")
                            xml_content = get_xml(alt_acc, timeout=timeout)
                            json_content = get_json(alt_acc, timeout=timeout)
                            if xml_content:
                                parsed_xml = parse_uniprot_xml(xml_content)
            except Exception as e:
                logger.debug(f"UniProt cross-reference search failed: {e}")

            if not parsed_xml.get("protein_name") or parsed_xml.get("protein_name") == "Unknown":
                parsed_xml["protein_name"] = ncbi_data.get("protein_name", accession)
            if not parsed_xml.get("organism") or parsed_xml.get("organism") == "Unknown":
                parsed_xml["organism"] = ncbi_data.get("organism", "Unknown")
            if not parsed_xml.get("taxid") and ncbi_data.get("taxid"):
                parsed_xml["taxid"] = ncbi_data.get("taxid")

    fasta_path = None
    if fasta_content:
        fasta_path = os.path.join(seq_dir, f"{accession}.fasta")
        raw_fasta_path = os.path.join(raw_dir, f"{accession}.fasta")
        with open(fasta_path, "w", encoding="utf-8") as f:
            f.write(fasta_content)
        with open(raw_fasta_path, "w", encoding="utf-8") as f:
            f.write(fasta_content)
        logger.info(f"Saved FASTA to {fasta_path}")

    xml_path = None
    if xml_content:
        xml_path = os.path.join(raw_dir, f"{accession}.xml")
        with open(xml_path, "w", encoding="utf-8") as f:
            f.write(xml_content)
        logger.info(f"Saved XML to {xml_path}")

    json_path = None
    if json_content:
        json_path = os.path.join(raw_dir, f"{accession}.json")
        with open(json_path, "w", encoding="utf-8") as f:
            json.dump(json_content, f, indent=2)
        logger.info(f"Saved JSON to {json_path}")

    if not parsed_xml and xml_content:
        parsed_xml = parse_uniprot_xml(xml_content)

    # Extract clean sequence from FASTA if available
    sequence = ""
    header = ""
    if fasta_content:
        lines = fasta_content.strip().split("\n")
        header = lines[0] if lines else ""
        sequence = "".join(lines[1:]).replace(" ", "").upper()

    return {
        "accession": accession,
        "header": header,
        "sequence": sequence,
        "sequence_length": len(sequence),
        "fasta_path": fasta_path,
        "xml_path": xml_path,
        "json_path": json_path,
        "organism": parsed_xml.get("organism", "Unknown"),
        "taxid": parsed_xml.get("taxid"),
        "protein_name": parsed_xml.get("protein_name", "Unknown"),
        "gene_names": parsed_xml.get("gene_names", []),
        "function_comments": parsed_xml.get("function_comments", []),
        "ptm_features": parsed_xml.get("ptm_features", []),
        "isoforms": parsed_xml.get("isoforms", []),
        "pdb_cross_references": parsed_xml.get("pdb_cross_references", []),
    }
