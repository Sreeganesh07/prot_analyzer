"""
Unit tests for Stage 1: Data Ingestion & Input Validation.
"""

import os
import shutil
import tempfile
import unittest
from unittest.mock import MagicMock, patch

from proteinscope.stages.ingest import (
    get_fasta,
    get_json,
    get_xml,
    ingest_uniprot_data,
    parse_uniprot_xml,
    validate_accession,
)


class TestIngest(unittest.TestCase):

    def setUp(self):
        self.test_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.test_dir, ignore_errors=True)

    def test_validate_accession(self):
        # Valid accessions
        self.assertTrue(validate_accession("P0DTC2"))
        self.assertTrue(validate_accession("P36334"))
        self.assertTrue(validate_accession("Q9BYF1"))
        self.assertTrue(validate_accession("A0A024RBG1"))

        # Invalid accessions
        self.assertFalse(validate_accession(""))
        self.assertFalse(validate_accession("INVALID_ACCESSION"))
        self.assertFalse(validate_accession("12345"))
        self.assertFalse(validate_accession("P123"))

    def test_parse_uniprot_xml(self):
        sample_xml = """<?xml version="1.0" encoding="UTF-8"?>
        <uniprot xmlns="http://uniprot.org/uniprot">
            <entry>
                <protein>
                    <recommendedName>
                        <fullName>Spike glycoprotein</fullName>
                    </recommendedName>
                </protein>
                <gene>
                    <name type="primary">S</name>
                </gene>
                <organism>
                    <name type="scientific">Severe acute respiratory syndrome coronavirus 2</name>
                    <dbReference type="NCBI Taxonomy" id="2697049"/>
                </organism>
                <comment type="function">
                    <text>Attaches the virion to the host cell membrane receptor ACE2.</text>
                </comment>
                <feature type="glycosylation site" description="N-linked (GlcNAc...)">
                    <location>
                        <position position="165"/>
                    </location>
                </feature>
                <dbReference type="PDB" id="6VXX">
                    <property type="method" value="EM"/>
                    <property type="resolution" value="2.80 A"/>
                </dbReference>
            </entry>
        </uniprot>"""

        parsed = parse_uniprot_xml(sample_xml)
        self.assertEqual(parsed["protein_name"], "Spike glycoprotein")
        self.assertEqual(parsed["organism"], "Severe acute respiratory syndrome coronavirus 2")
        self.assertEqual(parsed["taxid"], "2697049")
        self.assertIn("S", parsed["gene_names"])
        self.assertEqual(len(parsed["ptm_features"]), 1)
        self.assertEqual(parsed["ptm_features"][0]["position"], "165")
        self.assertEqual(len(parsed["pdb_cross_references"]), 1)
        self.assertEqual(parsed["pdb_cross_references"][0]["pdb_id"], "6VXX")

    @patch("requests.get")
    def test_ingest_uniprot_data_mocked(self, mock_get):
        mock_fasta = ">sp|P0DTC2|SPIKE_SARS2 Spike glycoprotein\nMFVFLVLLPLVSSQCVNLT\n"
        mock_json = {"primaryAccession": "P0DTC2", "uniProtkbId": "SPIKE_SARS2"}
        mock_xml = """<uniprot xmlns="http://uniprot.org/uniprot"><entry><protein><recommendedName><fullName>Spike</fullName></recommendedName></protein></entry></uniprot>"""

        def side_effect(url, **kwargs):
            mock_resp = MagicMock()
            mock_resp.status_code = 200
            if url.endswith(".fasta"):
                mock_resp.text = mock_fasta
            elif url.endswith(".json"):
                mock_resp.json.return_value = mock_json
            elif url.endswith(".xml"):
                mock_resp.text = mock_xml
            return mock_resp

        mock_get.side_effect = side_effect

        data = ingest_uniprot_data("P0DTC2", self.test_dir)
        self.assertEqual(data["accession"], "P0DTC2")
        self.assertEqual(data["sequence"], "MFVFLVLLPLVSSQCVNLT")
        self.assertEqual(data["sequence_length"], 19)
        self.assertTrue(os.path.exists(data["fasta_path"]))


if __name__ == "__main__":
    unittest.main()
