"""
Unit tests for Stage 3: BLAST Similarity Search.
"""

import os
import shutil
import tempfile
import unittest

from proteinscope.stages.blast import (
    interpret_blast_hit,
    parse_blast_xml_content,
    run_blast_search,
)


class TestBlast(unittest.TestCase):

    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.sample_blast_xml = """<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastp</BlastOutput_program>
  <BlastOutput_version>BLASTP 2.14.0+</BlastOutput_version>
  <BlastOutput_reference>Stephen F. Altschul et al.</BlastOutput_reference>
  <BlastOutput_db>nr</BlastOutput_db>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>test_query</BlastOutput_query-def>
  <BlastOutput_query-len>100</BlastOutput_query-len>
  <BlastOutput_param>
    <Parameters>
      <Parameters_matrix>BLOSUM62</Parameters_matrix>
      <Parameters_expect>10</Parameters_expect>
      <Parameters_gap-open>11</Parameters_gap-open>
      <Parameters_gap-extend>1</Parameters_gap-extend>
      <Parameters_filter>F</Parameters_filter>
    </Parameters>
  </BlastOutput_param>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_query-ID>Query_1</Iteration_query-ID>
      <Iteration_query-def>test_query</Iteration_query-def>
      <Iteration_query-len>100</Iteration_query-len>
      <Iteration_hits>
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>ref|YP_009724390.1|</Hit_id>
          <Hit_def>surface glycoprotein [Severe acute respiratory syndrome coronavirus 2]</Hit_def>
          <Hit_accession>YP_009724390</Hit_accession>
          <Hit_len>1273</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>205.5</Hsp_bit-score>
              <Hsp_score>520</Hsp_score>
              <Hsp_evalue>1.2e-65</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>100</Hsp_query-to>
              <Hsp_hit-from>1</Hsp_hit-from>
              <Hsp_hit-to>100</Hsp_hit-to>
              <Hsp_query-frame>0</Hsp_query-frame>
              <Hsp_hit-frame>0</Hsp_hit-frame>
              <Hsp_identity>100</Hsp_identity>
              <Hsp_positive>100</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>100</Hsp_align-len>
              <Hsp_qseq>MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSN</Hsp_qseq>
              <Hsp_hseq>MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSN</Hsp_hseq>
              <Hsp_midline>||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
      <Iteration_stat>
        <Statistics>
          <Statistics_db-num>1</Statistics_db-num>
          <Statistics_db-len>1273</Statistics_db-len>
          <Statistics_hsp-len>0</Statistics_hsp-len>
          <Statistics_eff-space>0</Statistics_eff-space>
          <Statistics_kappa>0.041</Statistics_kappa>
          <Statistics_lambda>0.267</Statistics_lambda>
          <Statistics_entropy>0.14</Statistics_entropy>
        </Statistics>
      </Iteration_stat>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>"""

    def tearDown(self):
        shutil.rmtree(self.test_dir, ignore_errors=True)

    def test_interpret_blast_hit(self):
        # High confidence & safe homology
        res1 = interpret_blast_hit(1e-50, 85.0)
        self.assertIn("High confidence", res1["confidence_level"])
        self.assertIn("Safe homology", res1["identity_zone"])

        # Moderate confidence & twilight zone
        res2 = interpret_blast_hit(1e-5, 25.0)
        self.assertIn("Moderate confidence", res2["confidence_level"])
        self.assertIn("Twilight zone", res2["identity_zone"])

        # Weak & midnight zone
        res3 = interpret_blast_hit(0.5, 12.0)
        self.assertIn("Weak", res3["confidence_level"])
        self.assertIn("Midnight zone", res3["identity_zone"])

    def test_parse_blast_xml_content(self):
        hits = parse_blast_xml_content(self.sample_blast_xml, query_len=100, max_hits=5)
        self.assertEqual(len(hits), 1)
        hit = hits[0]
        self.assertEqual(hit["rank"], 1)
        self.assertEqual(hit["accession"], "YP_009724390")
        self.assertEqual(hit["identity_percent"], 100.0)
        self.assertEqual(hit["query_coverage_percent"], 100.0)
        self.assertLess(hit["evalue"], 1e-50)

    def test_run_blast_search_with_cached_xml(self):
        analysis_dir = os.path.join(self.test_dir, "analysis")
        os.makedirs(analysis_dir, exist_ok=True)
        xml_file = os.path.join(analysis_dir, "blast_results.xml")
        with open(xml_file, "w", encoding="utf-8") as f:
            f.write(self.sample_blast_xml)

        result = run_blast_search(
            sequence="MFVFLVLLPL",
            output_dir=self.test_dir,
            max_hits=5,
            run_online=False
        )
        self.assertEqual(result["total_hits_found"], 1)
        self.assertTrue(os.path.exists(os.path.join(analysis_dir, "blast_results.json")))


if __name__ == "__main__":
    unittest.main()
