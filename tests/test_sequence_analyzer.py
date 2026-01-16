import pytest
import sys
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src import SequenceAnalyzer

@pytest.fixture
def sequence_analyzer():
  test_dir = os.path.dirname(os.path.abspath(__file__))
  data_path = os.path.join(test_dir, 'sample.fasta')
  return SequenceAnalyzer(data_path)

class TestValidateSequence:
  def test_validate_sequence(self, sequence_analyzer):
    assert all(sequence_analyzer.validate_sequence(r.seq) for r in sequence_analyzer.records)

class TestValidateBadSequence:
  def test_validate_sequence(self, sequence_analyzer):
    validated_seqs = []
    sequence_analyzer.records.append(SeqRecord(Seq("AGTTAX"), id="ABCD1234", name="Test record", description=("A Fake test strand to use for validating sequences")))
    for record in sequence_analyzer.records:
      validated_seqs.append(sequence_analyzer.validate_sequence(record.seq))
    assert not validated_seqs[-1]
    # "https://www.youtube.com/watch?v=WxMFCfFRY2w"

class TestGetCodons:
  def test_get_codons(self, sequence_analyzer):
    test_seq_id = 'MP007459.1'
    codons = sequence_analyzer.get_codons(test_seq_id)
    assert len(codons) == 4
    for codon in codons:
      assert len(codon) == 3

class TestIncompleteGetCodons:
  def test_get_codons(self, sequence_analyzer):
    sequence_analyzer.records.append(SeqRecord(Seq("AGTT"), id="ABCD1234", name="Test record", description=("A Fake test strand to use for validating sequences")))
    codons = sequence_analyzer.get_codons('ABCD1234')
    assert len(codons) == 1

class TestPointMutation:
  def test_point_mutation(self, sequence_analyzer):
    test_position = 0
    test_seq_id = 'MP007459.1'
    test_new_base = "T"
    new_mutation = sequence_analyzer.point_mutation(test_seq_id, test_position, test_new_base)
    assert new_mutation == 'TTGCGTAGCTAA'
    test_position = 11
    test_new_base = 'G'
    new_mutation = sequence_analyzer.point_mutation(test_seq_id, test_position, test_new_base)
    assert new_mutation == 'TTGCGTAGCTAG'

class TestRandomMutation:
  def test_random_mutation(self, sequence_analyzer):
    test_seq_id = 'MP007459.1'
    new_mutation = sequence_analyzer.simulate_random_mutation(test_seq_id)
    assert new_mutation != 'ATGCGTAGCTAA'

