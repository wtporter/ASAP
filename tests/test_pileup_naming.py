"""
Regression tests for de-novo SNP naming in _process_pileup.

Verifies that a SNP not present in the assay's configured snp_dict is named
using nucleotide notation ({reference}{position}{variant}) instead of the old
"unknown" placeholder, and that the position-0 "position of interest"
wildcard override (used by TB-style assays) is preserved.
"""
import os
import sys

import pysam
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from asap.newBamProcessor import _process_pileup
from asap.assayInfo import Amplicon, SNP, Significance

REF_LEN = 20
REF_NAME = "test_ref"
AMPLICON_SEQ = "A" * REF_LEN

HEADER = pysam.AlignmentHeader.from_dict({
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [{"LN": REF_LEN, "SN": REF_NAME}],
})


def _make_read(name, seq, start):
    """Return an AlignedSegment covering positions [start, start+len(seq))."""
    r = pysam.AlignedSegment(HEADER)
    r.query_name = name
    r.query_sequence = seq
    r.flag = 0
    r.reference_id = 0
    r.reference_start = start
    r.mapping_quality = 60
    r.cigar = [(0, len(seq))]  # all match
    r.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    return r


@pytest.fixture
def snp_pileup(tmp_path):
    """
    Write a small synthetic BAM (10 reads, all carrying T at amplicon position
    6 against an all-"A" reference) and yield a real pysam pileup over it.
    """
    reads = []
    for i in range(10):
        seq = list(AMPLICON_SEQ)
        seq[5] = "T"  # 0-based index 5 -> 1-based position 6
        reads.append(_make_read(f"r{i}", "".join(seq), 0))

    raw = str(tmp_path / "raw.bam")
    with pysam.AlignmentFile(raw, "wb", header=HEADER) as bam:
        for r in reads:
            bam.write(r)
    sorted_bam = raw + ".sorted.bam"
    pysam.sort("-o", sorted_bam, raw)
    pysam.index(sorted_bam)

    samdata = pysam.AlignmentFile(sorted_bam, "rb")
    pileup = samdata.pileup(REF_NAME, max_depth=10000000, ignore_orphans=False, ignore_overlaps=False)
    yield pileup
    samdata.close()


def _run_pileup(pileup, amplicon):
    return _process_pileup(
        pileup, amplicon,
        depth=1, proportion=0.05, mutdepth=1, offset=0, wholegenome=False,
        base_qual=20, con_prop=0.5, fill_gap_char="false", fill_del_char="false",
        n_read_array=[0] * REF_LEN,
    )


def test_de_novo_snp_named_with_nucleotide_notation(snp_pileup):
    """A SNP not in the assay's snp_dict gets {reference}{position}{variant}."""
    amplicon = Amplicon(AMPLICON_SEQ)

    pileup_dict = _run_pileup(snp_pileup, amplicon)

    snps = pileup_dict['SNPs']
    assert len(snps) == 1
    assert snps[0]['name'] == "A6T"
    assert snps[0]['reference'] == "A"
    assert snps[0]['variant'] == "T"
    assert snps[0]['position'] == "6"


def test_position_zero_wildcard_overrides_de_novo_name(snp_pileup):
    """A TB-style position=0 'catch-all' SNP config still names de-novo SNPs
    'position of interest' rather than nucleotide notation."""
    amplicon = Amplicon(AMPLICON_SEQ)
    amplicon.add_SNP(SNP(position=0, reference='A', variant=None,
                          significance=Significance("position of interest")))

    pileup_dict = _run_pileup(snp_pileup, amplicon)

    snps = pileup_dict['SNPs']
    assert len(snps) == 1
    assert snps[0]['name'] == "position of interest"
    assert snps[0]['significance'].message == "position of interest"
