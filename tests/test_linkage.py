"""
Unit tests for _tally_allele_linkage.

Creates synthetic BAM files using pysam so each case exercises the
read-pair merging logic without requiring real sequencing data.
"""
import os
import sys
import pytest
import pysam

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from asap.newBamProcessor import _build_fragment_allele_table, _tally_allele_linkage

REF_LEN = 1000
REF_NAME = "test_ref"

HEADER = pysam.AlignmentHeader.from_dict({
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [{"LN": REF_LEN, "SN": REF_NAME}],
})


def _make_read(name, seq, start, flag=0, is_read1=True):
    """Return an AlignedSegment covering positions [start, start+len(seq))."""
    r = pysam.AlignedSegment(HEADER)
    r.query_name = name
    r.query_sequence = seq
    r.flag = flag
    r.reference_id = 0
    r.reference_start = start
    r.mapping_quality = 60
    r.cigar = [(0, len(seq))]  # all match
    r.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    return r


def _write_bam(reads, path):
    with pysam.AlignmentFile(path, "wb", header=HEADER) as bam:
        for r in reads:
            bam.write(r)
    pysam.sort("-o", path + ".sorted.bam", path)
    pysam.index(path + ".sorted.bam")
    return path + ".sorted.bam"


@pytest.fixture
def tmp_bam(tmp_path):
    """Fixture helper: write a list of reads, sort+index, return open samdata."""
    def _factory(reads):
        raw = str(tmp_path / "raw.bam")
        sorted_bam = _write_bam(reads, raw)
        return pysam.AlignmentFile(sorted_bam, "rb")
    return _factory


# ---------------------------------------------------------------------------
# Case 1: two SNP positions on the same read — perfectly linked
# ---------------------------------------------------------------------------
def test_linked_on_same_read(tmp_bam):
    """Both variant positions (100, 110) are on every read → always linked."""
    reads = []
    for i in range(20):
        # 25-base read starting at 95, covers pos 100 (offset 5) and 110 (offset 15)
        seq = "A" * 5 + "T" + "A" * 9 + "C" + "A" * 10   # T at 100, C at 110
        reads.append(_make_read(f"r{i}", seq, 95))

    samdata = tmp_bam(reads)
    pos_table, _, _ = _build_fragment_allele_table(samdata, [100, 110])
    result = _tally_allele_linkage(pos_table, [100, 110])
    samdata.close()

    assert len(result) == 1
    combo = next(iter(result))
    assert frozenset([(100, "T"), (110, "C")]) == combo
    assert result[combo] == 20


# ---------------------------------------------------------------------------
# Case 2: two positions on different mates of a read pair — should be linked
# ---------------------------------------------------------------------------
def test_linked_across_mates(tmp_bam):
    """Position 50 is on read1, position 200 is on read2; same fragment → linked."""
    reads = []
    for i in range(15):
        # read1: covers 40–90 (pos 50 at offset 10) → A variant
        r1 = _make_read(f"pair{i}", "A" * 10 + "T" + "A" * 39, 40, flag=0x43)  # paired, read1
        # read2: covers 190–240 (pos 200 at offset 10) → G variant
        r2 = _make_read(f"pair{i}", "A" * 10 + "G" + "A" * 39, 190, flag=0x83)  # paired, read2
        reads.extend([r1, r2])

    samdata = tmp_bam(reads)
    pos_table, _, _ = _build_fragment_allele_table(samdata, [50, 200])
    result = _tally_allele_linkage(pos_table, [50, 200])
    samdata.close()

    assert len(result) == 1
    combo = next(iter(result))
    assert frozenset([(50, "T"), (200, "G")]) == combo
    assert result[combo] == 15


# ---------------------------------------------------------------------------
# Case 3: unlinked SNPs — two independent allele combinations expected
# ---------------------------------------------------------------------------
def test_unlinked_snps(tmp_bam):
    """Half the reads have A@100+C@110, other half have T@100+G@110; two combos."""
    reads = []
    for i in range(10):
        seq = "A" * 5 + "A" + "A" * 9 + "C" + "A" * 5   # A@100, C@110 (reference-like)
        reads.append(_make_read(f"ref{i}", seq, 95))
    for i in range(10):
        seq = "A" * 5 + "T" + "A" * 9 + "G" + "A" * 5   # T@100, G@110 (alt-like)
        reads.append(_make_read(f"alt{i}", seq, 95))

    samdata = tmp_bam(reads)
    pos_table, _, _ = _build_fragment_allele_table(samdata, [100, 110])
    result = _tally_allele_linkage(pos_table, [100, 110])
    samdata.close()

    assert len(result) == 2
    counts = {frozenset(k): v for k, v in result.items()}
    assert counts[frozenset([(100, "A"), (110, "C")])] == 10
    assert counts[frozenset([(100, "T"), (110, "G")])] == 10


# ---------------------------------------------------------------------------
# Case 4: reads that do not span all positions are excluded
# ---------------------------------------------------------------------------
def test_reads_not_spanning_all_excluded(tmp_bam):
    """Reads covering only pos 100 (not 300) must not contribute to the tally."""
    reads = []
    for i in range(10):
        # 20-base read: covers 95–115, reaches pos 100 but NOT pos 300
        seq = "A" * 5 + "T" + "A" * 14
        reads.append(_make_read(f"short{i}", seq, 95))

    samdata = tmp_bam(reads)
    pos_table, _, _ = _build_fragment_allele_table(samdata, [100, 300])
    result = _tally_allele_linkage(pos_table, [100, 300])
    samdata.close()

    assert len(result) == 0


# ---------------------------------------------------------------------------
# Case 5: supplementary / secondary reads are skipped
# ---------------------------------------------------------------------------
def test_supplementary_secondary_excluded(tmp_bam):
    """Supplementary (flag 0x800) and secondary (flag 0x100) reads are ignored."""
    reads = []
    seq = "A" * 5 + "T" + "A" * 9 + "C" + "A" * 10
    for i in range(5):
        reads.append(_make_read(f"sup{i}", seq, 95, flag=0x800))   # supplementary
        reads.append(_make_read(f"sec{i}", seq, 95, flag=0x100))   # secondary

    samdata = tmp_bam(reads)
    pos_table, _, _ = _build_fragment_allele_table(samdata, [100, 110])
    result = _tally_allele_linkage(pos_table, [100, 110])
    samdata.close()

    assert len(result) == 0


# ---------------------------------------------------------------------------
# Case 6: caller applies a min_reads threshold to the raw tally
# ---------------------------------------------------------------------------
def test_min_reads_threshold(tmp_bam):
    """A combo seen only 3 times can be filtered out by a min_reads=10 threshold."""
    reads = []
    seq_rare  = "A" * 5 + "T" + "A" * 9 + "C" + "A" * 10  # rare combo
    seq_common = "A" * 5 + "A" + "A" * 9 + "A" + "A" * 10  # common (ref)
    for i in range(3):
        reads.append(_make_read(f"rare{i}", seq_rare, 95))
    for i in range(20):
        reads.append(_make_read(f"common{i}", seq_common, 95))

    samdata = tmp_bam(reads)
    pos_table, _, _ = _build_fragment_allele_table(samdata, [100, 110])
    result = _tally_allele_linkage(pos_table, [100, 110])
    samdata.close()

    # The rare combo (3 reads) is below min_reads=10, only common survives
    surviving = {k: v for k, v in result.items() if v >= 10}
    assert len(surviving) == 1
    assert result[frozenset([(100, "A"), (110, "A")])] == 20


# ---------------------------------------------------------------------------
# Case 7: overlapping mates disagree on masking -- a real call wins
# ---------------------------------------------------------------------------
def test_masked_pos_table_disjoint_on_overlap(tmp_bam):
    """
    Mate1 covers position 100 with an 'N' (e.g. primer-masked); mate2 of the
    same fragment also covers position 100, but with a real base call. The
    fragment must end up in pos_table (real call wins) and NOT also in
    masked, so callers can't double-count it.
    """
    # mate1: 95-114, position 100 is offset 5 -> 'N' (masked)
    seq1 = "A" * 5 + "N" + "A" * 14
    r1 = _make_read("pair0", seq1, 95, flag=0x43)  # paired, read1
    # mate2: 90-109, position 100 is offset 10 -> 'T' (real call)
    seq2 = "A" * 10 + "T" + "A" * 9
    r2 = _make_read("pair0", seq2, 90, flag=0x83)  # paired, read2

    samdata = tmp_bam([r1, r2])
    pos_table, _, masked = _build_fragment_allele_table(samdata, [100])
    samdata.close()

    assert pos_table[100].get("pair0") == "T"
    assert "pair0" not in masked[100]
