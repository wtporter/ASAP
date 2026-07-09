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
from asap.allele_linkage import _build_fragment_allele_table, _tally_allele_linkage

REF_LEN = 1000
REF_NAME = "test_ref"

HEADER = pysam.AlignmentHeader.from_dict({
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [{"LN": REF_LEN, "SN": REF_NAME}],
})


def _make_read(name, seq, start, flag=0, is_read1=True, quals=None):
    """Return an AlignedSegment covering positions [start, start+len(seq)).

    `quals`, if given, is a list/str of per-base Phred scores overriding the
    default uniform Q40 ("I") -- used to reproduce maskPrimers.py's real
    output, which writes primer-masked 'N' bases with literal quality 0.
    """
    r = pysam.AlignedSegment(HEADER)
    r.query_name = name
    r.query_sequence = seq
    r.flag = flag
    r.reference_id = 0
    r.reference_start = start
    r.mapping_quality = 60
    r.cigar = [(0, len(seq))]  # all match
    if quals is not None:
        r.query_qualities = list(quals)
    else:
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
    """Purpose: verify that two SNP positions found on the same read are
    reported as perfectly linked.

    Function under test: _build_fragment_allele_table /
    _tally_allele_linkage -- per-fragment base-call extraction followed by a
    tally of (pos, base) combinations across fragments present at all
    requested positions.

    Test input: 20 synthetic 25-bp reads, each starting at ref position 95,
    with a fixed 'T' at offset 5 (genome pos 100) and 'C' at offset 15
    (genome pos 110); _build_fragment_allele_table/_tally_allele_linkage
    called for positions=[100, 110].

    Expected result: result has exactly 1 combo,
    frozenset({(100, 'T'), (110, 'C')}), with count == 20.
    """
    reads = []
    for i in range(20):
        # 25-base read starting at 95, covers pos 100 (offset 5) and 110 (offset 15)
        seq = "A" * 5 + "T" + "A" * 9 + "C" + "A" * 10   # T at 100, C at 110
        reads.append(_make_read(f"r{i}", seq, 95))

    samdata = tmp_bam(reads)
    pos_table, _, _, _ = _build_fragment_allele_table(samdata, [100, 110])
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
    """Purpose: verify that SNP positions located on different mates of the
    same read pair (fragment) are still linked together.

    Function under test: _build_fragment_allele_table -- fragment-level
    merge step, where read1 and read2 share a query_name so base calls from
    both mates are combined into one fragment entry before
    _tally_allele_linkage tallies combos.

    Test input: 15 read pairs; read1 (flag 0x43) covers ref 40-89 with 'T' at
    offset 10 (pos 50); read2 (flag 0x83) covers ref 190-239 with 'G' at
    offset 10 (pos 200); positions=[50, 200].

    Expected result: result has exactly 1 combo,
    frozenset({(50, 'T'), (200, 'G')}), with count == 15.
    """
    reads = []
    for i in range(15):
        # read1: covers 40–90 (pos 50 at offset 10) → A variant
        r1 = _make_read(f"pair{i}", "A" * 10 + "T" + "A" * 39, 40, flag=0x43)  # paired, read1
        # read2: covers 190–240 (pos 200 at offset 10) → G variant
        r2 = _make_read(f"pair{i}", "A" * 10 + "G" + "A" * 39, 190, flag=0x83)  # paired, read2
        reads.extend([r1, r2])

    samdata = tmp_bam(reads)
    pos_table, _, _, _ = _build_fragment_allele_table(samdata, [50, 200])
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
    """Purpose: verify that two distinct allele combinations present in the
    population are tallied as two separate linkage groups with correct
    counts.

    Function under test: _tally_allele_linkage -- the returned Counter
    accumulates one entry per distinct (pos, base) frozenset observed across
    fragments.

    Test input: 20 reads at ref 95-119 covering positions 100 & 110: 10
    "ref"-named reads carry A@100/C@110, 10 "alt"-named reads carry
    T@100/G@110; positions=[100, 110].

    Expected result: result has exactly 2 entries, each with count == 10:
    frozenset({(100, 'A'), (110, 'C')}) and
    frozenset({(100, 'T'), (110, 'G')}).
    """
    reads = []
    for i in range(10):
        seq = "A" * 5 + "A" + "A" * 9 + "C" + "A" * 5   # A@100, C@110 (reference-like)
        reads.append(_make_read(f"ref{i}", seq, 95))
    for i in range(10):
        seq = "A" * 5 + "T" + "A" * 9 + "G" + "A" * 5   # T@100, G@110 (alt-like)
        reads.append(_make_read(f"alt{i}", seq, 95))

    samdata = tmp_bam(reads)
    pos_table, _, _, _ = _build_fragment_allele_table(samdata, [100, 110])
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
    """Purpose: verify that fragments which don't cover every requested
    position contribute nothing to the linkage tally.

    Function under test: _tally_allele_linkage -- only considers fragments
    present in the position table for ALL requested positions (driven by the
    smallest position table); fragments missing from any requested position
    are skipped entirely.

    Test input: 10 reads, 20-bp starting at ref 95 (covers genome positions
    95-114, reaching pos 100 but not pos 300); positions=[100, 300].

    Expected result: result is empty (len(result) == 0), since no fragment
    spans both pos 100 and pos 300.
    """
    reads = []
    for i in range(10):
        # 20-base read: covers 95–115, reaches pos 100 but NOT pos 300
        seq = "A" * 5 + "T" + "A" * 14
        reads.append(_make_read(f"short{i}", seq, 95))

    samdata = tmp_bam(reads)
    pos_table, _, _, _ = _build_fragment_allele_table(samdata, [100, 300])
    result = _tally_allele_linkage(pos_table, [100, 300])
    samdata.close()

    assert len(result) == 0


# ---------------------------------------------------------------------------
# Case 5: supplementary / secondary reads are skipped
# ---------------------------------------------------------------------------
def test_supplementary_secondary_excluded(tmp_bam):
    """Purpose: verify that supplementary and secondary alignments are
    ignored, so they don't pollute the per-position allele table or the
    linkage tally.

    Function under test: _build_fragment_allele_table -- read-filtering step
    that skips reads with flag bits 0x800 (supplementary) or 0x100
    (secondary) before recording base calls.

    Test input: 5 reads flagged supplementary (0x800) and 5 flagged
    secondary (0x100), all otherwise carrying T@100/C@110; positions=[100, 110].

    Expected result: result is empty (len(result) == 0), since none of the
    flagged reads contribute to pos_table and so no fragment spans both
    positions.
    """
    reads = []
    seq = "A" * 5 + "T" + "A" * 9 + "C" + "A" * 10
    for i in range(5):
        reads.append(_make_read(f"sup{i}", seq, 95, flag=0x800))   # supplementary
        reads.append(_make_read(f"sec{i}", seq, 95, flag=0x100))   # secondary

    samdata = tmp_bam(reads)
    pos_table, _, _, _ = _build_fragment_allele_table(samdata, [100, 110])
    result = _tally_allele_linkage(pos_table, [100, 110])
    samdata.close()

    assert len(result) == 0


# ---------------------------------------------------------------------------
# Case 6: caller applies a min_reads threshold to the raw tally
# ---------------------------------------------------------------------------
def test_min_reads_threshold(tmp_bam):
    """Purpose: verify that the raw tally returned by _tally_allele_linkage
    can be downstream-filtered by a min_reads threshold to drop rare allele
    combinations.

    Function under test: _tally_allele_linkage -- performs no filtering
    itself; it returns the raw Counter of (pos, base) combo -> fragment
    count, leaving any min_reads filtering to the caller.

    Test input: 23 reads at positions=[100, 110]: 3 "rare" reads carrying
    T@100/C@110, 20 "common" reads carrying A@100/A@110 (reference-like).

    Expected result: the raw result contains both combos (counts 3 and 20);
    after filtering entries with count >= 10, only
    frozenset({(100, 'A'), (110, 'A')}) -> 20 survives.
    """
    reads = []
    seq_rare  = "A" * 5 + "T" + "A" * 9 + "C" + "A" * 10  # rare combo
    seq_common = "A" * 5 + "A" + "A" * 9 + "A" + "A" * 10  # common (ref)
    for i in range(3):
        reads.append(_make_read(f"rare{i}", seq_rare, 95))
    for i in range(20):
        reads.append(_make_read(f"common{i}", seq_common, 95))

    samdata = tmp_bam(reads)
    pos_table, _, _, _ = _build_fragment_allele_table(samdata, [100, 110])
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
    """Purpose: verify that when overlapping mates of the same fragment
    disagree on masking at a shared position, the real call wins and the
    fragment is not double-counted as masked.

    Function under test: _build_fragment_allele_table -- masking
    reconciliation logic: if one mate reports 'N' (e.g. primer-masked) and
    the other mate reports a real base for the same fragment/position, the
    real call is recorded in pos_table and the fragment is excluded from
    masked for that position.

    Test input: a single fragment "pair0" with two overlapping mates: mate1
    (read1, flag 0x43) covers ref 95-114 with 'N' at offset 5 (pos 100);
    mate2 (read2, flag 0x83) covers ref 90-109 with 'T' at offset 10
    (pos 100); _build_fragment_allele_table called for positions=[100].

    Expected result: pos_table[100]["pair0"] == "T" (real call wins), and
    "pair0" not in masked[100].
    """
    # mate1: 95-114, position 100 is offset 5 -> 'N' (masked)
    seq1 = "A" * 5 + "N" + "A" * 14
    r1 = _make_read("pair0", seq1, 95, flag=0x43)  # paired, read1
    # mate2: 90-109, position 100 is offset 10 -> 'T' (real call)
    seq2 = "A" * 10 + "T" + "A" * 9
    r2 = _make_read("pair0", seq2, 90, flag=0x83)  # paired, read2

    samdata = tmp_bam([r1, r2])
    pos_table, _, masked, _ = _build_fragment_allele_table(samdata, [100])
    samdata.close()

    assert pos_table[100].get("pair0") == "T"
    assert "pair0" not in masked[100]


# ---------------------------------------------------------------------------
# Case 8: a Q0 masked base must still land in `masked`, not vanish
# ---------------------------------------------------------------------------
def test_masked_base_zero_quality_still_recorded(tmp_bam):
    """Purpose: verify that a primer-masked 'N' base written with literal
    quality 0 (maskPrimers.py's real output -- see maskPrimers.py:137,161,
    which zeroes quality for the entire masked region) is still recorded in
    `masked`, not silently dropped.

    This guards against a real regression class: pysam's `pileup()` defaults
    to `min_base_quality=13`, under which a Q0 base is invisible to
    `pileupcolumn.pileups` entirely (not flagged, just absent). The
    fetch()-based `_build_fragment_allele_table` has no quality filtering at
    all, so it must see this fragment regardless of quality.

    Test input: a single read covering ref 95-114, with an 'N' at offset 5
    (pos 100) whose quality is 0; all other bases quality 40 ("I").

    Expected result: pos_table[100] is empty for this fragment, and
    masked[100] contains it.
    """
    seq = "A" * 5 + "N" + "A" * 14
    quals = [40] * 5 + [0] + [40] * 14
    r = _make_read("masked0", seq, 95, quals=quals)

    samdata = tmp_bam([r])
    pos_table, _, masked, _ = _build_fragment_allele_table(samdata, [100])
    samdata.close()

    assert "masked0" not in pos_table[100]
    assert "masked0" in masked[100]
