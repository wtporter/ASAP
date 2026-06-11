"""
Unit tests for identityFilter._identity_filter and the filter_pairs parameter.

Uses synthetic BAM files with MD/NM tags to drive _passes_identity (which
relies on get_aligned_pairs(with_seq=True), which works off the MD tag alone
when no reference fasta is supplied).
"""
import os
import sys
import pysam
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from asap import identityFilter

READ_LEN = 20
PERCID = 0.97

# 18/20 = 0.9 -- fails at PERCID=0.97
FAIL_MD, FAIL_NM = "5G5G8", 2
# 20/20 = 1.0 -- passes at PERCID=0.97
PASS_MD, PASS_NM = "20", 0

HEADER = pysam.AlignmentHeader.from_dict({
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [
        {"LN": 1000, "SN": "checked_ref"},
        {"LN": 1000, "SN": "other_ref"},
    ],
})


def _make_read(name, start, is_read1, md_tag, nm, ref_id=0):
    r = pysam.AlignedSegment(HEADER)
    r.query_name = name
    r.query_sequence = "A" * READ_LEN
    r.flag = 0x1 | 0x2 | (0x40 if is_read1 else 0x80)
    r.reference_id = ref_id
    r.reference_start = start
    r.mapping_quality = 60
    r.cigar = [(0, READ_LEN)]
    r.query_qualities = pysam.qualitystring_to_array("I" * READ_LEN)
    r.next_reference_id = ref_id
    r.next_reference_start = start
    r.set_tag("MD", md_tag)
    r.set_tag("NM", nm)
    return r


def _build_bam(tmp_path, reads, name):
    raw = str(tmp_path / f"{name}_raw.bam")
    with pysam.AlignmentFile(raw, "wb", header=HEADER) as bam:
        for r in reads:
            bam.write(r)
    sorted_bam = str(tmp_path / f"{name}_sorted.bam")
    pysam.sort("-o", sorted_bam, raw)
    pysam.index(sorted_bam)
    return sorted_bam


def _read_output(out_fp):
    with pysam.AlignmentFile(out_fp, "rb", check_sq=False) as out_bam:
        return {("R1" if r.is_read1 else "R2"): r for r in out_bam.fetch(until_eof=True)}


def test_filter_pairs_drops_both_mates_when_one_fails(tmp_path, monkeypatch):
    """Purpose: verify that with filter_pairs=True, if either mate of a read
    pair fails the percent-identity check, BOTH mates are marked unmapped in
    the output (pair-aware filtering).

    Function under test: identityFilter._identity_filter -- the two-pass
    filter_pairs=True logic: pass 1 collects query_names of reads on checked
    references that fail _passes_identity; pass 2 marks any read (or its
    mate) in that set as unmapped.

    Test input: one read pair "PAIR1" on "checked_ref": mate1 (read1) has
    FAIL_MD/FAIL_NM (18/20=0.9 identity, fails PERCID=0.97); mate2 (read2)
    has PASS_MD/PASS_NM (20/20=1.0, passes). _identity_filter called with
    ref_names=None, percid=0.97, filter_pairs=True.

    Expected result: discarded_reads == 2; both R1 and R2 in the output BAM
    have is_unmapped == True.
    """
    monkeypatch.chdir(tmp_path)
    reads = [
        _make_read("PAIR1", 0, True, FAIL_MD, FAIL_NM),
        _make_read("PAIR1", 100, False, PASS_MD, PASS_NM),
    ]
    bam_path = _build_bam(tmp_path, reads, "test1")
    out_fp = str(tmp_path / "out1.bam")
    with pysam.AlignmentFile(bam_path, "rb") as samdata:
        (_, discarded_reads) = identityFilter._identity_filter(samdata, None, PERCID, False, True, out_fp)

    assert discarded_reads == 2
    out_reads = _read_output(out_fp)
    assert out_reads["R1"].is_unmapped
    assert out_reads["R2"].is_unmapped


def test_filter_pairs_false_drops_only_failing_mate(tmp_path, monkeypatch):
    """Purpose: verify that with filter_pairs=False, only the individual
    mate that fails the identity check is marked unmapped, while its passing
    mate is left untouched (no pair-aware propagation).

    Function under test: identityFilter._identity_filter -- the single-pass
    filter_pairs=False logic, where each read is evaluated independently via
    _passes_identity and only failing reads are marked unmapped, regardless
    of mate status.

    Test input: same "PAIR1" setup as
    test_filter_pairs_drops_both_mates_when_one_fails (mate1 fails identity,
    mate2 passes); _identity_filter called with ref_names=None, percid=0.97,
    filter_pairs=False.

    Expected result: discarded_reads == 1; R1.is_unmapped == True and
    R2.is_unmapped == False.
    """
    monkeypatch.chdir(tmp_path)
    reads = [
        _make_read("PAIR1", 0, True, FAIL_MD, FAIL_NM),
        _make_read("PAIR1", 100, False, PASS_MD, PASS_NM),
    ]
    bam_path = _build_bam(tmp_path, reads, "test2")
    out_fp = str(tmp_path / "out2.bam")
    with pysam.AlignmentFile(bam_path, "rb") as samdata:
        (_, discarded_reads) = identityFilter._identity_filter(samdata, None, PERCID, False, False, out_fp)

    assert discarded_reads == 1
    out_reads = _read_output(out_fp)
    assert out_reads["R1"].is_unmapped
    assert not out_reads["R2"].is_unmapped


@pytest.mark.parametrize("filter_pairs", [True, False])
def test_both_mates_pass_are_kept(tmp_path, monkeypatch, filter_pairs):
    """Purpose: verify that when both mates of a pair pass the
    percent-identity check, neither is marked unmapped, regardless of the
    filter_pairs setting.

    Function under test: identityFilter._identity_filter /
    _passes_identity -- _passes_identity returns True for both mates
    (PASS_MD/PASS_NM = 20/20 = 1.0 identity >= PERCID), so neither the
    pair-aware (filter_pairs=True) nor the independent (filter_pairs=False)
    code path adds them to the failing set.

    Test input: "PAIR1" with both mates using PASS_MD/PASS_NM (1.0
    identity); _identity_filter called with ref_names=None, percid=0.97, and
    filter_pairs parametrized over [True, False].

    Expected result: discarded_reads == 0 for both parametrizations;
    R1.is_unmapped == False and R2.is_unmapped == False.
    """
    monkeypatch.chdir(tmp_path)
    reads = [
        _make_read("PAIR1", 0, True, PASS_MD, PASS_NM),
        _make_read("PAIR1", 100, False, PASS_MD, PASS_NM),
    ]
    bam_path = _build_bam(tmp_path, reads, f"test3_{filter_pairs}")
    out_fp = str(tmp_path / f"out3_{filter_pairs}.bam")
    with pysam.AlignmentFile(bam_path, "rb") as samdata:
        (_, discarded_reads) = identityFilter._identity_filter(samdata, None, PERCID, False, filter_pairs, out_fp)

    assert discarded_reads == 0
    out_reads = _read_output(out_fp)
    assert not out_reads["R1"].is_unmapped
    assert not out_reads["R2"].is_unmapped


def test_filter_pairs_scoped_to_checked_references(tmp_path, monkeypatch):
    """Purpose: verify that percent-identity filtering (and pair-aware
    dropping) is scoped only to the references listed in ref_names -- reads
    aligned to other (unchecked) references pass through untouched even if
    their mate fails on a checked reference.

    Function under test: identityFilter._identity_filter -- reference-scoping
    logic: when ref_names is provided, only reads whose reference_name is in
    ref_names are evaluated by _passes_identity / contribute to the failing
    query_name set; reads on other references are written through unmodified
    even under filter_pairs=True.

    Test input: "PAIR1" with mate1 on "checked_ref" (FAIL_MD/FAIL_NM, fails)
    and mate2 on "other_ref" (PASS_MD/PASS_NM, passes); _identity_filter
    called with ref_names=["checked_ref"], percid=0.97, filter_pairs=True.

    Expected result: discarded_reads == 1; R1.is_unmapped == True;
    R2.is_unmapped == False and R2.reference_name == "other_ref" (untouched
    despite filter_pairs=True, because "other_ref" isn't checked).
    """
    monkeypatch.chdir(tmp_path)
    reads = [
        _make_read("PAIR1", 0, True, FAIL_MD, FAIL_NM, ref_id=0),   # checked_ref, fails
        _make_read("PAIR1", 0, False, PASS_MD, PASS_NM, ref_id=1),  # other_ref, not checked
    ]
    bam_path = _build_bam(tmp_path, reads, "test4")
    out_fp = str(tmp_path / "out4.bam")
    with pysam.AlignmentFile(bam_path, "rb") as samdata:
        (_, discarded_reads) = identityFilter._identity_filter(samdata, ["checked_ref"], PERCID, False, True, out_fp)

    assert discarded_reads == 1
    out_reads = _read_output(out_fp)
    assert out_reads["R1"].is_unmapped
    assert not out_reads["R2"].is_unmapped
    assert out_reads["R2"].reference_name == "other_ref"
