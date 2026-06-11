"""
Unit tests for generateSMORbam._write_bam and generateSMORbam_correction._write_bam.

Creates synthetic BAM files using pysam to verify that grouping name-sorted
reads with itertools.groupby correctly accounts for every input read --
matched pairs, singletons (mate not aligned to this reference), and
non-overlapping pairs -- regardless of where a singleton falls in the
name-sorted order.
"""
import os
import sys
import csv
import pysam

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from asap import generateSMORbam
from asap import generateSMORbam_correction

REF_LEN = 200
REF_NAME = "test_ref"
READ_LEN = 50

# Deterministic 200bp "reference" sequence; reads are slices of this so that
# overlapping reads always agree at every shared position.
REF_SEQ = "ACGT" * (REF_LEN // 4)

HEADER = pysam.AlignmentHeader.from_dict({
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [{"LN": REF_LEN, "SN": REF_NAME}],
})


def _make_read(name, start, is_read1=True):
    r = pysam.AlignedSegment(HEADER)
    r.query_name = name
    seq = REF_SEQ[start:start + READ_LEN]
    r.query_sequence = seq
    r.flag = 0x1 | 0x2 | (0x40 if is_read1 else 0x80)
    r.reference_id = 0
    r.reference_start = start
    r.mapping_quality = 60
    r.cigar = [(0, len(seq))]
    r.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    return r


def _build_bam(tmp_path):
    """
    Builds a BAM with, in name-sorted order:
      AAA: a normal overlapping pair      -> 1 consensus read
      BBB: a singleton (mate not present) -> 1 singleton read
      CCC: a normal overlapping pair      -> 1 consensus read
      DDD: a non-overlapping pair         -> 2 dropped reads

    BBB sits between two real pairs so a naive "consecutive pairs of 2"
    grouping (the old bug) desyncs and silently loses CCC and DDD.
    """
    reads = [
        _make_read("AAA", 0, True),
        _make_read("AAA", 30, False),
        _make_read("BBB", 10, True),
        _make_read("CCC", 50, True),
        _make_read("CCC", 70, False),
        _make_read("DDD", 0, True),
        _make_read("DDD", 150, False),
    ]
    raw = str(tmp_path / "raw.bam")
    with pysam.AlignmentFile(raw, "wb", header=HEADER) as bam:
        for r in reads:
            bam.write(r)
    sorted_bam = str(tmp_path / "sorted.bam")
    pysam.sort("-o", sorted_bam, raw)
    pysam.index(sorted_bam)
    return sorted_bam, len(reads)


def _read_stats(stats_path):
    with open(stats_path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 1
    return {k: int(v) for k, v in rows[0].items() if k != "ref_name"}


def _assert_reconciles(stats, total_reads):
    assert stats["input_reads"] == total_reads
    assert stats["consensus_reads"] == 2   # AAA, CCC
    assert stats["singleton_reads"] == 1   # BBB
    assert stats["pairs_dropped"] == 2     # DDD
    assert (2 * stats["consensus_reads"] + stats["pairs_dropped"]
            + stats["singleton_reads"]) == stats["input_reads"]


def test_generateSMORbam_correction_accounts_for_every_read(tmp_path, monkeypatch):
    """Purpose: verify that itertools.groupby-based pair grouping accounts
    for every input read, including a singleton positioned between two valid
    pairs in name-sorted order (regression test for a grouping desync bug).

    Function under test: generateSMORbam_correction._write_bam -- groups
    name-sorted reads by query_name; classifies each group as a singleton
    (counted, skipped), a non-overlapping pair (both reads -> pairs_dropped),
    or an overlapping pair (merged via _get_consensus into one consensus
    read); writes per-reference stats to smor_stats.tsv.

    Test input: a 7-read BAM, name-sorted as AAA (overlapping pair), BBB
    (singleton), CCC (overlapping pair), DDD (non-overlapping pair);
    generateSMORbam_correction._write_bam(samdata, out_file, "N", 10).

    Expected result: smor_stats.tsv shows input_reads==7, consensus_reads==2
    (AAA, CCC), singleton_reads==1 (BBB), pairs_dropped==2 (DDD), and
    2*consensus_reads + pairs_dropped + singleton_reads == input_reads (all 7
    reads accounted for).
    """
    bam_path, total_reads = _build_bam(tmp_path)
    monkeypatch.chdir(tmp_path)

    out_file = str(tmp_path / "out_SMOR.bam")
    with pysam.AlignmentFile(bam_path, "rb") as samdata:
        generateSMORbam_correction._write_bam(samdata, out_file, "N", 10)

    stats = _read_stats(tmp_path / "smor_stats.tsv")
    _assert_reconciles(stats, total_reads)


def test_generateSMORbam_accounts_for_every_read(tmp_path, monkeypatch):
    """Purpose: verify the same read-accounting guarantee as
    test_generateSMORbam_correction_accounts_for_every_read, but for the
    non-quality-corrected generateSMORbam._write_bam implementation, ensuring
    its grouping logic doesn't desync around a mid-sequence singleton either.

    Function under test: generateSMORbam._write_bam -- same
    itertools.groupby grouping/classification logic as
    generateSMORbam_correction._write_bam, but consensus merging via
    generateSMORbam._get_consensus (overlap-region merge without
    quality-based correction).

    Test input: the same 7-read BAM as
    test_generateSMORbam_correction_accounts_for_every_read (AAA pair, BBB
    singleton, CCC pair, DDD non-overlapping pair);
    generateSMORbam._write_bam(samdata, out_file, "N", 0, False).

    Expected result: same reconciliation as the correction test:
    input_reads==7, consensus_reads==2, singleton_reads==1,
    pairs_dropped==2, and 2*consensus_reads + pairs_dropped +
    singleton_reads == input_reads.
    """
    bam_path, total_reads = _build_bam(tmp_path)
    monkeypatch.chdir(tmp_path)

    out_file = str(tmp_path / "out_SMOR.bam")
    with pysam.AlignmentFile(bam_path, "rb") as samdata:
        generateSMORbam._write_bam(samdata, out_file, "N", 0, False)

    stats = _read_stats(tmp_path / "smor_stats.tsv")
    _assert_reconciles(stats, total_reads)
