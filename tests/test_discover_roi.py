"""
Unit tests for _apply_discover_roi, focused on the discover_roi_min_snp_perc
pre-filter.

Creates synthetic BAM files using pysam so each case exercises the
read-level linkage logic without requiring real sequencing data.
"""
import os
import sys
from collections import Counter
from xml.etree import ElementTree

import pytest
import pysam

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from asap.newBamProcessor import (
    _apply_discover_roi,
    _build_fragment_allele_table,
    _add_linked_snps_node,
)

REF_LEN = 1000
REF_NAME = "test_ref"

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


def _linked_names(snp):
    return {entry['target_name'] for entry in snp.get('linked_snps', [])}


# ---------------------------------------------------------------------------
# Case 1/2: a low-frequency "noise" SNP (4%) is excluded by the default
# min_snp_perc=0.05, both as an anchor and as a linked candidate, but is
# considered when min_snp_perc=0.0.
# ---------------------------------------------------------------------------
def _make_anchor_linked_noise_bam(tmp_bam):
    """
    25 reads, 41bp starting at ref 95 -> covers amp positions 100, 110, 120.

    - 21 reads: T@100, C@110, A@120
    - 1 read:   T@100, C@110, G@120  (the lone "noise" variant carrier)
    - 3 reads:  A@100, A@110, A@120  (reference)
    """
    reads = []
    for i in range(21):
        seq = list("A" * 41)
        seq[5] = "T"   # amp 100
        seq[15] = "C"  # amp 110
        reads.append(_make_read(f"tc{i}", "".join(seq), 95))

    seq = list("A" * 41)
    seq[5] = "T"
    seq[15] = "C"
    seq[25] = "G"  # amp 120 -- the single "noise" carrier
    reads.append(_make_read("noise_carrier", "".join(seq), 95))

    for i in range(3):
        reads.append(_make_read(f"ref{i}", "A" * 41, 95))

    samdata = tmp_bam(reads)
    pos_table, reach, masked, _ = _build_fragment_allele_table(samdata, [100, 110, 120], REF_NAME)
    samdata.close()
    return pos_table, reach, masked


def _make_snp_list():
    anchor = {'name': 'anchor', 'position': '101', 'reference': 'A', 'variant': 'T',
              'depth': '25', 'basecalls': Counter({'T': 22, 'A': 3})}
    linked = {'name': 'linked', 'position': '111', 'reference': 'A', 'variant': 'C',
              'depth': '25', 'basecalls': Counter({'C': 22, 'A': 3})}
    noise = {'name': 'noise', 'position': '121', 'reference': 'A', 'variant': 'G',
             'depth': '25', 'basecalls': Counter({'G': 1, 'A': 24})}
    return [anchor, linked, noise]


def test_low_freq_snp_excluded_by_default_min_snp_perc(tmp_bam):
    """Purpose: verify that with the default min_snp_perc=0.05 threshold, a
    low-frequency (4%) "noise" SNP is excluded entirely from linkage
    discovery -- both as a candidate that gets its own linked_snps list, and
    as a linked candidate for other SNPs -- while a high-frequency (88%)
    anchor/linked pair still link to each other.

    Function under test: _apply_discover_roi -- the discover_roi_min_snp_perc
    pre-filter, which excludes SNPs with basecalls[variant]/depth <
    min_snp_perc from both sides of the pairwise linkage comparison.

    Test input: a 25-fragment synthetic BAM at amplicon positions
    100/110/120 (21 reads T@100/C@110/A@120, 1 "noise" read
    T@100/C@110/G@120, 3 reference reads); snp_list = [anchor (22/25=88%),
    linked (22/25=88%), noise (1/25=4%)]; _apply_discover_roi called with
    min_perc=0.0, min_reads=1, min_snp_perc=0.05 (default).

    Expected result: 'linked_snps' not in noise; 'noise' does not appear in
    anchor's or linked's linked_snps; 'linked' appears in anchor's
    linked_snps and 'anchor' appears in linked's linked_snps.
    """
    pos_table, reach, masked = _make_anchor_linked_noise_bam(tmp_bam)
    snp_list = _make_snp_list()
    anchor, linked, noise = snp_list

    _apply_discover_roi(snp_list, pos_table, reach, masked, offset=0,
                         min_perc=0.0, min_reads=1, min_snp_perc=0.05)

    assert 'linked_snps' not in noise
    assert 'noise' not in _linked_names(anchor)
    assert 'noise' not in _linked_names(linked)
    assert 'linked' in _linked_names(anchor)
    assert 'anchor' in _linked_names(linked)


def test_low_freq_snp_included_when_min_snp_perc_zero(tmp_bam):
    """Purpose: verify that lowering min_snp_perc to 0.0 allows the same 4%
    "noise" SNP to fully participate in linkage discovery -- it gets its own
    linked_snps entry and is reported as a candidate linked to 'anchor'.

    Function under test: _apply_discover_roi -- same
    discover_roi_min_snp_perc pre-filter as
    test_low_freq_snp_excluded_by_default_min_snp_perc, but with
    min_snp_perc=0.0 so noise's frequency (1/25=0.04 >= 0.0) passes the
    pre-filter and is eligible for pairwise linkage comparison.

    Test input: same BAM and snp_list (anchor, linked, noise) as
    test_low_freq_snp_excluded_by_default_min_snp_perc, but
    _apply_discover_roi called with min_perc=0.0, min_reads=1,
    min_snp_perc=0.0.

    Expected result: 'linked_snps' is present in noise and non-empty;
    'noise' appears in anchor's linked_snps.
    """
    pos_table, reach, masked = _make_anchor_linked_noise_bam(tmp_bam)
    snp_list = _make_snp_list()
    anchor, linked, noise = snp_list

    _apply_discover_roi(snp_list, pos_table, reach, masked, offset=0,
                         min_perc=0.0, min_reads=1, min_snp_perc=0.0)

    assert 'linked_snps' in noise
    assert _linked_names(noise)
    assert 'noise' in _linked_names(anchor)


# ---------------------------------------------------------------------------
# Case 3: a SNP at exactly min_snp_perc is included (inclusive minimum, >=)
# ---------------------------------------------------------------------------
def test_min_snp_perc_boundary_is_inclusive(tmp_bam):
    """Purpose: verify the discover_roi_min_snp_perc pre-filter uses an
    inclusive (>=) comparison, so a SNP whose frequency exactly equals the
    threshold (5%) is still included in linkage discovery.

    Function under test: _apply_discover_roi -- same
    discover_roi_min_snp_perc pre-filter, exercising the boundary case where
    the frequency comparison is `>= min_snp_perc`, not strictly `>`.

    Test input: 20 fragments at amplicon positions 100/130 (1 read carries
    T@100 + C@130 ("both"), 17 reads carry only T@100, 2 reference reads);
    snp_list = [anchor2 (18/20=90%), boundary (1/20=5% == min_snp_perc)];
    _apply_discover_roi called with min_perc=0.0, min_reads=1,
    min_snp_perc=0.05.

    Expected result: 'linked_snps' is present in boundary, and 'anchor2'
    appears in boundary's linked_snps -- the 5% SNP was not excluded by the
    pre-filter.
    """
    # 20 reads, 41bp starting at ref 95 -> covers amp positions 100, 130.
    reads = []
    # 1 read carrying both anchor2's variant (T@100) and boundary's variant (C@130)
    seq = list("A" * 41)
    seq[5] = "T"    # amp 100
    seq[35] = "C"   # amp 130
    reads.append(_make_read("both", "".join(seq), 95))

    # 17 reads carrying only anchor2's variant
    for i in range(17):
        seq = list("A" * 41)
        seq[5] = "T"
        reads.append(_make_read(f"anchor_only{i}", "".join(seq), 95))

    # 2 reference reads
    for i in range(2):
        reads.append(_make_read(f"ref{i}", "A" * 41, 95))

    samdata = tmp_bam(reads)
    pos_table, reach, masked, _ = _build_fragment_allele_table(samdata, [100, 130], REF_NAME)
    samdata.close()

    anchor2 = {'name': 'anchor2', 'position': '101', 'reference': 'A', 'variant': 'T',
               'depth': '20', 'basecalls': Counter({'T': 18, 'A': 2})}
    boundary = {'name': 'boundary', 'position': '131', 'reference': 'A', 'variant': 'C',
                'depth': '20', 'basecalls': Counter({'C': 1, 'A': 19})}
    snp_list = [anchor2, boundary]

    _apply_discover_roi(snp_list, pos_table, reach, masked, offset=0,
                         min_perc=0.0, min_reads=1, min_snp_perc=0.05)

    # boundary's frequency is exactly 1/20 == 0.05 == min_snp_perc -> included
    assert 'linked_snps' in boundary
    assert 'anchor2' in _linked_names(boundary)


# ---------------------------------------------------------------------------
# Case 4: linked_snps entry fields use the new semantic names
# ---------------------------------------------------------------------------
def test_linked_snps_entry_fields(tmp_bam):
    """Purpose: verify that for a partially-linked SNP pair the linked_snps
    entry carries the renamed fields with correct values.

    Function under test: _apply_discover_roi -- the linked_snps entry
    construction, specifically:
        linkage_pct          = co_count / (co_count + anchor_only) * 100
        sample_frequency_pct = co_count / shared_read_depth * 100
        target_only          = reads with s_b's variant but NOT s_a's
        neither_variant      = shared_read_depth - co - anchor_only - target_only

    Test input: amp positions 100 (s_a, A->T) and 110 (s_b, A->C); 41bp
    reads starting at ref 95 unless noted:
        - 6 reads: T@100, C@110              (both_variants)
        - 3 reads: T@100, A@110              (anchor_only: confident, unlinked)
        - 1 read (10bp): T@100, too short to reach 110 (read_too_short)
        - 5 reads: A@100, A@110              (neither_variant: reference,
          spans both but no variant at either position)
    snp_list = [s_a (depth=15, T count=10), s_b (depth=14, C count=6)];
    _apply_discover_roi called with min_perc=0.0, min_reads=1,
    min_snp_perc=0.0.

    Expected result: co_count=6, anchor_only=3, target_only=0,
    neither_variant=5, read_too_short=1, shared_read_depth=14,
    linkage_pct=66.7, sample_frequency_pct=42.9.
    """
    reads = []
    # 6 reads: T@100, C@110 (both variants)
    for i in range(6):
        seq = list("A" * 41)
        seq[5] = "T"
        seq[15] = "C"
        reads.append(_make_read(f"linked{i}", "".join(seq), 95))

    # 3 reads: T@100, A@110 (anchor only)
    for i in range(3):
        seq = list("A" * 41)
        seq[5] = "T"
        reads.append(_make_read(f"standalone{i}", "".join(seq), 95))

    # 1 read: T@100, too short to reach 110
    seq = list("A" * 10)
    seq[5] = "T"
    reads.append(_make_read("short0", "".join(seq), 95))

    # 5 reads: A@100, A@110 (reference -- spans both, neither variant)
    for i in range(5):
        reads.append(_make_read(f"ref{i}", "A" * 41, 95))

    samdata = tmp_bam(reads)
    pos_table, reach, masked, _ = _build_fragment_allele_table(samdata, [100, 110], REF_NAME)
    samdata.close()

    s_a = {'name': 's_a', 'position': '101', 'reference': 'A', 'variant': 'T',
           'depth': '15', 'basecalls': Counter({'T': 10, 'A': 5})}
    s_b = {'name': 's_b', 'position': '111', 'reference': 'A', 'variant': 'C',
           'depth': '14', 'basecalls': Counter({'C': 6, 'A': 8})}
    snp_list = [s_a, s_b]

    _apply_discover_roi(snp_list, pos_table, reach, masked, offset=0,
                         min_perc=0.0, min_reads=1, min_snp_perc=0.0)

    entry = next(e for e in s_a['linked_snps'] if e['target_name'] == 's_b')
    assert entry['co_count'] == 6
    assert entry['anchor_only'] == 3
    assert entry['target_only'] == 0
    assert entry['neither_variant'] == 5
    assert entry['read_too_short'] == 1
    assert entry['shared_read_depth'] == 14
    assert entry['linkage_pct'] == 66.7
    assert entry['sample_frequency_pct'] == 42.9


# ---------------------------------------------------------------------------
# Case 5: _add_linked_snps_node produces the expected XML structure
# ---------------------------------------------------------------------------
def test_add_linked_snps_node_xml_structure():
    """Purpose: verify _add_linked_snps_node serializes a linked_snps entry
    to a <linked_snp> element with anchor/target attributes and a
    <read_evidence> child element.

    Function under test: _add_linked_snps_node -- pure XML serialization of
    a single linked_snps entry as produced by _apply_discover_roi.

    Test input: a synthetic entry dict and anchor_name='s_a'.

    Expected result: <linked_snp anchor_variant="s_a" target_variant="s_b"
    shared_read_depth="14" co_occurring_count="6" linkage_pct="66.7"
    sample_frequency_pct="42.9"> with a <read_evidence> child carrying
    both_variants, anchor_only, target_only, neither_variant, masked, and
    read_too_short.
    """
    snp_node = ElementTree.Element('snp')
    entry = {
        'target_name': 's_b',
        'shared_read_depth': 14,
        'co_count': 6,
        'linkage_pct': 66.7,
        'sample_frequency_pct': 42.9,
        'both_variants': 6,
        'anchor_only': 3,
        'target_only': 0,
        'neither_variant': 5,
        'masked': 0,
        'read_too_short': 1,
    }

    _add_linked_snps_node(snp_node, [entry], anchor_name='s_a')

    link_node = snp_node.find('linked_snps/linked_snp')
    assert link_node.attrib['anchor_variant'] == 's_a'
    assert link_node.attrib['target_variant'] == 's_b'
    assert link_node.attrib['shared_read_depth'] == '14'
    assert link_node.attrib['co_occurring_count'] == '6'
    assert link_node.attrib['linkage_pct'] == '66.7'
    assert link_node.attrib['sample_frequency_pct'] == '42.9'

    ev = link_node.find('read_evidence')
    assert ev is not None
    assert ev.attrib['both_variants'] == '6'
    assert ev.attrib['anchor_only'] == '3'
    assert ev.attrib['target_only'] == '0'
    assert ev.attrib['neither_variant'] == '5'
    assert ev.attrib['masked'] == '0'
    assert ev.attrib['read_too_short'] == '1'
