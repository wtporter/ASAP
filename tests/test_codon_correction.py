"""
Unit tests for _apply_codon_correction and _add_codon_merges_node, and for
_apply_discover_roi's exclusion of codon_merge partners from linked_snps.

pos_table/reach/masked/deleted are constructed by hand (rather than via a
synthetic BAM) since _apply_codon_correction and _apply_discover_roi only
consume the plain {pos: {qname: base}} / {pos: (lo, hi)} / {pos: {qname, ...}}
dicts produced by _build_fragment_allele_table.
"""
import os
import sys
from collections import Counter
from xml.etree import ElementTree

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from asap.newBamProcessor import (
    CdsFeature,
    _apply_codon_correction,
    _add_codon_merges_node,
    _apply_discover_roi,
)


def _make_pos_table_3way(pos_a, pos_b, pos_c, rows):
    """Build a 3-position pos_table; rows is a list of (base_a, base_b,
    base_c, count) tuples, one per group of fragments."""
    pos_table = {pos_a: {}, pos_b: {}, pos_c: {}}
    i = 0
    for base_a, base_b, base_c, count in rows:
        for _ in range(count):
            qname = f"r{i}"
            pos_table[pos_a][qname] = base_a
            pos_table[pos_b][qname] = base_b
            pos_table[pos_c][qname] = base_c
            i += 1
    return pos_table


def _linked_names(snp):
    return {entry['target_name'] for entry in snp.get('linked_snps', [])}


# ---------------------------------------------------------------------------
# Case 1: two-SNP codon, basic codon distribution.
# ---------------------------------------------------------------------------
def test_complete_merge_codon_annotation():
    """Purpose: verify that a 2-SNP codon gets one codon_merges entry per SNP
    showing the full 3-bp codon distribution.

    Codon boundary (100,103,'AAA'): positions 100/101/102 (translated 101/102/103).
    17 reads: A at 100, A at 101, A at 102 -> 'AAA' (reference codon).
    3 reads: T at 100, T at 101, A at 102 -> 'TTA' (both SNPs mutated).
    codon_depth=20; dominant call 'AAA' at 85.0%.
    """
    pos_table = _make_pos_table_3way(100, 101, 102, [
        ('A', 'A', 'A', 17),
        ('T', 'T', 'A', 3),
    ])
    snp_a = {'name': 'A101T', 'position': '101', 'reference': 'A', 'variant': 'T',
             'depth': '20', 'basecalls': Counter({'A': 17, 'T': 3})}
    snp_b = {'name': 'A102T', 'position': '102', 'reference': 'A', 'variant': 'T',
             'depth': '20', 'basecalls': Counter({'A': 17, 'T': 3})}
    snp_list = [snp_a, snp_b]
    cds = CdsFeature('ORF1', '+', [(100, 103, 'AAA')])

    result = _apply_codon_correction(snp_list, pos_table, {}, {}, [cds], offset=0,
                                      error_threshold=0.05, min_reads=1)

    assert result is None
    assert len(snp_list) == 2
    assert len(snp_a['codon_merges']) == 1
    assert len(snp_b['codon_merges']) == 1

    cm_a = snp_a['codon_merges'][0]
    cm_b = snp_b['codon_merges'][0]
    assert cm_a is cm_b  # same entry object shared by both SNPs

    assert cm_a['name'] == 'ORF1_codon_1'
    assert cm_a['region'] == 'ORF1'
    assert cm_a['direction'] == 'forward'
    assert cm_a['position'] == '101-103'
    assert cm_a['reference'] == 'AAA'
    assert cm_a['codon_depth'] == 20
    assert cm_a['codon_call'] == 'AAA'
    assert cm_a['codon_call_count'] == 17
    assert cm_a['codon_call_percentage'] == 85.0
    assert cm_a['codon_distribution'] == {'AAA': 17, 'TTA': 3}

    # codon_partner_names set is populated for discover_roi exclusion
    assert snp_a.get('codon_partner_names') == {'A102T'}
    assert snp_b.get('codon_partner_names') == {'A101T'}


def test_add_codon_merges_node_xml():
    """Purpose: verify _add_codon_merges_node serializes a codon_merges entry
    to the new <codon_merge> XML format with <codon_call>, <codon_distribution>,
    and <excluded_reads> children.
    """
    snp_node = ElementTree.Element('snp')
    codon_merges = [{
        'name': 'ORF1_codon_1',
        'region': 'ORF1',
        'direction': 'forward',
        'position': '101-103',
        'codon_depth': 20,
        'reference': 'AAA',
        'codon_call': 'AAA',
        'codon_call_count': 17,
        'codon_call_percentage': 85.0,
        'codon_distribution': {'AAA': 17, 'TTA': 3},
        'excl_has_n': 2,
        'excl_no_span': 5,
    }]

    _add_codon_merges_node(snp_node, codon_merges)

    cm_node = snp_node.find('codon_merge')
    assert cm_node.attrib == {
        'name': 'ORF1_codon_1',
        'region': 'ORF1',
        'direction': 'forward',
        'position': '101-103',
        'codon_depth': '20',
        'reference': 'AAA',
    }

    codon_call = cm_node.find('codon_call')
    assert codon_call.text == 'AAA'
    assert codon_call.attrib == {'codon_count': '17', 'percentage': '85.0'}

    codon_dist = cm_node.find('codon_distribution')
    # sorted by descending count: AAA=17 before TTA=3
    assert codon_dist.attrib == {'AAA': '17', 'TTA': '3'}

    excluded = cm_node.find('excluded_reads')
    assert excluded.attrib == {'has_n': '2', 'no_span': '5'}


# ---------------------------------------------------------------------------
# Case 2: mixed codon distribution -- multiple distinct codon sequences.
# ---------------------------------------------------------------------------
def test_mixed_codon_distribution():
    """Purpose: verify that multiple distinct codon sequences each appear as
    separate entries in codon_distribution.

    Codon (100,103,'AAA'): positions 100/101/102 (translated 101/102/103).
    10 reads: A|A|A -> 'AAA'; 5 reads: T|C|A -> 'TCA'; 5 reads: T|A|A -> 'TAA'.
    """
    pos_table = _make_pos_table_3way(100, 101, 102, [
        ('A', 'A', 'A', 10),
        ('T', 'C', 'A', 5),
        ('T', 'A', 'A', 5),
    ])
    snp_a = {'name': 'A101T', 'position': '101', 'reference': 'A', 'variant': 'T',
             'depth': '20', 'basecalls': Counter({'A': 10, 'T': 10})}
    snp_b = {'name': 'A102C', 'position': '102', 'reference': 'A', 'variant': 'C',
             'depth': '20', 'basecalls': Counter({'A': 15, 'C': 5})}
    snp_list = [snp_a, snp_b]
    cds = CdsFeature('ORF1', '+', [(100, 103, 'AAA')])

    _apply_codon_correction(snp_list, pos_table, {}, {}, [cds], offset=0,
                             error_threshold=0.05, min_reads=1)

    assert len(snp_list) == 2
    cm = snp_a['codon_merges'][0]
    assert cm['codon_depth'] == 20
    assert cm['codon_call'] == 'AAA'
    assert cm['codon_distribution'] == {'AAA': 10, 'TCA': 5, 'TAA': 5}


# ---------------------------------------------------------------------------
# Case 3: duplicate CDS features with identical boundary each produce their
# own codon_merges entry (one per CDS, no deduplication).
# ---------------------------------------------------------------------------
def test_duplicate_cds_boundary_produces_one_entry_per_cds():
    """Purpose: verify that when two CDS features (e.g. "ORF1ab" and "ORF1a")
    share an identical codon boundary, each CDS produces its own codon_merges
    entry. SNPs at the shared boundary get 2 entries (one per CDS), while
    SNPs at the unique ORF1ab-only boundary get 1 entry.

    cds_features = [CdsFeature("ORF1ab", [(100,103,'AAA'),(200,203,'GGG')]),
                    CdsFeature("ORF1a",  [(100,103,'AAA')])]

    SNPs at translated 101/102 -> 2 entries each (ORF1ab codon 1 + ORF1a codon 1).
    SNPs at translated 201/202 -> 1 entry each (ORF1ab codon 2 only).
    """
    pos_table = {}
    pos_table.update(_make_pos_table_3way(100, 101, 102, [
        ('A', 'A', 'A', 17), ('T', 'T', 'A', 3),
    ]))
    pos_table.update(_make_pos_table_3way(200, 201, 202, [
        ('G', 'G', 'G', 17), ('C', 'C', 'G', 3),
    ]))

    snp_a = {'name': 'A101T', 'position': '101', 'reference': 'A', 'variant': 'T',
             'depth': '20', 'basecalls': Counter({'A': 17, 'T': 3})}
    snp_b = {'name': 'A102T', 'position': '102', 'reference': 'A', 'variant': 'T',
             'depth': '20', 'basecalls': Counter({'A': 17, 'T': 3})}
    snp_c = {'name': 'G201C', 'position': '201', 'reference': 'G', 'variant': 'C',
             'depth': '20', 'basecalls': Counter({'G': 17, 'C': 3})}
    snp_d = {'name': 'G202C', 'position': '202', 'reference': 'G', 'variant': 'C',
             'depth': '20', 'basecalls': Counter({'G': 17, 'C': 3})}
    snp_list = [snp_a, snp_b, snp_c, snp_d]

    cds_orf1ab = CdsFeature('ORF1ab', '+', [(100, 103, 'AAA'), (200, 203, 'GGG')])
    cds_orf1a = CdsFeature('ORF1a', '+', [(100, 103, 'AAA')])

    _apply_codon_correction(snp_list, pos_table, {}, {}, [cds_orf1ab, cds_orf1a],
                             offset=0, error_threshold=0.05, min_reads=1)

    # Shared boundary -> 2 entries each
    assert len(snp_a['codon_merges']) == 2
    assert len(snp_b['codon_merges']) == 2
    assert {cm['region'] for cm in snp_a['codon_merges']} == {'ORF1ab', 'ORF1a'}
    assert {cm['name'] for cm in snp_a['codon_merges']} == {
        'ORF1ab_codon_1', 'ORF1a_codon_1',
    }

    # Unique boundary -> 1 entry each
    assert len(snp_c['codon_merges']) == 1
    assert len(snp_d['codon_merges']) == 1
    assert snp_c['codon_merges'][0]['region'] == 'ORF1ab'
    assert snp_c['codon_merges'][0]['name'] == 'ORF1ab_codon_2'


# ---------------------------------------------------------------------------
# Case 4: a SNP belonging to two distinct (non-identical) codon boundaries
# (overlapping reading frames) accumulates one codon_merges entry per codon.
# ---------------------------------------------------------------------------
def test_multi_codon_membership():
    """Purpose: verify that a SNP participating in two different qualifying
    codons (overlapping reading frames) accumulates one codon_merges entry
    per codon, each with a distinct name and position.

    amp positions 100/101/102 (translated 101/102/103):
    codon1 = (100,102,'AT') pairs translated 101/102,
    codon2 = (101,103,'TG') pairs translated 102/103.
    The SNP at translated 102 (amp 101) is the shared middle position.
    """
    pos_table = _make_pos_table_3way(100, 101, 102, [
        ('A', 'T', 'G', 17),
        ('T', 'A', 'C', 3),
    ])

    snp_left = {'name': 'A101T', 'position': '101', 'reference': 'A', 'variant': 'T',
                'depth': '20', 'basecalls': Counter({'A': 17, 'T': 3})}
    snp_mid = {'name': 'T102A', 'position': '102', 'reference': 'T', 'variant': 'A',
               'depth': '20', 'basecalls': Counter({'T': 17, 'A': 3})}
    snp_right = {'name': 'G103C', 'position': '103', 'reference': 'G', 'variant': 'C',
                 'depth': '20', 'basecalls': Counter({'G': 17, 'C': 3})}
    snp_list = [snp_left, snp_mid, snp_right]

    cds = CdsFeature('ORFx', '+', [(100, 102, 'AT'), (101, 103, 'TG')])

    _apply_codon_correction(snp_list, pos_table, {}, {}, [cds], offset=0,
                             error_threshold=0.05, min_reads=1)

    assert len(snp_left['codon_merges']) == 1
    assert len(snp_right['codon_merges']) == 1
    assert len(snp_mid['codon_merges']) == 2
    assert {cm['name'] for cm in snp_mid['codon_merges']} == {
        'ORFx_codon_1', 'ORFx_codon_2',
    }


# ---------------------------------------------------------------------------
# Case 5: no-op / regression cases
# ---------------------------------------------------------------------------
def test_noop_empty_inputs():
    """Purpose: verify _apply_codon_correction is a no-op (returns None,
    mutates nothing) when cds_features or snp_list is empty.
    """
    assert _apply_codon_correction([], {}, {}, {}, [], offset=0,
                                    error_threshold=0.05, min_reads=1) is None

    snp = {'name': 'A101T', 'position': '101', 'reference': 'A', 'variant': 'T',
           'depth': '20', 'basecalls': Counter({'A': 17, 'T': 3})}
    cds = CdsFeature('ORF1', '+', [(100, 103, 'AAA')])
    assert _apply_codon_correction([snp], {}, {}, {}, [], offset=0,
                                    error_threshold=0.05, min_reads=1) is None
    assert _apply_codon_correction([snp], {}, {}, {}, [cds], offset=0,
                                    error_threshold=0.05, min_reads=1) is None
    assert 'codon_merges' not in snp


def test_three_snp_codon_not_annotated():
    """Purpose: verify a codon containing three polymorphic positions is left
    unannotated -- only len(codon_snps) == 2 codons are processed.
    """
    snp_a = {'name': 'A101T', 'position': '101', 'reference': 'A', 'variant': 'T',
             'depth': '20', 'basecalls': Counter({'A': 17, 'T': 3})}
    snp_b = {'name': 'A102T', 'position': '102', 'reference': 'A', 'variant': 'T',
             'depth': '20', 'basecalls': Counter({'A': 17, 'T': 3})}
    snp_c = {'name': 'A103T', 'position': '103', 'reference': 'A', 'variant': 'T',
             'depth': '20', 'basecalls': Counter({'A': 17, 'T': 3})}
    snp_list = [snp_a, snp_b, snp_c]
    cds = CdsFeature('ORF1', '+', [(100, 103, 'AAA')])

    _apply_codon_correction(snp_list, {}, {}, {}, [cds], offset=0,
                             error_threshold=0.05, min_reads=1)

    assert 'codon_merges' not in snp_a
    assert 'codon_merges' not in snp_b
    assert 'codon_merges' not in snp_c
    assert len(snp_list) == 3


def test_min_reads_gate_skips_codon():
    """Purpose: verify a codon is skipped when the dominant combination has
    fewer than min_reads supporting fragments.
    """
    pos_table = _make_pos_table_3way(100, 101, 102, [
        ('A', 'A', 'A', 2), ('T', 'T', 'A', 1),
    ])
    snp_a = {'name': 'A101T', 'position': '101', 'reference': 'A', 'variant': 'T',
             'depth': '3', 'basecalls': Counter({'A': 2, 'T': 1})}
    snp_b = {'name': 'A102T', 'position': '102', 'reference': 'A', 'variant': 'T',
             'depth': '3', 'basecalls': Counter({'A': 2, 'T': 1})}
    snp_list = [snp_a, snp_b]
    cds = CdsFeature('ORF1', '+', [(100, 103, 'AAA')])

    _apply_codon_correction(snp_list, pos_table, {}, {}, [cds], offset=0,
                             error_threshold=0.05, min_reads=5)

    assert 'codon_merges' not in snp_a
    assert 'codon_merges' not in snp_b


# ---------------------------------------------------------------------------
# Case 6: pairs covered by codon_merges are excluded from linked_snps.
# ---------------------------------------------------------------------------
def test_codon_merge_pairs_excluded_from_linked_snps():
    """Purpose: verify that a SNP pair annotated via codon_merges does NOT also
    appear in each other's linked_snps, while an unrelated SNP pair still links.

    Codon pair at amp 100/101/102 (ORF1 codon 1); unrelated pair at amp
    200/201 not covered by any CDS boundary.
    """
    pos_table = {}
    pos_table.update(_make_pos_table_3way(100, 101, 102,
                                          [('A', 'A', 'A', 17), ('T', 'T', 'A', 3)]))
    # Unrelated pair at 200/201 (2-position table, no 3rd position needed for discover_roi)
    pos_table.setdefault(200, {})
    pos_table.setdefault(201, {})
    for i in range(17):
        pos_table[200][f's{i}'] = 'G'
        pos_table[201][f's{i}'] = 'G'
    for i in range(3):
        pos_table[200][f'sv{i}'] = 'C'
        pos_table[201][f'sv{i}'] = 'C'

    reach = {
        100: (90, 150), 101: (90, 150), 102: (90, 150),
        200: (190, 250), 201: (190, 250),
    }
    masked = {p: set() for p in [100, 101, 102, 200, 201]}

    snp_a = {'name': 'A101T', 'position': '101', 'reference': 'A', 'variant': 'T',
             'depth': '20', 'basecalls': Counter({'A': 17, 'T': 3})}
    snp_b = {'name': 'A102T', 'position': '102', 'reference': 'A', 'variant': 'T',
             'depth': '20', 'basecalls': Counter({'A': 17, 'T': 3})}
    snp_e = {'name': 'G201C', 'position': '201', 'reference': 'G', 'variant': 'C',
             'depth': '20', 'basecalls': Counter({'G': 17, 'C': 3})}
    snp_f = {'name': 'G202C', 'position': '202', 'reference': 'G', 'variant': 'C',
             'depth': '20', 'basecalls': Counter({'G': 17, 'C': 3})}
    snp_list = [snp_a, snp_b, snp_e, snp_f]

    cds = CdsFeature('ORF1', '+', [(100, 103, 'AAA')])

    _apply_codon_correction(snp_list, pos_table, {}, masked, [cds], offset=0,
                             error_threshold=0.05, min_reads=1)
    _apply_discover_roi(snp_list, pos_table, reach, masked, offset=0,
                         min_perc=0.0, min_reads=1, min_snp_perc=0.0)

    # codon pair has codon_merges; neither appears in the other's linked_snps
    assert 'codon_merges' in snp_a
    assert 'codon_merges' in snp_b
    assert 'A102T' not in _linked_names(snp_a)
    assert 'A101T' not in _linked_names(snp_b)

    # unrelated pair has no codon_merges and appears in linked_snps
    assert 'codon_merges' not in snp_e
    assert 'codon_merges' not in snp_f
    assert 'G202C' in _linked_names(snp_e)
    assert 'G201C' in _linked_names(snp_f)


# ---------------------------------------------------------------------------
# Case 7: full-codon deletion (all 3 positions deleted) appears in distribution.
# ---------------------------------------------------------------------------
def test_full_codon_deletion_in_distribution():
    """Purpose: verify that reads with deletions at all 3 codon positions
    appear as '___' in codon_distribution and contribute to codon_depth.

    5 reads with reference A|A|A; 3 reads with deletion at all 3 positions.
    deleted = {100: {'d0','d1','d2'}, 101: {'d0','d1','d2'}, 102: {'d0','d1','d2'}}.
    codon_depth = 8; '___' count = 3.
    """
    pos_table = _make_pos_table_3way(100, 101, 102, [('A', 'A', 'A', 5)])
    deleted = {
        100: {'d0', 'd1', 'd2'},
        101: {'d0', 'd1', 'd2'},
        102: {'d0', 'd1', 'd2'},
    }

    snp_a = {'name': 'A101_', 'position': '101', 'reference': 'A', 'variant': '_',
             'depth': '8', 'basecalls': Counter({'_': 3, 'A': 5})}
    snp_b = {'name': 'A102_', 'position': '102', 'reference': 'A', 'variant': '_',
             'depth': '8', 'basecalls': Counter({'_': 3, 'A': 5})}
    snp_list = [snp_a, snp_b]
    cds = CdsFeature('ORF1', '+', [(100, 103, 'AAA')])

    _apply_codon_correction(snp_list, pos_table, deleted, {}, [cds], offset=0,
                             error_threshold=0.05, min_reads=1)

    assert len(snp_a['codon_merges']) == 1
    cm = snp_a['codon_merges'][0]
    assert cm['codon_depth'] == 8
    assert cm['codon_distribution'] == {'AAA': 5, '___': 3}
    assert cm['codon_call'] == 'AAA'
    assert cm['codon_call_count'] == 5


# ---------------------------------------------------------------------------
# Case 8: partial deletion (1-2 positions deleted) appears in distribution.
# ---------------------------------------------------------------------------
def test_partial_deletion_in_codon_distribution():
    """Purpose: verify that a read with a deletion at one codon position and
    real bases at the other two appears in codon_distribution as a mixed
    sequence like '_AA'.

    10 reads A|A|A; 3 reads deleted at position 100, reference at 101/102.
    deleted = {100: {'pd0','pd1','pd2'}}.
    Expected codon_distribution: {'AAA': 10, '_AA': 3}.
    """
    pos_table = _make_pos_table_3way(100, 101, 102, [('A', 'A', 'A', 10)])
    # 3 reads: deleted at 100, A at 101, A at 102
    for i in range(3):
        qname = f'pd{i}'
        pos_table[101][qname] = 'A'
        pos_table[102][qname] = 'A'
    deleted = {100: {'pd0', 'pd1', 'pd2'}}

    snp_a = {'name': 'A101T', 'position': '101', 'reference': 'A', 'variant': 'T',
             'depth': '13', 'basecalls': Counter({'A': 13})}
    snp_b = {'name': 'A102T', 'position': '102', 'reference': 'A', 'variant': 'T',
             'depth': '13', 'basecalls': Counter({'A': 13})}
    snp_list = [snp_a, snp_b]
    cds = CdsFeature('ORF1', '+', [(100, 103, 'AAA')])

    _apply_codon_correction(snp_list, pos_table, deleted, {}, [cds], offset=0,
                             error_threshold=0.05, min_reads=1)

    cm = snp_a['codon_merges'][0]
    assert cm['codon_depth'] == 13
    assert cm['codon_distribution'] == {'AAA': 10, '_AA': 3}


# ---------------------------------------------------------------------------
# Case 9: excluded_reads counts are reported correctly.
# ---------------------------------------------------------------------------
def test_excluded_reads_counts():
    """Purpose: verify excl_has_n and excl_no_span are correct.

    10 reads span all 3 positions (included in codon_depth).
    2 reads have N at position 101 (in masked[101]) -> excl_has_n >= 2.
    3 reads are present at position 100 only, not at 101/102 (short reads)
       -> they appear in any_definite but not in tally -> excl_no_span >= 3.
    """
    pos_table = _make_pos_table_3way(100, 101, 102, [
        ('A', 'A', 'A', 7),
        ('T', 'T', 'A', 3),
    ])
    # 3 short reads at position 100 only
    for i in range(3):
        pos_table[100][f'short{i}'] = 'A'

    masked = {101: {'n0', 'n1'}}  # 2 reads have N at position 101

    snp_a = {'name': 'A101T', 'position': '101', 'reference': 'A', 'variant': 'T',
             'depth': '15', 'basecalls': Counter({'A': 12, 'T': 3})}
    snp_b = {'name': 'A102T', 'position': '102', 'reference': 'A', 'variant': 'T',
             'depth': '13', 'basecalls': Counter({'A': 10, 'T': 3})}
    snp_list = [snp_a, snp_b]
    cds = CdsFeature('ORF1', '+', [(100, 103, 'AAA')])

    _apply_codon_correction(snp_list, pos_table, {}, masked, [cds], offset=0,
                             error_threshold=0.05, min_reads=1)

    cm = snp_a['codon_merges'][0]
    assert cm['codon_depth'] == 10
    assert cm['excl_has_n'] == 2     # 2 reads with N at position 101
    assert cm['excl_no_span'] == 3   # 3 short reads present only at position 100
