"""
Unit tests for _apply_codon_correction and _add_codon_merges_node, and for
_apply_discover_roi's exclusion of codon_merge partners from linked_snps.

pos_table/reach/masked are constructed by hand (rather than via a synthetic
BAM) since _apply_codon_correction and _apply_discover_roi only consume the
plain {pos: {qname: base}} / {pos: (lo, hi)} / {pos: {qname, ...}} dicts
produced by _build_fragment_allele_table.
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


def _make_pos_table(pos_a, pos_b, combo_counts):
    """Build a {pos_a: {qname: base}, pos_b: {qname: base}} table where
    combo_counts maps (base_a, base_b) -> number of fragments carrying that
    combination at (pos_a, pos_b)."""
    pos_table = {pos_a: {}, pos_b: {}}
    i = 0
    for (base_a, base_b), count in combo_counts.items():
        for _ in range(count):
            qname = f"r{i}"
            pos_table[pos_a][qname] = base_a
            pos_table[pos_b][qname] = base_b
            i += 1
    return pos_table


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
    return {entry['name'] for entry in snp.get('linked_snps', [])}


# ---------------------------------------------------------------------------
# Case 1: two-SNP codon, complete merge (mirrors the real T5118A/T5119A
# example: variant frequencies match, so linkage=="complete").
# ---------------------------------------------------------------------------
def test_complete_merge_codon_annotation():
    """Purpose: verify that a 2-SNP codon where both SNPs' variant
    frequencies match within error_threshold gets annotated with one
    codon_merges entry per SNP, cross-referencing each other, with
    linkage=="complete" and combos classified as reference/variant.

    Function under test: _apply_codon_correction.

    Test input: codon spanning amp positions 100/101 (offset=0 -> translated
    101/102); pos_table with 17 fragments A|A and 3 fragments T|T
    (spanning_depth=20); both SNPs A->T with basecalls {A:17, T:3} (freq=0.15
    each); error_threshold=0.05, min_reads=1.

    Expected result: snp_list unchanged in length; both SNPs get exactly one
    codon_merges entry, cross-referencing linked_snp by name,
    spanning_depth==20, linkage=="complete", combos==
    [{"A|A",17,85.0,"reference"}, {"T|T",3,15.0,"variant"}].
    """
    pos_table = _make_pos_table(100, 101, {('A', 'A'): 17, ('T', 'T'): 3})
    snp_a = {'name': 'A101T', 'position': '101', 'reference': 'A', 'variant': 'T',
             'depth': '20', 'basecalls': Counter({'A': 17, 'T': 3})}
    snp_b = {'name': 'A102T', 'position': '102', 'reference': 'A', 'variant': 'T',
             'depth': '20', 'basecalls': Counter({'A': 17, 'T': 3})}
    snp_list = [snp_a, snp_b]
    cds = CdsFeature('ORF1', '+', [(100, 103)])

    result = _apply_codon_correction(snp_list, pos_table, [cds], offset=0,
                                      error_threshold=0.05, min_reads=1)

    assert result is None
    assert len(snp_list) == 2
    assert len(snp_a['codon_merges']) == 1
    assert len(snp_b['codon_merges']) == 1

    cm_a = snp_a['codon_merges'][0]
    cm_b = snp_b['codon_merges'][0]
    assert cm_a['linked_snp'] == 'A102T'
    assert cm_b['linked_snp'] == 'A101T'
    assert cm_a['spanning_depth'] == 20
    assert cm_b['spanning_depth'] == 20
    assert cm_a['linkage'] == 'complete'
    assert cm_b['linkage'] == 'complete'
    assert cm_a['combos'] == [
        {'bases': 'A|A', 'count': 17, 'percent': 85.0, 'type': 'reference'},
        {'bases': 'T|T', 'count': 3, 'percent': 15.0, 'type': 'variant'},
    ]
    assert cm_b['combos'] == cm_a['combos']


def test_add_codon_merges_node_xml():
    """Purpose: verify _add_codon_merges_node serializes a codon_merges entry
    to a <codon_merge> element with <combo> children carrying the bases,
    count, percent and type attributes.

    Function under test: _add_codon_merges_node -- pure XML serialization.

    Test input: a synthetic codon_merges list with one entry
    (linked_snp="A102T", spanning_depth=20, linkage="complete") containing
    two combos (a "reference" and a "variant" combo).

    Expected result: the resulting <codon_merge> element has the expected
    attributes, and its two <combo> children have the expected
    bases/count/percent/type attributes in order.
    """
    snp_node = ElementTree.Element('snp')
    codon_merges = [{
        'linked_snp': 'A102T',
        'spanning_depth': 20,
        'linkage': 'complete',
        'combos': [
            {'bases': 'A|A', 'count': 17, 'percent': 85.0, 'type': 'reference'},
            {'bases': 'T|T', 'count': 3, 'percent': 15.0, 'type': 'variant'},
        ],
    }]

    _add_codon_merges_node(snp_node, codon_merges)

    cm_node = snp_node.find('codon_merge')
    assert cm_node.attrib == {'linked_snp': 'A102T', 'spanning_depth': '20', 'linkage': 'complete'}
    combos = cm_node.findall('combo')
    assert len(combos) == 2
    assert combos[0].attrib == {'bases': 'A|A', 'count': '17', 'percent': '85.0', 'type': 'reference'}
    assert combos[1].attrib == {'bases': 'T|T', 'count': '3', 'percent': '15.0', 'type': 'variant'}


# ---------------------------------------------------------------------------
# Case 2: two-SNP codon, partial merge -- mismatched variant frequencies
# produce linkage=="partial" and a "discordant" combo.
# ---------------------------------------------------------------------------
def test_partial_merge_codon_annotation():
    """Purpose: verify that a 2-SNP codon where the SNPs' variant frequencies
    differ by more than error_threshold gets linkage=="partial", and that a
    fragment carrying only one of the two variants is classified as
    "discordant" (not "reference" or "variant").

    Function under test: _apply_codon_correction.

    Test input: codon spanning amp positions 100/101 (translated 101/102);
    pos_table with 10 fragments A|A (reference/reference), 5 fragments T|C
    (both variants together), 5 fragments T|A (only the first SNP's
    variant); s_a is A->T (basecalls {A:10, T:10}, freq=0.5), s_b is A->C
    (basecalls {A:15, C:5}, freq=0.25); error_threshold=0.05, min_reads=1.

    Expected result: both SNPs remain in snp_list unchanged; linkage==
    "partial"; spanning_depth==20; combos include bases "A|A"->"reference",
    "T|C"->"variant", "T|A"->"discordant".
    """
    pos_table = _make_pos_table(100, 101, {
        ('A', 'A'): 10,
        ('T', 'C'): 5,
        ('T', 'A'): 5,
    })
    snp_a = {'name': 'A101T', 'position': '101', 'reference': 'A', 'variant': 'T',
             'depth': '20', 'basecalls': Counter({'A': 10, 'T': 10})}
    snp_b = {'name': 'A102C', 'position': '102', 'reference': 'A', 'variant': 'C',
             'depth': '20', 'basecalls': Counter({'A': 15, 'C': 5})}
    snp_list = [snp_a, snp_b]
    cds = CdsFeature('ORF1', '+', [(100, 103)])

    _apply_codon_correction(snp_list, pos_table, [cds], offset=0,
                             error_threshold=0.05, min_reads=1)

    assert len(snp_list) == 2
    assert snp_a is snp_list[0]
    assert snp_b is snp_list[1]

    cm_a = snp_a['codon_merges'][0]
    assert cm_a['linkage'] == 'partial'
    assert cm_a['spanning_depth'] == 20

    types = {c['bases']: c['type'] for c in cm_a['combos']}
    assert types['A|A'] == 'reference'
    assert types['T|C'] == 'variant'
    assert types['T|A'] == 'discordant'


# ---------------------------------------------------------------------------
# Case 3: duplicate CDS features with an identical codon boundary (e.g. SC2's
# ORF1ab/ORF1a) must be processed once, not once per CDS feature.
# ---------------------------------------------------------------------------
def test_duplicate_cds_boundary_deduped():
    """Purpose: verify that when two CDS features (e.g. "ORF1ab" and "ORF1a")
    share an identical codon boundary, the codon is only annotated once --
    each SNP gets a single codon_merges entry, not one per CDS feature --
    while a second, non-duplicated boundary is still annotated normally.

    Function under test: _apply_codon_correction -- the all_codon_boundaries
    dedup across cds_features.

    Test input: two codons (amp 100/101 and amp 200/201), each a complete
    merge as in case 1 (17 ref/ref + 3 variant/variant fragments,
    spanning_depth=20). cds_features = [CdsFeature("ORF1ab", boundaries=
    [(100,103), (200,203)]), CdsFeature("ORF1a", boundaries=[(100,103)])].

    Expected result: every SNP ends up with len(codon_merges) == 1.
    """
    pos_table = {}
    pos_table.update(_make_pos_table(100, 101, {('A', 'A'): 17, ('T', 'T'): 3}))
    pos_table.update(_make_pos_table(200, 201, {('G', 'G'): 17, ('C', 'C'): 3}))

    snp_a = {'name': 'A101T', 'position': '101', 'reference': 'A', 'variant': 'T',
             'depth': '20', 'basecalls': Counter({'A': 17, 'T': 3})}
    snp_b = {'name': 'A102T', 'position': '102', 'reference': 'A', 'variant': 'T',
             'depth': '20', 'basecalls': Counter({'A': 17, 'T': 3})}
    snp_c = {'name': 'G201C', 'position': '201', 'reference': 'G', 'variant': 'C',
             'depth': '20', 'basecalls': Counter({'G': 17, 'C': 3})}
    snp_d = {'name': 'G202C', 'position': '202', 'reference': 'G', 'variant': 'C',
             'depth': '20', 'basecalls': Counter({'G': 17, 'C': 3})}
    snp_list = [snp_a, snp_b, snp_c, snp_d]

    cds_orf1ab = CdsFeature('ORF1ab', '+', [(100, 103), (200, 203)])
    cds_orf1a = CdsFeature('ORF1a', '+', [(100, 103)])

    _apply_codon_correction(snp_list, pos_table, [cds_orf1ab, cds_orf1a], offset=0,
                             error_threshold=0.05, min_reads=1)

    assert len(snp_a['codon_merges']) == 1
    assert len(snp_b['codon_merges']) == 1
    assert len(snp_c['codon_merges']) == 1
    assert len(snp_d['codon_merges']) == 1


# ---------------------------------------------------------------------------
# Case 4: a SNP belonging to two distinct (non-identical) codon boundaries
# (overlapping reading frames) accumulates one codon_merges entry per codon.
# ---------------------------------------------------------------------------
def test_multi_codon_membership():
    """Purpose: verify that a SNP participating in two different qualifying
    codons (e.g. overlapping reading frames in different ORFs) accumulates
    one codon_merges entry per codon, each with the correct linked_snp.

    Function under test: _apply_codon_correction -- multiple codon_boundaries
    sharing a SNP position.

    Test input: amp positions 100/101/102 (translated 101/102/103); codon1 =
    (100,102) pairs translated 101/102, codon2 = (101,103) pairs translated
    102/103. 17 fragments carry A|T|G, 3 fragments carry T|A|C across
    positions 100/101/102. SNP at translated 102 (amp 101) is the shared
    "middle" SNP.

    Expected result: the middle SNP gets codon_merges of length 2, with
    linked_snp values {"A101T", "G103C"}; the other two SNPs each get
    codon_merges of length 1.
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

    cds = CdsFeature('ORFx', '+', [(100, 102), (101, 103)])

    _apply_codon_correction(snp_list, pos_table, [cds], offset=0,
                             error_threshold=0.05, min_reads=1)

    assert len(snp_left['codon_merges']) == 1
    assert len(snp_right['codon_merges']) == 1
    assert len(snp_mid['codon_merges']) == 2
    assert {cm['linked_snp'] for cm in snp_mid['codon_merges']} == {'A101T', 'G103C'}


# ---------------------------------------------------------------------------
# Case 5: no-op / regression cases
# ---------------------------------------------------------------------------
def test_noop_empty_inputs():
    """Purpose: verify _apply_codon_correction is a no-op (returns None,
    mutates nothing) when cds_features or snp_list is empty.

    Function under test: _apply_codon_correction -- early-return guard.
    """
    assert _apply_codon_correction([], {}, [], offset=0,
                                    error_threshold=0.05, min_reads=1) is None

    snp = {'name': 'A101T', 'position': '101', 'reference': 'A', 'variant': 'T',
           'depth': '20', 'basecalls': Counter({'A': 17, 'T': 3})}
    cds = CdsFeature('ORF1', '+', [(100, 103)])
    assert _apply_codon_correction([snp], {}, [], offset=0,
                                    error_threshold=0.05, min_reads=1) is None
    assert _apply_codon_correction([snp], {}, [cds], offset=0,
                                    error_threshold=0.05, min_reads=1) is None
    assert 'codon_merges' not in snp


def test_three_snp_codon_not_annotated():
    """Purpose: verify a codon containing three polymorphic positions is
    left unannotated -- only len(codon_snps) == 2 codons are processed.

    Function under test: _apply_codon_correction -- the
    `if len(codon_snps) != 2: continue` guard.

    Test input: codon (100,103) covers translated 101/102/103, and snp_list
    has a SNP at each of those three positions.

    Expected result: none of the three SNPs gets a 'codon_merges' key.
    """
    snp_a = {'name': 'A101T', 'position': '101', 'reference': 'A', 'variant': 'T',
             'depth': '20', 'basecalls': Counter({'A': 17, 'T': 3})}
    snp_b = {'name': 'A102T', 'position': '102', 'reference': 'A', 'variant': 'T',
             'depth': '20', 'basecalls': Counter({'A': 17, 'T': 3})}
    snp_c = {'name': 'A103T', 'position': '103', 'reference': 'A', 'variant': 'T',
             'depth': '20', 'basecalls': Counter({'A': 17, 'T': 3})}
    snp_list = [snp_a, snp_b, snp_c]
    cds = CdsFeature('ORF1', '+', [(100, 103)])

    _apply_codon_correction(snp_list, {}, [cds], offset=0,
                             error_threshold=0.05, min_reads=1)

    assert 'codon_merges' not in snp_a
    assert 'codon_merges' not in snp_b
    assert 'codon_merges' not in snp_c
    assert len(snp_list) == 3


def test_min_reads_gate_skips_codon():
    """Purpose: verify a codon is skipped entirely (no codon_merges added)
    when the most common allele combination has fewer than min_reads
    supporting fragments.

    Function under test: _apply_codon_correction -- the
    `if dominant_combo_count < min_reads: continue` guard.

    Test input: codon (100,103) over translated 101/102; pos_table with 2
    A|A fragments + 1 T|T fragment (dominant_combo_count=2); min_reads=5.

    Expected result: neither SNP gets a 'codon_merges' key.
    """
    pos_table = _make_pos_table(100, 101, {('A', 'A'): 2, ('T', 'T'): 1})
    snp_a = {'name': 'A101T', 'position': '101', 'reference': 'A', 'variant': 'T',
             'depth': '3', 'basecalls': Counter({'A': 2, 'T': 1})}
    snp_b = {'name': 'A102T', 'position': '102', 'reference': 'A', 'variant': 'T',
             'depth': '3', 'basecalls': Counter({'A': 2, 'T': 1})}
    snp_list = [snp_a, snp_b]
    cds = CdsFeature('ORF1', '+', [(100, 103)])

    _apply_codon_correction(snp_list, pos_table, [cds], offset=0,
                             error_threshold=0.05, min_reads=5)

    assert 'codon_merges' not in snp_a
    assert 'codon_merges' not in snp_b


# ---------------------------------------------------------------------------
# Case 6: pairs already covered by <codon_merge> are excluded from
# linked_snps, while unrelated pairs still get linked_snps as usual.
# ---------------------------------------------------------------------------
def test_codon_merge_pairs_excluded_from_linked_snps():
    """Purpose: verify that a SNP pair annotated via codon_merges (by
    _apply_codon_correction) does NOT also show up in each other's
    linked_snps (from _apply_discover_roi), while an unrelated SNP pair --
    not part of any codon -- still links normally.

    Functions under test: _apply_codon_correction + _apply_discover_roi run
    in sequence (matching the real call order), focused on the
    `codon_partners` exclusion in _apply_discover_roi.

    Test input: codon pair at amp 100/101 (translated 101/102, A->T, 17
    ref/ref + 3 variant/variant fragments -- same as case 1, would otherwise
    qualify for linked_snps under loose thresholds); unrelated pair at amp
    200/201 (translated 201/202, G->C, same 17/3 split) not covered by any
    CDS codon boundary. discover_roi called with min_perc=0.0, min_reads=1,
    min_snp_perc=0.0.

    Expected result: codon_merges populated for the amp 100/101 pair as in
    case 1; neither SNP appears in the other's linked_snps. The amp 200/201
    pair has no codon_merges and DOES appear in each other's linked_snps.
    """
    pos_table = {}
    pos_table.update(_make_pos_table(100, 101, {('A', 'A'): 17, ('T', 'T'): 3}))
    pos_table.update(_make_pos_table(200, 201, {('G', 'G'): 17, ('C', 'C'): 3}))

    reach = {100: (90, 150), 101: (90, 150), 200: (190, 250), 201: (190, 250)}
    masked = {100: set(), 101: set(), 200: set(), 201: set()}

    snp_a = {'name': 'A101T', 'position': '101', 'reference': 'A', 'variant': 'T',
             'depth': '20', 'basecalls': Counter({'A': 17, 'T': 3})}
    snp_b = {'name': 'A102T', 'position': '102', 'reference': 'A', 'variant': 'T',
             'depth': '20', 'basecalls': Counter({'A': 17, 'T': 3})}
    snp_e = {'name': 'G201C', 'position': '201', 'reference': 'G', 'variant': 'C',
             'depth': '20', 'basecalls': Counter({'G': 17, 'C': 3})}
    snp_f = {'name': 'G202C', 'position': '202', 'reference': 'G', 'variant': 'C',
             'depth': '20', 'basecalls': Counter({'G': 17, 'C': 3})}
    snp_list = [snp_a, snp_b, snp_e, snp_f]

    cds = CdsFeature('ORF1', '+', [(100, 103)])

    _apply_codon_correction(snp_list, pos_table, [cds], offset=0,
                             error_threshold=0.05, min_reads=1)
    _apply_discover_roi(snp_list, pos_table, reach, masked, offset=0,
                         min_perc=0.0, min_reads=1, min_snp_perc=0.0)

    assert snp_a['codon_merges'][0]['linked_snp'] == 'A102T'
    assert snp_b['codon_merges'][0]['linked_snp'] == 'A101T'
    assert 'A102T' not in _linked_names(snp_a)
    assert 'A101T' not in _linked_names(snp_b)

    assert 'codon_merges' not in snp_e
    assert 'codon_merges' not in snp_f
    assert 'G202C' in _linked_names(snp_e)
    assert 'G201C' in _linked_names(snp_f)
