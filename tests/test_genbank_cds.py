"""
Unit tests for _parse_genbank_cds.

Uses minimal synthetic GenBank files written to temp files so the tests
are fast and self-contained.  One integration-style test uses the real
H37Rv test fixture if it is available.
"""
import os
import sys
import textwrap
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from asap.newBamProcessor import _parse_genbank_cds, _load_genbank_records

H37RV_GB = os.path.join(
    os.path.dirname(__file__),
    "../nextflow/tests/preparejson/H37Rv_NC0009623.gb",
)


def _write_gb(tmp_path, name, seq, features_block):
    """Write a minimal GenBank file and return its path."""
    seq_clean = seq.upper().replace(" ", "").replace("\n", "")
    seq_len = len(seq_clean)

    # Format ORIGIN block (10 bases per group, 6 groups per line)
    origin_lines = []
    for i in range(0, seq_len, 60):
        chunk = seq_clean[i:i + 60]
        groups = " ".join(chunk[j:j + 10] for j in range(0, len(chunk), 10))
        origin_lines.append(f"       {i + 1:>9} {groups.lower()}")

    # Re-indent the (already dedented) features block under FEATURES: 5 spaces
    # for the feature key, 21 for qualifiers -- skbio requires this to treat
    # them as part of the FEATURES table rather than top-level section headers.
    indented_features = "\n".join(
        ("     " + ln if ln.strip() else ln)
        for ln in features_block.rstrip("\n").split("\n")
    )

    # Built line-by-line (not via textwrap.dedent on an f-string) so the
    # interpolated features_block/origin_lines -- which carry their own
    # indentation -- don't affect dedent's common-prefix calculation for the
    # LOCUS/DEFINITION/... header lines, which skbio requires to start at
    # column 0.
    lines = [
        f"LOCUS       {name:<16} {seq_len} bp    DNA     linear   VRL 01-JAN-2000",
        "DEFINITION  Synthetic test sequence.",
        f"ACCESSION   {name}",
        f"VERSION     {name}.1",
        "FEATURES             Location/Qualifiers",
        indented_features,
        "ORIGIN",
        *origin_lines,
        "//",
    ]
    content = "\n".join(lines) + "\n"
    gb_path = str(tmp_path / f"{name}.gb")
    with open(gb_path, "w") as fh:
        fh.write(content)
    return gb_path


# ---------------------------------------------------------------------------
# Case 1: simple forward-strand CDS — verify codon count and boundaries
# ---------------------------------------------------------------------------
def test_forward_cds_basic(tmp_path):
    """Purpose: verify in-frame forward-strand CDS overlap is detected and
    codon boundaries are computed correctly.

    Function under test: _parse_genbank_cds -- amplicon location, CDS-overlap
    detection, and codon-boundary computation when the amplicon starts
    exactly on a codon boundary (zero frame offset).

    Test input: 30-bp genome = 10 "ATG" codons, with a single forward CDS
    spanning 1..30 (gene="aaa"). Amplicon = genome[3:18] (positions 3-17,
    0-based) = 5 codons worth.

    Expected result: exactly 1 CdsFeature returned with name=="aaa",
    strand=="+", and exactly 5 codon_boundaries, each a 3-bp span fully
    within amplicon-local [0, 15).
    """
    seq = "ATG" * 10  # 30 bp
    features = textwrap.dedent("""\
             CDS             1..30
                             /gene="aaa"
                             /codon_start=1
                             /product="hypothetical"
    """)
    gb = _write_gb(tmp_path, "FWD", seq, features)

    amplicon = seq[3:18]  # 15 bp, 5 complete codons
    results = _parse_genbank_cds(gb, amplicon)

    assert len(results) == 1
    feat = results[0]
    assert feat.name == "aaa"
    assert feat.strand == "+"
    assert len(feat.codon_boundaries) == 5
    # All codon boundaries must be within amplicon local coords [0, 15)
    for start, end, ref_seq in feat.codon_boundaries:
        assert 0 <= start < 15
        assert 0 < end <= 15
        assert end - start == 3
        assert len(ref_seq) == 3


# ---------------------------------------------------------------------------
# Case 2: amplicon starts mid-codon — reading frame boundary skip
# ---------------------------------------------------------------------------
def test_frame_skip_mid_codon(tmp_path):
    """Purpose: verify that when the amplicon starts mid-codon (not aligned
    to the CDS reading frame), the leading partial codon is skipped.

    Function under test: _parse_genbank_cds -- frame-offset calculation that
    determines how many bases into a codon the amplicon start falls, and
    skips that many bases before reporting the first codon_boundary.

    Test input: same 30-bp "ATG"x10 genome with forward CDS 1..30
    (gene="bbb"). Amplicon = genome[1:16], which starts 1 nt into codon 1
    (frame_in_codon == 1).

    Expected result: the first codon_boundary starts at amplicon-local
    position 2 (a 2-nt skip accounts for the 1-nt frame offset).
    """
    seq = "ATG" * 10  # positions 0-29 (0-based)
    features = textwrap.dedent("""\
             CDS             1..30
                             /gene="bbb"
                             /codon_start=1
                             /product="hypothetical"
    """)
    gb = _write_gb(tmp_path, "MIDCODON", seq, features)

    # Start 1 nt into codon 1 (genome pos 1) → frame_in_codon = 1 → skip 2
    amplicon = seq[1:16]  # 15 bp starting mid-codon
    results = _parse_genbank_cds(gb, amplicon)

    assert len(results) == 1
    feat = results[0]
    # First complete codon should start at amplicon position 2 (after 2-nt skip)
    first_start = min(s for s, e, r in feat.codon_boundaries)
    assert first_start == 2


# ---------------------------------------------------------------------------
# Case 3: no CDS overlap — returns empty list
# ---------------------------------------------------------------------------
def test_no_cds_overlap(tmp_path):
    """Purpose: verify that an amplicon located outside any annotated CDS
    yields an empty result (not an error or a spurious match).

    Function under test: _parse_genbank_cds -- overlap-detection step, which
    excludes CDS features whose genomic range does not intersect the located
    amplicon coordinates.

    Test input: 80-bp non-repeating synthetic genome with one CDS at 1..30
    (gene="ccc"). Amplicon = genome[40:70], entirely outside the CDS.

    Expected result: results == [].
    """
    # Non-repeating 80 bp sequence (each 10 bp block distinct) so the
    # amplicon at [40:70] is only found at its true location, not via an
    # earlier identical-content match overlapping the CDS at [0:30].
    seq = "A" * 10 + "C" * 10 + "G" * 10 + "T" * 10 + "AC" * 5 + "GT" * 5 + "AG" * 5 + "CT" * 5  # 80 bp
    features = textwrap.dedent("""\
             CDS             1..30
                             /gene="ccc"
                             /codon_start=1
                             /product="hypothetical"
    """)
    gb = _write_gb(tmp_path, "NOOVERLAP", seq, features)

    amplicon = seq[40:70]  # completely outside CDS
    results = _parse_genbank_cds(gb, amplicon)

    assert results == []


# ---------------------------------------------------------------------------
# Case 4: amplicon not found in genome → empty list with no exception
# ---------------------------------------------------------------------------
def test_amplicon_not_in_genome(tmp_path):
    """Purpose: verify that an amplicon sequence which cannot be located
    anywhere in the genome returns an empty result without raising.

    Function under test: _parse_genbank_cds -- amplicon-location step (exact
    match, falling back to local alignment for short sequences), which finds
    no match and returns early.

    Test input: 32-bp all-"A" genome with CDS 1..30 (gene="ddd"). Amplicon =
    "CCCCCCCCCCCCCCCC" (16 C's), which does not occur anywhere in the genome.

    Expected result: results == [].
    """
    seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"  # 32 A's
    features = textwrap.dedent("""\
             CDS             1..30
                             /gene="ddd"
                             /codon_start=1
                             /product="hypothetical"
    """)
    gb = _write_gb(tmp_path, "NOTFOUND", seq, features)

    amplicon = "CCCCCCCCCCCCCCCC"  # 16 C's — not in the genome
    results = _parse_genbank_cds(gb, amplicon)

    assert results == []


# ---------------------------------------------------------------------------
# Case 5: reverse-complement amplicon is located correctly
# ---------------------------------------------------------------------------
def test_reverse_complement_amplicon(tmp_path):
    """Purpose: verify that an amplicon supplied as the reverse complement of
    a genomic region still resolves to the correct CDS and codon boundaries.

    Function under test: _parse_genbank_cds -- reverse-complement detection
    branch, which locates the RC match when the amplicon doesn't match the
    forward strand directly, tracks strand orientation, and computes codon
    boundaries the same way as the forward-strand case.

    Test input: same 30-bp "ATG"x10 genome with forward CDS 1..30
    (gene="eee"). Amplicon = reverse complement of genome[3:18] (the same
    5-codon region as test_forward_cds_basic, but RC'd).

    Expected result: exactly 1 CdsFeature returned with name=="eee" and
    exactly 5 codon_boundaries.
    """
    from skbio import DNA

    seq = "ATG" * 10  # 30 bp forward
    features = textwrap.dedent("""\
             CDS             1..30
                             /gene="eee"
                             /codon_start=1
                             /product="hypothetical"
    """)
    gb = _write_gb(tmp_path, "RCTEST", seq, features)

    # Amplicon is the reverse complement of genome positions 3-18
    fwd_sub = seq[3:18]
    rc_amplicon = str(DNA(fwd_sub).reverse_complement())
    results = _parse_genbank_cds(gb, rc_amplicon)

    # Should find the CDS (the overlap region is the same 5 codons)
    assert len(results) == 1
    assert results[0].name == "eee"
    assert len(results[0].codon_boundaries) == 5


# ---------------------------------------------------------------------------
# Case 6: multiple GenBank files — CDS found regardless of which file/order
# ---------------------------------------------------------------------------
def test_multi_file_genbank(tmp_path):
    """Purpose: verify _parse_genbank_cds accepts a list of GenBank files and
    finds the matching CDS regardless of which file contains it or the order
    the files are given in.

    Function under test: _parse_genbank_cds -- multi-file iteration, which
    searches each provided GenBank file (each loaded/cached via
    _load_genbank_records) for the amplicon until a match is found.

    Test input: two synthetic GenBank files -- "MULTI_OTHER" (32-bp all-A
    genome, CDS gene="zzz", does not contain the amplicon) and "MULTI_HIT"
    (30-bp "ATG"x10 genome, CDS gene="fff", amplicon = hit_seq[3:18]).
    _parse_genbank_cds is called with both file orderings:
    [other, hit] and [hit, other].

    Expected result: both orderings return exactly 1 CdsFeature with
    name=="fff" and 5 codon_boundaries.
    """
    # File with no overlapping CDS for the amplicon used below.
    other_seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"  # 32 A's
    other_features = textwrap.dedent("""\
             CDS             1..30
                             /gene="zzz"
                             /codon_start=1
                             /product="hypothetical"
    """)
    gb_other = _write_gb(tmp_path, "MULTI_OTHER", other_seq, other_features)

    # File containing the amplicon and its overlapping CDS.
    hit_seq = "ATG" * 10  # 30 bp
    hit_features = textwrap.dedent("""\
             CDS             1..30
                             /gene="fff"
                             /codon_start=1
                             /product="hypothetical"
    """)
    gb_hit = _write_gb(tmp_path, "MULTI_HIT", hit_seq, hit_features)

    amplicon = hit_seq[3:18]  # 5 complete codons

    for gb_files in ([gb_other, gb_hit], [gb_hit, gb_other]):
        results = _parse_genbank_cds(gb_files, amplicon)
        assert len(results) == 1
        assert results[0].name == "fff"
        assert len(results[0].codon_boundaries) == 5


# ---------------------------------------------------------------------------
# Case 7: GenBank records are cached across calls
# ---------------------------------------------------------------------------
def test_genbank_records_cached(tmp_path):
    """Purpose: verify repeated _parse_genbank_cds calls for the same file
    reuse the parsed-record cache instead of re-reading/re-parsing the file.

    Function under test: _load_genbank_records -- decorated with
    @functools.lru_cache, so the first call for a given path is a cache miss
    (parses the file) and subsequent calls with the same path are cache hits.

    Test input: a single synthetic GenBank file "CACHETEST" (30-bp "ATG"x10,
    CDS gene="ggg"). The cache is cleared via cache_clear(), then
    _parse_genbank_cds is called twice with the identical (gb, amplicon) args.

    Expected result: after the first call, cache_info() shows misses==1,
    hits==0; after the second call, misses==1 (unchanged), hits==1.
    """
    seq = "ATG" * 10
    features = textwrap.dedent("""\
             CDS             1..30
                             /gene="ggg"
                             /codon_start=1
                             /product="hypothetical"
    """)
    gb = _write_gb(tmp_path, "CACHETEST", seq, features)
    amplicon = seq[3:18]

    _load_genbank_records.cache_clear()

    _parse_genbank_cds(gb, amplicon)
    info = _load_genbank_records.cache_info()
    assert info.misses == 1
    assert info.hits == 0

    _parse_genbank_cds(gb, amplicon)
    info = _load_genbank_records.cache_info()
    assert info.misses == 1
    assert info.hits == 1


# ---------------------------------------------------------------------------
# Integration test: real H37Rv GenBank + rpoB amplicon excerpt
# ---------------------------------------------------------------------------
@pytest.mark.skipif(
    not os.path.exists(H37RV_GB),
    reason="H37Rv GenBank test fixture not available"
)
def test_rpob_in_h37rv(tmp_path):
    """Purpose: end-to-end sanity check that _parse_genbank_cds works against
    a real annotated genome and a real amplicon, finding the rpoB CDS and
    valid codon boundaries (skipped if the H37Rv fixture is unavailable).

    Function under test: _parse_genbank_cds -- the full pipeline (amplicon
    location, CDS-overlap detection, codon-boundary computation) exercised
    against real GenBank annotation data rather than synthetic minimal
    records.

    Test input: the real H37Rv GenBank file
    (../nextflow/tests/preparejson/H37Rv_NC0009623.gb) and an 80-bp excerpt of
    the rpoB amplicon sequence taken from TB.json.

    Expected result: results has >= 1 CdsFeature; at least one feature name
    contains "rpoB" or "Rv0667"; for that feature, every codon_boundary
    (start, end) lies within [0, len(amplicon)) and end-start == 3.
    """
    # Short excerpt from TB.json rpoB amplicon (known to overlap rpoB CDS)
    rpob_excerpt = "CCGAGCGGGGTGATGTCAACCCAGTGGGTGGCCTGGAAGAGGTGCTCTACGAGCTGTCTCCGATCGAGGACTTCTCCGGG"

    results = _parse_genbank_cds(H37RV_GB, rpob_excerpt)

    # Must find at least one CDS feature
    assert len(results) >= 1

    # The rpoB gene (Rv0667) should be among them
    names = [f.name for f in results]
    assert any("rpoB" in n or "Rv0667" in n for n in names), \
        f"Expected rpoB/Rv0667 in results, got: {names}"

    # All codon boundaries must be within the amplicon
    rpob_feat = next(f for f in results if "rpoB" in f.name or "Rv0667" in f.name)
    amp_len = len(rpob_excerpt)
    for start, end, ref_seq in rpob_feat.codon_boundaries:
        assert 0 <= start < amp_len, f"codon start {start} out of amplicon range"
        assert 0 < end <= amp_len, f"codon end {end} out of amplicon range"
        assert end - start == 3, f"codon span {end - start} != 3"
        assert len(ref_seq) == 3, f"ref_seq '{ref_seq}' is not 3 bp"


# ---------------------------------------------------------------------------
# Case 8: polyprotein product disambiguates duplicate gene names
# ---------------------------------------------------------------------------
def test_polyprotein_product_name_disambiguation(tmp_path):
    """Purpose: verify that when two CDS features share the same gene= value
    but differ in their product= qualifier (e.g. ORF1ab/ORF1a in SARS-CoV-2),
    the more specific name from product is used for each.

    Function under test: _parse_genbank_cds -- feat_name resolution using
    product qualifier when product ends in "polyprotein".

    Test input: 90-bp genome, two overlapping CDS features both with
    gene="ORF1ab" but product="ORF1ab polyprotein" and "ORF1a polyprotein".

    Expected result: two CdsFeatures named "ORF1ab" and "ORF1a" (not two
    identical "ORF1ab" entries).
    """
    seq = "ATG" * 30  # 90 bp
    features = textwrap.dedent("""\
             CDS             1..90
                             /gene="ORF1ab"
                             /codon_start=1
                             /product="ORF1ab polyprotein"
             CDS             1..60
                             /gene="ORF1ab"
                             /codon_start=1
                             /product="ORF1a polyprotein"
    """)
    gb = _write_gb(tmp_path, "SC2POLY", seq, features)

    amplicon = seq[0:30]  # 30 bp, 10 codons, overlaps both CDS features
    results = _parse_genbank_cds(gb, amplicon)

    names = [f.name for f in results]
    assert "ORF1ab" in names, f"Expected 'ORF1ab' in results, got: {names}"
    assert "ORF1a" in names, f"Expected 'ORF1a' in results, got: {names}"
    assert names.count("ORF1ab") == 1, f"Duplicate 'ORF1ab' entries: {names}"
    assert names.count("ORF1a") == 1, f"Duplicate 'ORF1a' entries: {names}"
