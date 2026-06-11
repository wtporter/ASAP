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
    """
    30-bp genome = 10 ATG codons.  Amplicon = genome[3:18] (positions 3-17,
    0-based) = 5 codons worth.  Expect exactly 5 codon_boundaries all in [0,15).
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
    for start, end in feat.codon_boundaries:
        assert 0 <= start < 15
        assert 0 < end <= 15
        assert end - start == 3


# ---------------------------------------------------------------------------
# Case 2: amplicon starts mid-codon — reading frame boundary skip
# ---------------------------------------------------------------------------
def test_frame_skip_mid_codon(tmp_path):
    """
    Amplicon starts 1 nt into a codon.  The first complete codon in the
    amplicon should start at amplicon position 2 (skipping 1 frame offset).
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
    first_start = min(s for s, e in feat.codon_boundaries)
    assert first_start == 2


# ---------------------------------------------------------------------------
# Case 3: no CDS overlap — returns empty list
# ---------------------------------------------------------------------------
def test_no_cds_overlap(tmp_path):
    """Amplicon is from a region with no CDS features → empty list."""
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
    """When the amplicon sequence doesn't appear in the genome, return []."""
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
    """
    The amplicon is the RC of a genomic region that overlaps a CDS.
    _parse_genbank_cds should still find the CDS and return codon boundaries.
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
    """
    _parse_genbank_cds should accept a list of GenBank files and search all of
    them for the amplicon, regardless of which file contains the match or the
    order the files are given in.
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
    """
    Repeated _parse_genbank_cds calls for the same file should hit the
    _load_genbank_records cache instead of re-reading/re-parsing it.
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
    """
    Use a known rpoB amplicon excerpt (from TB.json) and verify that
    _parse_genbank_cds finds the rpoB CDS and returns codon boundaries
    within the amplicon window.
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
    for start, end in rpob_feat.codon_boundaries:
        assert 0 <= start < amp_len, f"codon start {start} out of amplicon range"
        assert 0 < end <= amp_len, f"codon end {end} out of amplicon range"
        assert end - start == 3, f"codon span {end - start} != 3"
