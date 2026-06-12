#!/usr/bin/env python3
# encoding: utf-8
'''
asap.allele_linkage -- Read-level allele linkage and codon correction

@author:     TGen North
@copyright:  2025 TGen North. All rights reserved.
@license:    ACADEMIC AND RESEARCH LICENSE -- see ../LICENSE
'''

from collections import Counter


def _build_fragment_allele_table(samdata, positions, ref_name=None):
    """
    Single pass over the BAM, recording the base each DNA fragment (read pair
    merged by query_name) carries at each requested 0-based amplicon-local
    position. Supplementary / secondary / unmapped alignments are ignored.

    Also computes, per requested position, the (min, max) reference span
    reached by any fragment covering that position -- the fragment's own
    extent widened to its mate via template_length, falling back to
    reference_length for single-end/long reads. Callers use this to skip
    allele-linkage checks for SNP pairs no fragment could ever co-cover.

    Fragments with an 'N' base at a requested position (e.g. primer-masked
    bases left aligned but uncalled) are recorded separately in `masked`
    rather than `pos_table`/`reach`, since the allele there is unknown. If a
    fragment's mates disagree (one mate masked, the other called) at the
    same position, the real call wins: `pos_table` and `masked` are kept
    disjoint per position, so each fragment falls into exactly one of them.

    Fragments with a deletion (gap) at a requested position are recorded in
    `deleted` rather than `pos_table`. These are fragments whose alignment
    spans the position (the reference base is consumed) but the read itself
    has a gap there. `pos_table`, `masked`, and `deleted` are kept disjoint
    per position; a real base call (from either mate) takes priority over a
    deletion call.

    ref_name: restrict fetch to this contig (pass the amplicon reference name).
    Returns (pos_table, reach, masked, deleted):
        pos_table: {pos: {qname: base, ...}, ...} for each position in `positions`
        reach:     {pos: (min_frag_start, max_frag_end), ...}
        masked:    {pos: {qname, ...}, ...} fragments with an 'N' at pos, and
                    no real call at pos from either mate
        deleted:   {pos: {qname, ...}, ...} fragments with a deletion at pos,
                    no real call at pos from either mate
    """
    pos_set = set(positions)
    pos_table = {p: {} for p in pos_set}
    masked = {p: set() for p in pos_set}
    deleted = {p: set() for p in pos_set}
    reach = {}

    if not pos_set:
        return pos_table, reach, masked, deleted

    fetch_iter = samdata.fetch(contig=ref_name) if ref_name else samdata.fetch()
    for read in fetch_iter:
        if read.is_supplementary or read.is_secondary or read.is_unmapped:
            continue
        other_end = (read.reference_start + read.template_length) if read.template_length else read.reference_end
        frag_start = min(read.reference_start, other_end)
        frag_end = max(read.reference_end, other_end)

        qname = read.query_name
        for qpos, rpos in read.get_aligned_pairs():
            if rpos not in pos_set:
                continue
            if qpos is None:
                # Deletion in read relative to reference at rpos.
                deleted[rpos].add(qname)
                continue
            base = read.query_sequence[qpos].upper()
            if base == 'N':
                # Masked base (e.g. primer trimming) -- the allele here is
                # unknown, so this fragment can't inform linkage at rpos.
                masked[rpos].add(qname)
                continue
            pos_table[rpos][qname] = base
            lo, hi = reach.get(rpos, (frag_start, frag_end))
            reach[rpos] = (min(lo, frag_start), max(hi, frag_end))

    # A fragment's mates can disagree on base/masked/deleted state (e.g. one
    # mate spans a deletion, the other a normal base). If either mate produced
    # a real base call, that wins; otherwise a deletion call wins over masked.
    for p in pos_set:
        called = set(pos_table[p].keys())
        masked[p] -= called
        deleted[p] -= called
        deleted[p] -= masked[p]

    return pos_table, reach, masked, deleted


def _tally_allele_linkage(pos_table, positions):
    """
    Given a position -> {qname: base} table (see _build_fragment_allele_table),
    count how often each allele combination co-occurs across `positions` on
    the same DNA fragment.

    Only fragments present at ALL requested positions are counted.
    Returns a Counter mapping frozenset{(pos, base), ...} -> fragment count.
    """
    positions = list(positions)
    tables = [pos_table.get(p, {}) for p in positions]
    if not all(tables):
        return Counter()

    # Iterate the smallest table and check membership in the others, so cost
    # scales with the shallowest position's depth rather than total fragments.
    smallest_idx = min(range(len(tables)), key=lambda i: len(tables[i]))
    smallest_pos, smallest = positions[smallest_idx], tables[smallest_idx]
    others = [(p, t) for i, (p, t) in enumerate(zip(positions, tables)) if i != smallest_idx]

    tally = Counter()
    for qname, base in smallest.items():
        alleles = {smallest_pos: base}
        for p, t in others:
            if qname not in t:
                break
            alleles[p] = t[qname]
        else:
            tally[frozenset(alleles.items())] += 1

    return tally


def _tally_codon_linkage(pos_table, deleted, positions):
    """
    Like _tally_allele_linkage but combines real base calls (pos_table) and
    deletion events (deleted) so that partial and full codon deletions are
    included in the tally. A read is counted only if it has a definite call
    — real base OR deletion — at EVERY requested position.

    Returns a Counter mapping frozenset{(pos, base_or_'_'), ...} -> count.
    Reads absent from both pos_table and deleted at any position are excluded
    (either masked/N or short read that doesn't reach that position).
    """
    combined = {}
    for p in positions:
        combined[p] = {}
        for qname, base in pos_table.get(p, {}).items():
            combined[p][qname] = base
        for qname in deleted.get(p, set()):
            combined[p][qname] = '_'

    if not all(combined[p] for p in positions):
        return Counter()

    smallest = min(positions, key=lambda p: len(combined[p]))
    others = [(p, combined[p]) for p in positions if p != smallest]

    tally = Counter()
    for qname, base in combined[smallest].items():
        alleles = {smallest: base}
        for p, t in others:
            if qname not in t:
                break
            alleles[p] = t[qname]
        else:
            tally[frozenset(alleles.items())] += 1
    return tally


def _amp_to_translated(amp_pos, offset):
    """Convert 0-based amplicon-local position to gene-relative (translated) position."""
    if offset < 0 and (amp_pos + 1) >= abs(offset):
        return amp_pos + offset + 2
    return amp_pos + offset + 1


def _translated_to_amp(trans_pos, offset):
    """Convert gene-relative (translated) position to 0-based amplicon-local position."""
    if offset < 0 and trans_pos >= 1:
        return trans_pos - offset - 2
    return trans_pos - offset - 1


def _snp_variant_freq(snp):
    """Return (freq, count) of the variant allele for a SNP dict entry."""
    depth = int(snp.get('depth', 0))
    if depth == 0:
        return 0.0, 0
    basecalls = snp.get('basecalls') or Counter()
    count = basecalls.get(snp['variant'], 0)
    return count / depth, count


def _apply_codon_correction(snp_list, pos_table, deleted, masked, cds_features, offset, error_threshold, min_reads):
    """
    Annotate pairs of SNPs that fall in the same codon (per GenBank CDS
    features) with full 3-bp codon-level read counts.

    For each codon containing exactly two polymorphic positions, both SNPs are
    annotated with a 'codon_merges' entry describing the observed codon
    distribution across all reads that fully span (or fully delete) all 3
    positions. Each CDS annotation produces its own entry even if two CDS
    features share identical codon boundaries (e.g. ORF1ab + ORF1a).

    Each codon_merges entry has:
        name                 - "{cds_name}_codon_{N}" (1-based codon index)
        region               - CDS gene name
        direction            - "forward" or "reverse"
        position             - "start-end" in translated (gene-relative) coords
        codon_depth          - reads with a definite call (base OR deletion)
                               at all 3 codon positions
        reference            - 3-bp reference sequence (reading-direction order)
        codon_call           - dominant observed 3-bp sequence
        codon_call_count     - count of the dominant codon
        codon_call_percentage - dominant_count / codon_depth * 100
        codon_distribution   - dict mapping every observed 3-bp sequence
                               (using '_' for per-position deletions) to its
                               read count, sorted by descending count
        excl_has_n           - reads with N at ≥1 codon position (masked)
        excl_no_span         - reads present at ≥1 position but absent from
                               ≥1 other position (short / soft-clipped)

    Codons with fewer than min_reads supporting the dominant combination are
    skipped. Mutates snp_list entries in place. Returns None.
    """
    if not cds_features or not snp_list:
        return

    snp_by_trans = {int(s['position']): s for s in snp_list}

    for cds in cds_features:
        direction = 'forward' if cds.strand == '+' else 'reverse'
        for codon_idx, (codon_amp_start, codon_amp_end, ref_seq) in enumerate(cds.codon_boundaries):
            codon_num = codon_idx + 1
            codon_amp_positions = list(range(codon_amp_start, codon_amp_end))

            codon_trans = {_amp_to_translated(p, offset) for p in codon_amp_positions}
            codon_snps = [snp_by_trans[tp] for tp in codon_trans if tp in snp_by_trans]
            if len(codon_snps) != 2:
                continue
            codon_snps.sort(key=lambda s: int(s['position']))

            tally = _tally_codon_linkage(pos_table, deleted, codon_amp_positions)
            if not tally:
                continue

            dominant_combo_count = max(tally.values())
            if dominant_combo_count < min_reads:
                continue

            codon_depth = sum(tally.values())

            codon_distribution = Counter()
            for combo, count in tally.items():
                allele_by_amp = dict(combo)
                codon_seq = ''.join(allele_by_amp[p] for p in codon_amp_positions)
                codon_distribution[codon_seq] += count

            dominant_seq = max(codon_distribution, key=lambda k: (codon_distribution[k], k))
            dominant_count = codon_distribution[dominant_seq]
            dominant_pct = round(dominant_count / max(codon_depth, 1) * 100, 1)

            # Exclusion counts
            any_definite = set().union(
                *[set(pos_table.get(p, {})) | deleted.get(p, set())
                  for p in codon_amp_positions]
            )
            excl_has_n = set().union(
                *[masked.get(p, set()) for p in codon_amp_positions]
            )
            excl_no_span_count = len(any_definite) - codon_depth
            excl_has_n_count = len(excl_has_n)

            pos_lo = _amp_to_translated(codon_amp_start, offset)
            pos_hi = _amp_to_translated(codon_amp_end - 1, offset)
            position_str = f"{min(pos_lo, pos_hi)}-{max(pos_lo, pos_hi)}"

            entry = {
                'name': f"{cds.name}_codon_{codon_num}",
                'region': cds.name,
                'direction': direction,
                'position': position_str,
                'codon_depth': codon_depth,
                'reference': ref_seq,
                'codon_call': dominant_seq,
                'codon_call_count': dominant_count,
                'codon_call_percentage': dominant_pct,
                'codon_distribution': dict(codon_distribution),
                'excl_has_n': excl_has_n_count,
                'excl_no_span': excl_no_span_count,
            }
            partner_names = {s['name'] for s in codon_snps}
            for snp in codon_snps:
                snp.setdefault('codon_merges', []).append(entry)
                # Track partner names separately for _apply_discover_roi exclusion
                snp.setdefault('codon_partner_names', set()).update(
                    partner_names - {snp['name']}
                )


def _apply_discover_roi(snp_list, pos_table, reach, masked, offset, min_perc, min_reads, min_snp_perc):
    """
    For each SNP in snp_list, find other SNPs from the same list that
    co-occur on the same reads and annotate them with a 'linked_snps' list.

    For each ordered pair (s_a, s_b), the count_a reads carrying s_a's variant
    are split into four buckets relative to s_b's position:
        linked          - fragment spans both positions and also carries
                          s_b's variant allele
        standalone      - fragment spans both positions but carries a
                          different (non-N) allele at s_b's position
        masked          - fragment spans both positions but s_b's position
                          is an uncalled 'N' (e.g. primer-masked)
        non-overlapping - fragment/read never reaches s_b's position

    `reach[amp_a]` gives the (min, max) reference span covered by any
    fragment touching amp_a; SNP B positions outside that span are skipped
    outright, since no fragment could ever cover both. `masked[amp_b]` is
    the set of fragment names with an 'N' at amp_b (see
    _build_fragment_allele_table).

    Each entry in 'linked_snps' is a dict with:
        target_name, shared_read_depth, co_count,
        linkage_pct, sample_frequency_pct,
        both_variants, anchor_only, target_only,
        neither_variant, masked, read_too_short

    `linkage_pct` is co_count / (co_count + anchor_only) -- of reads
    carrying s_a's variant that confidently reach s_b, the percentage also
    carrying s_b's variant. `sample_frequency_pct` is co_count /
    shared_read_depth -- of all reads spanning both positions (any allele),
    the percentage carrying both variants simultaneously.
    `shared_read_depth` counts all reads with confident calls at both
    positions regardless of allele.

    `min_snp_perc` is a pre-filter on `callable_snps`: a SNP's own variant
    frequency must be >= this value to be considered at all, whether as an
    anchor (s_a) or as a linked candidate (s_b). This keeps low-frequency
    "position of interest" SNPs (which are tracked regardless of the global
    --proportion threshold) from inflating the O(n^2) search below.

    SNP pairs already annotated via `codon_merges` (same-codon pairs handled
    by _apply_codon_correction) are excluded from `linked_snps` -- the
    `codon_merge` element already gives the richer reference/variant/discordant
    breakdown for those pairs, so `linked_snps` is reserved for cross-codon
    or otherwise unrelated co-occurrences.
    """
    if len(snp_list) < 2:
        return

    # Only consider SNPs with actual depth and a variant frequency at or
    # above min_snp_perc (skip dummy/no-coverage and low-frequency entries)
    callable_snps = [
        s for s in snp_list
        if int(s.get('depth', 0)) > 0 and _snp_variant_freq(s)[0] >= min_snp_perc
    ]
    if len(callable_snps) < 2:
        return

    for i, s_a in enumerate(callable_snps):
        amp_a = _translated_to_amp(int(s_a['position']), offset)
        var_a = s_a['variant']
        _, count_a = _snp_variant_freq(s_a)
        linked = []

        reach_a = reach.get(amp_a)
        if reach_a is None:
            continue
        lo, hi = reach_a

        # Fragments carrying s_a's variant allele at amp_a, used below to
        # count how many of them are masked ('N') at amp_b.
        qnames_a = {q for q, b in pos_table.get(amp_a, {}).items() if b == var_a}

        codon_partners = s_a.get('codon_partner_names', set())

        for j, s_b in enumerate(callable_snps):
            if i == j:
                continue
            if s_b['name'] in codon_partners:
                continue
            amp_b = _translated_to_amp(int(s_b['position']), offset)
            if amp_b < lo or amp_b > hi:
                continue
            var_b = s_b['variant']

            tally = _tally_allele_linkage(pos_table, [amp_a, amp_b])
            if not tally:
                continue

            linked_combo = frozenset([(amp_a, var_a), (amp_b, var_b)])
            linked_count = tally.get(linked_combo, 0)
            if linked_count < min_reads:
                continue

            linked_pct = linked_count / max(count_a, 1) * 100
            if linked_pct / 100 < min_perc:
                continue

            shared_read_depth = sum(tally.values())
            anchor_only = sum(
                v for combo, v in tally.items()
                if (amp_a, var_a) in combo and (amp_b, var_b) not in combo
            )
            target_only = sum(
                v for combo, v in tally.items()
                if (amp_b, var_b) in combo and (amp_a, var_a) not in combo
            )
            neither_variant = shared_read_depth - linked_count - anchor_only - target_only
            masked_count = len(qnames_a & masked.get(amp_b, set()))
            read_too_short = max(count_a - linked_count - anchor_only - masked_count, 0)

            # linkage_pct: of anchor-variant reads that confidently reach
            # the target position, what fraction also carry the target variant
            confident_count = linked_count + anchor_only
            linkage_pct = round(linked_count / max(confident_count, 1) * 100, 1)

            # sample_frequency_pct: frequency of the co-occurrence across
            # all reads that span both positions (any allele)
            sample_frequency_pct = round(linked_count / max(shared_read_depth, 1) * 100, 1)

            linked.append({
                'target_name': s_b['name'],
                'shared_read_depth': shared_read_depth,
                'co_count': linked_count,
                'linkage_pct': linkage_pct,
                'sample_frequency_pct': sample_frequency_pct,
                'both_variants': linked_count,
                'anchor_only': anchor_only,
                'target_only': target_only,
                'neither_variant': neither_variant,
                'masked': masked_count,
                'read_too_short': read_too_short,
            })

        if linked:
            s_a['linked_snps'] = linked
