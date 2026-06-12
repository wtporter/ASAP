#!/usr/bin/env python3
# encoding: utf-8
'''
asap.genbank_cds -- GenBank CDS feature parsing and codon boundary derivation

@author:     TGen North
@copyright:  2025 TGen North. All rights reserved.
@license:    ACADEMIC AND RESEARCH LICENSE -- see ../LICENSE
'''

import os
import logging
import functools
from collections import namedtuple

import skbio.io
from skbio import DNA
from skbio.alignment import local_pairwise_align_nucleotide


CdsFeature = namedtuple('CdsFeature', ['name', 'strand', 'codon_boundaries'])
# codon_boundaries: list of (start, end, ref_seq) 0-based amplicon-local coords per codon


@functools.lru_cache(maxsize=None)
def _load_genbank_records(gb_file):
    """
    Read and cache all DNA records (with interval metadata) from a GenBank
    file. Cached because _parse_genbank_cds is called once per amplicon, and
    would otherwise re-read/re-parse the same file(s) for every amplicon.
    """
    return list(skbio.io.registry.read(gb_file, format='genbank', constructor=DNA))


def _parse_genbank_cds(gb_files, amplicon_sequence):
    """
    Extract CDS features from one or more GenBank files that overlap the
    amplicon and return codon boundaries in amplicon-local (0-based)
    coordinates.

    `gb_files` may be a single path or an iterable of paths; records from all
    files are searched. Locates the amplicon in each GenBank record via exact
    string match (fast even on full genomes). Falls back to
    local_pairwise_align_nucleotide only for records < 100 kb (avoids
    O(n*m) on whole-genome files).

    Returns a list of CdsFeature; returns [] with a warning when the amplicon
    cannot be located or no CDS features overlap.
    """
    if isinstance(gb_files, (str, bytes, os.PathLike)):
        gb_files = [gb_files]

    amp_seq = amplicon_sequence.upper().replace('\n', '').replace(' ', '')
    amp_len = len(amp_seq)
    results = []

    for gb_file in gb_files:
        for gb_seq in _load_genbank_records(gb_file):
            genome_str = str(gb_seq).upper()
            genome_len = len(genome_str)

            # --- locate the amplicon in the genome ---
            genome_offset = genome_str.find(amp_seq)
            amp_is_rc = False

            if genome_offset == -1:
                rc_seq = str(DNA(amp_seq).reverse_complement()).upper()
                rc_offset = genome_str.find(rc_seq)
                if rc_offset != -1:
                    # Amplicon is the reverse complement of a genome region.
                    # amplicon[0] == RC of genome[rc_offset + amp_len - 1]
                    genome_offset = rc_offset
                    amp_is_rc = True

            if genome_offset == -1:
                if genome_len > 100_000:
                    logging.warning(
                        "Amplicon not found by exact match in large GenBank record "
                        f"({gb_file}, {genome_len} bp); skipping CDS derivation"
                    )
                    continue
                try:
                    alignment, score, start_end = local_pairwise_align_nucleotide(
                        DNA(amp_seq), gb_seq
                    )
                    aligned_amp_len = start_end[0][1] - start_end[0][0]
                    if aligned_amp_len < amp_len * 0.9:
                        logging.warning(
                            f"Amplicon aligns at <90% to {gb_file}; skipping CDS derivation"
                        )
                        continue
                    genome_offset = start_end[1][0]
                except Exception as exc:
                    logging.warning(f"Alignment failed for {gb_file}: {exc}")
                    continue

            genome_amp_end = genome_offset + amp_len

            # --- iterate CDS features that overlap the amplicon window ---
            if not hasattr(gb_seq, 'interval_metadata'):
                continue

            for feature in gb_seq.interval_metadata.query(metadata={'type': 'CDS'}):
                gene = (feature.metadata.get('gene') or '').strip('"')
                product = (feature.metadata.get('product') or '').strip('"')
                locus_tag = (feature.metadata.get('locus_tag') or '').strip('"')
                # When product ends in "polyprotein" it carries a more specific
                # gene symbol than the parent `gene` qualifier (e.g. "ORF1a"
                # vs the parent "ORF1ab" that covers both ORF1a and ORF1ab).
                if product.lower().endswith('polyprotein'):
                    feat_name = product[:-len('polyprotein')].strip() or gene or locus_tag or 'unknown'
                else:
                    feat_name = gene or locus_tag or 'unknown'

                # skbio stores strand as int: 1 (forward) or -1 (reverse)
                raw_strand = feature.metadata.get('strand', 1)
                feat_strand = '+' if raw_strand in (1, '+', '1') else '-'

                # codon_start is 1-based in GenBank; convert to 0-based reading-frame offset
                codon_start_offset = int(feature.metadata.get('codon_start', 1)) - 1

                # Collect all genomic positions in CDS order (handles compound locations)
                cds_genome_positions = []
                for (feat_start, feat_end) in feature.bounds:
                    cds_genome_positions.extend(range(feat_start, feat_end))
                if feat_strand == '-':
                    cds_genome_positions = list(reversed(cds_genome_positions))

                # Find which CDS positions fall within the amplicon window
                overlapping_cds_indices = [
                    i for i, gp in enumerate(cds_genome_positions)
                    if genome_offset <= gp < genome_amp_end
                ]
                if len(overlapping_cds_indices) < 3:
                    continue  # need at least one full codon

                first_cds_idx = overlapping_cds_indices[0]
                # Frame: how many positions into a codon is the first overlapping base?
                frame_in_codon = (first_cds_idx + codon_start_offset) % 3
                # Skip to the next codon boundary
                skip = (3 - frame_in_codon) % 3
                coding_cds_indices = overlapping_cds_indices[skip:]

                # Convert to amplicon-local coords
                def _genome_to_amp(gp):
                    if amp_is_rc:
                        return amp_len - 1 - (gp - genome_offset)
                    return gp - genome_offset

                codon_boundaries = []
                for i in range(0, len(coding_cds_indices) - 2, 3):
                    triplet_cds = coding_cds_indices[i:i + 3]
                    if len(triplet_cds) < 3:
                        break
                    amp_positions = [_genome_to_amp(cds_genome_positions[ci])
                                     for ci in triplet_cds]
                    if feat_strand == '+':
                        ref_seq = ''.join(amp_seq[p] for p in sorted(amp_positions))
                    else:
                        ref_seq = str(DNA(
                            ''.join(amp_seq[p] for p in sorted(amp_positions))
                        ).reverse_complement())
                    codon_boundaries.append((min(amp_positions), max(amp_positions) + 1, ref_seq))

                if codon_boundaries:
                    results.append(CdsFeature(
                        name=feat_name,
                        strand=feat_strand,
                        codon_boundaries=codon_boundaries,
                    ))

    if not results:
        logging.debug(f"No overlapping CDS features found in {gb_files} for amplicon")
    return results
