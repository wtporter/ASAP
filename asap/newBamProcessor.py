#!/usr/bin/env python3
# encoding: utf-8
'''
asap.bamProcessor -- Process BAM alignment files with an AssayInfo JSON file and generate XML for the results

asap.bamProcessor

@author:     Darrin Lemmer

@copyright:  2015,2025 TGen North. All rights reserved.

@license:    ACADEMIC AND RESEARCH LICENSE -- see ../LICENSE

@contact:    dlemmer@tgen.org
'''


import sys
import os
import re
import argparse
import functools
import logging

import pysam
from collections import Counter, namedtuple
from xml.etree import ElementTree
import skbio.io
from skbio import DNA
from skbio.alignment import local_pairwise_align_nucleotide

from asap import assayInfo
from asap import __version__
# https://github.com/martinblech/xmltodict
import json
import xmltodict


__all__ = []
__updated__ = '2025-05-16'
__date__ = '2015-07-16'

DEBUG = 1
INFO = 0
TESTRUN = 0
PROFILE = 0
REMOVE_TEMP = True

low_level_cutoff = 0.01
high_level_cutoff = 0.50
proportion = 0.1

def pairwise(iterable):
    from itertools import tee
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)

def _write_parameters(node, data):
    for k, v in data.items():
        subnode = ElementTree.SubElement(node, k)
        subnode.text = str(v)
    return node

def _get_n_counts(pileup_iterator, amplicon_length):
    """
    Calculates and returns an array of N-base counts for each position
    in the amplicon, ensuring each unique alignment is only processed once.
    
    This function uses a hybrid approach to count N-bases in both aligned and 
    soft-clipped primer regions.
    
    Args:
        pileup_iterator: An iterator for the pileup data.
        amplicon_length: The length of the amplicon sequence.

    Returns:
        A list of integers representing the count of N-bases at each position.
    """
    n_read_array = [0] * amplicon_length
    processed_alignments = set()

    for pileupcolumn in pileup_iterator:
        for pileupread in pileupcolumn.pileups:
            try:
                alignment = pileupread.alignment
                alignment_id = (alignment.query_name, alignment.reference_start, alignment.query_length)

                if alignment_id in processed_alignments:
                    continue
                
                current_read_sequence = alignment.query_sequence.upper()
                
                # --- HYBRID APPROACH START ---
                
                # 1. Count Ns in soft-clipped/unaligned regions at the start
                first_aligned_pos_in_query = alignment.get_aligned_pairs()[0][0]
                if first_aligned_pos_in_query is not None and first_aligned_pos_in_query > 0:
                    for i in range(first_aligned_pos_in_query):
                        if current_read_sequence[i] == 'N':
                            ref_pos_with_n = alignment.reference_start - (first_aligned_pos_in_query - i)
                            if 0 <= ref_pos_with_n < amplicon_length:
                                n_read_array[ref_pos_with_n] += 1

                # 2. Count Ns in soft-clipped/unaligned regions at the end
                last_aligned_pos_in_query = alignment.get_aligned_pairs()[-1][0]
                read_length = len(current_read_sequence)
                if last_aligned_pos_in_query is not None and last_aligned_pos_in_query < read_length - 1:
                    for i in range(last_aligned_pos_in_query + 1, read_length):
                        if current_read_sequence[i] == 'N':
                            ref_pos_with_n = alignment.reference_end + (i - last_aligned_pos_in_query -1)
                            if 0 <= ref_pos_with_n < amplicon_length:
                                n_read_array[ref_pos_with_n] += 1
                                
                # 3. Count Ns in the aligned regions
                for query_pos, ref_pos in alignment.get_aligned_pairs():
                    if query_pos is not None and ref_pos is not None:
                        if current_read_sequence[query_pos] == 'N' and 0 <= ref_pos < amplicon_length:
                            n_read_array[ref_pos] += 1
                
                # --- HYBRID APPROACH END ---
                
                processed_alignments.add(alignment_id)

            except Exception:
                pass

    return n_read_array

def _process_pileup(pileup, amplicon, depth, proportion, mutdepth, offset, wholegenome, base_qual, con_prop, fill_gap_char, fill_del_char, n_read_array):
    global low_level_cutoff, high_level_cutoff
    pileup_dict = {}
    snp_dict = _create_snp_dict(amplicon)
    deletion_counter = Counter() #keep track of deletions by read name
    consensus_seq = ""
    if fill_gap_char != "false": # If the flag was actually passed (i.e., fill_gap_char is 'n' or a custom char)
        gapfilled_consensus_seq = "" # For calculating the actual sequence
    else:
        # If the flag was NOT passed, initialize to the final string we want in the XML
        gapfilled_consensus_seq = "fill_gaps was not provided... consider using --fill_gaps n"
    snp_list = []
    breadth_positions = 0
    avg_depth_total = avg_depth_positions = 0
    amplicon_length = len(amplicon.sequence)
    depth_array = [0] * amplicon_length
    quality_discard_array = [0] * amplicon_length
    prop_array = ["0"] * amplicon_length
    previous_position = 0
    # for each position in alignment/pileup
    for pileupcolumn in pileup:
        base_counter = Counter()
        position = pileupcolumn.pos+1
        # This fills gaps in the alignment with n's or user defined char
        if fill_gap_char != "false":
            if previous_position+1 < position: #We've skipped some positions in the alignment
                #print("%i, %i" % (previous_position, position))
                for i in range(previous_position+1, position):
                    gapfilled_consensus_seq += fill_gap_char #Fill in the gap
        previous_position = position
        depth_array[pileupcolumn.pos] = pileupcolumn.n
        depth_passed = False
        passed_Qual_filter = 0
        for pileupread in pileupcolumn.pileups:
            #print("processing read, qual=%i" % pileupread.alignment.query_qualities[pileupread.query_position])
            try:
                if pileupread.is_del:
                    #This position in the alignment is a deletion in the query sequence, therefore it has no quality score
                    # Let's use the average of the quality scores of the two aligned bases flanking the deletion
                    qscore = (pileupread.alignment.query_qualities[pileupread.query_position_or_next] +
                              pileupread.alignment.query_qualities[pileupread.query_position_or_next - 1]) / 2
                    if qscore >= base_qual:
                        passed_Qual_filter += 1
                        base_counter.update({"_" : 1})
                    else:
                        quality_discard_array[pileupcolumn.pos] += 1
                elif pileupread.alignment.query_qualities[pileupread.query_position] >= base_qual: # check here
                    passed_Qual_filter += 1
                    if pileupread.indel < 0: #This means the next position is a deletion, we'll process later
                        for d in range(1, abs(pileupread.indel)+1):
                            deletion_counter.update({str(position + d)})
                    if pileupread.indel > 0: #This means the next position is an insertion, unlike with deletions, this we can process now
                        start = pileupread.query_position
                        end = pileupread.query_position + pileupread.indel + 1
                        base_counter.update({pileupread.alignment.query_sequence[start:end]: 1})
                    else:
                        base_counter.update(pileupread.alignment.query_sequence[pileupread.query_position])
                else:
                    quality_discard_array[pileupcolumn.pos] += 1
            except Exception as e:
                if str(e.__class__.__name__) != "TypeError":
                    print("Unexpected error:", sys.exc_info()[0])
                    pass
                quality_discard_array[pileupcolumn.pos] += 1 #check here
                pass

        column_depth = passed_Qual_filter #check this, this will count bases that have been filtered out by quality?
        depth_array[pileupcolumn.pos] = passed_Qual_filter #reset to depth that passed qual filter

        if column_depth > 0: #TODO: This is going to end up being specific to these TB assays (with flanking sequence), maybe have a clever way to make this line optional
            avg_depth_positions += 1
            avg_depth_total += column_depth
        if column_depth >= depth:
            breadth_positions += 1
            depth_passed = True
        ordered_list = base_counter.most_common()
        if not ordered_list: #No coverage, should only happen here if all reads were thrown out because of quality
            consensus_seq += "N"
            if fill_gap_char != "false":
                gapfilled_consensus_seq += "N"
            continue
        alignment_call = ordered_list[0][0]
        alignment_call_proportion = ordered_list[0][1] / column_depth
        prop_array[pileupcolumn.pos] = "%.3f" % alignment_call_proportion
        reference_call = amplicon.sequence[pileupcolumn.pos]
        if reference_call == '-':
            reference_call = '_' #Need to use '_' instead of '-' for gaps because of XSLT
        #if alignment_call != reference_call:
        #    snp_call = alignment_call
        #    snp_count = ordered_list[0][1]
        #    snp_call_proportion = alignment_call_proportion
        #elif len(ordered_list) > 1:
        #    snp_call = ordered_list[1][0]
        #    snp_count = ordered_list[1][1]
        #    snp_call_proportion = ordered_list[1][1] / column_depth
        # Initialize SNP variables TP added 2025 ########
        snp_call = None
        snp_count = None
        snp_call_proportion = None
        # Find the first valid SNP candidate (a non-reference, non-ambiguous base)
        for base, count in ordered_list:
            if base != reference_call and base.upper() != "N":
                snp_call = base
                snp_count = count
                snp_call_proportion = count / column_depth
                break # Exit the loop once a valid SNP is found

        #Generate consensus call at this pos
        #consensus_seq += alignment_call if alignment_call_proportion >= consensus_proportion else "N"
        # unless the alignment_call is a deletion, and > consensus proportion (Previously 50%) -- don't ever replace deletions with Ns
        # or if coverage is less than the depth threshold, then always call N
        if not depth_passed: # N's if we don't have enough coverage
            consensus_seq += "N"
            if fill_gap_char != "false":
                gapfilled_consensus_seq += "N"
        elif alignment_call != "_":
            if alignment_call_proportion >= con_prop: #Add call to consensus and gap_filled...
                consensus_seq += alignment_call
                if fill_gap_char != "false":
                    gapfilled_consensus_seq += alignment_call
            else: #Consensus proportion not high enough
                consensus_seq += "N"
                if fill_gap_char != "false":
                    gapfilled_consensus_seq += "N"
        else:
            if alignment_call_proportion < con_prop: #TP Changed from <= 0.5 to be less than con_prop, this matches above calling.
                consensus_seq += "N"
                if fill_gap_char != "false":
                    gapfilled_consensus_seq += "N"
            else:
                if fill_del_char != "false": #Put in gaps if user requested them
                    consensus_seq += fill_del_char
                    if fill_gap_char != "false":
                        gapfilled_consensus_seq += fill_del_char

        if position >= abs(offset) and offset < 0: #if the offset is negative, ie. amplicon starts before beginning of the gene, then when converting to gene-based coordinates need to make offset 1 unit more positive to account for there being no 0-base in gene-coordinates
            translated = position + (offset + 1)
        else:
            translated = position + offset #normal case where gene encompasses the amplicon
        if position in snp_dict:
            for (name, reference, variant, significance) in snp_dict[position]:
                snp = {'name':name, 'position':str(translated), 'depth':str(column_depth), 'reference':reference, 'variant':variant, 'basecalls':base_counter}
                variant_proportion = base_counter[variant]/column_depth
                variant_count = base_counter[variant]
                if variant_proportion >= proportion and variant_count >= mutdepth:
                    snp['significance'] = significance
                    if variant_proportion <= low_level_cutoff:
                        snp['level'] = "low"
                    elif variant_proportion >= high_level_cutoff:
                        snp['level'] = "high"
                if not depth_passed:
                    snp['flag'] = "low coverage"
                snp_list.append(snp)
            # We've covered it, now remove it from the dict so we can see what we might have missed
            del snp_dict[position]
        elif depth_passed and snp_call and snp_count >= mutdepth and snp_call_proportion >= proportion:
            snp = {'name':f"{reference_call}{translated}{snp_call}", 'position':str(translated), 'depth':str(column_depth), 'reference':reference_call, 'variant':snp_call, 'basecalls':base_counter}
            if 0 in snp_dict:
                (name, *rest, significance) = snp_dict[0][0]
                snp['name'] = name
                snp['significance'] = significance
            snp_list.append(snp)
    #Check for any positions_of_interest that weren't covered
    snp_dict.pop(0, None)
    for position in snp_dict.keys():
        for (name, reference, variant, significance) in snp_dict[position]:
            snp = {'name':name, 'position':str(position), 'depth':str(0), 'reference':reference, 'variant':variant}
            snp_list.append(snp)
    if not wholegenome: #If reference is whole genome, none of these are going to make sense, and they will make the output too large
        pileup_dict['consensus_sequence'] = consensus_seq
        pileup_dict['gapfilled_consensus_sequence'] = gapfilled_consensus_seq #TP added
        # if fill_gap_char: #TP removed...
        #     pileup_dict['gapfilled_consensus_sequence'] = gapfilled_consensus_seq
        pileup_dict['depths'] = ",".join(str(n) for n in depth_array)
        pileup_dict['proportions'] = ",".join(prop_array)
        pileup_dict['n_reads'] = ",".join(str(n) for n in n_read_array)
    pileup_dict['breadth'] = str(breadth_positions/amplicon_length * 100)
    pileup_dict['quality_discards'] = ",".join(str(n) for n in quality_discard_array)
    pileup_dict['SNPs'] = snp_list
    pileup_dict['average_depth'] = str(avg_depth_total/avg_depth_positions) if avg_depth_positions else "0"
    return pileup_dict

def _create_snp_dict(amplicon):
    snp_dict = {}
    for snp in amplicon.SNPs:
        name = snp.name if snp.name else "position of interest"
        if not snp.variant or snp.variant == "any":
            for v in {'A', 'G', 'C', 'T'}:
                if v != snp.reference:
                    if snp.position in snp_dict:
                        snp_dict[snp.position].append((name, snp.reference, v, snp.significance))
                    else:
                        snp_dict[snp.position] = [(name, snp.reference, v, snp.significance)]
        else:
            if snp.position in snp_dict:
                snp_dict[snp.position].append((name, snp.reference, snp.variant, snp.significance))
            else:
                snp_dict[snp.position] = [(name, snp.reference, snp.variant, snp.significance)]
    return snp_dict

def _add_snp_node(parent, snp):
    snp_attributes = {k:snp[k] for k in ('name', 'position', 'depth', 'reference')}
    snp_node = ElementTree.SubElement(parent, 'snp', snp_attributes)
    base_counter = snp.get('basecalls')
    snpcall = snp['variant']
    depth = int(snp['depth'])
    snpcount = base_counter[snpcall] if base_counter else 0
    percent = snpcount/depth*100 if depth else 0
    snpcall_node = ElementTree.SubElement(snp_node, 'snp_call', {'count':str(snpcount), 'percent':str(percent)})
    snpcall_node.text = snpcall
    if 'significance' in snp or 'flag' in snp:
        significance_node = ElementTree.SubElement(snp_node, 'significance')
        if 'significance' in snp:
            significance_node.text = snp['significance'].message
            if snp['significance'].resistance:
                significance_node.set("resistance", snp['significance'].resistance)
            if 'level' in snp:
                significance_node.set("level", snp['level'])
        if 'flag' in snp:
            significance_node.set('flag', snp['flag'])
    if base_counter:
        ElementTree.SubElement(snp_node, 'base_distribution', {k:str(v) for k,v in base_counter.items()})
    return snp_node


# ---------------------------------------------------------------------------
# GenBank CDS / coordinate layer
# ---------------------------------------------------------------------------

CdsFeature = namedtuple('CdsFeature', ['name', 'strand', 'codon_boundaries'])
# codon_boundaries: list of (start, end) 0-based amplicon-local coords per codon


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
                feat_name = (feature.metadata.get('gene') or
                             feature.metadata.get('locus_tag') or
                             'unknown').strip('"')

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
                    codon_boundaries.append((min(amp_positions), max(amp_positions) + 1))

                if codon_boundaries:
                    results.append(CdsFeature(
                        name=feat_name,
                        strand=feat_strand,
                        codon_boundaries=codon_boundaries,
                    ))

    if not results:
        logging.debug(f"No overlapping CDS features found in {gb_files} for amplicon")
    return results


# ---------------------------------------------------------------------------
# Read-level allele linkage helper
# ---------------------------------------------------------------------------

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


def _apply_codon_correction(snp_list, pos_table, deleted, cds_features, offset, error_threshold, min_reads):
    """
    Annotate pairs of SNPs that fall in the same codon (per GenBank CDS
    features) with read-level allele-linkage information.

    For each codon containing exactly two polymorphic positions (codons with
    three or more polymorphic positions are left unannotated), each SNP is
    annotated with a 'codon_merges' list entry describing how its variant
    allele co-occurs on the same DNA fragments with the other SNP's allele.
    Both SNP entries are kept unchanged in snp_list -- nothing is merged,
    replaced, or removed.

    Each codon_merges entry has:
        linked_snp      - name of the other SNP in the codon
        spanning_depth  - total fragments with confident calls at both positions
        linkage         - "complete" if abs(freq_a - freq_b) <= error_threshold,
                          else "partial" (carried over from the old
                          complete-merge / codon_link distinction)
        snp_percentage_linked        - linked_count / (spanning_depth +
                                        both_deleted) * 100.  Of all reads
                                        with a definite allele-pair
                                        determination at both positions (either
                                        both called bases, or both deleted),
                                        the percentage showing the linked
                                        variant pair.  Same value on both
                                        SNPs' entries for a pair.
        total_percentage_depth_linked - linked_count / this SNP's own depth
                                        * 100.  Can differ slightly between
                                        the two SNPs in a pair.
        combos          - list of {bases, count, percent, type}, one per
                          observed allele combination (from
                          _tally_allele_linkage), ordered by descending count.
                          'bases' is "X|Y" with X/Y in ascending-translated-
                          position order, so both SNPs in a pair report
                          combos in the same order. 'type' classifies the
                          combo as "reference" (matches both SNPs'
                          reference alleles), "variant" (matches both SNPs'
                          called variant alleles), or "discordant" (anything
                          else, e.g. a fragment carrying only one of the two
                          variants).
                          linked_count = count of the "variant" combo from
                          tally + both_deleted (fragments with a deletion at
                          BOTH positions, from the `deleted` table).  Exactly
                          one of these is non-zero per pair: tally can never
                          yield a '_|_' combo (deleted positions have no
                          callable base in pos_table), so deletion-variant
                          pairs use both_deleted exclusively.

    A SNP can accumulate multiple codon_merges entries if it participates in
    more than one qualifying codon (e.g. overlapping CDS features in
    different reading frames).

    min_reads gates whether a codon is annotated at all: if the most common
    allele combination spanning both positions has fewer than min_reads
    supporting fragments, the codon is skipped.

    GenBank annotation sets (e.g. SARS-CoV-2 ORF1ab + ORF1a) can contain
    multiple CDS features whose codon_boundaries describe the SAME reading
    frame over an overlapping region; codon boundary tuples are deduplicated
    across all cds_features before processing so each distinct (start, end)
    codon is handled once.

    Mutates snp_list entries in place. Returns None.
    """
    if not cds_features or not snp_list:
        return

    snp_by_trans = {int(s['position']): s for s in snp_list}

    all_codon_boundaries = {
        boundary
        for cds in cds_features
        for boundary in cds.codon_boundaries
    }

    for codon_amp_start, codon_amp_end in all_codon_boundaries:
        codon_trans = {_amp_to_translated(p, offset)
                      for p in range(codon_amp_start, codon_amp_end)}
        codon_snps = [snp_by_trans[tp] for tp in codon_trans if tp in snp_by_trans]
        if len(codon_snps) != 2:
            continue

        codon_snps.sort(key=lambda s: int(s['position']))
        s_a, s_b = codon_snps

        amp_a = _translated_to_amp(int(s_a['position']), offset)
        amp_b = _translated_to_amp(int(s_b['position']), offset)
        tally = _tally_allele_linkage(pos_table, [amp_a, amp_b])
        if not tally:
            continue

        spanning_depth = sum(tally.values())
        dominant_combo_count = max(tally.values())
        if dominant_combo_count < min_reads:
            continue

        freq_a, _ = _snp_variant_freq(s_a)
        freq_b, _ = _snp_variant_freq(s_b)
        linkage = "complete" if abs(freq_a - freq_b) <= error_threshold else "partial"

        ref_pair = (s_a['reference'], s_b['reference'])
        var_pair = (s_a['variant'], s_b['variant'])

        combos = []
        for combo, count in tally.items():
            allele_by_amp = dict(combo)
            base_a, base_b = allele_by_amp[amp_a], allele_by_amp[amp_b]
            if (base_a, base_b) == ref_pair:
                combo_type = "reference"
            elif (base_a, base_b) == var_pair:
                combo_type = "variant"
            else:
                combo_type = "discordant"
            combos.append({
                'bases': f"{base_a}|{base_b}",
                'count': count,
                'percent': count / spanning_depth * 100,
                'type': combo_type,
            })
        combos.sort(key=lambda c: (-c['count'], c['bases']))

        linked_combo_count = next((c['count'] for c in combos if c['type'] == 'variant'), 0)
        both_deleted = len(deleted.get(amp_a, set()) & deleted.get(amp_b, set()))
        linked_count = linked_combo_count + both_deleted
        total_valid = spanning_depth + both_deleted
        snp_percentage_linked = linked_count / max(total_valid, 1) * 100
        depth_a = int(s_a.get('depth', 0))
        depth_b = int(s_b.get('depth', 0))
        total_pct_depth_linked_a = linked_count / max(depth_a, 1) * 100
        total_pct_depth_linked_b = linked_count / max(depth_b, 1) * 100

        s_a.setdefault('codon_merges', []).append({
            'linked_snp': s_b['name'],
            'spanning_depth': spanning_depth,
            'linkage': linkage,
            'snp_percentage_linked': snp_percentage_linked,
            'total_percentage_depth_linked': total_pct_depth_linked_a,
            'combos': combos,
        })
        s_b.setdefault('codon_merges', []).append({
            'linked_snp': s_a['name'],
            'spanning_depth': spanning_depth,
            'linkage': linkage,
            'snp_percentage_linked': snp_percentage_linked,
            'total_percentage_depth_linked': total_pct_depth_linked_b,
            'combos': combos,
        })


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
        name, ref_pos_allele, linked_pct, linked_count,
        standalone_pct, standalone_count, masked_pct, masked_count,
        nonoverlap_pct, nonoverlap_count, percentage_linked,
        comparable_count, spanning_depth, snp_depth

    `linked_pct` (== "snp_percentage_linked" in the XML) is linked_count /
    count_a -- of ALL of s_a's variant-supporting reads, the percentage also
    carrying s_b's variant. `percentage_linked` is linked_count /
    comparable_count -- of s_a's variant-supporting reads with a confident
    call at s_b's position (comparable_count = linked_count +
    standalone_count), the percentage carrying s_b's variant.
    `spanning_depth` is broader still: the total reads with confident calls
    at both positions regardless of allele at s_a's position.

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

        codon_partners = {cm['linked_snp'] for cm in s_a.get('codon_merges', [])}

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

            spanning_depth = sum(tally.values())
            standalone_count = sum(
                v for combo, v in tally.items()
                if combo != linked_combo and (amp_a, var_a) in combo
            )
            masked_count = len(qnames_a & masked.get(amp_b, set()))
            nonoverlap_count = max(count_a - linked_count - standalone_count - masked_count, 0)

            standalone_pct = standalone_count / max(count_a, 1) * 100
            masked_pct = masked_count / max(count_a, 1) * 100
            nonoverlap_pct = nonoverlap_count / max(count_a, 1) * 100

            # Of s_a's reads with a confident call at s_b's position (i.e.
            # excluding non-overlapping and masked reads), the percentage
            # that carry s_b's linked allele. A perfect link is 100%
            # regardless of how many reads don't overlap or are masked.
            confident_count = linked_count + standalone_count
            percentage_linked = linked_count / max(confident_count, 1) * 100

            linked.append({
                'name': s_b['name'],
                'ref_pos_allele': f"{s_b['reference']}{s_b['position']}",
                'linked_pct': round(linked_pct, 1),
                'linked_count': linked_count,
                'standalone_pct': round(standalone_pct, 1),
                'standalone_count': standalone_count,
                'masked_pct': round(masked_pct, 1),
                'masked_count': masked_count,
                'nonoverlap_pct': round(nonoverlap_pct, 1),
                'nonoverlap_count': nonoverlap_count,
                'percentage_linked': round(percentage_linked, 1),
                'comparable_count': confident_count,
                'spanning_depth': spanning_depth,
                'snp_depth': int(s_a.get('depth', 0)),
            })

        if linked:
            s_a['linked_snps'] = linked


def _add_linked_snps_node(snp_node, linked_snps):
    """Add a <linked_snps> sub-element to a SNP XML node."""
    ls_node = ElementTree.SubElement(snp_node, 'linked_snps')
    for entry in linked_snps:
        link_node = ElementTree.SubElement(ls_node, 'linked_snp', {
            'name': entry['name'],
            'ref_pos': entry['ref_pos_allele'],
            'spanning_depth': str(entry['spanning_depth']),
            'linked_depth': str(entry['linked_count']),
            'comparable_depth': str(entry['comparable_count']),
            'percentage_linked': str(entry['percentage_linked']),
            'snp_percentage_linked': str(entry['linked_pct']),
        })
        full_dist_node = ElementTree.SubElement(link_node, 'full_distribution')
        full_dist_node.text = (
            f"linked {entry['linked_pct']}% (N={entry['linked_count']}), "
            f"standalone {entry['standalone_pct']}% (N={entry['standalone_count']}), "
            f"masked {entry['masked_pct']}% (N={entry['masked_count']}), "
            f"non-overlapping {entry['nonoverlap_pct']}% (N={entry['nonoverlap_count']})"
        )


def _add_codon_merges_node(snp_node, codon_merges):
    """Add one <codon_merge> child element per entry in codon_merges."""
    for entry in codon_merges:
        cm_node = ElementTree.SubElement(snp_node, 'codon_merge', {
            'linked_snp': entry['linked_snp'],
            'spanning_depth': str(entry['spanning_depth']),
            'linkage': entry['linkage'],
            'snp_percentage_linked': f"{entry['snp_percentage_linked']:.1f}",
            'total_percentage_depth_linked': f"{entry['total_percentage_depth_linked']:.1f}",
        })
        for combo in entry['combos']:
            ElementTree.SubElement(cm_node, 'combo', {
                'bases': combo['bases'],
                'count': str(combo['count']),
                'percent': f"{combo['percent']:.1f}",
                'type': combo['type'],
            })


class CLIError(Exception):
    '''Generic exception to raise and log different fatal errors.'''
    def __init__(self, msg):
        super(CLIError).__init__(type(self))
        self.msg = "E: %s" % msg
    def __str__(self):
        return self.msg
    def __unicode__(self):
        return self.msg


def evaluateItem(item, sample_node):
    foundNode = findItem(item, sample_node)
    if not foundNode:
        print("Could not find node for the item in sample node, valued as False")
        return False
        pass
    print("Found Node: ", foundNode)

    if item.evaluation == "depth_greater_than":
        if "depth" in foundNode.attrib:
            return float(foundNode.attrib["depth"]) > float(item.value)
            pass
        else:
            return False
        pass
    if item.evaluation == "depth_less_than":
        if "depth" in foundNode.attrib:
            return float(foundNode.attrib["depth"]) < float(item.value)
            pass
        else:
            return False
        pass

    if item.item_type == "snp":
        if "percentage" in item.evaluation or "count" in item.evaluation:
            for call in foundNode:
                if call.tag == "snp_call":
                    foundNode = call
                    print("Looking at snp_call: ", foundNode)
                    break # This should ensure the topmost snp_call
                    pass

    if item.evaluation == "percentage_greater_than":
        if "percentage" in foundNode.attrib:
            return float(foundNode.attrib["percentage"]) > float(item.value)
            pass
        else:
            return False
        pass
    if item.evaluation == "percentage_less_than":
        if "percentage" in foundNode.attrib:
            return float(foundNode.attrib["percentage"]) < float(item.value)
            pass
        else:
            return False
        pass

    if item.evaluation == "count_greater_than":
        if "count" in foundNode.attrib:
            return float(foundNode.attrib["count"]) > float(item.value)
            pass
        else:
            return False
        pass
    if item.evaluation == "count_less_than":
        if "count" in foundNode.attrib:
            return float(foundNode.attrib["count"]) < float(item.value)
            pass
        else:
            return False
        pass
    # Evaluation empty is simply true
    if not item.evaluation or item.evaluation == "" or str(item.evaluation).lower() == "none":
        return True
        pass
    print("Item contains invalid evaluation: ", item.evaluation)
    return False

def findItem(item, sample_node):
    foundNode = None
    for childEle in sample_node:
        if childEle.tag == "assay":
            if item.item_type.lower() == "assay":
                if item.identity_key not in childEle.attrib:
                    continue
                    pass
                if str(childEle.attrib[item.identity_key]).lower() == item.identity_value.lower():
                    return childEle
                    pass
                pass
        for possibleAmpEle in childEle:
            if possibleAmpEle.tag == "amplicon":
                if item.item_type.lower() == "amplicon":
                    if item.identity_key not in possibleAmpEle.attrib:
                        continue
                        pass
                    if str(possibleAmpEle.attrib[item.identity_key]).lower() == item.identity_value.lower():
                        return possibleAmpEle
                        pass
                    pass
                for possibleSubEle in possibleAmpEle:
                    if possibleSubEle.tag == "roi":
                        if item.item_type.lower() == "roi":
                            if item.identity_key not in possibleSubEle.attrib:
                                continue
                                pass
                            if str(possibleSubEle.attrib[item.identity_key]).lower() == item.identity_value.lower():
                                foundNode = possibleSubEle
                                pass
                            pass
                        if item.item_type.lower() == "mutation":
                            for possibleMutEle in possibleSubEle:
                                if possibleMutEle.tag == "mutation":
                                    if item.identity_key == "base" or item.identity_key == "bases" or item.identity_key == "text" or item.identity_key == "mutation":
                                        if str(possibleMutEle.text).lower() == item.identity_value.lower():
                                            foundNode = possibleMutEle
                                            pass
                                        pass
                                    if item.identity_key not in possibleMutEle.attrib:
                                        continue
                                        pass
                                    if str(possibleMutEle.attrib[item.identity_key]).lower() == item.identity_value.lower():
                                        foundNode = possibleMutEle
                                        pass
                                    pass
                                pass
                            pass
                        pass
                    if possibleSubEle.tag == "snp":
                        if item.item_type.lower() == "snp":
                            if item.identity_key not in possibleSubEle.attrib:
                                continue
                                pass
                            if str(possibleSubEle.attrib[item.identity_key]).lower() == item.identity_value.lower():
                                foundNode = possibleSubEle
                                pass
                            pass
                        pass
                    pass
                pass
            pass
        pass
    #End of findItem func
    return foundNode

def evaluateOperation(operation, sample_node):
    print("Evaluating Operation: ", operation)
    print("On sample: ", sample_node)
    boolList = []
    for child in operation.children:
        if isinstance(child, assayInfo.ITEM):
            print("Evaluating Truthness of Item: ", child)
            boolList.append(evaluateItem(child, sample_node))
            pass
        if isinstance(child, assayInfo.Operation):
            print("Nested Operation: ", child)
            boolList.append(evaluateOperation(child, sample_node))
            pass
        pass
    print("Evaled Bools: ", boolList)
    if operation.operation_type == "AND":
        endBool = True
        for evaldBool in boolList:
            endBool = evaldBool and endBool
            pass
        return endBool
        pass
    if operation.operation_type == "OR":
        endBool = False
        for evaldBool in boolList:
            endBool = evaldBool or endBool
            pass
        return endBool
        pass
    if operation.operation_type == "NOT":
        # Should always have only 1 child
        return not boolList[0]
        pass

def main(argv=None): # IGNORE:C0111
    '''Command line options.'''

    global proportion

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    if __name__ == '__main__':
        program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    else:
        program_shortdesc = __doc__.split("\n")[1]
    #program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_license = '''%s

  Created by TGen North on %s.
  Copyright 2015 TGen North. All rights reserved.

  Available for academic and research use only under a license
  from The Translational Genomics Research Institute (TGen)
  that is free for non-commercial use.

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
''' % (program_shortdesc, str(__date__))

    try:
        # Setup argument parser
        parser = argparse.ArgumentParser(prog='TGen-ASAP', description=program_license, formatter_class=argparse.RawTextHelpFormatter)
        required_group = parser.add_argument_group("required arguments")
        required_group.add_argument("-j", "--json", metavar="FILE", required=True, type=argparse.FileType('r'), help="JSON file of assay descriptions. [REQUIRED]")
        required_group.add_argument("-b", "--bam", metavar="FILE", required=True, type=argparse.FileType('rb'), default=sys.stdin, help="BAM file to analyze. [REQUIRED]")
        parser.add_argument("-d", "--depth", default=100, type=int, help="minimum read depth required to consider a position covered. [default: 100]")
        parser.add_argument("--breadth", default=0.8, type=float, help="minimum breadth of coverage required to consider an amplicon as present. [default: 0.8]")
        parser.add_argument("-p", "--proportion", type=float, help="minimum proportion required to call a mutation at a given locus. [default: 0.1]") #Don't explicitly set default because I need to be certain whether user set the value
        parser.add_argument("-m", "--mutation-depth", dest="mutdepth", default=5, type=int, help="minimum number of reads required to call a mutation at a given locus. [default: 5]")
        parser.add_argument("-V", "--version", action="version", version=program_version_message)
        parser.add_argument("-D", "--debug", action="store_true", default=False, help="write <sample_name>.log file with debugging information")
        parser.add_argument("-w", "--whole-genome", action="store_true", dest="wholegenome", default=False, help="JSON file uses a whole genome reference, so don't write out the consensus, depth, and proportion arrays for each sample")
        parser.add_argument("--allele-output-threshold", dest="allele_min_reads", default=8, type=int, help="cutoff of # of reads below which allels for amino acids and nucleotide alleles will not be output [default: 8]")
        parser.add_argument('-o', '--out', metavar="FILE", type=argparse.FileType('w'), default=sys.stdout, help="output filename [default: stdout]")
        parser.add_argument("--output-format", type=str.lower, choices=('xml', 'json'), default='xml', help="output format [default: xml]")
        parser.add_argument("--min-base-qual", dest="bqual", default=5, type=int, help="What is the minimum base quality score to use a position (phred scale, i.e. 10=90, 20=99, 30=99.9 accuracy) [default: 5]")
        parser.add_argument("--consensus-proportion", default=0.8, type=float, help="minimum proportion required to call at base at that position, else 'N'. [default: 0.8]")
        parser.add_argument("--fill-gaps", nargs="?", const="n", default="n", dest="gap_char", help="fill no coverage gaps in the consensus sequence [default: n], optional parameter is either the character to use for filling [defaut: n] or `false` for no gap filled array")
        #parser.add_argument("--fill-gaps", nargs="?", const="n", default=None, dest="gap_char", help="fill no coverage gaps in the consensus sequence [default: False], optional parameter is the character to use for filling [defaut: n]") # TP edited
        parser.add_argument("--mark-deletions", nargs="?", const="_", dest="del_char", help="fill deletions in the consensus sequence [default: _] or `false` for consensus without deletions.")
        #parser.add_argument("--mark-deletions", nargs="?", const="_", dest="del_char", help="fill deletions in the consensus sequence [default: False], optional parameter is the character to use for filling [defaut: _]") # TP edited
        parser.add_argument("--original-bam", metavar="FILE", dest="original_bam", type=argparse.FileType('rb'), default=None, help="Original aligned BAM (pre-ASAP filtering) used to report mapped_reads. [default: use --bam]")
        parser.add_argument("--fastp-json", metavar="FILE", dest="fastp_json", type=argparse.FileType('r'), default=None, help="fastp/fastplong JSON file; provides total_reads and trimmed_reads attributes in output. [default: none]")
        parser.add_argument("--primer-stats", metavar="FILE", dest="primer_stats", type=argparse.FileType('r'), default=None, help="primer_masking_stats.tsv from MASK_PRIMERS step. [default: none]")
        parser.add_argument("--identity-stats", metavar="FILE", dest="identity_stats", type=argparse.FileType('r'), default=None, help="identity_filter_stats.tsv from IDENTITY_FILTER step. [default: none]")
        parser.add_argument("--smor-stats", metavar="FILE", dest="smor_stats", type=argparse.FileType('r'), default=None, help="smor_stats.tsv from SMOR/SMOR_CORRECTION step. [default: none]")
        parser.add_argument("--codon-correction", action="store_true", dest="codon_correction", default=False, help="annotate SNP pairs in the same codon with read-level allele-linkage info (<codon_merge>/<combo>) using GenBank CDS annotations. Requires --codon-correction-genbank. [default: False]")
        parser.add_argument("--codon-correction-genbank", metavar="FILE", dest="codon_correction_genbank", nargs='+', default=None, help="One or more GenBank files for --codon-correction codon boundary derivation. [default: none]")
        parser.add_argument("--codon-correction-error", dest="codon_correction_error", type=float, default=0.05, help="frequency tolerance (0-1) for 'complete' vs. 'partial' codon_merge linkage classification. [default: 0.05]")
        parser.add_argument("--codon-correction-min-reads", dest="codon_correction_min_reads", type=int, default=10, help="minimum spanning reads to confirm codon linkage. [default: 10]")
        parser.add_argument("--discover-roi", action="store_true", dest="discover_roi", default=False, help="add <linked_snps> field to SNP XML nodes showing read-level co-occurring variants. [default: False]")
        parser.add_argument("--discover-roi-min-perc", dest="discover_roi_min_perc", type=float, default=0.1, help="minimum co-occurrence proportion (0-1) to report a linked SNP. [default: 0.1]")
        parser.add_argument("--discover-roi-min-reads", dest="discover_roi_min_reads", type=int, default=10, help="minimum co-occurring read count to report a linked SNP. [default: 10]")
        parser.add_argument("--discover-roi-min-snp-perc", dest="discover_roi_min_snp_perc", type=float, default=0.05, help="minimum variant frequency (0-1) for a SNP to be considered (as anchor or candidate) in discover-roi linkage analysis. [default: 0.05]")

        # Process arguments
        args = parser.parse_args()

        json_fp = args.json
        bam_fp = args.bam
        bam_file = bam_fp.name
        depth = args.depth
        breadth = args.breadth
        proportion = args.proportion
        mutdepth = args.mutdepth
        debug = args.debug
        allele_min_reads = args.allele_min_reads
        wholegenome = args.wholegenome
        base_qual = args.bqual
        con_prop = args.consensus_proportion
        fill_gap_char = args.gap_char
        fill_del_char = args.del_char
        codon_correction = args.codon_correction
        codon_correction_genbank = args.codon_correction_genbank
        codon_correction_error = args.codon_correction_error
        codon_correction_min_reads = args.codon_correction_min_reads
        discover_roi = args.discover_roi
        discover_roi_min_perc = args.discover_roi_min_perc
        discover_roi_min_reads = args.discover_roi_min_reads
        discover_roi_min_snp_perc = args.discover_roi_min_snp_perc

        #out_dir = args.odir
        #if not out_dir:
        #    out_dir = os.getcwd()

        #out_dir = dispatcher.expandPath(out_dir)
        #if not os.path.exists(out_dir):
        #    os.makedirs(out_dir)

        operation_list = []
        operation_err = ""
        try:
            operation_list = assayInfo.parseOperation(args.json.name)
            pass
        except Exception as e:
            operation_err = str(e)
        assay_list = assayInfo.parseJSON(args.json.name)

        def _primary_mapped(bam_path):
            import re
            flagstat = pysam.flagstat(bam_path)
            m = re.search(r'^(\d+) \+ \d+ primary mapped', flagstat, re.MULTILINE)
            return m.group(1) if m else None

        def _load_ref_stats(fp):
            """Load a ref_name-keyed TSV stats file. Returns {} if fp is None or null sentinel."""
            import csv
            if fp is None:
                return {}
            try:
                if os.path.basename(fp.name) == 'null' or os.path.getsize(fp.name) == 0:
                    return {}
            except OSError:
                return {}
            stats = {}
            fp.seek(0)
            for row in csv.DictReader(fp, delimiter='\t'):
                stats[row['ref_name']] = row
            return stats

        samdata = pysam.AlignmentFile(bam_fp.name, "rb")
        sample_dict = {}
        if 'RG' in samdata.header.to_dict() :
            sample_dict['name'] = samdata.header.to_dict()['RG'][0]['ID']
        else:
            sample_dict['name'] = os.path.splitext(os.path.basename(bam_fp.name))[0]
        # Use original pre-filter BAM for mapped_reads (primary alignments only)
        bam_for_count = args.original_bam.name if args.original_bam else bam_fp.name
        sample_dict['mapped_reads'] = _primary_mapped(bam_for_count) or str(samdata.mapped)
        sample_dict['unmapped_reads'] = str(samdata.unmapped)
        sample_dict['unassigned_reads'] = str(samdata.nocoordinate)
        # Add pre-QC and post-QC read counts from fastp/fastplong JSON when available
        if args.fastp_json:
            import json as _json
            fastp_data = _json.load(args.fastp_json)
            sample_dict['total_reads'] = str(fastp_data['summary']['before_filtering']['total_reads'])
            sample_dict['trimmed_reads'] = str(fastp_data['summary']['after_filtering']['total_reads'])
        # Restore SMOR flag when the input BAM was produced by generateSMORbam
        if '_SMOR' in os.path.basename(bam_fp.name):
            sample_dict['SMOR'] = 'True'

        # Load per-reference stats from upstream filter steps
        primer_stats  = _load_ref_stats(args.primer_stats)
        identity_stats = _load_ref_stats(args.identity_stats)
        smor_stats    = _load_ref_stats(args.smor_stats)

        # Open original BAM for per-amplicon aligned_reads counts
        orig_samdata = pysam.AlignmentFile(args.original_bam.name, "rb") if args.original_bam else None

        sample_dict['depth_filter'] = str(depth)
        sample_dict['proportion_filter'] = str(proportion)
        sample_dict['breadth_filter'] = str(breadth)
        sample_dict['mutation_depth_filter'] = str(mutdepth)
        # minidom.parseString will raise xml.parsers.expat.ExpatError: not well-formed (invalid token)
        # if json_file or bam_file contain the python string representation of a file-like object.
        # e.g. bam_file="<_io.BufferedReader name=\'/shared/Targeted_sequence_fastqs/ASAP/COD-10-24_S302.bam\'>"
        sample_dict['json_file'] = json_fp.name
        sample_dict['bam_file'] = bam_fp.name
        sample_node = ElementTree.Element("sample", sample_dict)

        if INFO or debug:
            if not os.path.isdir("./bamProcessorLogs"):
                os.mkdir("./bamProcessorLogs")
                pass
            logfile = "./bamProcessorLogs/%s.log" % sample_dict['name']
            pass

        if INFO:
            logging.basicConfig(level=logging.INFO,
                                format='%(asctime)s %(levelname)-8s %(message)s',
                                datefmt='%m/%d/%Y %H:%M:%S',
                                filename=logfile,
                                filemode='w')
            pass
        if debug:
            logging.basicConfig(level=logging.DEBUG,
                                format='%(asctime)s %(levelname)-8s %(message)s',
                                datefmt='%m/%d/%Y %H:%M:%S',
                                filename=logfile,
                                filemode='w')
        if operation_err != "":
            logging.info("Operation fetch did not succeed: "+operation_err)
            pass
        logging.info("----------------------bamProcessor STARTED---------------------------")
        logging.info("JSON: "+str(json_fp))
        logging.info("BAM: "+str(bam_fp.name))
        for assay in assay_list:
            assay_dict = {}
            assay_dict['name'] = assay.name
            assay_dict['type'] = assay.assay_type
            assay_dict['function'] = assay.target.function or ""
            assay_dict['gene'] = assay.target.gene_name or ""
            assay_dict['start'] = assay.target.start_position or ""
            assay_dict['end'] = assay.target.end_position or ""
            logging.info("Assay: "+str(assay.name))
            #offset is where the amplicon sits relative to a reference, subtract 1 to make 0-based position
            #report the lesser of start and end in case amplicon is on reverse strand
            try:
                offset = min(int(assay.target.start_position), int(assay.target.end_position))-1
                ref_positions = list(range(min(int(assay.target.start_position), int(assay.target.end_position)),max(int(assay.target.start_position), int(assay.target.end_position))+1))
                #remove the zero if present
                if 0 in ref_positions:
                    ref_positions.remove(0)
            except:
                offset = 0
                ref_positions = None
            assay_node = ElementTree.SubElement(sample_node, "assay", assay_dict)
            ref_name = assay.name
            reverse_comp = assay.target.reverse_comp
            for amplicon in assay.target.amplicons:
                logging.info("+++AMPLICON+++")
                if not ref_positions:
                    ref_positions = list(range(1, len(amplicon.sequence)+1))
                temp_file = None
                ref_name = assay.name + "_%s" % amplicon.variant_name if amplicon.variant_name else assay.name
                logging.info("Now Checking Amplicon: "+str(ref_name))
                amplicon_dict = {}
                seq_counter = None
                if samdata.closed:
                    samdata = pysam.AlignmentFile(bam_file, "rb")
                amplicon_dict['reads'] = str(samdata.count(ref_name))
                if amplicon.variant_name:
                    amplicon_dict['variant'] = amplicon.variant_name
                # Per-amplicon read funnel from upstream filter steps
                if orig_samdata:
                    amplicon_dict['aligned_reads'] = str(orig_samdata.count(ref_name))
                if ref_name in primer_stats:
                    amplicon_dict['primer_reads']    = primer_stats[ref_name]['primer_reads']
                    amplicon_dict['no_primer_reads'] = primer_stats[ref_name]['no_primer_reads']
                if ref_name in identity_stats:
                    amplicon_dict['identity_input']     = identity_stats[ref_name]['input_reads']
                    amplicon_dict['identity_discarded'] = identity_stats[ref_name]['discarded_reads']
                if ref_name in smor_stats:
                    amplicon_dict['smor_input']           = smor_stats[ref_name]['input_reads']
                    amplicon_dict['smor_pairs_dropped']   = smor_stats[ref_name]['pairs_dropped']
                    amplicon_dict['smor_consensus_reads'] = smor_stats[ref_name]['consensus_reads']
                    amplicon_dict['smor_singleton_reads'] = smor_stats[ref_name].get('singleton_reads', '0')
                amplicon_node = ElementTree.SubElement(assay_node, "amplicon", amplicon_dict)
                if seq_counter:
                    ElementTree.SubElement(amplicon_node, "sequence_distribution", {k:str(v) for k,v in seq_counter.items()})
                if samdata.count(ref_name) == 0:
                    significance_node = ElementTree.SubElement(amplicon_node, "significance", {"flag":"no coverage"})
                    #Check for indeterminate resistances
                    resistances = set()
                    if amplicon.significance and amplicon.significance.resistance:
                        resistances.add(amplicon.significance.resistance)
                    for snp in amplicon.SNPs:
                        name = snp.name if snp.name else "position of interest"
                        dummy_snp = {'name':name, 'position':str(snp.position), 'depth':"0", 'reference':snp.reference, 'variant':snp.variant, 'basecalls':None}
                        _add_snp_node(amplicon_node, dummy_snp)
                        if snp.significance.resistance:
                            resistances.add(snp.significance.resistance)
                    if resistances:
                        significance_node.set("resistance", ",".join(resistances))
                else:
                    if amplicon.significance or samdata.count(ref_name) < depth:
                        significance_node = ElementTree.SubElement(amplicon_node, "significance")
                        if amplicon.significance:
                            significance_node.text = amplicon.significance.message
                            if amplicon.significance.resistance:
                                significance_node.set("resistance", amplicon.significance.resistance)
                        if samdata.count(ref_name) < depth:
                            significance_node.set("flag", "low coverage")
                            #Check for indeterminate resistances
                            resistances = set()
                            if amplicon.significance and amplicon.significance.resistance:
                                resistances.add(amplicon.significance.resistance)
                            for snp in amplicon.SNPs:
                                if snp.significance.resistance:
                                    resistances.add(snp.significance.resistance)
                            if resistances:
                                significance_node.set("resistance", ",".join(resistances))
                    # Warning: not designed to handle greater than 10 million X coverage
                    pileup_for_n_counting = samdata.pileup(ref_name, max_depth=10000000, ignore_orphans=False, ignore_overlaps=False)
                    amplicon_length = len(amplicon.sequence)

                    # First, run the N-counting function to get the N-read array.
                    n_read_array = _get_n_counts(pileup_for_n_counting, amplicon_length)
                    pileup = samdata.pileup(ref_name, max_depth=10000000, ignore_orphans=False, ignore_overlaps=False)
                    amplicon_data = _process_pileup(pileup, amplicon, depth, proportion, mutdepth, offset, wholegenome, base_qual, con_prop, fill_gap_char, fill_del_char, n_read_array)
                    if float(amplicon_data['breadth']) < breadth*100:
                        significance_node = amplicon_node.find("significance")
                        if significance_node is None:
                            significance_node = ElementTree.SubElement(amplicon_node, "significance")
                        if not significance_node.get("flag"):
                            significance_node.set("flag", "insufficient breadth of coverage")
                    # Build a fragment-allele table once per amplicon (single BAM
                    # pass), shared by codon correction and discover-roi to avoid
                    # re-scanning the BAM for every SNP/codon pair.
                    pos_table, reach, masked, deleted = {}, {}, {}, {}
                    if (codon_correction and codon_correction_genbank) or discover_roi:
                        amp_positions = {
                            _translated_to_amp(int(s['position']), offset)
                            for s in amplicon_data['SNPs'] if int(s.get('depth', 0)) > 0
                        }
                        if len(amp_positions) >= 2:
                            pos_table, reach, masked, deleted = _build_fragment_allele_table(samdata, amp_positions, ref_name)
                    # Apply codon correction (annotate same-codon SNPs via GenBank CDS)
                    if codon_correction and codon_correction_genbank:
                        cds_features = _parse_genbank_cds(codon_correction_genbank, amplicon.sequence)
                        _apply_codon_correction(
                            amplicon_data['SNPs'], pos_table, deleted,
                            cds_features, offset, codon_correction_error, codon_correction_min_reads
                        )
                    # Annotate co-occurring SNPs (discover-roi)
                    if discover_roi:
                        _apply_discover_roi(
                            amplicon_data['SNPs'], pos_table, reach, masked, offset,
                            discover_roi_min_perc, discover_roi_min_reads, discover_roi_min_snp_perc
                        )
                    # Handle SNPs
                    for snp in amplicon_data['SNPs']:
                        snp_node = _add_snp_node(amplicon_node, snp)
                        if 'codon_merges' in snp:
                            _add_codon_merges_node(snp_node, snp['codon_merges'])
                        if 'linked_snps' in snp:
                            _add_linked_snps_node(snp_node, snp['linked_snps'])
                    del amplicon_data['SNPs']
                    _write_parameters(amplicon_node, amplicon_data)
                    if not wholegenome:
                        ref_positions_node = ElementTree.SubElement(amplicon_node, "ref_positions")
                        ref_positions_node.text = ",".join(str(n) for n in ref_positions)
                if temp_file and REMOVE_TEMP:
                    samdata.close()
                    os.remove(temp_file)
                    os.remove(temp_file+".bai")
                    samdata = pysam.AlignmentFile(bam_file, "rb")

        # Close File
        if samdata.is_open():
            samdata.close()
        if orig_samdata and orig_samdata.is_open():
            orig_samdata.close()

        # Handle Operations
        # For each true operation add a 'significance' into the output detailing the operation
        for operation in operation_list:
            if evaluateOperation(operation, sample_node):
                significance_node = ElementTree.SubElement(sample_node, "operation")
                significance_node.set("flag", str(operation))
                significance_node.set("message", str(operation.message))
                pass
            pass


        logging.info("Writing output: "+str(args.out))
        _write_output(args.out, sample_node, args.output_format)

    except KeyboardInterrupt:
        logging.info("Process ended via KeyboardInterrupt.")
        pass
    except Exception as e:
        if DEBUG or TESTRUN:
            raise(e)
        indent = len("newBamProcessor") * " "
        sys.stderr.write("newBamProcessor: " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        logging.info("An Exception Occured! "+str(e))
        return 2
    logging.info("-------------------------bamProcessor FINISHED EXIT(0)---------------------------")
    return 0

def _write_output(file_obj, xml_element, output_format='xml'):
    if output_format == 'xml':
        from xml.dom import minidom
        dom = minidom.parseString(ElementTree.tostring(xml_element))
        file_obj.write(dom.toprettyxml(indent="  "))
    elif output_format == 'json':
        xml_str = ElementTree.tostring(xml_element)
        # The 'sample' root node is discarded as an unnecessary layer for the JSON object.
        xml_obj = xmltodict.parse(xml_str)['sample']
        # FIXME: The output is en/decoded multiple times because it seemed
        # easier to use the json object_hook to ensure each key had a
        # a consistent type than to write a nested loop with type checks
        # and conversions modifying the object as it was traversed.
        #
        # Ideally the output should start as a python object that is
        # encoded to XML or JSON once.
        json_encoded_xml = json.loads(json.dumps(xml_obj), object_hook=cast_json_output_types)
        json.dump(json_encoded_xml, file_obj, separators=(',', ':'))
    else:
        raise Exception('unsupported output format: %s' % output_format)

# cast_json_output_types is a json decoder object_hook intended to be used on
# a ASAP output decoded from XML:
# - casts numbers from strings to float/int
# - keys that are expected to contain 0 to n elements are lists (or undefined)
#   eliminating the 1 element object case.
# - As a special addition, values that were stored in the XML as strings of
#   comma separated values are converted to an array of an appropriate type.
def cast_json_output_types(e):
    ## Sample
    if '@breadth_filter' in e:
        e['@breadth_filter'] = float(e['@breadth_filter'])
    if '@depth_filter' in e:
        e['@depth_filter'] = int(e['@depth_filter'])
    # @json_file
    if '@mapped_reads' in e:
        e['@mapped_reads'] = int(e['@mapped_reads'])
    # @name
    if '@proportion_filter' in e:
        e['@proportion_filter'] = float(e['@proportion_filter'])
    if '@unassigned_reads' in e:
        e['@unassigned_reads'] = int(e['@unassigned_reads'])
    if '@unmapped_reads' in e:
        e['@unmapped_reads'] = int(e['@unmapped_reads'])
    if 'assay' in e and not isinstance(e['assay'], list):
        e['assay'] = [(e['assay'])]

    ## Assay
    # @function
    # @gene
    # @name
    # @type
    # amplicon
    if 'amplicon' in e and not isinstance(e['amplicon'], list):
        e['amplicon'] = [e['amplicon']]

    ## Amplicon
    if '@reads' in e:
        e['@reads'] = int(e['@reads'])
    #if 'significance' in e and not isinstance(e['significance'], dict):
    # consensus_sequence
    if 'breadth' in e:
        e['breadth'] = float(e['breadth'])
    if 'depths' in e:
        e['depths'] = [int(v) for v in e['depths'].split(',')]
    if 'proportions' in e:
        e['proportions'] = [float(v) for v in e['proportions'].split(',')]
    if 'average_depth' in e:
        e['average_depth'] = float(e['average_depth'])
    if 'snp' in e and not isinstance(e['snp'], list):
        e['snp'] = [e['snp']]
    ## SNP
    if '@depth' in e:
        e['@depth'] = int(e['@depth'])
    # @name
    if '@position' in e:
        e['@position'] = int(e['@position'])
    # @reference
    # snp_call
    if 'base_distribution' in e:
        e['base_distribution'] = {k: int(v) for k, v in e['base_distribution'].items()}

    ## CodonMerge
    if '@spanning_depth' in e:
        e['@spanning_depth'] = int(e['@spanning_depth'])
    if 'combo' in e and not isinstance(e['combo'], list):
        e['combo'] = [e['combo']]
    if 'codon_merge' in e and not isinstance(e['codon_merge'], list):
        e['codon_merge'] = [e['codon_merge']]

    ## SnpCall
    if '@count' in e:
        e['@count'] = int(e['@count'])
    if '@percent' in e:
        e['@percent'] = float(e['@percent'])
    # #text: "T"

    if '@changes' in e:
        e['@changes'] = int(e['@changes'])
    if 'mutation' in e and not isinstance(e['mutation'], list):
        e['mutation'] = [e['mutation']]

    ## Significance
    if 'significance' in e and isinstance(e['significance'], str):
        e['significance'] = {'#text': e['significance']}
    if '@resistance' in e:
        e['@resistance'] = e['@resistance'].split(',')

    return e

if __name__ == "__main__":
    if DEBUG:
        pass
    if TESTRUN:
        import doctest
        doctest.testmod()
    if PROFILE:
        import cProfile
        import pstats
        profile_filename = 'asap.newBamProcessor_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())
