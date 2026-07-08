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
import logging

import pysam
import statistics
from collections import Counter, defaultdict
from xml.etree import ElementTree

from asap import assayInfo
from asap import __version__
from asap.genbank_cds import _parse_genbank_cds
from asap.allele_linkage import (
    _build_fragment_allele_table, _translated_to_amp,
    _apply_codon_correction, _apply_discover_roi,
)
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

def _process_pileup(pileup, amplicon, depth, proportion, mutdepth, offset, wholegenome, base_qual, con_prop, fill_gap_char, fill_del_char):
    """
    Single-pass replacement for the old _get_n_counts + _process_pileup pair:
    consumes one pileup iterator and produces both the N-read array (per
    unique alignment, quality-independent) and the consensus/depth/SNP-calling
    output (per pileupcolumn, quality-filtered), instead of scanning the BAM
    twice for the same amplicon.
    """
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
    n_read_array = [0] * amplicon_length
    processed_alignments = set()
    depth_array = [0] * amplicon_length
    quality_discard_array = [0] * amplicon_length
    prop_array = ["0"] * amplicon_length
    previous_position = 0
    # for each position in alignment/pileup
    for pileupcolumn in pileup:
        base_counter = Counter()
        base_quality_scores = defaultdict(list)
        base_R1_counter = Counter()
        base_R2_counter = Counter()
        base_SE_counter = Counter()
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
            # --- N-counting (from the old _get_n_counts pass): runs once per
            # unique alignment regardless of base quality, using the same
            # hybrid soft-clip + aligned-region scan and the same dedup key. ---
            try:
                alignment = pileupread.alignment
                alignment_id = (alignment.query_name, alignment.reference_start, alignment.query_length)

                if alignment_id not in processed_alignments:
                    current_read_sequence = alignment.query_sequence.upper()

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

                    processed_alignments.add(alignment_id)

            except Exception as e:
                logging.warning(f"Skipped read in N-counting: {e}")
            #print("processing read, qual=%i" % pileupread.alignment.query_qualities[pileupread.query_position])
            try:
                if pileupread.is_del:
                    #This position in the alignment is a deletion in the query sequence, therefore it has no quality score
                    # Let's use the average of the quality scores of the two aligned bases flanking the deletion
                    quals = pileupread.alignment.query_qualities
                    q_next = pileupread.query_position_or_next
                    q_prev = q_next - 1
                    q_next = min(q_next, len(quals) - 1)
                    q_prev = max(q_prev, 0)
                    qscore = (quals[q_next] + quals[q_prev]) / 2
                    if qscore >= base_qual:
                        passed_Qual_filter += 1
                        base_counter.update({"_": 1})
                        base_quality_scores["_"].append(qscore)
                        if pileupread.alignment.is_read1:
                            base_R1_counter.update({"_": 1})
                        elif pileupread.alignment.is_read2:
                            base_R2_counter.update({"_": 1})
                        else:
                            base_SE_counter.update({"_": 1})
                    else:
                        quality_discard_array[pileupcolumn.pos] += 1
                elif pileupread.alignment.query_qualities[pileupread.query_position] >= base_qual: # check here
                    passed_Qual_filter += 1
                    qual_score = pileupread.alignment.query_qualities[pileupread.query_position]
                    if pileupread.indel < 0: #This means the next position is a deletion, we'll process later
                        for d in range(1, abs(pileupread.indel)+1):
                            deletion_counter.update({str(position + d)})
                    if pileupread.indel > 0: #This means the next position is an insertion, unlike with deletions, this we can process now
                        start = pileupread.query_position
                        end = pileupread.query_position + pileupread.indel + 1
                        base_key = pileupread.alignment.query_sequence[start:end]
                        base_counter.update({base_key: 1})
                        base_quality_scores[base_key].append(qual_score)
                        if pileupread.alignment.is_read1:
                            base_R1_counter.update({base_key: 1})
                        elif pileupread.alignment.is_read2:
                            base_R2_counter.update({base_key: 1})
                        else:
                            base_SE_counter.update({base_key: 1})
                    else:
                        base_key = pileupread.alignment.query_sequence[pileupread.query_position]
                        base_counter.update(base_key)
                        base_quality_scores[base_key].append(qual_score)
                        if pileupread.alignment.is_read1:
                            base_R1_counter.update({base_key: 1})
                        elif pileupread.alignment.is_read2:
                            base_R2_counter.update({base_key: 1})
                        else:
                            base_SE_counter.update({base_key: 1})
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
                snp = {'name':name, 'position':str(translated), 'depth':str(column_depth), 'reference':reference, 'variant':variant, 'basecalls':base_counter, 'base_qualities':base_quality_scores, 'base_R1':base_R1_counter, 'base_R2':base_R2_counter, 'base_SE':base_SE_counter}
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
            snp = {'name':f"{reference_call}{translated}{snp_call}", 'position':str(translated), 'depth':str(column_depth), 'reference':reference_call, 'variant':snp_call, 'basecalls':base_counter, 'base_qualities':base_quality_scores, 'base_R1':base_R1_counter, 'base_R2':base_R2_counter, 'base_SE':base_SE_counter}
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
        pileup_dict['quality_discards'] = ",".join(str(n) for n in quality_discard_array)
    pileup_dict['breadth'] = str(breadth_positions/amplicon_length * 100)
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
        base_R1 = snp.get('base_R1')
        base_R2 = snp.get('base_R2')
        base_SE = snp.get('base_SE')
        if base_R1 or base_R2:
            strand_node = ElementTree.SubElement(snp_node, 'base_strand_distribution')
            for base in base_counter:
                attrs = {
                    'base': base,
                    'R1':   str(base_R1.get(base, 0) if base_R1 else 0),
                    'R2':   str(base_R2.get(base, 0) if base_R2 else 0),
                }
                if base_SE and base_SE.get(base, 0):
                    attrs['SE'] = str(base_SE[base])
                ElementTree.SubElement(strand_node, 'strand', attrs)
        base_qualities = snp.get('base_qualities')
        if base_qualities:
            qual_node = ElementTree.SubElement(snp_node, 'base_quality')
            for base in base_counter:
                if base in base_qualities and base_qualities[base]:
                    q = base_qualities[base]
                    ElementTree.SubElement(qual_node, 'qual', {
                        'base':   base,
                        'mean':   f"{statistics.mean(q):.1f}",
                        'median': f"{statistics.median(q):.1f}",
                        'min':    str(min(q)),
                        'max':    str(max(q)),
                    })
    return snp_node



def _add_linked_snps_node(snp_node, linked_snps, anchor_name):
    """Add a <linked_snps> sub-element to a SNP XML node."""
    ls_node = ElementTree.SubElement(snp_node, 'linked_snps')
    for entry in linked_snps:
        link_node = ElementTree.SubElement(ls_node, 'linked_snp', {
            'anchor_variant': anchor_name,
            'target_variant': entry['target_name'],
            'shared_read_depth': str(entry['shared_read_depth']),
            'co_occurring_count': str(entry['co_count']),
            'linkage_pct': str(entry['linkage_pct']),
            'sample_frequency_pct': str(entry['sample_frequency_pct']),
        })
        ElementTree.SubElement(link_node, 'read_evidence', {
            'both_variants': str(entry['both_variants']),
            'anchor_only': str(entry['anchor_only']),
            'target_only': str(entry['target_only']),
            'neither_variant': str(entry['neither_variant']),
            'masked': str(entry['masked']),
            'read_too_short': str(entry['read_too_short']),
        })


def _add_codon_merges_node(snp_node, codon_merges):
    """Add one <codon_merge> child element per entry in codon_merges."""
    for entry in codon_merges:
        cm_node = ElementTree.SubElement(snp_node, 'codon_merge', {
            'name': entry['name'],
            'region': entry['region'],
            'direction': entry['direction'],
            'position': entry['position'],
            'codon_depth': str(entry['codon_depth']),
            'reference': entry['reference'],
        })
        ElementTree.SubElement(cm_node, 'codon_call', {
            'codon_count': str(entry['codon_call_count']),
            'percentage': f"{entry['codon_call_percentage']:.1f}",
        }).text = entry['codon_call']
        dist_attribs = {
            seq: str(cnt)
            for seq, cnt in sorted(entry['codon_distribution'].items(),
                                   key=lambda kv: -kv[1])
        }
        ElementTree.SubElement(cm_node, 'codon_distribution', dist_attribs)
        ElementTree.SubElement(cm_node, 'excluded_reads', {
            'has_n': str(entry['excl_has_n']),
            'no_span': str(entry['excl_no_span']),
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
        # Report mapped_reads/unmapped_reads/unassigned_reads from the original,
        # pre-ASAP-filter BAM only, so all three come from one consistent snapshot
        # (not mixed with primer-masking/identity-filtering/SMOR collapsing applied
        # to bam_fp downstream). Amplicon-specific funnel fields (aligned_reads,
        # amplicon_reads, etc.) still track the filtered BAM separately below.
        orig_samdata = pysam.AlignmentFile(args.original_bam.name, "rb") if args.original_bam else None
        count_samdata = orig_samdata or samdata
        bam_for_count = args.original_bam.name if args.original_bam else bam_fp.name
        sample_dict['mapped_reads'] = _primary_mapped(bam_for_count) or str(count_samdata.mapped)
        sample_dict['unmapped_reads'] = str(count_samdata.unmapped)
        sample_dict['unassigned_reads'] = str(count_samdata.nocoordinate)
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
                # Cache this amplicon's read count instead of recomputing it
                # (via a fresh index scan) on every use below.
                amplicon_read_count = samdata.count(ref_name)
                amplicon_dict['reads'] = str(amplicon_read_count)
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
                if amplicon_read_count == 0:
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
                    if amplicon.significance or amplicon_read_count < depth:
                        significance_node = ElementTree.SubElement(amplicon_node, "significance")
                        if amplicon.significance:
                            significance_node.text = amplicon.significance.message
                            if amplicon.significance.resistance:
                                significance_node.set("resistance", amplicon.significance.resistance)
                        if amplicon_read_count < depth:
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
                    # Single pileup pass now produces both the N-read array and
                    # the consensus/depth/SNP data (was two separate full passes).
                    pileup = samdata.pileup(ref_name, max_depth=10000000, ignore_orphans=False, ignore_overlaps=False)
                    amplicon_data = _process_pileup(pileup, amplicon, depth, proportion, mutdepth, offset, wholegenome, base_qual, con_prop, fill_gap_char, fill_del_char)
                    if float(amplicon_data['breadth']) < breadth*100:
                        significance_node = amplicon_node.find("significance")
                        if significance_node is None:
                            significance_node = ElementTree.SubElement(amplicon_node, "significance")
                        if not significance_node.get("flag"):
                            significance_node.set("flag", "insufficient breadth of coverage")
                    # Parse CDS features first so all codon positions can be
                    # included in the fragment-allele table (the 3rd codon
                    # position may not be a SNP and would otherwise be absent).
                    cds_features = []
                    if codon_correction and codon_correction_genbank:
                        cds_features = _parse_genbank_cds(codon_correction_genbank, amplicon.sequence)
                    # Build a fragment-allele table once per amplicon (single BAM
                    # pass), shared by codon correction and discover-roi to avoid
                    # re-scanning the BAM for every SNP/codon pair.
                    pos_table, reach, masked, deleted = {}, {}, {}, {}
                    if cds_features or discover_roi:
                        amp_positions = {
                            _translated_to_amp(int(s['position']), offset)
                            for s in amplicon_data['SNPs'] if int(s.get('depth', 0)) > 0
                        }
                        # Include all 3 positions of every codon so non-SNP
                        # positions are available for the full-codon tally.
                        # Only for codons that actually contain a SNP --
                        # _apply_codon_correction skips any codon without
                        # exactly 2 SNPs, so codons with none are never used,
                        # and pulling in every codon of every CDS (i.e. nearly
                        # the whole genome for densely-coding references)
                        # blows up the fragment-allele table at high coverage.
                        snp_amp_positions = frozenset(amp_positions)
                        for cds in cds_features:
                            for start, end, _ref in cds.codon_boundaries:
                                if snp_amp_positions.intersection(range(start, end)):
                                    amp_positions.update(range(start, end))
                        if len(amp_positions) >= 2:
                            pos_table, reach, masked, deleted = _build_fragment_allele_table(samdata, amp_positions, ref_name)
                    # Apply codon correction (annotate same-codon SNPs via GenBank CDS)
                    if cds_features:
                        _apply_codon_correction(
                            amplicon_data['SNPs'], pos_table, deleted, masked,
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
                            _add_linked_snps_node(snp_node, snp['linked_snps'], snp['name'])
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
