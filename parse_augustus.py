#!/usr/bin/env python

import argparse
import re
from os.path import basename
from pvh.gff_utils import parse_gff_attributes, gff_string_from_list
from bx.intervals.intersection import IntervalTree, Interval
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import json

attribute_re = re.compile(' ?(?P<key>\S+) "(?P<value>[^"]+)";')
def parse_gtf_attributes(attribute_str):
    """Parse GTF attributes and return attributes dict

    GTF attributes are key value pairs with the value in double quotes and attributes terminated by ;
    @type attribute_str: str
    @rtype: dict
    """
    attributes = dict()
    for match in attribute_re.finditer(attribute_str):
        attributes[match.group('key')] = match.group('value')
    return attributes

hint_line_re = re.compile('#\s+(?P<hint_type>[^:]+):( +)?\d+ \((?P<hint_str>[^)]+)')
def parse_hint_group(line):
    match = hint_line_re.match(line)
    assert match != None, "could not match hint line regexp against: {}\n".format(line)
    return(match.group('hint_type'), match.group('hint_str'))

def parse_augustus(input_file, contig_name=None):
    """Given a augustus output file (with --gff=off), parse the GTF and proteins and return
    GFF3 lines list corresponding to the GTF and a list of BioPython SeqRecords corresponding to the proteins

    @type input_file: file
    @type contig_name: str
    @rtype list of list and dict of set
    """
    states = ()
    in_hints = False
    in_header = True
    in_gene = False
    in_protein = False
    in_cds = False
    in_gene_gtf = False
    blank_hints = [0.0, '', '']
    gene_name = None
    current_seq_str = None
    if contig_name == None:
        # contig name is the first part of filename: e.g. contig1.augustus.out
        filename = basename(input_file.name)
        fields = filename.split('.')
        contig_name = fields[0]
    gff_lines = []
    proteins = []
    CDSs = []
    evidence = dict()
    part_counter = dict()
    name_from_type = dict(start_codon='start', stop_codon='stop')
    for line in input_file:
        if in_header and line.startswith(':# Predicted genes for sequence number'):
            fields = line.rstrip().split()
            sequence_num = int(fields[-1])
            assert sequence_num == 1, "don't support augustus output with more than one sequence (contig), got sequence number {}".format(sequence_num)
            in_header = False
            # if not in_hints and line.startswith('# Constraints'):
            #     in_hints = True
            # elif in_hints and line.startswith('# Predicted genes for sequence number 1'):
            #     in_hints = False
            #     in_header = False
            # elif in_hints:
            #     # we're in the hints section
            #     fields = line.rstrip().split('\t')
            #     assert len(fields) == 9, "invalid hints line:\n{}".format(line)
            #     start = int(fields[3])
            #     end = int(fields[4])
            #     strand = fields[6]
            #     assert strand == '+' or strand == '-', 'unknown strand in line:\n{}'.format(line)
            #     if start == end:
            #         # fudge this by adding 1 to end: QUERY: is this wise? pvh
            #         end += 1
            #     elif start > end:
            #         # bad hint, ignore : TODO - fix upstream hint generator. pvh
            #         continue
            #     attributes = parse_gff_attributes(fields[8])
            #     transcript_id = attributes['grp']
            #     if strand == '+':
            #         forward_tree.insert(start, end, transcript_id)
            #     else:
            #         reverse_tree.insert(start, end, transcript_id)
        elif not in_gene and line.startswith('# start gene '):
            fields = line.rstrip().split()
            gene_id = fields[-1]
            gene_name = '_'.join([contig_name, gene_id])
            in_gene = True
            in_gene_gtf = True
        elif in_gene:
            if line.startswith('# end gene'):
                in_gene = False
                part_counter = dict()
            elif line.startswith('# hint groups fully obeyed:') or line.startswith('# incompatible hint groups:'):
                fields = line.rstrip().split()
                num_hint_groups = int(fields[-1])
                if num_hint_groups > 0:
                    in_obeyed_hints = in_incompatible_hints = False
                    if line.startswith('# hint groups fully obeyed:'):
                        in_obeyed_hints = True
                    elif line.startswith('# incompatible hint groups:'):
                        in_incompatible_hints = True
                    line = input_file.next()
                    (hint_type, hints_str) = parse_hint_group(line.rstrip())
                    # fix this below
                    assert hint_type == 'XNT', 'right now we can only deal with a single XNT hint type'
                    current_evidence = evidence.get(gene_name, blank_hints[:])
                    if in_obeyed_hints:
                        current_evidence[1] = hints_str
                    elif in_incompatible_hints:
                        current_evidence[2] = hints_str
                    evidence[gene_name] = current_evidence
            elif line.startswith('# % of transcript supported by hints (any source):'):
                fields = line.rstrip().split()
                hint_support = float(fields[-1])
                current_evidence = evidence.get(gene_name, blank_hints[:])
                current_evidence[0] = hint_support
                evidence[gene_name] = current_evidence
            elif line.startswith('# protein sequence') or line.startswith('# coding sequence'):
                if line.startswith('# protein sequence'):
                    in_protein = True
                    in_cds = False
                elif line.startswith('# coding sequence'):
                    in_cds = True
                    in_protein = False
                assert (in_cds == True or in_protein == True), "found state that shouldn't exist (in_cds: {} in_protein: {})at line:\n{}".format(in_cds, in_protein, line)
                in_gene_gtf = False
                line = line.rstrip()
                fields = line.split()
                # protein is last word starts with a '[' that we skip
                if line.endswith(']'):
                    # protein starts and ends on this line, grab it and don't set in_protein = True
                    current_seq_str = fields[-1][1:-1]
                    if in_protein:
                        alphabet = IUPAC.protein
                    elif in_cds:
                        alphabet = IUPAC.ambiguous_dna
                    seq = SeqRecord(Seq(current_seq_str, alphabet), id=gene_name, description='')
                    # ok we've got the whole protein or CDS, not in_protein or in_cds anymore
                    if in_protein:
                        proteins.append(seq)
                        in_protein = False
                    elif in_cds:
                        CDSs.append(seq)
                        in_cds = False
                else:
                    current_seq_str = fields[-1][1:]
            elif in_gene_gtf:
                line = line.rstrip()
                fields = line.split('\t')
                assert len(fields) == 9, "invalid GTF line in gene {}:\n{}\n".format(gene_name, line)
                (chrom, tool, seq_type, start_str, end_str, score_str, strand, frame_str, attribute_str) = fields
                start = int(start_str)
                end = int(end_str)
                attributes = dict()
                if seq_type == 'gene':
                    attributes = dict(ID=gene_name)
                elif seq_type == 'transcript':
                    # transcripts are given IDs like g1.t1 by augustus
                    transcript_id = '_'.join([gene_name,attribute_str.split('.')[1]])
                    attributes = dict(ID=transcript_id)
                else:
                    gtf_attributes = parse_gtf_attributes(attribute_str)
                    parent_id = '_'.join([gene_name,gtf_attributes['transcript_id'].split('.')[1]])
                    attributes['Parent'] = parent_id
                    if seq_type == 'CDS' or seq_type == 'intron':
                        index = part_counter.get(seq_type + parent_id, 1)
                        attributes['ID'] = '_'.join([parent_id,seq_type,str(index)])
                        part_counter[seq_type + parent_id] = index + 1
                    elif seq_type == 'stop_codon' or seq_type == 'start_codon':
                        name = name_from_type[seq_type]
                        attributes['ID'] = '_'.join([parent_id, name])
                gff_lines.append(gff_string_from_list([chrom, tool, seq_type, start_str, end_str, 
                                                       score_str, strand, frame_str, attributes]))
            elif in_protein or in_cds:
                line = line.rstrip()
                if line.startswith('# Evidence for and against this transcript'):
                    in_protein = False
                    in_cds = False
                elif line.endswith(']'):
                    current_seq_str += line[2:-1]
                    if in_protein:
                        alphabet = IUPAC.protein
                    elif in_cds:
                        alphabet = IUPAC.ambiguous_dna
                    seq = SeqRecord(Seq(current_seq_str, alphabet), id=gene_name, description='')
                    if in_protein:
                        proteins.append(seq)
                        in_protein = False
                    elif in_cds:
                        CDSs.append(seq)
                        in_cds = False
                else:
                    current_seq_str += line[2:]
    return (gff_lines, proteins, CDSs, evidence)



parser = argparse.ArgumentParser(description="Parse Augustus output and generate GFF3 of genes and predicted protein FASTA")
parser.add_argument('augustus_file', type=argparse.FileType(), help='Augustus file as generated with the --gff3=off flag')
parser.add_argument('gff3_output_file', type=argparse.FileType('w'))
parser.add_argument('protein_fasta_output_file', type=argparse.FileType('w'))
parser.add_argument('CDS_fasta_output_file', type=argparse.FileType('w'))
parser.add_argument('evidence_output_file', nargs='?', type=argparse.FileType('w'))

args = parser.parse_args()

(gff_lines, proteins, CDSs, evidence) = parse_augustus(args.augustus_file)
if len(gff_lines) > 0:
    args.gff3_output_file.write('##gff-version 3\n'+''.join(gff_lines)+'\n')
args.gff3_output_file.close()

json.dump(evidence, args.evidence_output_file)

SeqIO.write(proteins, args.protein_fasta_output_file, 'fasta')
args.protein_fasta_output_file.close()
SeqIO.write(CDSs, args.CDS_fasta_output_file, 'fasta')
args.CDS_fasta_output_file.close()
args.evidence_output_file.close()