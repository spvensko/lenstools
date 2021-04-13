#!/usr/bin/env python

import argparse
import vcf
from Bio import SeqIO
from Bio.SeqUtils import seq1
from Bio.SeqUtils import Seq
import re
import sys
from pprint import pprint
import numpy as np
from scipy import stats
import csv
import pysam
import hashlib
import os

from pysam.libcalignmentfile import IteratorColumnRegion

def get_args():
    """
    """
    parser = argparse.ArgumentParser(prog="lenstools",
                                     description="lenstools")

    subparsers = parser.add_subparsers(dest='command')

    # Subparser for generating peptides
    parser_make_peptides = subparsers.add_parser('make-peptides',
                                                 help="Make peptide FASTA file.")
    parser_make_peptides.add_argument('-t', '--tx-aa-fasta',
                                      help="Transcript amino acid sequences.",
                                      required=True)
    parser_make_peptides.add_argument('-v', '--vcf',
                                      help="Annotated (by snpeff/?) VCF).",
                                      required=True)
    parser_make_peptides.add_argument('-l', '--length',
                                      help="Total peptide length.",
                                      default=8)
    parser_make_peptides.add_argument('-o', '--mt-output',
                                      help="Mutant output file name.",
                                      required=True)
    parser_make_peptides.add_argument('-w', '--wt-output',
                                      help="Wildtype peptides output file name.",
                                      required=True)

    # Subparser for generating genomic context
    parser_make_context = subparsers.add_parser('make-genomic-context',
                                                 help="Make nucleotide FASTA file.")
    parser_make_context.add_argument('-c', '--tx-cds-fasta',
                                      help="Transcript CDS.",
                                      required=True)
    parser_make_context.add_argument('-v', '--vcf',
                                      help="Annotated (by snpeff/?) VCF).",
                                      required=True)
    parser_make_context.add_argument('-l', '--length',
                                      help="Total nucleotide length.",
                                      default=39)
    parser_make_context.add_argument('-o', '--output',
                                      help="output file name.",
                                      required=True)

    # Subparser for generating indel peptides
    parser_make_indel_peptides = subparsers.add_parser('make-indel-peptides',
                                                 help="Make indel peptide FASTA file.")
    parser_make_indel_peptides.add_argument('-t', '--tx-aa-fasta',
                                      help="Transcript amino acid sequences.",
                                      required=True)
    parser_make_indel_peptides.add_argument('-v', '--vcf',
                                      help="Annotated (by snpeff/?) VCF).",
                                      required=True)
    parser_make_indel_peptides.add_argument('-l', '--length',
                                      help="Total peptide length.",
                                      default=8)
    parser_make_indel_peptides.add_argument('-o', '--output',
                                      help="Output file name.",
                                      required=True)


    # Subparser for filtering variants for expressin
    parser_expressed_variants = subparsers.add_parser('expressed-variants',
                                                 help="Filter expressed variants.")
    parser_expressed_variants.add_argument('-a', '--abundances',
                                           help="Transcript abundance file.",
                                           required=True)
    parser_expressed_variants.add_argument('-m', '--metric',
                                           help="Column for expression from abundance file.",
                                           default="TPM")
    parser_expressed_variants.add_argument('-r', '--exclude-zeros',
                                           help="Exclude zero counts when calculating percentile.",
                                           action='store_true')
    parser_expressed_variants.add_argument('-p', '--percentile',
                                           help="Percentile for determining expression.",
                                           default=50)
    parser_expressed_variants.add_argument('-t', '--abundance-threshold',
                                           help="abundance threshold.",
                                           default=0)
    parser_expressed_variants.add_argument('-v', '--vcf',
                                      help="Annotated VCF).",
                                      required=True)
    parser_expressed_variants.add_argument('-o', '--output',
                                      help="Output file name.",
                                      required=True)

    # Subparser for filtering variants for having sufficient RNA coverage
    parser_covered_variants = subparsers.add_parser('rna-covered-variants',
                                                 help="Filter variants for RNA coverage.")
    parser_covered_variants.add_argument('-b', '--rna-bam',
                                           help="Transcript BAM file.",
                                           required=True)
    parser_covered_variants.add_argument('-c', '--required-coverage',
                                           help="Required coverage for variant",
                                           default=1)
    parser_covered_variants.add_argument('-v', '--vcf',
                                      help="Annotated VCF).",
                                      required=True)
    parser_covered_variants.add_argument('-o', '--output',
                                      help="Output file name.",
                                      required=True)

    # Subparser for filtering variants for expression
    parser_isolated_variants = subparsers.add_parser('isolated-variants',
                                                 help="Filter isolated variants.")
    parser_isolated_variants.add_argument('-g', '--germline-vcf',
                                           help="Germline VCF.",
                                           required=True)
    parser_isolated_variants.add_argument('-s', '--somatic-vcf',
                                           help="Somatic VCF.",
                                           required=True)
    parser_isolated_variants.add_argument('-p', '--proximity',
                                           help="Required clearance around variant (in bp).",
                                           default=30)
    parser_isolated_variants.add_argument('-a', '--allow-silent-and-homozygous',
                                           help="Allow silent and homozygous germline variants near variant.",
                                           action="store_true")
    parser_isolated_variants.add_argument('-o', '--output',
                                          help="Output file name.",
                                          required=True)

    # Subparser for calculalting agretopicity (mut BA/wt BA)
    parser_calculate_agretopicity = subparsers.add_parser('calculate-agretopicity',
                                                          help="Calcuate agreotopicity (mut BA/wt BA).")
    parser_calculate_agretopicity.add_argument('-w', '--wt-fasta',
                                               help="Wildtype peptide FASTA.",
                                               required=True)
    parser_calculate_agretopicity.add_argument('-m', '--mt-fasta',
                                               help="Mutant peptide FASTA.",
                                               required=True)
    parser_calculate_agretopicity.add_argument('-o', '--output',
                                               help="Output file name.",
                                               required=True)

    # Subparser for creating PyClone-VI inputs
    parser_make_pvi_inputs = subparsers.add_parser('make-pyclone-vi-inputs',
                                                          help="Calcuate agreotopicity (mut BA/wt BA).")
    parser_make_pvi_inputs.add_argument('-c', '--candidate-vcf',
                                        help="VCF containing candidate variants.",
                                        required=True)
    parser_make_pvi_inputs.add_argument('-m', '--mutect-vcf',
                                        help="Mutect VCF for variant depth information.",
                                        required=True)
    parser_make_pvi_inputs.add_argument('-s', '--sequenza-segments',
                                        help="Sequenza segments file.",
                                        required=True)
    parser_make_pvi_inputs.add_argument('--sequenza-solutions',
                                        help="Sequenza solutions file.",
                                        required=True)
    parser_make_pvi_inputs.add_argument('--samp-id', required=True)
    parser_make_pvi_inputs.add_argument('-o', '--output', required=True)

    # Subparser for adding SNV metadata
    parser_add_snv_metadata = subparsers.add_parser('add-snv-metadata',
                                                   help="Add metadata to snv netmhcpan output")
    parser_add_snv_metadata.add_argument('-m', '--mutant-peptides',
                                        help="FASTA file with mutant peptides.",
                                        required=True)
    parser_add_snv_metadata.add_argument('-u', '--mutant-nucs',
                                        help="FASTA file with mutant nucleotide context.",
                                        required=True)
    parser_add_snv_metadata.add_argument('-q', '--quants',
                                        help="Quant file.",
                                        required=True)
    parser_add_snv_metadata.add_argument('-c', '--cancer-cell-fraction',
                                        help="Cancer cell fraction file.",
                                        required=True)
    parser_add_snv_metadata.add_argument('-g', '--gtf')
    parser_add_snv_metadata.add_argument('-n', '--netmhcpan')
    parser_add_snv_metadata.add_argument('-o', '--output', required=True)

    # Subparser for adding SNV metadata
    parser_add_indel_metadata = subparsers.add_parser('add-indel-metadata',
                                                   help="Add metadata to indel netmhcpan output")
    parser_add_indel_metadata.add_argument('-m', '--mutant-peptides',
                                        help="FASTA file with mutant peptides.",
                                        required=True)
    parser_add_indel_metadata.add_argument('-q', '--quants',
                                        help="Quant file.",
                                        required=True)
    parser_add_indel_metadata.add_argument('-c', '--cancer-cell-fraction',
                                        help="Cancer cell fraction file.",
                                        required=True)
    parser_add_indel_metadata.add_argument('-g', '--gtf')
    parser_add_indel_metadata.add_argument('-n', '--netmhcpan')
    parser_add_indel_metadata.add_argument('-o', '--output', required=True)


    # Subparser for filtering expressed hERVs
    parser_expressed_hervs = subparsers.add_parser('expressed-hervs',
                                                   help="Filter expressed hervs.")
    parser_expressed_hervs.add_argument('-a', '--abundances',
                                        help="Transcript abundance file.",
                                        required=True)
    parser_expressed_hervs.add_argument('-m', '--metric',
                                        help="Column for expression from abundance file.",
                                        default="TPM")
    parser_expressed_hervs.add_argument('-r', '--exclude-zeros',
                                        help="Exclude zero counts when calculating percentile.",
                                        action='store_true')
    parser_expressed_hervs.add_argument('-p', '--percentile',
                                        help="Percentile for determining expression.",
                                        default=50)
    parser_expressed_hervs.add_argument('-t', '--abundance-threshold',
                                        help="abundance threshold.",
                                        default=0)
    parser_expressed_hervs.add_argument('-o', '--output',
                                        help="Output file name.",
                                        required=True)

    # Subparser for making hERV peptides
    parser_make_herv_peptides = subparsers.add_parser('make-herv-peptides',
                                                   help="Create hERV peptides.")
    parser_make_herv_peptides.add_argument('-e', '--expressed-hervs',
                                        help="Expressed hERVs",
                                        required=True)
    parser_make_herv_peptides.add_argument('-r', '--herv-ref',
                                        help="hERV reference FASTA.",
                                        required=True)
    parser_make_herv_peptides.add_argument('-o', '--output',
                                        help="Output file name.",
                                        required=True)


    parser_get_herv_metadata = subparsers.add_parser('get-herv-metadata',
                                                   help="Create hERV peptides.")
    parser_get_herv_metadata.add_argument('-a', '--abundances',
                                        help="Expressed hERVs",
                                        required=True)
    parser_get_herv_metadata.add_argument('-n', '--netmhcpan',
                                        help="Expressed hERVs",
                                        required=True)
    parser_get_herv_metadata.add_argument('-o', '--output',
                                        help="Output file name.",
                                        required=True)
    # Subparser for adding SNV metadata
    parser_consol_mqc_stats = subparsers.add_parser('consolidate-multiqc-stats',
                                                   help="Add metadata to indel netmhcpan output")
    parser_consol_mqc_stats.add_argument('-d', '--multiqc-data',
                                        help="Directory containing multiqc stats files.",
                                        required=True)
    parser_consol_mqc_stats.add_argument('-o', '--output', required=True)

    # Subparser for filtering expression from a provided list
    parser_expressed_self_genes = subparsers.add_parser('expressed-self-genes',
                                                 help="Filter expressed self_genes.")
    parser_expressed_self_genes.add_argument('-a', '--abundances',
                                           help="Transcript abundance file.",
                                           required=True)
    parser_expressed_self_genes.add_argument('-m', '--metric',
                                           help="Column for expression from abundance file.",
                                           default="TPM")
    parser_expressed_self_genes.add_argument('-r', '--exclude-zeros',
                                           help="Exclude zero counts when calculating percentile.",
                                           action='store_true')
    parser_expressed_self_genes.add_argument('-p', '--percentile',
                                           help="Percentile for determining expression.",
                                           default=50)
    parser_expressed_self_genes.add_argument('-t', '--abundance-threshold',
                                           help="abundance threshold.",
                                           default=0)
    parser_expressed_self_genes.add_argument('-g', '--gene-list',
                                             help="Genes list (e.g. self-antigens, CTAs, etc.)",
                                      required=True)
    parser_expressed_self_genes.add_argument('-f', '--gtf',
                                      help="GTF",
                                      required=True)
    parser_expressed_self_genes.add_argument('-o', '--output',
                                      help="Output file name.",
                                      required=True)


    # Subparser for making self-antigen peptides
    parser_make_self_peptides = subparsers.add_parser('make-self-antigen-peptides',
                                           help="Make self-antigen peptides.")
    parser_make_self_peptides.add_argument('-e', '--expressed-selfs',
                                      required=True)
    parser_make_self_peptides.add_argument('-r', '--tx-aa-fasta',
                                      help="GTF",
                                      required=True)
    parser_make_self_peptides.add_argument('-s', '--somatic-vcf',
                                      help="Somatic VCF",
                                      required=True)
    parser_make_self_peptides.add_argument('-g', '--germline-vcf',
                                      help="Somatic VCF",
                                      required=True)
    parser_make_self_peptides.add_argument('-o', '--output',
                                      help="Output file",
                                      required=True)


    # Subparser for getting self-antigen metadata
    parser_add_self_metadata = subparsers.add_parser('add-self-antigen-metadata',
                                                     help="Add self-antigen metadata.")
    parser_add_self_metadata.add_argument('-q', '--quants',
                                          required=True)
    parser_add_self_metadata.add_argument('-b', '--binding-affinities',
                                          help="Binding affinity data (netMHCpan format)",
                                          required=True)
    parser_add_self_metadata.add_argument('-g', '--gtf',
                                          help="GTF",
                                          required=True)
    parser_add_self_metadata.add_argument('-o', '--output',
                                          help="Output file",
                                          required=True)


   # Subparser for adding SNV metadata
    parser_filter_virdetect_by_counts = subparsers.add_parser('filter-virdetect-by-counts',
                                                   help="Generate list of expressed viruses")
    parser_filter_virdetect_by_counts.add_argument('-q', '--viral-quants',
                                        help="Viral counts file from VirDetect.",
                                        required=True)
    parser_filter_virdetect_by_counts.add_argument('-r', '--viral-ref',
                                        help="Viral reference FASTA (used for VirDetect).",
                                        required=True)
    parser_filter_virdetect_by_counts.add_argument('-m', '--min-threshold',
                                        help="Minimal count threshold.",
                                        default='1')
    parser_filter_virdetect_by_counts.add_argument('-o', '--output', required=True)

    # Subparser for making viral peptides
    parser_make_viral_peptides = subparsers.add_parser('make-viral-peptides',
                                                       help="Make viral peptides")
    parser_make_viral_peptides.add_argument('-e', '--expressed-viruses',
                                        help="Viral counts file from VirDetect.",
                                        required=True)
#    parser_make_viral_peptides.add_argument('-p', '--viral-pep-ref',
#                                        help="Viral peptide reference FASTA (used for VirDetect).",
#                                        required=True)
#    parser_make_viral_peptides.add_argument('-c', '--viral-cds-ref',
#                                        help="Viral CDS reference FASTA (used for VirDetect).",
#                                        required=True)
    parser_make_viral_peptides.add_argument('-c', '--sample-viral-ref',
                                        help="Viral CDS reference FASTA (used for VirDetect).",
                                        required=True)
    parser_make_viral_peptides.add_argument('-o', '--output', required=True)


    # Subparser for getting viral metadata
    parser_add_viral_metadata = subparsers.add_parser('add-viral-metadata',
                                                   help="Add viral metadata.")
    parser_add_viral_metadata.add_argument('-b', '--binding-affinities',
                                        help="Binding affinities (netMHCpan format).",
                                        required=True)
    parser_add_viral_metadata.add_argument('-q', '--viral-quants',
                                        help="Viral counts file from VirDetect.",
                                        required=True)
    parser_add_viral_metadata.add_argument('-r', '--viral-cds-ref',
                                        help="Viral CDS reference FASTA (used for VirDetect).",
                                        required=True)
    parser_add_viral_metadata.add_argument('-o', '--output', required=True)


    # Subparser for making fusion peptides
    parser_make_fusion_peptides = subparsers.add_parser('make-fusion-peptides',
                                                        help="Make fusion peptides FASTA.")
    parser_make_fusion_peptides.add_argument('-f', '--fusions',
                                            help="Predicted fusions (STARFusion format).",
                                            required=True)
    parser_make_fusion_peptides.add_argument('-o', '--output',
                                             required=True)

    # Subparser for adding fusion metadata
    parser_add_fusion_metadata = subparsers.add_parser('add-fusion-metadata',
                                                      help="Add metadata to fusion binding data.")
    parser_add_fusion_metadata.add_argument('-b', '--binding-affinities',
                                           help="Binding affinity data (netMHCpan format).",
                                           required=True)
    parser_add_fusion_metadata.add_argument('-f', '--fusions',
                                           help="Predicted fusions (STARFusion format).",
                                            required=True)
    parser_add_fusion_metadata.add_argument('-o', '--output',
                                            required=True)


    return parser.parse_args()


def load_tx_aas(args):
    """
    """
    tx_to_aa = {}
    for seq_record in SeqIO.parse(args.tx_aa_fasta, "fasta"):
        #tx = re.search('transcript:\S*', seq_record.description).group(0).split(':')[1].split('.')[0]
        tx = re.search('transcript:\S*', seq_record.description).group(0).split(':')[1]
        tx_to_aa[tx] = seq_record.seq
    return tx_to_aa


def load_tx_cds(args):
    """
    """
    tx_to_cds = {}
    for seq_record in SeqIO.parse(args.tx_cds_fasta, "fasta"):
        #tx = seq_record.id.split('.')[0]
        tx = seq_record.id
        tx_to_cds[tx] = seq_record.seq
    return tx_to_cds


def extract_potential_somatic_vars(args):
    """
    This will need to be cleaned heavily, but good for now.
    """
    filtered_records = []
    #vcf_reader = vcf.Reader(open(args.vcf), 'r', compressed=True)
    vcf_reader = vcf.Reader(open(args.vcf), 'r')
    for record in vcf_reader:
        print(record)
        possible_variants = [x for x in record.INFO['ANN']]
        for possible_variant in possible_variants:
            split_possible = possible_variant.split('|')
            if split_possible[1] == 'missense_variant' and split_possible[13] and split_possible[-1] == '':
                #transcript = split_possible[6].split('.')[0]
                transcript = split_possible[6]
                change = split_possible[10].lstrip('p.')
                #print(change)
                change = re.split('(\d+)', change)
                #print(change)
                pre = change[0]
                post = change[2]
                pos = change[1]
                pre1 = seq1(pre)
                post1 = seq1(post)
                chr_pos = record
                #print("{} {} {} {} {} {} {}".format(transcript, change, pre, pre1, post, post1, pos))
                pos_total = split_possible[13]
                print(pos_total)
                snapshot = pos_total.split('/')[0]
                tlen = pos_total.split('/')[1]
                if pos != snapshot:
                    sys.exit(1)
                filtered_records.append([transcript, pos, tlen, pre1, post1, chr_pos])

    return filtered_records

def extract_potential_somatic_nuc_changes(args):
    """
    This will need to be cleaned heavily, but good for now.
    """
    filtered_records = []
    #vcf_reader = vcf.Reader(open(args.vcf), 'r', compressed=True)
    vcf_reader = vcf.Reader(open(args.vcf), 'r')
    for record in vcf_reader:
        print(record)
        possible_variants = [x for x in record.INFO['ANN']]
        for possible_variant in possible_variants:
            split_possible = possible_variant.split('|')
            if split_possible[1] == 'missense_variant' and split_possible[13] and split_possible[-1] == '':
                #transcript = split_possible[6].split('.')[0]
                transcript = split_possible[6]
                change = split_possible[9].lstrip('c.')
                print(change)
                pos= re.split('(\d+)', change)[1]
                change = re.split('(\d+)', change)[2].split('>')
                print(change)
                pre = change[0]
                post = change[1]
                chr_pos = record
                print("{} {} {} {}".format(transcript, change, pre,post, pos))
                pos_total = split_possible[12]
                codon_pos = int((int(split_possible[13].split('/')[0])*3) - 2)
                snapshot = pos_total.split('/')[0]
                tlen = pos_total.split('/')[1]
                if pos != snapshot:
                    sys.exit(1)
                filtered_records.append([transcript, pos, tlen, pre, post, codon_pos, chr_pos])

    return filtered_records

def extract_potential_somatic_indels(args):
    """
    This will need to be cleaned heavily, but good for now.

    Focusing on deletions now, need to find a good example of a conservative_inframe_insertion
    """
    filtered_records = []
    #vcf_reader = vcf.Reader(open(args.vcf), 'r', compressed=True)
    vcf_reader = vcf.Reader(open(args.vcf), 'r')
    for record in vcf_reader:
        possible_variants = [x for x in record.INFO['ANN']]
        for possible_variant in possible_variants:
            print(possible_variant)
#            if re.search("\|conserative_inframe_insertion\|", possible_variant) or re.search("\|conserative_inframe_deletion\|", possible_variant):
            split_possible = possible_variant.split('|')
            if re.search('^conservative_inframe_[a-z]+$', split_possible[1]) and split_possible[13]:# and split_possible[-1] == '':
                print(split_possible)
#                transcript = split_possible[6].split('.')[0]
                transcript = split_possible[6]
                change = split_possible[10].lstrip('p.')
                print(change)
#                change = re.split('(\d+)', change)
#                print(change)
#                pre = change[0]
#                post = change[2]
#                pos = change[1]
#                pre1 = seq1(pre)
#                post1 = seq1(post)
#                chr_pos = record
#                print(chr_pos)
#                print("{} {} {} {} {} {} {}".format(transcript, change, pre, pre1, post, post1, pos))
                print("{} {}".format(transcript, change))
#                pos_total = split_possible[13]
#                snapshot = pos_total.split('/')[0]
#                tlen = pos_total.split('/')[1]
#                if pos != snapshot:
#                    sys.exit(1)
                filtered_records.append([transcript, change, record])

    return filtered_records

def make_peptides(args):
    """
    """
    tx_to_aa = load_tx_aas(args)
    filtered_records = extract_potential_somatic_vars(args)
    emitted_peptides = {}
    wildtype_peptides = {}
    for record in filtered_records:
        print(record)
        if record[0] in tx_to_aa.keys():
            print("In the dict!")
            tx_ref_seq = list(tx_to_aa[record[0]])
            print(tx_ref_seq)
            print(len(tx_ref_seq))
            #print(record[2])
            if len(tx_ref_seq) != int(record[2]):
                print("transcript {} shows different lengths between peptide fasta ({}) and snpeff! ({})".format(record[0], record[2], len(tx_ref_seq)))
            mut_seq = tx_ref_seq[:]
            pos = int(record[1]) - 1
            #print(mut_seq)
            #print(pos)
            #print(record[3])
            if mut_seq[pos] != record[3]:
                print("Reference amino acid doesn't match! Something has gone horribly wrong.")
            else:
                mut_seq[pos] = record[4]
                ref_peptide = ''.join(tx_ref_seq[max(pos-8+1,0):min(pos+8, len(mut_seq))])
                mut_peptide = ''.join(mut_seq[max(pos-8+1, 0):min(pos+8, len(mut_seq))])
                header_mut_peptide = ''.join(mut_seq[max(pos-13, 0):min(pos+14, len(mut_seq))])
                print("{}\n{}\n{}".format(record[0],ref_peptide, mut_peptide))
                var_md5sum = hashlib.md5("{}".format(':'.join([str(record[-1].CHROM), str(record[-1].POS), str(record[0]), str(record[-1].REF), str(record[-1].ALT)]))).hexdigest()[:16]
                emitted_peptides["{} {}:{} {} {} {} {}".format(var_md5sum, record[-1].CHROM, record[-1].POS, record[0], record[-1].REF, record[-1].ALT, header_mut_peptide)] = mut_peptide
                wildtype_peptides["{} {}:{} {} {} {}".format(var_md5sum, record[-1].CHROM, record[-1].POS, record[0], record[-1].REF, record[-1].ALT)] = ref_peptide

    with open(args.mt_output, 'w') as ofo:
        for k, v in emitted_peptides.items():
            ofo.write('>{}\n{}\n'.format(k, v))
    with open(args.wt_output, 'w') as ofo:
        for k, v in wildtype_peptides.items():
            ofo.write('>{}\n{}\n'.format(k, v))


def make_genomic_context(args):
    """
    """
    tx_to_cds = load_tx_cds(args)
    filtered_records = extract_potential_somatic_nuc_changes(args)
    emitted_nucs = {}
    for record in filtered_records:
        print(record)
        if record[0] in tx_to_cds.keys():
            print("In the dict!")
            tx_ref_seq = list(tx_to_cds[record[0]])
#            print(tx_ref_seq)
            print(len(tx_ref_seq))
#            #print(record[2])
            if len(tx_ref_seq) != int(record[2]):
                print("transcript {} shows different lengths between peptide fasta ({}) and snpeff! ({})".format(record[0], record[2], len(tx_ref_seq)))
            mut_seq = tx_ref_seq[:]
            pos = int(record[1]) - 1
            if mut_seq[pos] != record[3]:
                print("Reference amino acid doesn't match! Something has gone horribly wrong.")
            else:
                # Have to be careful here... we want to start at the first base of the affected codon
                mut_seq[pos] = record[4]
                codon_pos = record[5] - 1
                ref_nuc = ''.join(tx_ref_seq[max(codon_pos-args.length-1, 0):min(codon_pos+args.length, len(tx_ref_seq)) + 3])
                mut_nuc = ''.join(mut_seq[max(codon_pos-args.length, 0):min(codon_pos+args.length, len(mut_seq)) + 3])
                print("{}\n{}\n{}".format(record[0],ref_nuc, mut_nuc))
                print(':'.join([str(record[-1].CHROM), str(record[-1].POS), str(record[0]), str(record[-1].REF), str(record[-1].ALT)]))
                var_md5sum = hashlib.md5("{}".format(':'.join([str(record[-1].CHROM), str(record[-1].POS), str(record[0]), str(record[-1].REF), str(record[-1].ALT)])).encode('utf-8')).hexdigest()[:16]
                emitted_nucs["{}\t{}:{}\t{}\t{}\t{}".format(var_md5sum, record[-1].CHROM, record[-1].POS, record[0], record[-1].REF, record[-1].ALT)] = mut_nuc

    with open(args.output, 'w') as ofo:
        for k, v in emitted_nucs.items():
            ofo.write('>{}\t{}\n'.format(k, v))

def make_indel_peptides(args):
    """
    """
    tx_to_aa = load_tx_aas(args)
#    print(tx_to_aa['ENST00000369529.1'])
    filtered_records = extract_potential_somatic_indels(args)
    emitted_peptides = {}
    for record in filtered_records:
        print(record)
        if record[0] in tx_to_aa.keys():
            print("In the dict!")
            tx_ref_seq = list(tx_to_aa[record[0]])
            #print(tx_ref_seq)
            print(len(tx_ref_seq))
#            if len(tx_ref_seq) != int(record[2]):
#                print("transcript {} shows different lengths between peptide fasta ({}) and snpeff! ({})".format(record[0], record[2], len(tx_ref_seq)))
            mut_seq = tx_ref_seq[:]
            if re.search("del$", record[1]):
                del_rec = record[1].strip('del')
                start_pos = int(del_rec.split('_')[0][3:]) - 1
                start_aa = seq1(del_rec.split('_')[0][:3])
                stop_pos = int(del_rec.split('_')[1][3:]) - 1
                stop_aa = seq1(del_rec.split('_')[1][:3])
                print("{} {}".format(start_pos, start_aa))
                print("{} {}".format(stop_pos, stop_aa))
                if mut_seq[start_pos] != start_aa or mut_seq[stop_pos] != stop_aa:
                    print("Reference amino acid doesn't match! Something has gone horribly wrong.")
                else:
                    #mut_seq[pos] = record[4]
                    new_mut_seq = ''.join(mut_seq[:start_pos] + mut_seq[stop_pos+1:])
#                    print(new_mut_seq)
                    ref_peptide = ''.join(tx_ref_seq[max(start_pos-8,0):min(stop_pos+8, len(mut_seq))])
                    mut_peptide = ''.join(new_mut_seq[max(start_pos-8, 0):min(start_pos+8, len(new_mut_seq))])
                    print("{}\n{}\n{}".format(record[0],ref_peptide, mut_peptide))
                  #  emitted_peptides["{} {} {} {} {}".format(record[0], record[-1].CHROM, record[-1].POS, record[-1].REF, record[-1].ALT)] = mut_peptide
                    var_md5sum = hashlib.md5("{}".format(':'.join([str(record[-1].CHROM), str(record[-1].POS), str(record[0]), str(record[-1].REF), str(record[-1].ALT)]))).hexdigest()[:16]
                    emitted_peptides["{} {}:{} {} {} {}".format(var_md5sum, record[-1].CHROM, record[-1].POS, record[0], record[-1].REF, record[-1].ALT)] = mut_peptide

    with open(args.output, 'w') as ofo:
        for k, v in emitted_peptides.items():
            ofo.write('>{}\n{}\n'.format(k, v))


def rna_covered_variants(args):
    """
    """
    filtered_chrom_pos = []
    filtered_vcf = []
    rna_bam = pysam.AlignmentFile(args.rna_bam, "rb")
    vcf_reader = vcf.Reader(filename=args.vcf, compressed=False)
    for record in vcf_reader:
        possible_variants = [x for x in record.INFO['ANN']]
        for possible_variant_idx, possible_variant in enumerate(possible_variants):
            split_possible = possible_variant.split('|')
            #Currently filtering for only missense. Need to fix this.
            if split_possible[1] == 'missense_variant' and split_possible[13]:# and split_possible[-1] == '':
                transcript = split_possible[6].partition('.')[0]
                var_chr_pos = "{}:{}".format(record.CHROM, record.POS)
                coverage = len([str(i) for i in pileup_truncated(rna_bam, record.CHROM, record.POS - 1, record.POS)])
                if coverage > 0:
                    pos, total_depth, var_depth = [[i.reference_pos + 1, i.get_num_aligned(), ''.join(i.get_query_sequences()).upper().count(str(record.ALT[0]))] for i in pileup_truncated(rna_bam, record.CHROM, record.POS - 1, record.POS)][0]
                    if int(var_depth) > int(args.required_coverage):
                        #possible_variant + "|VAR_READ_DEPTH={}".format(var_depth)
                        filtered_chrom_pos.append([record.CHROM, record.POS])
            elif re.search('conservative_inframe_deletion', split_possible[1]):
                print(possible_variant)
                transcript = split_possible[6].partition('.')[0]
                var_chr_pos = "{}:{}".format(record.CHROM, record.POS)
                print(record.ALT[0])
                reads = [i.get_overlap(record.POS - 1, record.POS + len(record.REF) - 1) for i in rna_bam.fetch(record.CHROM, record.POS - 1, record.POS + len(record.REF) - 1)]
                overlap = [i.get_overlap(record.POS - 1, record.POS + len(record.REF) - 1) for i in rna_bam.fetch(record.CHROM, record.POS - 1, record.POS + len(record.REF) - 1)]
                total_depth = len(reads)
                # This isn't specific enough, need to be searching for specific start and length.
#                var_depth = [x for x in reads if not(re.search(str(record.REF), x))].count(str(len(record.ALT[0])))
#                ref_depth = [x for x in reads if re.search(str(record.REF), x)].count(str(len(record.REF)))
                var_depth = overlap.count(len(record.ALT[0]))
                ref_depth = overlap.count(len(record.REF))
                print(overlap)
                if int(var_depth) > int(args.required_coverage):
                    #possible_variant + "|VAR_READ_DEPTH={}".format(var_depth)
                    filtered_chrom_pos.append([record.CHROM, record.POS])
            elif re.search('conservative_inframe_insertion', split_possible[1]):
                print(possible_variant)
                transcript = split_possible[6].partition('.')[0]
                var_chr_pos = "{}:{}".format(record.CHROM, record.POS)
                print(record.ALT[0])
                reads = [i.get_forward_sequence() for i in rna_bam.fetch(record.CHROM, record.POS - 1, record.POS + len(record.REF) - 1)]
                print(reads)
                var_depth = len([x for x in reads if re.search(str(record.ALT[0]), x)])
                ref_depth = len([x for x in reads if not(re.search(str(record.ALT[0]), x))])
                print("{}\t{}\t{}".format(len(reads), ref_depth, var_depth))
                if int(var_depth) > int(args.required_coverage):
                    #possible_variant + "|VAR_READ_DEPTH={}".format(var_depth)
                    filtered_chrom_pos.append([record.CHROM, record.POS])

    vcf_reader = vcf.Reader(filename=args.vcf, compressed=False)
    for record in vcf_reader:
        if [record.CHROM, record.POS] in filtered_chrom_pos:
            filtered_vcf.append(record)

    print(len(filtered_vcf))

    vcf_writer = vcf.Writer(open(args.output, 'w'), vcf_reader)
    for filtered_record in filtered_vcf:
        vcf_writer.write_record(filtered_record)



def pileup_truncated(bam,contig, start, stop):
    """
    Taken from: https://github.com/pysam-developers/pysam/issues/851#issuecomment-585990105
    Obtain Pysam columns only at selected region
    """
    has_coord, rtid, rstart, rstop = bam.parse_region(contig, start, stop)
    for i in IteratorColumnRegion(bam, tid=rtid, start=rstart, stop=rstop,truncate=True):
        yield i
#    yield from IteratorColumnRegion(bam,
#                                    tid=rtid,
#                                    start=rstart,
#                                    stop=rstop,truncate=True)

def expressed_variants(args):
    """
    """
    tx_abundances = load_tx_abundances(args)
    tx_threshold = 0
    if not(args.abundance_threshold):
        tx_threshold = get_tx_threshold(args, tx_abundances)
    else:
        tx_threshold = float(args.abundance_threshold)
    expressed_txids = get_expressed_txs(args, tx_threshold)
    print("tx_thresshold: {}".format(tx_threshold))
    print("# of expressed transcripts: {}".format(len(expressed_txids)))
    print("Some expressed transcripts: {}".format(expressed_txids[:10]))
    filtered_records = filter_vcf_by_expression(args, expressed_txids)
    print("# of filtered transcripts: {}".format(len(filtered_records)))
    write_expressed_vcf(args, filtered_records)

def load_tx_abundances(args):
    """
    """
    count_column = ''
    counts = np.array([])
    print("okay!")
    with open(args.abundances) as fo:
        header = []
        cfo = csv.reader(fo, delimiter='\t')
        for line_idx, line in enumerate(cfo):
            if line_idx == 0:
                header = line
                for i,j in enumerate(header):
                    if j ==  args.metric:
                        count_column = i
            elif line_idx:
                count = np.log2(float(line[count_column]) + 1)
                if args.exclude_zeros:
                    if count> 0:
                       counts = np.append(counts, count)
                else:
                    counts = np.append(counts, count)
    print("Max count: {}".format(np.max(counts)))

    return counts

def get_tx_threshold(args, counts):
    """
    """
    threshold = np.percentile(counts, int(args.percentile))
    print(args.percentile)
    print(threshold)
    return threshold

def get_expressed_txs(args, threshold):
    """
    """
    count_column = ''
    txid_column = ''
    expressed_txids = []
    print("okay!")
    with open(args.abundances) as fo:
        header = []
        cfo = csv.reader(fo, delimiter='\t')
        for line_idx, line in enumerate(cfo):
            if line_idx == 0:
                header = line
                for i,j in enumerate(header):
                    if j ==  args.metric:
                        count_column = i
                    elif j == 'Name':
                        txid_column = i
            elif line_idx:
                count = float(line[count_column])
                print("{} {}".format(count, threshold))
                if float(count) > float(threshold):
                    expressed_txids.append(line[txid_column].split('.')[0])
    return expressed_txids

def filter_vcf_by_expression(args, expressed_txids):
    """
    """
    filtered_records = []
    vcf_reader = vcf.Reader(filename=args.vcf, compressed=False)
    for record in vcf_reader:
        possible_variants = [x for x in record.INFO['ANN']]
        for possible_variant in possible_variants:
            split_possible = possible_variant.split('|')
            if split_possible[1] == 'missense_variant' and split_possible[13] and split_possible[-1] == '':
                #Having to strip version off transcript for this.
                transcript = split_possible[6].partition('.')[0]
                print("Transcript: {}".format(transcript))
                if transcript in expressed_txids:
                    filtered_records.append(record)
    return filtered_records

def write_expressed_vcf(args, filtered_records):
    """
    """
    vcf_reader = vcf.Reader(open(args.vcf), 'r')
    vcf_writer = vcf.Writer(open(args.output, 'w'), vcf_reader)
    for filtered_record in filtered_records:
        vcf_writer.write_record(filtered_record)

def write_isolated_vcf(args, filtered_records):
    """
    """
    vcf_reader = vcf.Reader(filename=args.somatic_vcf, compressed=True)
    vcf_writer = vcf.Writer(open(args.output, 'w'), vcf_reader)
    for filtered_record in filtered_records:
        vcf_writer.write_record(filtered_record)


def isolated_variants(args):
    """
    """
    cand_vars = get_candidate_variants(args)
    print("# of candidate variants: {}".format(len(cand_vars)))
    isolated_vars = get_all_isolated_variants(args)
    print("# of isolated variants: {}".format(len(isolated_vars)))
    isolated_cand_vars = []
    for record in cand_vars:
        if record in isolated_vars:
            isolated_cand_vars.append(record)
    print("# of isolated candidate variants: {}".format(len(isolated_cand_vars)))
    write_isolated_vcf(args, isolated_cand_vars)


def get_candidate_variants(args):
    """
    """
    cand_vars = []
#    vcf_reader = vcf.Reader(filename=args.candidate_vcf, compressed=True)
    vcf_reader = vcf.Reader(filename=args.somatic_vcf, compressed=True)
    for record in vcf_reader:
        cand_vars.append(record)
    return cand_vars

def get_all_isolated_variants(args):
    """Assumes variants are sorted for now.
    """
    isolated_vars = []
    germline_vcf_reader = vcf.Reader(filename=args.germline_vcf, compressed=True)
    somatic_vcf_reader = vcf.Reader(filename=args.somatic_vcf, compressed=True)
    all_records = []
    for record_idx, record in enumerate(germline_vcf_reader):
        all_records.append(record)
    for record_idx, record in enumerate(somatic_vcf_reader):
        all_records.append(record)
    print("# of total variants: {}".format(len(all_records)))
    for record_idx, record in enumerate(all_records[:-1]):
        upstream_record = all_records[record_idx - 1]
        downstream_record = all_records[record_idx + 1]
        isolated_variant = False

        upstream_isolated = False
        downstream_isolated = False

        if upstream_record.CHROM == record.CHROM:
            if np.absolute(upstream_record.POS - record.POS) > int(args.proximity):
                upstream_isolated = True
        else:
            upstream_isolated = True

        if downstream_record.CHROM == record.CHROM:
            if np.absolute(downstream_record.POS - record.POS) > int(args.proximity):
                downstream_isolated = True
        else:
            downstream_isolated = True


        if upstream_isolated and downstream_isolated:
            isolated_vars.append(record)

    return isolated_vars

def calculate_agretopicity(args):
    """
    """

    wt_nms = {}
    wt_peptides = {}
    with open(args.wt_fasta) as wto:
        for line in wto.readlines():
            line = line.rstrip('\n').split(' ')
            line = [x for x in line if x]
            if len(line) > 14 and line[0] not in ['Pos', 'Protein']:
#                print(line)
                k = "{}_{}_{}_{}".format(line[1], line[10], line[0], len(line[9]))
#                print(k)
                wt_nms[k] = line[15]
                wt_peptides[k] = line[9]
    #print(wt_nms)

    mts_w_agreto = []

    header = ''

    with open(args.mt_fasta) as mto:
        for line in mto.readlines():
            line = line.rstrip('\n').split(' ')
            line = [x for x in line if x]
            if len(line) > 14 and line[0] not in ['Pos', 'Protein']:
                m = "{}_{}_{}_{}".format(line[1], line[10], line[0], len(line[9]))
                result = line[:16]
                if m in wt_nms.keys():
                    print("{}\t{}\t{}\t{}".format(m, line[15], wt_nms[m], float(line[15])/float(wt_nms[m])))
                    result.extend([wt_nms[m], wt_peptides[m], str(float(line[15])/float(wt_nms[m]))])
                else:
                    result.extend(["NA", "NA", "NA"])
                mts_w_agreto.append(','.join(result))
            elif len(line) > 14 and line[0] == 'Pos':
                print(header)
                line.extend(['Wildtype Aff(nM)', 'Wildtype Peptide', 'Agretopocity'])
                header = ','.join(line)

    with open(args.output, 'w') as outf:
        outf.write("{}\n".format(header))
        for line in mts_w_agreto:
            outf.write("{}\n".format(line))



def create_lens_report(args):
    """
    """
    pass


def add_snv_metadata(args):
    """
    """
    print(args)
    checksum_to_meta_map = {}
    with open(args.mutant_peptides) as mto:
        for line in mto.readlines():
            if line.startswith('>'):
                line = line.rstrip().split(' ')
                print(line)
                checksum = line[0].lstrip('>')[:-1]
                print(checksum)
                var_pos = line[1]
                tx_id = line[2].split('.')[0]
                ref = line[3]
                alt = line[4]
                aa_context = line[5]
                checksum_to_meta_map[checksum] = {}
                checksum_to_meta_map[checksum]["var_pos"] = var_pos
                checksum_to_meta_map[checksum]["tx_id"] = tx_id
                checksum_to_meta_map[checksum]["ref"] = ref
                checksum_to_meta_map[checksum]["alt"] = alt
                checksum_to_meta_map[checksum]["aa_context"] = aa_context

    with open(args.mutant_nucs) as mno:
        for line in mno.readlines():
            if line.startswith('>'):
                line = line.rstrip().split('\t')
                checksum = line[0].lstrip('>')[:-1]
                nuc_context = line[-1]
                checksum_to_meta_map[checksum]["nuc_context"] = nuc_context

    tx_to_tpm = {}
    with open(args.quants) as qfo:
        tpm_col_idx = ''
        tx_col_idx = ''
        for line_idx, line in enumerate(qfo.readlines()):
            if line_idx == 0:
                print(line)
                tpm_col_idx = line.split('\t').index('TPM')
                tx_col_idx = line.split('\t').index('Name')
            else:
                line = line.split('\t')
                tx_to_tpm[line[tx_col_idx].split('.')[0]] = line[tpm_col_idx]

    tx_to_gene = {}
    with open(args.gtf) as gtfo:
        for line in gtfo.readlines():
            line = line.split('\t')
            if len(line) > 3 and line[2] == 'transcript':
                gene_name = str(line).split('gene_name "')[1].split('"')[0]
                tx_id = str(line).split('transcript_id "')[1].split('"')[0]
                tx_to_gene[tx_id] = gene_name

    print(tx_to_gene.keys())

    var_to_ccf = {}
    with open(args.cancer_cell_fraction) as ccfo:
        for line_idx, line in enumerate(ccfo.readlines()):
            if line_idx != 0:
                line = line.rstrip().replace("b'", '').replace("'", '').split()
                var = line[0]
                ccf = line[3]
                var_to_ccf[var] = ccf

    print(var_to_ccf.keys()[:10])

    header = []

    new_lines = []
    with open(args.netmhcpan) as mno:
        for line_idx, line in enumerate(mno.readlines()):
            line = line.rstrip()
            line = line.split(',')
            if line_idx == 0:
                header = line
                try:
                    header.remove('BindLevel')
                except:
                    pass
                header.extend(['gene_name', 'tx_id', 'variant_position', 'reference_allele', 'alternate_allele', 'tpm', 'ccf', 'amino_acid_context', 'nucleotide_context'])
            if len(line) < 16:
                continue
            if line[10] not in checksum_to_meta_map.keys():
                continue
            csum = line[10]
            tx_id = checksum_to_meta_map[csum]['tx_id']
            gene_name = tx_to_gene[tx_id]
            var_pos = checksum_to_meta_map[csum]['var_pos']
            ref = checksum_to_meta_map[csum]['ref']
            alt = checksum_to_meta_map[csum]['alt'].replace('[','').replace(']','')
            tpm = tx_to_tpm[tx_id]
            aa_context = checksum_to_meta_map[csum]['aa_context']
            nuc_context = checksum_to_meta_map[csum]['nuc_context']
            translated_nuc = str(Seq(nuc_context).translate()).replace('*', '')
            if translated_nuc != aa_context:
                print("One of your antigen's nucleotide context does not match the amino acid context. Check this!")
                print(nuc_context)
                print(translated_nuc)
                print(aa_context)
                sys.exit()
            ccf = 'NA'
            if var_pos in var_to_ccf.keys():
                ccf = var_to_ccf[var_pos]
            if 'SB' in line:
                line.remove('<=')
                line.remove('SB')
            if 'WB' in line:
                line.remove('WB')
            line.extend([gene_name, tx_id, var_pos, ref, alt, tpm, ccf, aa_context, nuc_context])
            new_lines.append(line)
    with open(args.output, 'w') as ofo:
        ofo.write("{}\n".format('\t'.join(header)))
        for new_line in new_lines:
            ofo.write("{}\n".format('\t'.join(new_line)))

def add_indel_metadata(args):
    """
    """
    print(args)
    checksum_to_meta_map = {}
    with open(args.mutant_peptides) as mto:
        for line in mto.readlines():
            if line.startswith('>'):
                line = line.rstrip().split(' ')
                checksum = line[0].lstrip('>')[:-1]
                var_pos = line[1]
                tx_id = line[2].split('.')[0]
                ref = line[3]
                alt = line[4]
                checksum_to_meta_map[checksum] = {}
                checksum_to_meta_map[checksum]["var_pos"] = var_pos
                checksum_to_meta_map[checksum]["tx_id"] = tx_id
                checksum_to_meta_map[checksum]["ref"] = ref
                checksum_to_meta_map[checksum]["alt"] = alt

    tx_to_tpm = {}
    with open(args.quants) as qfo:
        tpm_col_idx = ''
        tx_col_idx = ''
        for line_idx, line in enumerate(qfo.readlines()):
            if line_idx == 0:
                tpm_col_idx = line.split('\t').index('TPM')
                tx_col_idx = line.split('\t').index('Name')
            else:
                line = line.split('\t')
                tx_to_tpm[line[tx_col_idx].split('.')[0]] = line[tpm_col_idx]

    tx_to_gene = {}
    with open(args.gtf) as gtfo:
        for line in gtfo.readlines():
            line = line.split('\t')
            if len(line) > 3 and line[2] == 'transcript':
                gene_name = str(line).split('gene_name "')[1].split('"')[0]
                tx_id = str(line).split('transcript_id "')[1].split('"')[0]
                tx_to_gene[tx_id] = gene_name

    var_to_ccf = {}
    with open(args.cancer_cell_fraction) as ccfo:
        for line_idx, line in enumerate(ccfo.readlines()):
            if line_idx != 0:
                line = line.rstrip().replace("b'", '').replace("'", '').split()
                var = line[0]
                ccf = line[3]
                var_to_ccf[var] = ccf

    header = []

    new_lines = []
    with open(args.netmhcpan) as mno:
        for line_idx, line in enumerate(mno.readlines()):
            print(line)
            line = line.rstrip()
            line = line.split()
            if not(line):
                continue
            if line[0] == 'Pos':
                header = line
                header.remove('BindLevel')
                header.extend(['gene_name', 'tx_id', 'variant_position', 'reference_allele', 'alternate allele', 'TPM', 'CCF'])
            if len(line) < 16:
                continue
            if line[10] not in checksum_to_meta_map.keys():
                continue
            csum = line[10]
            tx_id = checksum_to_meta_map[csum]['tx_id']
            gene_name = tx_to_gene[tx_id]
            var_pos = checksum_to_meta_map[csum]['var_pos']
            ref = checksum_to_meta_map[csum]['ref']
            alt = checksum_to_meta_map[csum]['alt'].replace('[','').replace(']','')
            tpm = tx_to_tpm[tx_id]
            ccf = 'NA'
            if var_pos in var_to_ccf.keys():
                ccf = var_to_ccf[var_pos]
            if 'SB' in line:
                line.remove('<=')
                line.remove('SB')
            if 'WB' in line:
                line.remove('WB')
            line.extend([gene_name, tx_id, var_pos, ref, alt, tpm, ccf])
            new_lines.append(line)


    with open(args.output, 'w') as ofo:
        ofo.write("{}\n".format('\t'.join(header)))
        for new_line in new_lines:
            ofo.write("{}\n".format('\t'.join(new_line)))




def make_pyclone_vi_inputs(args):
    """
    This requires an isec vcf, a proper vcf (with depth info), and sequenza results info.
    """
    # This dictionary will be populated with normal and tumor depth information
    # from the proper VCFs. It'll initially have keys populated by the isec VCF.
    # vars[var]['var_depth'] =, vars[var]['ref_depth'] =, vars[var]['cn'] =
    vars = {}

    cellularity = 0
    with open(args.sequenza_solutions) as sso:
        for line_idx, line in enumerate(sso.readlines()):
            if line_idx == 1:
                line = line.split('\t')
                cellularity = line[0]


    with open(args.candidate_vcf) as cvo:
        for line in cvo.readlines():
            if line.startswith('#'):
                pass
            else:
                line = line.rstrip().split('\t')
#                if len(line[3]) == 1 and len(line[4]) == 1:
                if len(line[3]) > 0 and len(line[4]) > 0:
                    line = [i for i in line if i != '.']
                    vars['{}'.format('_'.join(line[:2]))] = {}
                    vars['{}'.format('_'.join(line[:2]))]['ref_depth'] = 0
                    vars['{}'.format('_'.join(line[:2]))]['alt_depth'] = 0
                    vars['{}'.format('_'.join(line[:2]))]['major_cn'] = 0
                    vars['{}'.format('_'.join(line[:2]))]['minor_cn'] = 0
    print("Candidate variant count: {}".format(len(vars)))
    #Note: not all variants may not be in depth VCF. This can be fixed later.


    mutect_vcf_reader = vcf.Reader(open(args.mutect_vcf), compressed=False)
    depths = {}
    for record in mutect_vcf_reader:
        call = record.genotype('tum')
        var_key = "{}_{}".format(record.CHROM, record.POS)
        depths[var_key] = {}
        depths[var_key]["ref_depth"] = call.data.AD[0]
        depths[var_key]["alt_depth"] = call.data.AD[1]

    copy_no = {}
    with open(args.sequenza_segments) as sfo:
        for line_idx, line in enumerate(sfo.readlines()):
            if line_idx == 0:
                continue
            line = line.split('\t')
            chr, start, stop, maj_count, min_count = [line[0], line[1], line[2], line[10], line[11]]
            print("{}".format(','.join([chr, start, stop, maj_count, min_count])))
            if chr not in copy_no.keys():
                copy_no[chr] = {}
            if "{}-{}".format(start, stop) not in copy_no[chr].keys():
                copy_no[chr]["{}-{}".format(start, stop)] = {}
            copy_no[chr]["{}-{}".format(start, stop)]["maj_count"] = maj_count
            copy_no[chr]["{}-{}".format(start, stop)]["min_count"] = min_count

    easily_parsable = []

    for var in vars.keys():
        print(var)
        chr, pos = var.split('_')
        if var in depths.keys():
            easily_parsable.append([chr, pos])

#    print(easily_parsable)

    pyclone_inp = []

    for chr in copy_no.keys():
        print(chr)
        chr_vars = sorted([x for x in easily_parsable if x[0] == chr], key=lambda x: x[1])
        print(len(chr_vars))
        for segment in sorted(copy_no[chr].keys(), key=lambda segment: int(segment.split('-')[0])):
            print("SEGMENT: {}".format(segment))
            start, stop = segment.split('-')
            segment_vars = sorted([x for x in chr_vars if int(x[1]) > int(start) and int(x[1]) < int(stop)], key=lambda x: int(x[1]))
            print(segment_vars)
            for segment_var in segment_vars:
                var = "{}_{}".format(segment_var[0], segment_var[1])
                ref_depth = str(depths[var]['ref_depth'])
                alt_depth = str(depths[var]['alt_depth'])
#                depths[var]['major_cn'] = str(copy_no[chr][segment]["maj_count"])
#                depths[var]['minor_cn'] = str(copy_no[chr][segment]["min_count"])
                ref_cn = str(copy_no[chr][segment]["maj_count"])
                alt_cn = str(copy_no[chr][segment]["min_count"])
                var = var.replace('_', ':')
                pyclone_inp.append("{}\n".format('\t'.join([var, args.samp_id, ref_depth, alt_depth, ref_cn, alt_cn, '2', '0.001', cellularity])))




#    print(pyclone_inp)
    with open(args.output, 'w') as ofo:
        ofo.write("mutation_id\tsample_id\tref_counts\talt_counts\tmajor_cn\tminor_cn\tnormal_cn\terror_rate\ttumour_content\n")
        for i in pyclone_inp:
            ofo.write(i)


def expressed_hervs(args):
    """
    """
    tx_abundances = load_tx_abundances(args)
    tx_threshold = 0
    if not(args.abundance_threshold):
        print("No abundance threshold!")
        tx_threshold = get_tx_threshold(args, tx_abundances)
    else:
        print("Abundance threshold: {}".format(args.abundance_threshold))
#        tx = re.search('transcript:\S*', seq_record.description).group(0).split(':')[1]
        tx_threshold = float(args.abundance_threshold)
    expressed_hervs = get_expressed_txs(args, tx_threshold)
#    print("tx_thresshold: {}".format(tx_threshold))
#    print("# of expressed transcripts: {}".format(len(expressed_txids)))
#    print("Some expressed transcripts: {}".format(expressed_txids[:10]))
#    print("# of filtered transcripts: {}".format(len(filtered_records)))
    with open(args.output, 'w') as ofo:
        for herv in expressed_hervs:
            ofo.write("{}\n".format(herv))


def consolidate_multiqc_stats(args):
    """
    """

    sample_level_stats = {}

    stats_files = ["multiqc_fastqc.txt", "multiqc_picard_RnaSeqMetrics.txt", "multiqc_samtools_stats.txt"]

    for stats_file in stats_files:
        with open(os.path.join(os.getcwd(), 'multiqc_data', stats_file)) as sfo:
            metric_to_idx = {}
            for line_idx, line in enumerate(sfo.readlines()):
                line = line.rstrip().split('\t')
                if line_idx == 0:
                    for metric_idx, metric in enumerate(line[1:]):
                        metric_to_idx[metric] = metric_idx + 1
                else:
                    print(line)
                    sample = line[0]
                    print("Sample: {}".format(line[0]))
                    suffix = ''
                    if re.search("_[12]\.", sample) or re.search("_[12]$", sample):
                        print("Reformating sample name")
                        sample = '_'.join(line[0].split('_')[:-1])
                        suffix = ' (file {})'.format(line[metric_to_idx['Filename']].split('_')[-1])
                    print("New Sample: {}".format(sample))
                    print("Suffix Sample: {}".format(suffix))

                    if sample not in sample_level_stats.keys():
                        sample_level_stats[sample] = {}
                    for metric, metric_idx in metric_to_idx.items():
                        sample_level_stats[sample]["{}{}".format(metric, suffix)] = line[metric_idx]

    print(sample_level_stats)


    with open(args.output, 'w') as ofo:
        header = []
        for sample in sample_level_stats.keys():
            for metric in sample_level_stats[sample]:
                if metric not in header:
                    header.append(metric)

        print(header)
        header = sorted(header)

        header.insert(0, 'sample')

        ofo.write("{}\n".format('\t'.join(header)))

        for sample in sample_level_stats.keys():
            print(sample_level_stats[sample])
            sample_line = []
            sample_line.append(sample)
            for metric in header[1:]:
                if metric in sample_level_stats[sample].keys():
                    sample_line.append(sample_level_stats[sample][metric])
                else:
                    sample_line.append('NA')
            ofo.write("{}\n".format('\t'.join(sample_line)))

def make_herv_peptides(args):
    """
    """
    expressed_hervs = []
    expressed_hervs_seqs = {}
    expressed_hervs_aas = {}
    with open(args.expressed_hervs) as eho:
        for line in eho.readlines():
            line = line.rstrip()
            expressed_hervs.append(line)

    for seq_record in SeqIO.parse(args.herv_ref, "fasta"):
        if seq_record.description in expressed_hervs:
            expressed_hervs_seqs[seq_record.description] = seq_record.seq

    print(expressed_hervs_seqs)

    for i in [0, 1, 2]:
        for id, seq in expressed_hervs_seqs.items():
            aa_seq = seq[i:].translate(to_stop=False)
            expressed_hervs_aas["{}_{}".format(id, i)] = aa_seq


    with open(args.output, 'w') as ofo:
        for k, v in expressed_hervs_aas.items():
            if len(v) > 8:
                ofo.write(">{}\n{}\n".format(k, v))



def get_herv_metadata(args):
    """
    """
    #tx_abundances = load_tx_abundances(args)

    #print(tx_abundances)

    tx_to_tpm = {}
    with open(args.abundances) as qfo:
        tpm_col_idx = ''
        tx_col_idx = ''
        for line_idx, line in enumerate(qfo.readlines()):
            if line_idx == 0:
                print(line)
                tpm_col_idx = line.split('\t').index('TPM')
                tx_col_idx = line.split('\t').index('Name')
            else:
                line = line.split('\t')
                tx_to_tpm[line[tx_col_idx].split('.')[0][:15].replace(':','_')] = line[tpm_col_idx]

    print(tx_to_tpm)

    herv_expression = {}

    output_lines = []

    with open(args.netmhcpan) as mno:
        for line_idx, line in enumerate(mno.readlines()):
            line = line.rstrip()
            line = line.split()
            if len(line) >= 16 and line[0] == 'Pos':
                header = line
                try:
                    header.remove('BindLevel')
                except:
                    pass
                header.extend(['read_count'])
            if len(line) < 16 or line[0] in ['Protein', 'Pos']:
                continue
            else:
                print(line)
                tpm = tx_to_tpm[line[10]]
                line.extend([tpm])
                output_lines.append(line)

    with open(args.output, 'w') as ofo:
        ofo.write("{}\n".format('\t'.join(header)))
        for line in output_lines:
            ofo.write("{}\n".format('\t'.join(line)))

def expressed_self_genes(args):
    """
    """
    expressed_selfs = []
    gene_list = []
    with open(args.gene_list) as glo:
        for line in glo.readlines():
            gene_list.append(line.rstrip())
    tx_abundances = load_tx_abundances(args)
    tx_threhsold = 0
    if not(args.abundance_threshold):
        tx_threshold = get_tx_threshold(args, tx_abundances)
    else:
        tx_threshold = float(args.abundance_threshold)
    expressed_txids = get_expressed_txs(args, tx_threshold)
    print(expressed_txids[:10])
    #print("tx_thresshold: {}".format(tx_threshold))
    #print("# of expressed transcripts: {}".format(len(expressed_txids)))
    #print("Some expressed transcripts: {}".format(expressed_txids[:10]))
    #print("# of filtered transcripts: {}".format(len(filtered_records)))


    tx_to_gene = {}
    with open(args.gtf) as gtfo:
        for line in gtfo.readlines():
            line = line.split('\t')
            if len(line) > 3 and line[2] == 'transcript':
                gene_name = str(line).split('gene_name "')[1].split('"')[0]
                tx_id = str(line).split('transcript_id "')[1].split('"')[0].split('.')[0]
                tx_to_gene[tx_id] = gene_name



    for tx, gene in tx_to_gene.items():
        if gene in gene_list and tx in tx in expressed_txids:
            print(gene, tx)
            expressed_selfs.append("{}:{}".format(gene, tx))

    with open(args.output, 'w') as ofo:
        for expressed_self in expressed_selfs:
            ofo.write("{}\n".format(expressed_self))


def make_self_antigen_peptides(args):
    """
    """
    print("Okay!")
    expressed_selfs = []
    expressed_self_seqs = {}
    with open(args.expressed_selfs) as eso:
        for line in eso.readlines():
            expressed_selfs.append(line.rstrip())

    tx_to_aa = load_tx_aas(args)
    tx_to_aa_trunc = {k.partition('.')[0]:v for k, v in tx_to_aa.items()}




    # Need to incorporate the germline variants here.

    for expressed_self in expressed_selfs:
        tx = expressed_self.partition(':')[2]
        if tx in tx_to_aa_trunc.keys():
            aa_seq = tx_to_aa_trunc[tx]
            expressed_self_seqs[tx] = aa_seq

    relevant_germline_vars = {}
    #vcf_reader = vcf.Reader(open(args.vcf), 'r', compressed=True)
    vcf_reader = vcf.Reader(open(args.germline_vcf), 'r')
    for record in vcf_reader:
        #print(record)
        possible_variants = [x for x in record.INFO['ANN']]
        for possible_variant in possible_variants:
            split_possible = possible_variant.split('|')
            if split_possible[1] == 'missense_variant' and split_possible[13] and split_possible[-1] == '':
                transcript = split_possible[6].split('.')[0]
                #transcript = split_possible[6]
                #print(record)
                change = split_possible[10].lstrip('p.')
                change = re.split('(\d+)', change)
                pre = change[0]
                post = change[2]
                pos = change[1]
                pre1 = seq1(pre)
                post1 = seq1(post)
                chr_pos = record
                #print("{} {} {} {} {} {} {}".format(transcript, change, pre, pre1, post, post1, pos))
                pos_total = split_possible[13]
                #print(pos_total)
                snapshot = pos_total.split('/')[0]
                tlen = pos_total.split('/')[1]
                if pos != snapshot:
                    sys.exit(1)
                if transcript in expressed_self_seqs.keys():
                    print(record)
                    print(possible_variant)
                    if transcript not in relevant_germline_vars.keys():
                        relevant_germline_vars[transcript] = [[pos, tlen, pre1, post1, chr_pos, record.genotype('norm')['GT']]]
                    else:
                        relevant_germline_vars[transcript].append([pos, tlen, pre1, post1, chr_pos, record.genotype('norm')['GT']])

    tx_to_peps = {}
    for expressed_self, expressed_self_seq in expressed_self_seqs.items():
        print(expressed_self)
        if expressed_self in relevant_germline_vars.keys():
            relevant_vars = relevant_germline_vars[expressed_self]
            print("{}\n{}".format(expressed_self, relevant_germline_vars[expressed_self]))
            # Should probably be a named tuple
            seq_windows = []
            window = 21
            step = 8
            capture_index = 0
            print(expressed_self_seq)
            for pos_idx in range(len(expressed_self_seq)):
                if pos_idx == capture_index:
                    seq_windows.append([str(expressed_self_seq[pos_idx:(pos_idx + window)]), pos_idx+ 1, pos_idx + window + 1])
                    capture_index += step
            for seq_window in seq_windows:
                print("Seq window: {}".format(seq_window))
                seq_window_start = seq_window[1]
                seq_window_stop = seq_window[2]
                found_vars = 0
                for relevant_var in relevant_vars:
                    relevant_var_pos = relevant_var[0]
                    print("Relevant var pos: {}".format(relevant_var_pos))
                    if int(relevant_var_pos) > int(seq_window_start) and int(relevant_var_pos) < int(seq_window_stop):
                        found_vars = 1
                        print("Var is within window.")
                        print(seq_window)
                        print(relevant_var)
                        relative_var_pos = int(relevant_var_pos) - int(seq_window_start)
                        print(relative_var_pos)
                        new_seq = list(seq_window[0])
                        print(new_seq)
                        if new_seq[relative_var_pos] != relevant_var[2]:
                            print("The reference amino acid isn't what's expected.")
                        new_seq[relative_var_pos] = relevant_var[3]
                        if relevant_var[-1] in ['0/1', '0|1']:
                            print("Variant is heterozygote, grabbing both genotypes.")
                            tx_to_peps["{}_{}_het_alt0".format(expressed_self, seq_window_start)] = seq_window[0]
                            tx_to_peps["{}_{}_het_alt1".format(expressed_self, seq_window_start)] = ''.join(new_seq)
                        if relevant_var[-1] in ['1/1', '1|1']:
                            tx_to_peps["{}_{}_hom_alt".format(expressed_self, seq_window_start)] = ''.join(new_seq)
                        print("old seq: {}\nnew_seq: {}".format(seq_window[0], ''.join(new_seq)))
                if found_vars == 0:
                    print("{}\t{}".format(expressed_self, expressed_self_seq))
                    tx_to_peps["{}_{}".format(expressed_self, seq_window_start)] = seq_window[0]
        else:
            print("No variants")
            # If there are no applicable variants, simply emit the reference sequence.
            tx_to_peps[expressed_self] = expressed_self_seq

    with open(args.output, 'w') as ofo:
        for expressed_self, expressed_self_seq in sorted(tx_to_peps.items()):
            print(expressed_self)
            ofo.write(">{}\n{}\n".format(expressed_self, expressed_self_seq))


def add_self_antigen_metadata(args):
    """
    """
    tx_to_tpm = {}
    with open(args.quants) as qfo:
        tpm_col_idx = ''
        tx_col_idx = ''
        for line_idx, line in enumerate(qfo.readlines()):
            if line_idx == 0:
                print(line)
                tpm_col_idx = line.split('\t').index('TPM')
                tx_col_idx = line.split('\t').index('Name')
            else:
                line = line.split('\t')
                tx_to_tpm[line[tx_col_idx].split('.')[0]] = line[tpm_col_idx]

    tx_to_gene = {}
    with open(args.gtf) as gtfo:
        for line in gtfo.readlines():
            line = line.split('\t')
            if len(line) > 3 and line[2] == 'transcript':
                gene_name = str(line).split('gene_name "')[1].split('"')[0]
                tx_id = str(line).split('transcript_id "')[1].split('"')[0].partition('.')[0]
                tx_to_gene[tx_id] = gene_name

    output_lines = []

    with open(args.binding_affinities) as mno:
        for line_idx, line in enumerate(mno.readlines()):
            line = line.rstrip()
            line = line.split()
            if len(line) >= 16 and line[0] == 'Pos':
                header = line
                try:
                    header.remove('BindLevel')
                except:
                    pass
                header.extend(['gene_name', 'read_count'])
            if len(line) < 16 or line[0] in ['Protein', 'Pos']:
                continue
            else:
                print(line)
                gene_name = tx_to_gene[line[10]]
                tpm = tx_to_tpm[line[10]]
                line.extend([gene_name, tpm])
                output_lines.append(line)

    with open(args.output, 'w') as ofo:
        ofo.write("{}\n".format('\t'.join(header)))
        for line in output_lines:
            ofo.write("{}\n".format('\t'.join(line)))

def filter_virdetect_by_counts(args):
    """
    Need a check here to ensure the raw_counts and raw_ref lists contain the
    same number of elements.
    """
    raw_counts = []
    raw_refs = []
    with open(args.viral_quants) as vco:
        raw_counts = vco.readlines()[0].rstrip().split()[1:]
    with open(args.viral_ref) as vro:
        for line in vro.readlines():
            if line.startswith('>'):
                line = line.rstrip('\n').lstrip('>')
                raw_refs.append(line.split('|')[3])
    viral_counts = {raw_refs[i]:raw_counts[i] for i in range(len(raw_refs)) if int(raw_counts[i]) > int(args.min_threshold)}
    print(viral_counts)
    with open(args.output, 'w') as ofo:
        ofo.write("virus_ref\tread_count\n")
        for k, v in viral_counts.items():
            ofo.write("{}\t{}\n".format(k, v))


def make_viral_peptides(args):
    """
    """
    expressed_viruses = {}
    virus_to_cds = {}
    viral_peptide_seqs = {}
    expressed_viral_peptides = {}

    with open(args.expressed_viruses) as fo:
        for line_idx, line in enumerate(fo.readlines()):
            if line_idx != 0:
                line = line.rstrip().split('\t')
                expressed_viruses[line[0]] = line[1]

    viral_seqs = {}

    for seq_record in SeqIO.parse(args.sample_viral_ref, "fasta"):
        map_info = seq_record.description.split('|')[1].split(' ')[0]
        virus_id, unneeded, protein_id = map_info.partition('_cds_')
        protein_id = re.sub(r'_1$', '', protein_id)
        if virus_id not in virus_to_cds.keys():
            virus_to_cds[virus_id] = [protein_id]
        else:
            virus_to_cds[virus_id].append(protein_id)

        viral_seqs[protein_id] = seq_record.seq.translate()



    with open(args.sample_viral_ref) as fo:
        for line in fo.readlines():
            if line.startswith('>'):
                line = line.split('|')[1].split(' ')[0]
                virus_id, unneeded, protein_id = line.partition('_cds_')
                protein_id = re.sub(r'_1$', '', protein_id)
                if virus_id not in virus_to_cds.keys():
                    virus_to_cds[virus_id] = [protein_id]
                else:
                    virus_to_cds[virus_id].append(protein_id)

#    for seq_record in SeqIO.parse(args.viral_pep_ref, "fasta"):
#        peptide_id = seq_record.description.split(' ')[0]
#        viral_peptide_seqs[peptide_id] = seq_record.seq

    for expressed_virus in expressed_viruses.keys():
        if expressed_virus in virus_to_cds.keys():
            for protein in virus_to_cds[expressed_virus]:
                expressed_viral_peptides[protein] = viral_seqs[protein]
        else:
            print("Cannot find virus {} in cds data.".format(expressed_virus))

    with open(args.output, 'w') as ofo:
        for expressed_viral_peptide_id, expressed_viral_peptide_seq in expressed_viral_peptides.items():
            ofo.write(">{}\n{}\n".format(expressed_viral_peptide_id, expressed_viral_peptide_seq))


def add_viral_metadata(args):
    """
    """
    expressed_viruses = {}
    cds_to_virus = {}

    with open(args.viral_cds_ref) as fo:
        for line in fo.readlines():
            if line.startswith('>'):
                line = line.split('|')[1].split(' ')[0]
                virus_id, trash, protein_id = line.partition('_cds_')
                protein_id = re.sub(r'_1$', '', protein_id)
                cds_to_virus[protein_id] = virus_id

    with open(args.viral_quants) as fo:
        for line_idx, line in enumerate(fo.readlines()):
            if line_idx != 0:
                line = line.rstrip().split('\t')
                expressed_viruses[line[0]] = line[1]

    output_lines = []
    header = []

    with open(args.binding_affinities) as fo:
        for line_idx, line in enumerate(fo.readlines()):
            line = line.rstrip().split()
            if len(line) >= 16 and line[0] == 'Pos':
                header = line
                try:
                    header.remove('BindLevel')
                except:
                    pass
                header.extend(['virus_id', 'read_count'])
            if len(line) < 16 or line[0] in ['Protein', 'Pos']:
                continue
            else:
                virus_id = cds_to_virus[line[10].replace('_1', '.1')]
                read_counts = expressed_viruses[virus_id]
                print(read_counts)
                line.extend([virus_id, read_counts])
                output_lines.append(line)

    with open(args.output, 'w') as ofo:
        ofo.write("{}\n".format('\t'.join(header)))
        for line in output_lines:
            ofo.write("{}\n".format('\t'.join(line)))

def make_fusion_peptides(args):
     """
     """
     valid_fusions = {}
    with open(args.fusions) as fo:
         header = []
         header_map = {}
        for line_idx, line in enumerate(fo.readlines()):
            line = line.rstrip().split('\t')
             if line_idx  == 0:
                header_map = {j:i for i, j in enumerate(line)}
             else:
                 if line[header_map['FUSION_TRANSL']] != '.':
                     start_prot = (int(line[header_map['CDS_LEFT_RANGE']].split('-')[1]) - 1)/3
                     peptide = line[header_map['FUSION_TRANSL']][start_prot - 8:start_prot + 8]
                     if not(re.search('\*', peptide)):
                         valid_fusions[line[header_map['#FusionName']]] = peptide

     with open(args.output, 'w') as ofo:
         for valid_fusion_id, valid_fusion_seq in valid_fusions.items():
             ofo.write(">{}\n{}\n".format(valid_fusion_id, valid_fusion_seq))


def add_fusion_metadata(args):
    """
    """
    fusion_metadata = {}
    with open(args.fusions) as fo:
        col_to_idx = {}
        for line_idx, line in enumerate(fo.readlines()):
            line = line.rstrip().split()
            if line_idx ==  0:
                for col_idx, col in enumerate(line):
                    col_to_idx[col_idx] = col
            else:
                fusion_metadata[line[0]] = {}
                for elem_idx, elem in enumerate(line[1:]):
                    fusion_metadata[line[0]][col_to_idx[elem_idx + 1]] = elem

    extensions = ['JunctionReadCount', 'SpanningFragCount', 'SpliceType', 'LeftGene',
                  'LeftBreakpoint', 'RightGene', 'RightBreakpoint', 'LargeAnchorSupport',
                  'FFPM', 'annots']
    output_lines = []
    header = []

    with open(args.binding_affinities) as fo:
        for line_idx, line in enumerate(fo.readlines()):
            line = line.rstrip().split()
            if len(line) >=16 and line[0] == 'Pos':
                header = line
                try:
                    header.remove('BindLevel')
                except:
                    pass
                header.extend(extensions)
            if len(line) < 16 or line[0] in ['Protein', 'Pos']:
                continue
            else:
                relevant_metadata = fusion_metadata[line[10]]
                line_extension = []
                for extension in extensions:
                    line_extension.append(relevant_metadata[extension])
                line.extend(line_extension)
                output_lines.append(line)

    with open(args.output, 'w') as ofo:
        ofo.write("{}\n".format('\t'.join(header)))
        for line in output_lines:
            ofo.write("{}\n".format('\t'.join(line)))



def main():
    args = get_args()
    #print(args)
    if args.command == 'make-peptides':
        make_peptides(args)
    if args.command == 'make-indel-peptides':
        make_indel_peptides(args)
    if args.command == 'expressed-variants':
        expressed_variants(args)
    if args.command == 'rna-covered-variants':
        rna_covered_variants(args)
    if args.command == 'isolated-variants':
        isolated_variants(args)
    if args.command == 'calculate-agretopicity':
        calculate_agretopicity(args)
    if args.command == 'make-pyclone-vi-inputs':
        make_pyclone_vi_inputs(args)
    if args.command == 'add-snv-metadata':
        add_snv_metadata(args)
    if args.command == 'add-indel-metadata':
        add_indel_metadata(args)
    if args.command == 'make-genomic-context':
        make_genomic_context(args)
    if args.command == 'expressed-hervs':
        expressed_hervs(args)
    if args.command == 'make-herv-peptides':
        make_herv_peptides(args)
    if args.command == 'get-herv-metadata':
        get_herv_metadata(args)
    if args.command == 'expressed-self-genes':
        expressed_self_genes(args)
    if args.command == 'make-self-antigen-peptides':
        make_self_antigen_peptides(args)
    if args.command == 'add-self-antigen-metadata':
        add_self_antigen_metadata(args)
    if args.command == 'filter-virdetect-by-counts':
        filter_virdetect_by_counts(args)
    if args.command == 'make-viral-peptides':
        make_viral_peptides(args)
    if args.command == 'add-viral-metadata':
        add_viral_metadata(args)
     if args.command == 'make-fusion-peptides':
         make_fusion_peptides(args)
    if args.command == 'add-fusion-metadata':
        add_fusion_metadata(args)
    if args.command == 'consolidate-multiqc-stats':
        consolidate_multiqc_stats(args)

if __name__=='__main__':
    main()
