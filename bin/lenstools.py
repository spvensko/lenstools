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

    # Subparser for generating SNV peptides
    parser_make_snv_peptides = subparsers.add_parser('make-snv-peptides',
                                                     help="Make SNV peptide FASTA file.")
    parser_make_snv_peptides.add_argument('-t', '--tx-aa-fasta',
                                          help="Transcript (CDS) amino acid sequences.",
                                          required=True)
    parser_make_snv_peptides.add_argument('-v', '--somatic-vcf',
                                          help="Annotated (snpEff) somatic VCF).",
                                          required=True)
    #Revisit the peptide length argument. We probably want to emit 15-mers
    #since that allows for an 8-mer sliding window in which all sequences contain
    #the mutant peptide.
    parser_make_snv_peptides.add_argument('-l', '--length',
                                          help="Emitted peptide length.",
                                          default=8)
    parser_make_snv_peptides.add_argument('-o', '--mt-output',
                                          help="Mutant peptides output file.",
                                          required=True)
    parser_make_snv_peptides.add_argument('-w', '--wt-output',
                                          #Note: Wildtype means REFERENCE! Need
                                          #to incorporate patient variants for
                                          #agretopicity calculations.
                                          help="Wildtype peptides output file.",
                                          required=True)

    #This mode needs more information. How is it to be used? Does one simply
    #pass it the VCF? 
    # Subparser for generating genomic context
    parser_get_context = subparsers.add_parser('get-snv-genomic-context',
                                               help="Make SNV neoAg nt context FASTA file.")
    parser_get_context.add_argument('-c', '--tx-cds-fasta',
                                      help="Transcript (CDS) nucleotide sequenes.",
                                      required=True)
    parser_get_context.add_argument('-v', '--vcf',
                                      help="Annotated (snpEff) VCF).",
                                      required=True)
    parser_get_context.add_argument('-l', '--length',
                                      help="Emitted nucleotide sequence length (default: 39)",
                                      default=39)
    parser_get_context.add_argument('-o', '--output',
                                      help="output file.",
                                      required=True)

    # Subparser for generating indel peptides
    #Revisit the peptide length argument. We probably want to emit 15-mers
    #since that allows for an 8-mer sliding window in which all sequences contain
    #the mutant peptide.
    parser_make_indel_peptides = subparsers.add_parser('make-indel-peptides',
                                                       help="Make InDel peptide FASTA file.")
    parser_make_indel_peptides.add_argument('-t', '--tx-aa-fasta',
                                            help="Transcript (CDS) amino acid sequences.",
                                            required=True)
    parser_make_indel_peptides.add_argument('-c', '--tx-cds-fasta',
                                            help="Transcript (CDS) nucleotide sequences.",
                                            required=True)
    parser_make_indel_peptides.add_argument('-v', '--somatic-vcf',
                                            help="Annotated (snpEff) VCF).",
                                            required=True)
    parser_make_indel_peptides.add_argument('-l', '--length',
                                            help="Total peptide length.",
                                            default=8)
    parser_make_indel_peptides.add_argument('-o', '--output',
                                            help="Output file.",
                                            required=True)

    # Subparser for filtering variants for expression
    parser_expressed_variants = subparsers.add_parser('expressed-variants',
                                                 help="Filter variants for expression.")
    parser_expressed_variants.add_argument('-q', '--quants',
                                           help="Transcript abundance file (.quant).",
                                           required=True)
    parser_expressed_variants.add_argument('-m', '--metric',
                                           help="Expression metric (default: TPM).",
                                           default="TPM")
    parser_expressed_variants.add_argument('-r', '--exclude-zeros',
                                           help="Exclude zeros in expression percentile.",
                                           action='store_true')
    parser_expressed_variants.add_argument('-p', '--percentile',
                                           help="Expression percentile for filtering (default: 50).",
                                           default=50)
    parser_expressed_variants.add_argument('-t', '--abundance-threshold',
                                           help="Expression threshold for filtering (default: 0).",
                                           default=0)
    parser_expressed_variants.add_argument('-v', '--vcf',
                                           help="Annotated (snpEff) VCF).",
                                           required=True)
    parser_expressed_variants.add_argument('-o', '--output',
                                           help="Output file.",
                                           required=True)

    # Subparser for filtering variants for having sufficient RNA coverage
    # This primarily applies to SNVs/InDels now since other variant types are
    # derived from RNA-Seq reads (e.g. virdetect, hERVQuant, STARFusion,
    # NeoSplice).
    parser_covered_variants = subparsers.add_parser('rna-covered-variants',
                                                 help="Filter variants for evidence in RNA reads.")
    parser_covered_variants.add_argument('-b', '--tumor-rna-bam',
                                         help="Tumor RNA BAM file.",
                                         required=True)
    parser_covered_variants.add_argument('-c', '--required-coverage',
                                         help="Required coverage for variants (default: 1)",
                                         default=1)
    parser_covered_variants.add_argument('-v', '--vcf',
                                         help="Annotated (snpEff) VCF).",
                                         required=True)
    parser_covered_variants.add_argument('-o', '--output',
                                         help="Output file.",
                                         required=True)

    # Subparser for filtering isolated variants (e.g. no proximal germline or somatic variants).
    parser_isolated_variants = subparsers.add_parser('isolated-variants',
                                                     help="Filter isolated variants.")
    parser_isolated_variants.add_argument('-g', '--germline-vcf',
                                           help="Germline VCF.",
                                           required=True)
    parser_isolated_variants.add_argument('-s', '--somatic-vcf',
                                           help="Somatic VCF.",
                                           required=True)
    parser_isolated_variants.add_argument('-p', '--proximity',
                                           help="Required clearance (in bp) around variants (default: 30).",
                                           default=30)
    #Allowing proximal silent (het/hom) and homozygous missense variants still needs to be implemented.
    parser_isolated_variants.add_argument('-a', '--allow-silent-and-homozygous',
                                           help="Allow silent and homozygous germline variants near variants.",
                                           action="store_true")
    parser_isolated_variants.add_argument('-o', '--output',
                                          help="Output file.",
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
                                               help="Output file.",
                                               required=True)

    #This currently uses MuTect for getting depth information. Ideally, this
    #would be using other sources of depth information too (e.g. Strelka2),
    #averaging the variant depths when appropriate.
    # Subparser for creating PyClone-VI inputs
    parser_make_pvi_inputs = subparsers.add_parser('make-pyclone-vi-inputs',
                                                   help="Make input file required for PyClone-VI")
    parser_make_pvi_inputs.add_argument('-c', '--candidate-vcf',
                                        help="VCF containing candidate variants.",
                                        required=True)
    #MuTect VCF is specifically referenced here since the parser is designed
    #around its outputs. This should expanded to include, at least, Strelka2 as
    #well. 
    parser_make_pvi_inputs.add_argument('-m', '--mutect-vcf',
                                        help="Mutect VCF for variant depth information.",
                                        required=True)
    parser_make_pvi_inputs.add_argument('-s', '--sequenza-segments',
                                        help="Sequenza segments file.",
                                        required=True)
    parser_make_pvi_inputs.add_argument('--sequenza-solutions',
                                        help="Sequenza solutions file.",
                                        required=True)
    parser_make_pvi_inputs.add_argument('--samp-id',
                                        help="Sample identifier.",
                                        required=True)
    parser_make_pvi_inputs.add_argument('-o', '--output', 
                                        help="Output file.",
                                        required=True)

    # Subparser for adding SNV metadata
    parser_add_snv_metadata = subparsers.add_parser('add-snv-metadata',
                                                    help="Add metadata to SNV binding affinity data.")
    parser_add_snv_metadata.add_argument('-m', '--mutant-peptides',
                                         help="FASTA file with mutant peptides.",
                                         required=True)
    parser_add_snv_metadata.add_argument('-u', '--mutant-nucs',
                                         help="FASTA file with mutant nucleotide context.",
                                         required=True)
    parser_add_snv_metadata.add_argument('-q', '--quants',
                                         help="Transcript abundance file (Salmon format).",
                                         required=True)
    parser_add_snv_metadata.add_argument('-c', '--cancer-cell-fraction',
                                         help="Cancer cell fraction file (PyClone-VI format).",
                                         required=True)
    parser_add_snv_metadata.add_argument('-g', '--gtf',
                                         help="GTF file.",
                                         required=True)
    parser_add_snv_metadata.add_argument('-b', '--binding-affinities',
                                         help="Binding affinities file (netMHCpan format).",
                                         required=True)
    parser_add_snv_metadata.add_argument('-o', '--output', 
                                         help="Output file.",
                                         required=True)

    # Subparser for adding InDel metadata
    parser_add_indel_metadata = subparsers.add_parser('add-indel-metadata',
                                                      help="Add metadata to InDel binding affinity data.")
    parser_add_indel_metadata.add_argument('-m', '--mutant-peptides',
                                           help="FASTA file with mutant peptides.",
                                           required=True)
    parser_add_indel_metadata.add_argument('-q', '--quants',
                                           help="Transcript abundance file (Salmon format).",
                                           required=True)
    parser_add_indel_metadata.add_argument('-c', '--cancer-cell-fraction',
                                           help="Cancer cell fraction file (PyClone-VI format).",
                                           required=True)
    parser_add_indel_metadata.add_argument('-g', '--gtf',
                                           help="GTF file.",
                                           required=True)
    parser_add_indel_metadata.add_argument('-b', '--binding-affinities',
                                           help="Binding affinities file (netMHCpan format).",
                                           required=True)
    parser_add_indel_metadata.add_argument('-o', '--output',
                                           required=True)


    # Subparser for filtering expressed hERVs
    parser_expressed_hervs = subparsers.add_parser('expressed-hervs',
                                                   help="Filter expressed hervs.")
    parser_expressed_hervs.add_argument('-q', '--quants',
                                        help="Transcript abundance file (Salmon format).",
                                        required=True)
    parser_expressed_hervs.add_argument('-m', '--metric',
                                        help="Column for expression from abundance file.",
                                        default="TPM")
    parser_expressed_hervs.add_argument('-r', '--exclude-zeros',
                                        help="Exclude zeros for expression percentile.",
                                        action='store_true')
    parser_expressed_hervs.add_argument('-p', '--percentile',
                                        help="Expression percentile for filtering expression. (default: 50)",
                                        default=50)
    parser_expressed_hervs.add_argument('-t', '--abundance-threshold',
                                        help="Expression threshold for filtering. (default: 0)",
                                        default=0)
    parser_expressed_hervs.add_argument('-o', '--output',
                                        help="Output file.",
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
                                           help="Output file.",
                                           required=True)

    # Subparser for adding hERV metadata
    parser_add_herv_metadata = subparsers.add_parser('add-herv-metadata',
                                                     help="Add hERV peptide metadata.")
    parser_add_herv_metadata.add_argument('-q', '--quants',
                                          help="Transcript abundance file (Salmon format).",
                                          required=True)
    parser_add_herv_metadata.add_argument('-b', '--binding-affinities',
                                          help="Binding affinities file (netMHCpan format).",
                                          required=True)
    parser_add_herv_metadata.add_argument('-o', '--output',
                                          help="Output file.",
                                          required=True)

    # Subparser for consolidating multiqc statistics files
    parser_consol_mqc_stats = subparsers.add_parser('consolidate-multiqc-stats',
                                                   help="Consolidate multiqc stats files into a single file.")
    parser_consol_mqc_stats.add_argument('-d', '--multiqc-data',
                                         help="Directory containing multiqc stats files.",
                                         required=True)
    parser_consol_mqc_stats.add_argument('-o', '--output',
                                         help="Output file.",
                                         required=True)

    # Subparser for filtering self-antigens for expression
    parser_expressed_self_genes = subparsers.add_parser('expressed-self-genes',
                                                        help="Filter expressed self-antigen genes.")
    parser_expressed_self_genes.add_argument('-q', '--quants',
                                             help="Transcript abundance file (Salmon format).",
                                             required=True)
    parser_expressed_self_genes.add_argument('-m', '--metric',
                                             help="Expression metric.",
                                             default="TPM")
    parser_expressed_self_genes.add_argument('-r', '--exclude-zeros',
                                             help="Exclude zeros for expression percentile.",
                                             action='store_true')
    parser_expressed_self_genes.add_argument('-p', '--percentile',
                                             help="Expression percentile for filtering expression (default: 50).",
                                             default=50)
    parser_expressed_self_genes.add_argument('-t', '--abundance-threshold',
                                             help="Expression abundance threshold for filtering expression (default: 0)",
                                             default=0)
    parser_expressed_self_genes.add_argument('-g', '--gene-list',
                                             help="File containing genes of interest (e.g. self-antigens, CTAs, etc.)",
                                             required=True)
    parser_expressed_self_genes.add_argument('-f', '--gtf',
                                             help="GTF file.",
                                             required=True)
    parser_expressed_self_genes.add_argument('-o', '--output',
                                             help="Output file.",
                                             required=True)


    # Subparser for making self-antigen peptides
    parser_make_self_peptides = subparsers.add_parser('make-self-antigen-peptides',
                                                      help="Make self-antigen peptides.")
    parser_make_self_peptides.add_argument('-e', '--expressed-selfs',
                                           help="File containing list of expressed self-antigen genes.",
                                           required=True)
    parser_make_self_peptides.add_argument('-r', '--tx-aa-fasta',
                                           help="Transcript (CDS) amino acid FASTA",
                                           required=True)
    parser_make_self_peptides.add_argument('-s', '--somatic-vcf',
                                           help="Somatic VCF",
                                           required=True)
    parser_make_self_peptides.add_argument('-g', '--germline-vcf',
                                           help="Somatic VCF",
                                           required=True)
    parser_make_self_peptides.add_argument('-o', '--output',
                                           help="Output file.",
                                           required=True)


    # Subparser for getting self-antigen metadata
    parser_add_self_metadata = subparsers.add_parser('add-self-antigen-metadata',
                                                     help="Add self-antigen metadata.")
    parser_add_self_metadata.add_argument('-q', '--quants',
                                          help="Transcript abundance file (Salmon format).",
                                          required=True)
    parser_add_self_metadata.add_argument('-b', '--binding-affinities',
                                          help="Binding affinitities data (netMHCpan format)",
                                          required=True)
    parser_add_self_metadata.add_argument('-g', '--gtf',
                                          help="GTF file",
                                          required=True)
    parser_add_self_metadata.add_argument('-o', '--output',
                                          help="Output file.",
                                          required=True)


    # Subparser for filtering virdetect outputs for expressed viruses
    parser_filter_virdetect_by_counts = subparsers.add_parser('filter-virdetect-by-counts',
                                                   help="Generate list of expressed viruses")
    parser_filter_virdetect_by_counts.add_argument('-q', '--viral-quants',
                                                   help="Viral counts file from VirDetect.",
                                                   required=True)
    parser_filter_virdetect_by_counts.add_argument('-r', '--viral-ref',
                                                   help="Viral reference FASTA (used for VirDetect).",
                                                   required=True)
    parser_filter_virdetect_by_counts.add_argument('-m', '--min-threshold',
                                                   help="Minimal count threshold for filtering (default: 1).",
                                                   default='1')
    parser_filter_virdetect_by_counts.add_argument('-o', '--output', 
                                                   help="Output file.",
                                                   required=True)


    # Subparser for making viral peptides
    parser_make_viral_peptides = subparsers.add_parser('make-viral-peptides',
                                                       help="Make viral peptides")
    parser_make_viral_peptides.add_argument('-e', '--expressed-viruses',
                                            help="File containing list of expressed viruses.",
                                            required=True)
    parser_make_viral_peptides.add_argument('-c', '--sample-viral-ref',
                                            help="Viral CDS reference FASTA (specific to patient).",
                                            required=True)
    parser_make_viral_peptides.add_argument('-o', '--output',
                                            help='Output file.',
                                            required=True)


    # Subparser for getting viral metadata
    parser_add_viral_metadata = subparsers.add_parser('add-viral-metadata',
                                                      help="Add viral metadata.")
    parser_add_viral_metadata.add_argument('-b', '--binding-affinities',
                                           help="Binding affinities data (netMHCpan format).",
                                           required=True)
    parser_add_viral_metadata.add_argument('-q', '--viral-quants',
                                           help="Viral counts file from VirDetect.",
                                           required=True)
    parser_add_viral_metadata.add_argument('-r', '--viral-cds-ref',
                                           help="Viral CDS reference FASTA (used for VirDetect).",
                                           required=True)
    parser_add_viral_metadata.add_argument('-o', '--output', 
                                           help="Output file.",
                                           required=True)


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
                                                       help="Add metadata to fusion binding affinity data.")
    parser_add_fusion_metadata.add_argument('-b', '--binding-affinities',
                                            help="Binding affinities data (netMHCpan format).",
                                            required=True)
    parser_add_fusion_metadata.add_argument('-f', '--fusions',
                                            help="Predicted fusions (STARFusion format).",
                                            required=True)
    parser_add_fusion_metadata.add_argument('-o', '--output',
                                            help="Output file.",
                                            required=True)


    return parser.parse_args()


def load_tx_aas(args):
    """
    Loads a transcript -> amino acid sequence FASTA into a dictionary.
    """
    tx_to_aa = {}
    for seq_record in SeqIO.parse(args.tx_aa_fasta, "fasta"):
        tx = re.search('transcript:\S*', seq_record.description).group(0).split(':')[1]
        tx_to_aa[tx] = seq_record.seq
    return tx_to_aa


def load_tx_cds(args):
    """
    Loads a transcript -> nucleotide sequence FASTA into a dictionary.
    """
    tx_to_cds = {}
    for seq_record in SeqIO.parse(args.tx_cds_fasta, "fasta"):
        tx = seq_record.id
        tx_to_cds[tx] = seq_record.seq
    return tx_to_cds


def extract_missense_snvs(input_vcf):
    """
    Extracts missense SNVs from annotated somatic VCF file.
    """
    filtered_records = []

    vcf_reader = ''
    if input_vcf.endswith('gz'):
        vcf_reader = vcf.Reader(open(input_vcf), 'r', compressed=True)
    else:
        vcf_reader = vcf.Reader(open(input_vcf), 'r')

    for record in vcf_reader:
        annotations = [x for x in record.INFO['ANN']]
        for annotation in annotations:
            effects = annotation.split('|')
            variant_effect = effects[1]
            if variant_effect == 'missense_variant' and effects[13]:
                transcript = effects[6]
                seq3_change = effects[10].lstrip('p.')
                seq3_change = re.split('(\d+)', seq3_change)
                orig_seq3 = seq3_change[0]
                alt_seq3 = seq3_change[2]
                pos = seq3_change[1]
                orig_aa = seq1(orig_seq3)
                alt_aa = seq1(alt_seq3)
                tlen = effects[13].split('/')[1]
                codon_pos = int((int(effects[13].split('/')[0])*3) - 2)
                filtered_records.append([transcript, pos, tlen, orig_aa, alt_aa, codon_pos, record])
    return filtered_records

#def extract_potential_somatic_nuc_changes(args):
#    """
#    This will need to be cleaned heavily, but good for now.
#    """
#    filtered_records = []
#    #vcf_reader = vcf.Reader(open(args.vcf), 'r', compressed=True)
#    vcf_reader = vcf.Reader(open(args.vcf), 'r')
#    for record in vcf_reader:
#        print(record)
#        possible_variants = [x for x in record.INFO['ANN']]
#        for possible_variant in possible_variants:
#            split_possible = possible_variant.split('|')
#            if split_possible[1] == 'missense_variant' and split_possible[13] and split_possible[-1] == '':
#                #transcript = split_possible[6].split('.')[0]
#                transcript = split_possible[6]
#                change = split_possible[9].lstrip('c.')
#                print(change)
#                pos= re.split('(\d+)', change)[1]
#                change = re.split('(\d+)', change)[2].split('>')
#                print(change)
#                pre = change[0]
#                post = change[1]
#                chr_pos = record
#                print("{} {} {} {}".format(transcript, change, pre,post, pos))
#                pos_total = split_possible[12]
#                codon_pos = int((int(split_possible[13].split('/')[0])*3) - 2)
#                snapshot = pos_total.split('/')[0]
#                tlen = pos_total.split('/')[1]
#                if pos != snapshot:
#                    sys.exit(1)
#                filtered_records.append([transcript, pos, tlen, pre, post, codon_pos, chr_pos])
#
#    return filtered_records

def extract_conservative_inframe_indels(args):
    """
    This will need to be cleaned heavily, but good for now.

    Focusing on deletions now, need to find a good example of a conservative_inframe_insertion
    """
    filtered_records = []
    if args.somatic_vcf.endswith('gz'):
        vcf_reader = vcf.Reader(open(args.somatic_vcf), 'r', compressed=True)
    else:
        vcf_reader = vcf.Reader(open(args.somatic_vcf), 'r')
    
    for record in vcf_reader:
        annotations = [x for x in record.INFO['ANN']]
        for annotation in annotations:
            effects = annotation.split('|')
            if re.search('^conservative_inframe_[a-z]+$', effects[1]) and effects[13]:
                transcript = effects[6]
                seq3_change = effects[10].lstrip('p.')
                filtered_records.append([transcript, seq3_change, record])
    return filtered_records

def extract_disruptive_inframe_indels(args):
    """
    This will need to be cleaned heavily, but good for now.

    Focusing on deletions now, need to find a good example of a conservative_inframe_insertion
    """
    filtered_records = []
    if args.somatic_vcf.endswith('gz'):
        vcf_reader = vcf.Reader(open(args.somatic_vcf), 'r', compressed=True)
    else:
        vcf_reader = vcf.Reader(open(args.somatic_vcf), 'r')
    
    for record in vcf_reader:
        annotations = [x for x in record.INFO['ANN']]
        for annotation in annotations:
            effects = annotation.split('|')
            if re.search('^disruptive_inframe_[a-z]+$', effects[1]) and effects[13]:
                transcript = effects[6]
                seq3_change = effects[10].lstrip('p.')
                filtered_records.append([transcript, seq3_change, record])
    return filtered_records

def extract_frameshift_indels(args):
    """
    This will need to be cleaned heavily, but good for now.

    Focusing on deletions now, need to find a good example of a conservative_inframe_insertion
    """
    filtered_records = []
    if args.somatic_vcf.endswith('gz'):
        vcf_reader = vcf.Reader(open(args.somatic_vcf), 'r', compressed=True)
    else:
        vcf_reader = vcf.Reader(open(args.somatic_vcf), 'r')
    
    for record in vcf_reader:
        annotations = [x for x in record.INFO['ANN']]
        for annotation in annotations:
            effects = annotation.split('|')
            if re.search('^frameshift_variant$', effects[1]) and effects[13] and effects[7] == 'protein_coding' and effects[-1] == '':
                transcript = effects[6]
                dna_change = effects[9].lstrip('c.')
                filtered_records.append([transcript, dna_change, record])
    return filtered_records

def make_snv_peptides(args):
    """
    """
    tx_to_aa = load_tx_aas(args)
    missense_snvs = extract_missense_snvs(args.somatic_vcf)
    mutant_peptides = {}
    reference_peptides = {}
    for record in missense_snvs:
        if record[0] not in tx_to_aa.keys():
            continue

        tx_ref_seq = list(tx_to_aa[record[0]])
        if len(tx_ref_seq) != int(record[2]):
            print("transcript {} shows different lengths between peptide fasta ({}) and snpeff! ({})".format(record[0], record[2], len(tx_ref_seq)))

        mut_seq = tx_ref_seq[:]
        pos = int(record[1]) - 1
        if mut_seq[pos] != record[3]:
            print("Reference amino acid doesn't match! Something has gone horribly wrong.")
            continue

        mut_seq[pos] = record[4]
        ref_peptide = ''.join(tx_ref_seq[max(pos-8+1,0):min(pos+8, len(mut_seq))])
        mut_peptide = ''.join(mut_seq[max(pos-8+1, 0):min(pos+8, len(mut_seq))])
        header_mut_peptide = ''.join(mut_seq[max(pos-13, 0):min(pos+14, len(mut_seq))])
        #var_md5 is being used to create a unique identifier for the resulting mutant peptide to overcome netMHCpan's length limitations.
        var_md5 = hashlib.md5("{}".format(':'.join([str(record[-1].CHROM), str(record[-1].POS), str(record[0]), str(record[-1].REF), str(record[-1].ALT)]))).hexdigest()[:16]
        mutant_peptides["{} {}:{} {} {} {} {}".format(var_md5, record[-1].CHROM, record[-1].POS, record[0], record[-1].REF, record[-1].ALT, header_mut_peptide)] = mut_peptide
        reference_peptides["{} {}:{} {} {} {}".format(var_md5, record[-1].CHROM, record[-1].POS, record[0], record[-1].REF, record[-1].ALT)] = ref_peptide

    with open(args.mt_output, 'w') as ofo:
        for k, v in mutant_peptides.items():
            ofo.write('>{}\n{}\n'.format(k, v))
    with open(args.wt_output, 'w') as ofo:
        for k, v in reference_peptides.items():
            ofo.write('>{}\n{}\n'.format(k, v))


def get_snv_genomic_context(args):
    """
    """
    tx_to_cds = load_tx_cds(args)
    missense_snvs = extract_missense_snvs(args.vcf)
    emitted_nucs = {}
    for record in missense_snvs:
        if record[0] not in tx_to_cds.keys():
            continue
        tx_ref_seq = list(tx_to_cds[record[0]])
#        if len(tx_ref_seq) != int(record[2]):
#            print("transcript {} shows different lengths between peptide fasta ({}) and snpeff! ({})".format(record[0], record[2], len(tx_ref_seq)))

        mut_seq = tx_ref_seq[:]
        pos = int(record[1]) - 1
#        if mut_seq[pos] != record[3]:
#            print("Reference amino acid doesn't match! Something has gone horribly wrong.")
#            continue
        # Have to be careful here... we want to start at the first base of the affected codon
        mut_seq[pos] = record[4]
        codon_pos = record[5] - 1
        ref_nuc = ''.join(tx_ref_seq[max(codon_pos-args.length-1, 0):min(codon_pos+args.length, len(tx_ref_seq)) + 3])
        mut_nuc = ''.join(mut_seq[max(codon_pos-args.length, 0):min(codon_pos+args.length, len(mut_seq)) + 3])
        print("{}\n{}\n{}".format(record[0],ref_nuc, mut_nuc))
        print(':'.join([str(record[-1].CHROM), str(record[-1].POS), str(record[0]), str(record[-1].REF), str(record[-1].ALT)]))
        var_md5 = hashlib.md5("{}".format(':'.join([str(record[-1].CHROM), str(record[-1].POS), str(record[0]), str(record[-1].REF), str(record[-1].ALT)])).encode('utf-8')).hexdigest()[:16]
        emitted_nucs["{}\t{}:{}\t{}\t{}\t{}".format(var_md5, record[-1].CHROM, record[-1].POS, record[0], record[-1].REF, record[-1].ALT)] = mut_nuc

    with open(args.output, 'w') as ofo:
        for k, v in emitted_nucs.items():
            ofo.write('>{}\t{}\n'.format(k, v))

def make_indel_peptides(args):
    """
    """
    tx_to_aa = load_tx_aas(args)
    tx_nv_to_aa = {}
    for k, v in tx_to_aa.items():
        k_no_version = k.split('.')[0]
        tx_nv_to_aa[k_no_version] =  v
    tx_to_cds = load_tx_cds(args)
    tx_nv_to_cds = {}
    for k, v in tx_to_cds.items():
        k_no_version = k.split('.')[0]
        tx_nv_to_cds[k_no_version] =  v
    print("Getting indels")
    conserv_inframe_indels = extract_conservative_inframe_indels(args)
    disrupt_inframe_indels = extract_disruptive_inframe_indels(args)
    frameshift_indels = extract_frameshift_indels(args)
    print("Got indels")
#    print(conserv_inframe_indels)
#    print(disrupt_inframe_indels)
#    print(frameshift_indels)

    mutant_peptides = {}

    all_indels = conserv_inframe_indels[:]
    all_indels.extend(disrupt_inframe_indels)
#    print(all_indels)

    for record in all_indels:
        print(record[0])
        tx_no_version = record[0].split('.')[0]
        if record[0] not in tx_to_aa.keys() and tx_no_version not in tx_nv_to_aa.keys():
            continue
        tx_ref_seq = []
        if record[0] in tx_to_aa.keys():
            tx_ref_seq = list(tx_to_aa[record[0]])
        elif tx_no_version in tx_nv_to_aa.keys():
            tx_ref_seq = list(tx_nv_to_aa[tx_no_version])
        mut_seq = tx_ref_seq[:]
        print(record[1])
        if re.search("del$", record[1]):
            del_rec = record[1].strip('del')
            start_pos = int(del_rec.split('_')[0][3:]) - 1
            start_aa = seq1(del_rec.split('_')[0][:3])
            stop_pos = 0
            stop_aa = 'foo'
            # Dealing with single amino acid deletions
            if re.search('_', del_rec):
                stop_pos = int(del_rec.split('_')[1][3:]) - 1
                stop_aa = seq1(del_rec.split('_')[1][:3])
            else:
                stop_pos = start_pos
                stop_aa = start_aa
            print("{} {}".format(start_pos, start_aa))
            print("{} {}".format(stop_pos, stop_aa))
            if start_pos > len(mut_seq) or stop_pos > len(mut_seq):
                print("The mutation occurs outside of the peptide sequence. Check transcript versions.")
                continue
            if mut_seq[start_pos] != start_aa or mut_seq[stop_pos] != stop_aa:
                print("Reference amino acid doesn't match! Something has gone horribly wrong.")
                continue
            new_mut_seq = ''.join(mut_seq[:start_pos] + mut_seq[stop_pos+1:])
            ref_peptide = ''.join(tx_ref_seq[max(start_pos-7,0):min(stop_pos+7, len(mut_seq))])
            mut_peptide = ''.join(new_mut_seq[max(start_pos-7, 0):min(start_pos+7, len(new_mut_seq))])
            print("Starting peptide: {}".format(ref_peptide))
            print("mutated peptide: {}".format(mut_peptide))
            md5able_str = "{}".format(':'.join([str(record[-1].CHROM), str(record[-1].POS), str(record[0]), str(record[-1].REF), str(record[-1].ALT)])).encode('utf-8')
            var_md5 = hashlib.md5(md5able_str).hexdigest()[:16]
            mutant_peptides["{} {}:{} {} {} {} {}".format(var_md5, record[-1].CHROM, record[-1].POS, record[0], record[-1].REF, record[-1].ALT, 'inframe_deletion')] = mut_peptide
        elif re.search("[0-9]ins", record[1]):
            ins_rec, spacer, ins_seq = record[1].partition('ins')
            ins_seq1 = seq1(ins_seq)
            print(ins_seq1)
            insert_pos = int(ins_rec.split('_')[0][3:]) - 1
            insert_aa = seq1(ins_rec.split('_')[0][:3])
            insert_plus_one_pos = int(ins_rec.split('_')[1][3:]) - 1
            insert_plus_one_aa = seq1(ins_rec.split('_')[1][:3])
            if mut_seq[insert_pos] != insert_aa or mut_seq[insert_plus_one_pos] != insert_plus_one_aa:
                print("Reference amino acid doesn't match! Something has gone horribly wrong.")
                continue
            new_mut_seq = ''.join(mut_seq[:insert_pos+1] + [ins_seq1] + mut_seq[insert_plus_one_pos:])
            mut_peptide = ''.join(new_mut_seq[max(insert_pos-6, 0):min(insert_plus_one_pos+len(ins_seq1)+7, len(new_mut_seq))])
            ref_peptide = ''.join(mut_seq[max(insert_pos-6,0):min(insert_plus_one_pos+7, len(mut_seq))])
            print("Starting peptide: {}".format(ref_peptide))
            print("mutated peptide: {}".format(mut_peptide))
            md5able_str = "{}".format(':'.join([str(record[-1].CHROM), str(record[-1].POS), str(record[0]), str(record[-1].REF), str(record[-1].ALT)])).encode('utf-8')
            var_md5 = hashlib.md5(md5able_str).hexdigest()[:16]
            mutant_peptides["{} {}:{} {} {} {} {}".format(var_md5, record[-1].CHROM, record[-1].POS, record[0], record[-1].REF, record[-1].ALT, 'inframe_insertion')] = mut_peptide
        elif re.search("delins", record[1]):
            del_rec, buffer, ins_rec = record[1].partition('delins')
            ins_seq1 = seq1(ins_rec)
            print(ins_seq1)
            start_pos = int(del_rec.split('_')[0][3:]) - 1
            start_aa = seq1(del_rec.split('_')[0][:3])
            stop_pos = 0
            stop_aa = 'foo'
            # Dealing with single amino acid deletions
            if re.search('_', del_rec):
                stop_pos = int(del_rec.split('_')[1][3:]) - 1
                stop_aa = seq1(del_rec.split('_')[1][:3])
            else:
                stop_pos = start_pos
                stop_aa = start_aa
            print("{} {}".format(start_pos, start_aa))
            print("{} {}".format(stop_pos, stop_aa))
            if mut_seq[start_pos] != start_aa or mut_seq[stop_pos] != stop_aa:
                print("Reference amino acid doesn't match! Something has gone horribly wrong.")
                continue
            del_mut_seq = ''.join(mut_seq[:start_pos] + mut_seq[stop_pos+1:])
            ins_mut_seq = ''.join(del_mut_seq[:start_pos] + ins_seq1 + del_mut_seq[start_pos:])
            mut_peptide = ''.join(ins_mut_seq[max(start_pos-7, 0):min(start_pos+len(ins_seq1)+7, len(ins_mut_seq))])
            ref_peptide = ''.join(mut_seq[max(start_pos-7,0):min(stop_pos+len(ins_seq1)+7, len(mut_seq))])
            print("Starting peptide: {}".format(ref_peptide))
            print("mutated peptide: {}".format(mut_peptide))
            md5able_str = "{}".format(':'.join([str(record[-1].CHROM), str(record[-1].POS), str(record[0]), str(record[-1].REF), str(record[-1].ALT)])).encode('utf-8')
            var_md5 = hashlib.md5(md5able_str).hexdigest()[:16]
            mutant_peptides["{} {}:{} {} {} {} {}".format(var_md5, record[-1].CHROM, record[-1].POS, record[0], record[-1].REF, record[-1].ALT, 'inframe_indel')] = mut_peptide
        elif re.search("dup", record[1]):
            dup_rec = record[1].strip('dup')
            start_pos = int(dup_rec.split('_')[0][3:]) - 1
            start_aa = seq1(dup_rec.split('_')[0][:3])
            stop_pos = 0
            stop_aa = 'foo'
            # Dealing with single amino acid deletions
            if re.search('_', dup_rec):
                stop_pos = int(dup_rec.split('_')[1][3:]) - 1
                stop_aa = seq1(dup_rec.split('_')[1][:3])
            else:
                stop_pos = start_pos
                stop_aa = start_aa
            print("{} {}".format(start_pos, start_aa))
            print("{} {}".format(stop_pos, stop_aa))
            if start_pos > len(mut_seq) or stop_pos > len(mut_seq):
                print("The mutation occurs outside of the peptide sequence. Check transcript versions.")
                continue
            if mut_seq[start_pos] != start_aa or mut_seq[stop_pos] != stop_aa:
                print("Reference amino acid doesn't match! Something has gone horribly wrong.")
                continue
            duped_seq = mut_seq[start_pos:stop_pos+1]
            new_mut_seq = mut_seq[:start_pos] + duped_seq + mut_seq[start_pos:]
            mut_peptide = ''.join(new_mut_seq[max(stop_pos-6, 0):min(start_pos+len(duped_seq)+7, len(new_mut_seq))])
            print("Duped seq: {}".format(''.join(duped_seq)))
            print("Ref seq: {}".format(''.join(mut_seq[start_pos-5:stop_pos+5])))
            print("Mut seq: {}".format(''.join(new_mut_seq[start_pos-5:stop_pos+5])))
            print("Mut peptide: {}".format(mut_peptide))
            md5able_str = "{}".format(':'.join([str(record[-1].CHROM), str(record[-1].POS), str(record[0]), str(record[-1].REF), str(record[-1].ALT)])).encode('utf-8')
            var_md5 = hashlib.md5(md5able_str).hexdigest()[:16]
            mutant_peptides["{} {}:{} {} {} {} {}".format(var_md5, record[-1].CHROM, record[-1].POS, record[0], record[-1].REF, record[-1].ALT, 'inframe_dup')] = mut_peptide

    for record in frameshift_indels:
        tx_no_version = record[0].split('.')[0]
        if record[0] not in tx_to_aa.keys() and tx_no_version not in tx_nv_to_aa.keys():
            continue
        tx_ref_seq = []
        if record[0] in tx_to_cds.keys():
            tx_ref_seq = list(tx_to_cds[record[0]])
        elif tx_no_version in tx_nv_to_cds.keys():
            tx_ref_seq = list(tx_nv_to_cds[tx_no_version])
        mut_seq = tx_ref_seq[:]
        if re.search("del", record[1]):
            print("Frameshift deletion {}".format(record))
            del_rec, buffer, del_seq = record[1].partition('del')
            print(del_rec)
            print(del_seq)
            start_pos = 0
            stop_pos = 0
            start_base = 'A'
            stop_base = 'A'
            if re.search('_', del_rec):
                 start_pos = int(del_rec.split('_')[0]) - 1
                 start_base = del_seq[0]
                 stop_pos = int(del_rec.split('_')[1]) - 1
                 stop_base = del_seq[-1]
            else:
                start_pos = int(del_rec) - 1
                start_base = del_seq[0]
                stop_pos = start_pos
                stop_base = start_base
            print("{} {}".format(start_pos, start_base))
            print("{} {}".format(stop_pos, stop_base))
            if stop_pos > len(mut_seq):
                continue
            new_mut_seq = mut_seq[:start_pos] + mut_seq[stop_pos+1:]
            translated_ref = str(Seq(''.join(mut_seq)).translate(to_stop=True, cds=False))
            translated_mut = str(Seq(''.join(new_mut_seq)).translate(to_stop=True, cds=False))
            aa_start = int(start_pos/3)
            aa_stop = int(stop_pos/3)
            print("{}\t{}".format(aa_start, aa_stop))
            mut_peptide = translated_mut[aa_start-7:]
            md5able_str = "{}".format(':'.join([str(record[-1].CHROM), str(record[-1].POS), str(record[0]), str(record[-1].REF), str(record[-1].ALT)])).encode('utf-8')
            var_md5 = hashlib.md5(md5able_str).hexdigest()[:16]
            mutant_peptides["{} {}:{} {} {} {} {}".format(var_md5, record[-1].CHROM, record[-1].POS, record[0], record[-1].REF, record[-1].ALT, 'frameshift_deletion')] = mut_peptide
        if re.search("ins", record[1]):
            print("Frameshift insertion {}".format(record))
            ins_rec, buffer, ins_seq = record[1].partition('ins')
            print(ins_rec)
            print(ins_seq)
            start_pos = 0
            stop_pos = 0
            start_base = 'A'
            stop_base = 'A'
            if re.search('_', ins_rec):
                 start_pos = int(ins_rec.split('_')[0]) - 1
                 start_base = ins_seq[0]
                 stop_pos = int(ins_rec.split('_')[1]) - 1
                 stop_base = ins_seq[-1]
            else:
                start_pos = int(ins_rec) - 1
                start_base = ins_seq[0]
                stop_pos = start_pos
                stop_base = start_base
            print("{} {}".format(start_pos, start_base))
            print("{} {}".format(stop_pos, stop_base))
            if stop_pos > len(mut_seq):
                continue
            new_mut_seq = mut_seq[:start_pos] + [ins_seq] +  mut_seq[stop_pos:]
            translated_ref = str(Seq(''.join(mut_seq)).translate(to_stop=True, cds=False))
            translated_mut = str(Seq(''.join(new_mut_seq)).translate(to_stop=True, cds=False))
            aa_start = int(start_pos/3)
            aa_stop = int(stop_pos/3)
            print("{}\t{}".format(aa_start, aa_stop))
            print(translated_ref[aa_start-5:aa_stop+5])
            print(translated_mut[aa_start-5:aa_stop+5])
            mut_peptide = translated_mut[aa_start-7:]
            md5able_str = "{}".format(':'.join([str(record[-1].CHROM), str(record[-1].POS), str(record[0]), str(record[-1].REF), str(record[-1].ALT)])).encode('utf-8')
            var_md5 = hashlib.md5(md5able_str).hexdigest()[:16]
            mutant_peptides["{} {}:{} {} {} {} {}".format(var_md5, record[-1].CHROM, record[-1].POS, record[0], record[-1].REF, record[-1].ALT, 'frameshift_insertion')] = mut_peptide
        if re.search("dup", record[1]):
            print("Frameshift duplication {}".format(record))
            dup_rec, buffer, dup_seq = record[1].partition('dup')
            print(dup_rec)
            print(dup_seq)
            start_pos = 0
            stop_pos = 0
            start_base = 'A'
            stop_base = 'A'
            if re.search('_', dup_rec):
                 start_pos = int(dup_rec.split('_')[0]) - 1
                 start_base = dup_seq[0]
                 stop_pos = int(dup_rec.split('_')[1]) - 1
                 stop_base = dup_seq[-1]
            else:
                start_pos = int(dup_rec) - 1
                start_base = dup_seq[0]
                stop_pos = start_pos
                stop_base = start_base
            print("{} {}".format(start_pos, start_base))
            print("{} {}".format(stop_pos, stop_base))
            if stop_pos > len(mut_seq):
                continue
            new_mut_seq = mut_seq[:start_pos] + [dup_seq] +  mut_seq[stop_pos:]
            translated_ref = str(Seq(''.join(mut_seq)).translate(to_stop=True, cds=False))
            translated_mut = str(Seq(''.join(new_mut_seq)).translate(to_stop=True, cds=False))
            aa_start = int(start_pos/3)
            aa_stop = int(stop_pos/3)
            print("{}\t{}".format(aa_start, aa_stop))
            print(translated_ref[aa_start-5:aa_stop+5])
            print(translated_mut[aa_start-5:aa_stop+5])
            mut_peptide = translated_mut[aa_start-7:]
            md5able_str = "{}".format(':'.join([str(record[-1].CHROM), str(record[-1].POS), str(record[0]), str(record[-1].REF), str(record[-1].ALT)])).encode('utf-8')
            var_md5 = hashlib.md5(md5able_str).hexdigest()[:16]
            mutant_peptides["{} {}:{} {} {} {} {}".format(var_md5, record[-1].CHROM, record[-1].POS, record[0], record[-1].REF, record[-1].ALT, 'frameshift_duplication')] = mut_peptide



    with open(args.output, 'w') as ofo:
        for k, v in mutant_peptides.items():
            ofo.write('>{}\n{}\n'.format(k, v))


def rna_covered_variants(args):
    """
    """
    filtered_chrom_pos = []
    filtered_vcf = []
    rna_bam = pysam.AlignmentFile(args.rna_bam, "rb")
    vcf_reader = vcf.Reader(filename=args.vcf, compressed=False)
    for record in vcf_reader:
        annotations = [x for x in record.INFO['ANN']]
        for annotation_idx, annotation in enumerate(annotations):
            effects = annotation.split('|')
            #Currently filtering for only missense. Need to fix this.
            if effects[1] == 'missense_variant' and effects[13]:
                transcript = split_possible[6].partition('.')[0]
                var_chr_pos = "{}:{}".format(record.CHROM, record.POS)
                coverage = len([str(i) for i in pileup_truncated(rna_bam, record.CHROM, record.POS - 1, record.POS)])
                if coverage > 0:
                    pos, total_depth, var_depth = [[i.reference_pos + 1, i.get_num_aligned(), ''.join(i.get_query_sequences()).upper().count(str(record.ALT[0]))] for i in pileup_truncated(rna_bam, record.CHROM, record.POS - 1, record.POS)][0]
                    if int(var_depth) > int(args.required_coverage):
                        filtered_chrom_pos.append([record.CHROM, record.POS])
            elif re.search('conservative_inframe_deletion', effects[1]):
                transcript = split_possible[6].partition('.')[0]
                var_chr_pos = "{}:{}".format(record.CHROM, record.POS)
                reads = [i.get_overlap(record.POS - 1, record.POS + len(record.REF) - 1) for i in rna_bam.fetch(record.CHROM, record.POS - 1, record.POS + len(record.REF) - 1)]
                overlap = [i.get_overlap(record.POS - 1, record.POS + len(record.REF) - 1) for i in rna_bam.fetch(record.CHROM, record.POS - 1, record.POS + len(record.REF) - 1)]
                total_depth = len(reads)
                var_depth = overlap.count(len(record.ALT[0]))
                ref_depth = overlap.count(len(record.REF))
                if int(var_depth) > int(args.required_coverage):
                    filtered_chrom_pos.append([record.CHROM, record.POS])
            elif re.search('conservative_inframe_insertion', effects[1]):
                transcript = split_possible[6].partition('.')[0]
                var_chr_pos = "{}:{}".format(record.CHROM, record.POS)
                reads = [i.get_forward_sequence() for i in rna_bam.fetch(record.CHROM, record.POS - 1, record.POS + len(record.REF) - 1)]
                var_depth = len([x for x in reads if re.search(str(record.ALT[0]), x)])
                ref_depth = len([x for x in reads if not(re.search(str(record.ALT[0]), x))])
                print("{}\t{}\t{}".format(len(reads), ref_depth, var_depth))
                if int(var_depth) > int(args.required_coverage):
                    filtered_chrom_pos.append([record.CHROM, record.POS])
           #else: Need to consider frameshift indels here, too.

    vcf_reader = vcf.Reader(filename=args.vcf, compressed=False)
    for record in vcf_reader:
        if [record.CHROM, record.POS] in filtered_chrom_pos:
            filtered_vcf.append(record)

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
#    print("tx_thresshold: {}".format(tx_threshold))
#    print("# of expressed transcripts: {}".format(len(expressed_txids)))
#    print("Some expressed transcripts: {}".format(expressed_txids[:10]))
    filtered_records = filter_vcf_by_expression(args, expressed_txids)
#    print("# of filtered transcripts: {}".format(len(filtered_records)))
    write_expressed_vcf(args, filtered_records)

def load_tx_abundances(args):
    """
    """
    count_column = ''
    counts = np.array([])
    print("okay!")
    with open(args.quants) as fo:
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
    return threshold


def get_expressed_txs(args, threshold):
    """
    """
    count_column = ''
    txid_column = ''
    expressed_txids = []
    print("okay!")
    with open(args.quants) as fo:
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
        annotations = [x for x in record.INFO['ANN']]
        for annotation in annotations:
            effects = annotation.split('|')
            transcript = effects[6].partition('.')[0]
            #print("Transcript: {}".format(transcript))
            if transcript in expressed_txids:
                filtered_records.append(record)
    return tuple(set(filtered_records))


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

    #A few thoughts here - should really only be checking for germline variants
    #around somatic variants. It may be worth checking for proximal somatic
    #variants, but should probably result in a warning, if anything. This is also
    #where the het phased and hom/het synonymous checks can occur. The incorporation
    #of het phased variants will have to be in the make_*_peptides step.
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
                k = "{}_{}_{}_{}".format(line[1], line[10], line[0], len(line[9]))
                wt_nms[k] = line[15]
                wt_peptides[k] = line[9]

    mts_w_agreto = []

    header = ''

    with open(args.mt_fasta) as mto:
        for line in mto.readlines():
            line = line.rstrip('\n').split(' ')
            line = [x for x in line if x]
            if len(line) > 14 and line[0] == 'Pos':
                #Note: This is REFERENCE, _NOT_ GERMLINE!
                line.extend(['Reference Aff(nM)', 'Reference Peptide', 'Agretopocity'])
                header = ','.join(line)
            elif len(line) > 14 and line[0] not in ['Pos', 'Protein']:
                m = "{}_{}_{}_{}".format(line[1], line[10], line[0], len(line[9]))
                result = line[:16]
                if m in wt_nms.keys():
                    print("{}\t{}\t{}\t{}".format(m, line[15], wt_nms[m], float(line[15])/float(wt_nms[m])))
                    result.extend([wt_nms[m], wt_peptides[m], str(float(line[15])/float(wt_nms[m]))])
                else:
                    result.extend(["NA", "NA", "NA"])
                mts_w_agreto.append(','.join(result))

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
                aa_context = line[5]
                checksum_to_meta_map[checksum] = {}
                checksum_to_meta_map[checksum]["var_pos"] = var_pos
                checksum_to_meta_map[checksum]["tx_id"] = tx_id
                checksum_to_meta_map[checksum]["ref"] = ref
                checksum_to_meta_map[checksum]["alt"] = alt
                checksum_to_meta_map[checksum]["aa_context"] = aa_context

#    with open(args.mutant_nucs) as mno:
#        for line in mno.readlines():
#            if line.startswith('>'):
#                line = line.rstrip().split('\t')
#                checksum = line[0].lstrip('>')[:-1]
#                nuc_context = line[-1]
#                checksum_to_meta_map[checksum]["nuc_context"] = nuc_context

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
                tx_id_no_version = str(line).split('transcript_id "')[1].split('"')[0].split('.')[0]
                tx_to_gene[tx_id] = gene_name
                tx_to_gene[tx_id_no_version] = gene_name


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
    with open(args.binding_affinities) as mno:
        for line_idx, line in enumerate(mno.readlines()):
            line = line.rstrip().split(',')
            if line_idx == 0:
                header = line
                try:
                    header.remove('BindLevel')
                except:
                    pass
#                header.extend(['gene_name', 'tx_id', 'variant_position', 'reference_allele', 'alternate_allele', 'tpm', 'ccf', 'amino_acid_context', 'nucleotide_context'])
                header.extend(['gene_name', 'tx_id', 'variant_position', 'reference_allele', 'alternate_allele', 'tpm', 'ccf', 'amino_acid_context'])
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
            print(checksum_to_meta_map[csum])
            aa_context = checksum_to_meta_map[csum]['aa_context']
#            nuc_context = checksum_to_meta_map[csum]['nuc_context']
#            translated_nuc = str(Seq(nuc_context).translate()).replace('*', '')
#            if translated_nuc != aa_context:
#                print("One of your antigen's nucleotide context does not match the amino acid context. Check this!")
#                print(nuc_context)
#                print(translated_nuc)
#                print(aa_context)
#                sys.exit()
            ccf = 'NA'
            if var_pos in var_to_ccf.keys():
                ccf = var_to_ccf[var_pos]
            if 'SB' in line:
                line.remove('<=')
                line.remove('SB')
            if 'WB' in line:
                line.remove('WB')
#            line.extend([gene_name, tx_id, var_pos, ref, alt, tpm, ccf, aa_context, nuc_context])
            line.extend([gene_name, tx_id, var_pos, ref, alt, tpm, ccf, aa_context])
            new_lines.append(line)

    with open(args.output, 'w') as ofo:
        ofo.write("{}\n".format('\t'.join(header)))
        for new_line in new_lines:
            ofo.write("{}\n".format('\t'.join(new_line)))

def add_indel_metadata(args):
    """
    """
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
                tx_id_no_version = str(line).split('transcript_id "')[1].split('"')[0].split('.')[0]
                tx_to_gene[tx_id_no_version] = gene_name
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
    with open(args.binding_affinities) as mno:
        for line_idx, line in enumerate(mno.readlines()):
            print(line)
            line = line.rstrip()
            line = line.split()
            if not(line):
                continue
            if line[0] == 'Pos':
                header = line
                header.remove('BindLevel')
                header.extend(['gene_name', 'tx_id', 'variant_position', 'reference_allele', 'alternate allele', 'tpm', 'ccf'])
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
                ref_cn = str(copy_no[chr][segment]["maj_count"])
                alt_cn = str(copy_no[chr][segment]["min_count"])
                var = var.replace('_', ':')
                pyclone_inp.append("{}\n".format('\t'.join([var, args.samp_id, ref_depth, alt_depth, ref_cn, alt_cn, '2', '0.001', cellularity])))

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



def add_herv_metadata(args):
    """
    """
    #tx_abundances = load_tx_abundances(args)

    #print(tx_abundances)

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
                tx_to_tpm[line[tx_col_idx].split('.')[0][:15].replace(':','_')] = line[tpm_col_idx]

    print(tx_to_tpm)

    herv_expression = {}

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
                header.extend(['read_count'])
            if len(line) < 16 or line[0] in ['Protein', 'Pos']:
                continue
            else:
                if 'SB' in line:
                    line.remove('<=')
                    line.remove('SB')
                if 'WB' in line:
                    line.remove('WB')
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
    print("tx_thresshold: {}".format(tx_threshold))
    print("# of expressed transcripts: {}".format(len(expressed_txids)))
    print("Some expressed transcripts: {}".format(expressed_txids[:10]))
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
                if 'SB' in line:
                    line.remove('<=')
                    line.remove('SB')
                if 'WB' in line:
                    line.remove('WB')
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
    if args.command == 'make-snv-peptides':
        make_snv_peptides(args)
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
    if args.command == 'get-snv-genomic-context':
        get_snv_genomic_context(args)
    if args.command == 'expressed-hervs':
        expressed_hervs(args)
    if args.command == 'make-herv-peptides':
        make_herv_peptides(args)
    if args.command == 'add-herv-metadata':
        add_herv_metadata(args)
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
