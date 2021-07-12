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

#import pandas as pd
from glob import glob

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
    parser_make_snv_peptides.add_argument('-c', '--tx-cds-fasta',
                                          help="Transcript (CDS) nucleotide sequenes.",
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
    parser_make_snv_peptides.add_argument('-n', '--context-nt-length',
                                          help="Emitted contextual nucleotide sequence length (default: 39)",
                                          default=39)
    parser_make_snv_peptides.add_argument('-o', '--mt-output',
                                          help="Mutant peptides output file.",
                                          required=True)
    parser_make_snv_peptides.add_argument('-w', '--wt-output',
                                          #Note: Wildtype means REFERENCE! Need
                                          #to incorporate patient variants for
                                          #agretopicity calculations.
                                          help="Wildtype peptides output file.",
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
                                            default=7)
    parser_make_indel_peptides.add_argument('-o', '--output',
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


    # Subparser for making viral peptides
    parser_make_viral_peptides = subparsers.add_parser('make-viral-peptides',
                                                       help="Make viral peptides")
    parser_make_viral_peptides.add_argument('-e', '--expressed-viruses',
                                            help="File containing list of expressed viruses.",
                                            required=True)
    parser_make_viral_peptides.add_argument('-f', '--fasta',
                                            help="Viral CDS reference FASTA.",
                                            required=True)
    parser_make_viral_peptides.add_argument('-o', '--output',
                                            help='Output file.',
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

    # Subparser for making fusion peptides
    parser_make_fusion_peptides = subparsers.add_parser('make-fusion-peptides',
                                                        help="Make fusion peptides FASTA.")
    parser_make_fusion_peptides.add_argument('-f', '--fusions',
                                            help="Predicted fusions (STARFusion format).",
                                            required=True)
    parser_make_fusion_peptides.add_argument('-o', '--output',
                                             required=True)


    # Subparser for adding SNV metadata
    parser_add_snv_metadata = subparsers.add_parser('add-snv-metadata',
                                                    help="Add metadata to SNV binding affinity data.")
    parser_add_snv_metadata.add_argument('-m', '--mutant-peptides',
                                         help="FASTA file with mutant peptides.",
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

    # Subparser for filtering variants for expression
    parser_expressed_variants = subparsers.add_parser('filter-expressed-variants',
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
                                           help="Expression percentile for filtering (default: 90).",
                                           default=90)
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
    parser_covered_variants = subparsers.add_parser('check-rna-coverage',
                                                 help="Filter variants for evidence in RNA reads.")
    parser_covered_variants.add_argument('-b', '--tumor-rna-bam',
                                         help="Tumor RNA BAM file.",
                                         required=True)
    parser_covered_variants.add_argument('-c', '--required-coverage',
                                         help="Required coverage for variants (default: 1)",
                                         default=1)
    parser_covered_variants.add_argument('-f', '--peptide-fasta',
                                         help="LensTools-derived Neoantigen peptide fasta..",
                                         required=True)
    parser_covered_variants.add_argument('-t', '--tx-cds-fasta',
                                         help="Transcript (CDS) nucleotide sequences.",
                                         required=True)
    parser_covered_variants.add_argument('-o', '--output',
                                         help="Output file.",
                                         required=True)

    # Subparser for filtering isolated variants (e.g. no proximal germline or somatic variants).
    parser_isolated_variants = subparsers.add_parser('filter-isolated-variants',
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

    # Subparser for filtering expressed hERVs
    parser_expressed_hervs = subparsers.add_parser('filter-expressed-hervs',
                                                   help="Filter expressed hERVs.")
    parser_expressed_hervs.add_argument('-q', '--quants',
                                        help="Transcript abundance file (Salmon format).",
                                        required=True)
    parser_expressed_hervs.add_argument('-m', '--metric',
                                        help="Column for expression from abundance file.",
                                        default="TPM")
    parser_expressed_hervs.add_argument('-r', '--exclude-zeros',
                                        help="Exclude zeros for expression percentile.",
                                        action='store_false')
#    parser_expressed_hervs.add_argument('-p', '--percentile',
#                                        help="Expression percentile for filtering expression. (default: 50)",
#                                        default=50)
    parser_expressed_hervs.add_argument('-t', '--abundance-threshold',
                                        help="Expression threshold for filtering. (default: 10)")
    parser_expressed_hervs.add_argument('-o', '--output',
                                        help="Output file.",
                                        required=True)


    # Subparser for filtering virdetect outputs for expressed viruses
    parser_filter_virdetect_by_counts = subparsers.add_parser('filter-expressed-viruses',
                                                   help="Generate list of expressed viruses")
    parser_filter_virdetect_by_counts.add_argument('-q', '--viral-quants',
                                                   help="Viral counts file from VirDetect.",
                                                   required=True)
    parser_filter_virdetect_by_counts.add_argument('-r', '--viral-ref',
                                                   help="Viral reference FASTA (used for VirDetect).",
                                                   required=True)
    parser_filter_virdetect_by_counts.add_argument('-m', '--min-threshold',
                                                   help="Minimal count threshold for filtering (default: 1).")
    parser_filter_virdetect_by_counts.add_argument('-o', '--output',
                                                   help="Output file.",
                                                   required=True)

    # Subparser for filtering self-antigens for expression
    parser_expressed_self_genes = subparsers.add_parser('filter-expressed-self-genes',
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
                                             help="Expression percentile for filtering expression (default: 90).",
                                             default=90)
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




#    #This mode needs more information. How is it to be used? Does one simply
#    #pass it the VCF?
#    # Subparser for generating genomic context
#    parser_get_context = subparsers.add_parser('get-snv-genomic-context',
#                                               help="Make SNV neoAg nt context FASTA file.")
#    parser_get_context.add_argument('-c', '--tx-cds-fasta',
#                                      help="Transcript (CDS) nucleotide sequenes.",
#                                      required=True)
#    parser_get_context.add_argument('-v', '--vcf',
#                                      help="Annotated (snpEff) VCF).",
#                                      required=True)
#    parser_get_context.add_argument('-l', '--length',
#                                      help="Emitted nucleotide sequence length (default: 39)",
#                                      default=39)
#    parser_get_context.add_argument('-o', '--output',
#                                      help="output file.",
#                                      required=True)

    # Subparser for consolidating multiqc statistics files
    parser_consol_mqc_stats = subparsers.add_parser('consolidate-multiqc-stats',
                                                   help="Consolidate multiqc stats files into a single file.")
    parser_consol_mqc_stats.add_argument('-d', '--multiqc-data',
                                         help="Directory containing multiqc stats files.",
                                         required=True)
    parser_consol_mqc_stats.add_argument('-o', '--output',
                                         help="Output file.",
                                         required=True)


    # Subparser for adding RNA normals
    parser_add_rna_norms = subparsers.add_parser('add-rna-normals',
                                                 help="Add RNA normal runs to patient data.")
    parser_add_rna_norms.add_argument('-m', '--manifest',
                                      required=True)
    parser_add_rna_norms.add_argument('-p', '--prefix',
                                      required=True)
    parser_add_rna_norms.add_argument('-d', '--dataset',
                                      required=True)
    parser_add_rna_norms.add_argument('-n', '--pat_name',
                                      required=True)
    parser_add_rna_norms.add_argument('-o', '--output',
                                      required=True)

    # Subparser for adding splice metadata
    parser_add_rna_norms = subparsers.add_parser('add-splice-metadata',
                                                 help="Add splice metadata.")
    parser_add_rna_norms.add_argument('-s', '--splice-summary',
                                      required=True)
    parser_add_rna_norms.add_argument('-o', '--output',
                                      required=True)

    # Subparser for making LENS report
    parser_make_lens_report = subparsers.add_parser('make-lens-report')
    parser_make_lens_report.add_argument('--metadata-dir', '-d',
                                         help="Path to metadata directory.")
    parser_make_lens_report.add_argument('--output', '-o',
                                         help="Output file.")

    # Subparser for making LENS report
    parser_make_lens_report = subparsers.add_parser('make-antigens-barplot')
    parser_make_lens_report.add_argument('--reports-dir', '-d',
                                         help="Path containing reports.")
    parser_make_lens_report.add_argument('--output', '-o',
                                         help="Output file.")


    # Subparser for adding TCGA data
    parser_add_tcga_data = subparsers.add_parser('add-tcga-data')
    parser_add_tcga_data.add_argument('--report', '-r',
                                         help="LENS Report",
                                         required=True)
    parser_add_tcga_data.add_argument('--tumor-type', '-t',
                                         help="Tumor type (either full name or abbreviation).",
                                         required=True)
    parser_add_tcga_data.add_argument('--tcga-transcript-summary', '-s',
                                         help="TCGA transcript summary file.",
                                         required=True)
    parser_add_tcga_data.add_argument('--output', '-o',
                                         help="Output file.",
                                         required=True)


    return parser.parse_args()



def load_tcga_dicts():

    abbrev_to_type = {
    "LAML": "Acute Myeloid Leukemia",
    "ACC": "Adrenocortical carcinoma",
    "BLCA": "Bladder Urothelial Carcinoma",
    "LGG": "Brain Lower Grade Glioma",
    "BRCA": "Breast invasive carcinoma",
    "CESC": "Cervical squamous cell carcinoma and endocervical adenocarcinoma",
    "CHOL": "Cholangiocarcinoma",
    "LCML": "Chronic Myelogenous Leukemia",
    "COAD": "Colon adenocarcinoma",
    "CNTL": "Controls",
    "ESCA": "Esophageal carcinoma",
    "FPPP": "FFPE Pilot Phase II",
    "GBM": "Glioblastoma multiforme",
    "HNSC": "Head and Neck squamous cell carcinoma",
    "KICH": "Kidney Chromophobe",
    "KIRC": "Kidney renal clear cell carcinoma",
    "KIRP": "Kidney renal papillary cell carcinoma",
    "LIHC": "Liver hepatocellular carcinoma",
    "LUAD": "Lung adenocarcinoma",
    "LUSC": "Lung squamous cell carcinoma",
    "DLBC": "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma",
    "MESO": "Mesothelioma",
    "MISC": "Miscellaneous",
    "OV": "Ovarian serous cystadenocarcinoma",
    "PAAD": "Pancreatic adenocarcinoma",
    "PCPG": "Pheochromocytoma and Paraganglioma",
    "PRAD": "Prostate adenocarcinoma",
    "READ": "Rectum adenocarcinoma",
    "SARC": "Sarcoma",
    "SKCM": "Skin Cutaneous Melanoma",
    "STAD": "Stomach adenocarcinoma",
    "TGCT": "Testicular Germ Cell Tumors",
    "THYM": "Thymoma",
    "THCA": "Thyroid carcinoma",
    "UCS": "Uterine Carcinosarcoma",
    "UCEC": "Uterine Corpus Endometrial Carcinoma",
    "UVM": "Uveal Melanoma"}

    return abbrev_to_type

def load_tx_aas(args):
    """
    Loads a transcript -> amino acid sequence FASTA into a dictionary.
    """
    tx_to_aa = {}
    for seq_record in SeqIO.parse(args.tx_aa_fasta, "fasta"):
        tx = re.search('transcript:\S*', seq_record.description).group(0).split(':')[1]
        tx_no_version = tx.split('.')[0]
        tx_to_aa[tx] = seq_record.seq
        tx_to_aa[tx_no_version] = seq_record.seq
    return tx_to_aa


def load_tx_cds(args):
    """
    Loads a transcript -> nucleotide sequence FASTA into a dictionary.
    """
    tx_to_cds = {}
    for seq_record in SeqIO.parse(args.tx_cds_fasta, "fasta"):
        tx = seq_record.id
        tx_no_version = tx.split('.')[0]
        tx_to_cds[tx] = seq_record.seq
        tx_to_cds[tx_no_version] = seq_record.seq
    return tx_to_cds


def extract_missense_snvs(input_vcf):
    """
    Extracts missense SNVs from annotated somatic VCF file.
    """
    missense_records = {}

    print(input_vcf)

    vcf_reader = ''
    if input_vcf.endswith('gz'):
        vcf_reader = vcf.Reader(open(input_vcf), 'r', compressed=True)
    else:
        vcf_reader = vcf.Reader(open(input_vcf), 'r')

    for record in vcf_reader:
        print(record)
        annotations = [x for x in record.INFO['ANN']]
        for annotation in annotations:
            effects = annotation.split('|')
            variant_effect = effects[1]
            if variant_effect == 'missense_variant' and effects[13]:
                transcript = effects[6]

                nt_change = effects[9].lstrip('c.')
                nt_change = nt_change[-3:].partition('>')

                ref_nt = nt_change[0]
                alt_nt = nt_change[2]

                nt_pos = effects[12].split('/')[0]
                nt_len = effects[12].split('/')[1]

                aa3_change = effects[10].lstrip('p.')
                aa3_change = re.split('(\d+)', aa3_change)

                ref_aa3 = aa3_change[0]
                alt_aa3 = aa3_change[2]

                aa_pos = aa3_change[1]

                ref_aa = seq1(ref_aa3)
                alt_aa = seq1(alt_aa3)
                aa_len = effects[13].split('/')[1]

                codon_pos = int((int(effects[13].split('/')[0])*3) - 2)

                #filtered_records.append([transcript, pos, tlen, orig_aa, alt_aa, codon_pos, record])
                missense_records[str(record)] = {'transcript': transcript,
                                            'nt_pos': nt_pos,
                                            'nt_len': nt_len,
                                            'ref_nt': ref_nt,
                                            'alt_nt': alt_nt,
                                            'aa_pos': aa_pos,
                                            'aa_len': aa_len,
                                            'ref_aa': ref_aa,
                                            'alt_aa': alt_aa,
                                            'codon_pos': codon_pos,
                                            'meta': record}
    return missense_records

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

def extract_conservative_inframe_indels(input_vcf):
    """
    This will need to be cleaned heavily, but good for now.

    Focusing on deletions now, need to find a good example of a conservative_inframe_insertion
    """
    conserv_inframe_indels = {}
    vcf_reader = ''
    if input_vcf.endswith('gz'):
        vcf_reader = vcf.Reader(open(input_vcf), 'r', compressed=True)
    else:
        vcf_reader = vcf.Reader(open(input_vcf), 'r')

    for record in vcf_reader:
        annotations = [x for x in record.INFO['ANN']]
        for annotation in annotations:
            effects = annotation.split('|')
            if re.search('^conservative_inframe_[a-z]+$', effects[1]) and effects[13]:
                transcript = effects[6]
                aa3_change = effects[10].lstrip('p.')
                aa_len = effects[13].split('/')[1]
                nt_change = effects[9].lstrip('c.')
                nt_len = effects[12].split('/')[1]
                conserv_inframe_indels[str(record)] = {'transcript': transcript,
                                                       'aa3_change': aa3_change,
                                                       'nt_change': nt_change,
                                                       'aa_len': aa_len,
                                                       'nt_len': nt_len,
                                                       'meta': record}
    return conserv_inframe_indels


def extract_disruptive_inframe_indels(input_vcf):
    """
    This will need to be cleaned heavily, but good for now.

    Focusing on deletions now, need to find a good example of a conservative_inframe_insertion
    """
    disrupt_inframe_indels = {}
    if input_vcf.endswith('gz'):
        vcf_reader = vcf.Reader(open(input_vcf), 'r', compressed=True)
    else:
        vcf_reader = vcf.Reader(open(input_vcf), 'r')

    for record in vcf_reader:
        annotations = [x for x in record.INFO['ANN']]
        for annotation in annotations:
            effects = annotation.split('|')
            if re.search('^disruptive_inframe_[a-z]+$', effects[1]) and effects[13]:
                transcript = effects[6]
                aa3_change = effects[10].lstrip('p.')
                aa_len = effects[13].split('/')[1]
                nt_change = effects[9].lstrip('c.')
                nt_len = effects[12].split('/')[1]
                disrupt_inframe_indels[str(record)] = {'transcript': transcript,
                                                       'aa3_change': aa3_change,
                                                       'nt_change': nt_change,
                                                       'aa_len': aa_len,
                                                       'nt_len': nt_len,
                                                       'meta': record}
    return disrupt_inframe_indels

def extract_frameshift_indels(input_vcf):
    """
    This will need to be cleaned heavily, but good for now.

    Focusing on deletions now, need to find a good example of a conservative_inframe_insertion
    """
    frameshift_indels = {}
    if input_vcf.endswith('gz'):
        vcf_reader = vcf.Reader(open(input_vcf), 'r', compressed=True)
    else:
        vcf_reader = vcf.Reader(open(input_vcf), 'r')

    for record in vcf_reader:
        annotations = [x for x in record.INFO['ANN']]
        for annotation in annotations:
            effects = annotation.split('|')
            if re.search('^frameshift_variant$', effects[1]) and effects[13] and effects[7] == 'protein_coding' and effects[-1] == '':
                transcript = effects[6]
                aa_len = effects[13].split('/')[1]
                nt_change = effects[9].lstrip('c.')
                nt_len = effects[12].split('/')[1]
                frameshift_indels[str(record)] = {'transcript': transcript,
                                                          'nt_change': nt_change,
                                                          'aa_len': aa_len,
                                                          'nt_len': nt_len,
                                                          'meta': record}
    return frameshift_indels


def make_snv_peptides(args):
    """
    Make peptide sequences derived from mutant (somatic) and normal (germline) sequences.
    """
    tx_to_aa = load_tx_aas(args)
    tx_to_cds = load_tx_cds(args)

    missense_snvs = extract_missense_snvs(args.somatic_vcf)
    print(missense_snvs)

    mutant_peptides = {}
    reference_peptides = {}

    for entry in missense_snvs.keys():
        record = missense_snvs[entry]
        print(record)
        tx = record['transcript']
        tx_no_version = tx.split('.')[0]
        if tx not in tx_to_aa.keys() and tx_no_version not in tx_to_aa.keys():
            print("{} cannot be found in the transcript-to-amino acid dictionary. Continuing to next record...".format(record['transcript']))
            continue
        if tx not in tx_to_cds.keys() and tx_no_version not in tx_to_cds.keys():
            print("{} cannot be found in the transcript-to-coding sequence dictionary. Continuing to next record...".format(record['transcript']))
            continue

        try:
            bufr_aa = list(tx_to_aa[tx])
        except:
            bufr_aa = list(tx_to_aa[tx_no_version])
        try:
            bufr_nt = list(tx_to_cds[tx])
        except:
            bufr_nt = list(tx_to_cds[tx_no_version])

        if len(bufr_aa) != int(record['aa_len']):
            print("transcript {} shows different lengths between amino acid fasta ({}) and snpEff annotations! ({})".format(tx, record['aa_len'], len(bufr_aa)))
        if len(bufr_nt) != int(record['nt_len']):
            print("transcript {} shows different lengths between coding sequence fasta ({}) and snpEff annotations! ({})".format(tx, record['nt_len'], len(bufr_nt)))

        # Deep copy to mut_aa...
        # This is the step where the annotated germline reference should be
        # utilized...
        mut_aa = bufr_aa[:]
        mut_nt = bufr_nt[:]
        aa_pos = int(record['aa_pos']) - 1
        nt_pos = int(record['nt_pos']) - 1
        # This check should be performed prior to incorporating germline variants.

        if bufr_aa[aa_pos] != record['ref_aa']:
            print("Reference amino acid from snpEff annotation ({}) doesn't match amino acid from amino acid FASTA ({})! Continuing to next record...".format(record['ref_aa'], bufr_aa[aa_pos]))
            continue
        if bufr_nt[nt_pos] != record['ref_nt']:
            print("Reference nucleotide from snpEff annotation ({}) doesn't match nucleotide from amino acid FASTA ({})! Continuing to next record...".format(record['ref_nt'], bufr_nt[nt_pos]))
            continue

        # Applying the mutated amino acid
        mut_aa[aa_pos] = record['alt_aa']
        mut_nt[nt_pos] = record['alt_nt']

        ref_peptide = ''.join(bufr_aa[max(aa_pos-args.length+1,0):min(aa_pos+args.length, len(bufr_aa))])
        mut_peptide = ''.join(mut_aa[max(aa_pos-args.length+1, 0):min(aa_pos+args.length, len(mut_aa))])

        print("Reference peptide: {}".format(ref_peptide))
        print("Mutant peptide:    {}".format(mut_peptide))

        codon_pos = record['codon_pos'] - 1
        ref_nuc = ''.join(bufr_nt[max(codon_pos-args.context_nt_length, 0):min(codon_pos+args.context_nt_length, len(bufr_nt)) + 3])
        mut_nuc = ''.join(mut_nt[max(codon_pos-args.context_nt_length, 0):min(codon_pos+args.context_nt_length, len(mut_nt)) + 3])

        print("Reference genomic context: {}".format(ref_nuc))
        print("Mutant genomic context:    {}".format(mut_nuc))

        header_mut_peptide = ''.join(mut_aa[max(aa_pos-13, 0):min(aa_pos+14, len(mut_aa))])

        #var_md5 is being used to create a unique identifier for the resulting mutant peptide to overcome netMHCpan's length limitations.
        var_md5 = hashlib.md5("{}".format(':'.join([str(record['meta'].CHROM), str(record['meta'].POS), str(tx), str(record['meta'].REF), str(record['meta'].ALT)])).encode('utf-8')).hexdigest()[:16]
        mutant_peptides["MD5:{} VARIANT_POS:{}_{} TRANSCRIPT:{} REF:{} ALT:{} SNV_TYPE:{} PROTEIN_CONTEXT:{} GENOMIC_CONTEXT:{}".format(var_md5, record['meta'].CHROM, record['meta'].POS, tx, record['meta'].REF, record['meta'].ALT, 'missense', header_mut_peptide, mut_nuc)] = mut_peptide
        reference_peptides["MD5:{} VARIANT_POS:{}_{} TRANSCRIPT:{} REF:{} ALT:{}".format(var_md5, record['meta'].CHROM, record['meta'].POS, tx, record['meta'].REF, record['meta'].ALT)] = ref_peptide

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



def get_indel_peptide_header_and_aa(args):
    """
    """
    pass



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


    # We will need to retrieve the nucleotide sequence records as well.
    print("Getting indels...")
    conserv_inframe_indels = extract_conservative_inframe_indels(args.somatic_vcf)
    print("Retrieved {} conservative inframe indels...".format(len(conserv_inframe_indels)))
    disrupt_inframe_indels = extract_disruptive_inframe_indels(args.somatic_vcf)
    print("Retrieved {} disruptive inframe indels...".format(len(disrupt_inframe_indels)))
    frameshift_indels = extract_frameshift_indels(args.somatic_vcf)
    print("Retrieved {} frameshifte indels...".format(len(frameshift_indels)))

    mutant_peptides = {}

    conserv_disrupt_indels = dict(conserv_inframe_indels, **disrupt_inframe_indels)


    # CONSERVATIVE/DISRUPTIVE #
    for entry in conserv_disrupt_indels:
        record = conserv_disrupt_indels[entry]
        print("Current indel: {}".format(record))

        tx_no_version = record['transcript'].split('.')[0]
        if record['transcript'] not in tx_to_aa.keys() and tx_no_version not in tx_nv_to_aa.keys():
            continue

        tx_ref_aa = []
        if record['transcript'] in tx_to_aa.keys():
            tx_ref_aa = list(tx_to_aa[record['transcript']])
            tx_ref_nt = list(tx_to_cds[record['transcript']])
        elif tx_no_version in tx_nv_to_aa.keys():
            tx_ref_aa = list(tx_nv_to_aa[tx_no_version])
            tx_ref_nt = list(tx_nv_to_cds[tx_no_version])

        bufr_aa = tx_ref_aa[:]
        bufr_nt = tx_ref_nt[:]

        # DELETIONS #
        if re.search("del$", record['aa3_change']):
            del_rec_aa = record['aa3_change'].strip('del')
            del_rec_nt = record['nt_change'].strip('del')

            start_pos_aa = int(del_rec_aa.split('_')[0][3:]) - 1
            start_aa3 = del_rec_aa.split('_')[0][:3]
            start_aa = seq1(start_aa3)
            stop_pos_aa = 0
            stop_aa3 = 'foo'
            stop_aa = 'foo'

            print(del_rec_nt)
            bufr1, start_pos_nt, bufr2, stop_pos_nt, seq = re.split('(\d+)', del_rec_nt)

            # Post processing
            del_nt = seq.replace('del', '')
            start_pos_nt = int(start_pos_nt) - 1
            stop_pos_nt = int(stop_pos_nt) - 1

            if len(del_nt) > 0:
                start_nt = del_nt[0]
                stop_nt = del_nt[-1]
            else:
                start_nt = False
                stop_nt = False

            print("{} {} {} {} {}".format(start_pos_nt, stop_pos_nt, seq, start_nt, stop_nt))


            # Handling single amino acid deletions
            if re.search('_', del_rec_aa):
                stop_pos_aa = int(del_rec_aa.split('_')[1][3:]) - 1
                stop_aa3 = del_rec_aa.split('_')[1][:3]
                stop_aa = seq1(stop_aa3)
            else:
                stop_pos_aa = start_pos_aa
                stop_aa3 = start_aa3
                stop_aa = start_aa

            print("{} {} {}".format(start_pos_aa, start_aa, start_aa3))
            print("{} {} {}".format(stop_pos_aa, stop_aa, stop_aa3))

            if start_pos_aa > len(bufr_aa) or stop_pos_aa > len(bufr_aa):
                print("The mutation occurs outside of the peptide sequence. Check transcript versions. Continuing with next record...")
                continue
            if start_pos_nt > len(bufr_nt) or stop_pos_nt > len(bufr_nt):
                print("The mutation occurs outside of the coding sequence. Check transcript versions. Continuing with next record...")
                continue

            if bufr_aa[start_pos_aa] != start_aa or bufr_aa[stop_pos_aa] != stop_aa:
                print("Either first amino acid doesn't match (FASTA: {}, snpEff: {}) or last amino acid doesn't match (FASTA: {}, snpEff: {}).".format(bufr_aa[start_pos_aa], start_aa, bufr_aa[stop_pos_aa], stop_aa))
                print("Continuing with next record...")
                continue

            if (start_nt and stop_nt) and (bufr_nt[start_pos_nt] != start_nt or bufr_nt[stop_pos_nt] != stop_nt):
                print("Either first nucleotide doesn't match (FASTA: {}, snpEff: {}) or last nucleotide doesn't match (FASTA: {}, snpEff: {}).".format(bufr_nt[start_pos_nt], start_nt, bufr_nt[stop_pos_nt], stop_nt))
                print("Continuing with next record...")
                continue

            # Applying deletion
            mut_aa = ''.join(bufr_aa[:start_pos_aa] + bufr_aa[stop_pos_aa+1:])
            mut_nt = ''.join(bufr_nt[:start_pos_nt] + bufr_nt[stop_pos_nt+1:])

            # Extraction sequence of interest from reference and mutated sequences.
            # Should we getting the germline sequence here? May not matter as we're not emitting it.
            ref_peptide = ''.join(bufr_aa[max(start_pos_aa-args.length,0):min(stop_pos_aa+args.length, len(bufr_aa))])
            mut_peptide = ''.join(mut_aa[max(start_pos_aa-args.length, 0):min(start_pos_aa+args.length, len(mut_aa))])

            aa_context = ''.join(mut_aa[max(start_pos_aa-13, 0):min(start_pos_aa+14, len(mut_aa))])
            nt_context = ''.join(mut_nt[max(start_pos_nt-39, 0):min(start_pos_nt+40, len(mut_nt))])

            print("Reference peptide: {}".format(ref_peptide))
            print("Mutant peptide:    {}".format(mut_peptide))
            print("Protein context:   {}".format(aa_context))
            print("Genomic context:   {}".format(nt_context))

            genomic_context_lower_pos = record['meta'].POS - 39
            genomic_context_upper_pos = record['meta'].POS + 40

            if check_mut_pep_in_nt_context(mut_peptide, nt_context):
                md5able_str = "{}".format(':'.join([str(record['meta'].CHROM), str(record['meta'].POS), str(record['transcript']), str(record['meta'].REF), str(record['meta'].ALT)])).encode('utf-8')
                var_md5 = hashlib.md5(md5able_str).hexdigest()[:16]
                mutant_peptides["MD5:{} VARIANT_POS:{}_{} TRANSCRIPT:{} REF:{} ALT:{} INDEL_TYPE:{} PROTEIN_CONTEXT:{} GENOMIC_CONTEXT:{} GENOMIC_CONTEXT_RANGE:{}_{}-{}".format(var_md5, record['meta'].CHROM, record['meta'].POS, record['transcript'], record['meta'].REF, record['meta'].ALT, 'inframe_deletion', aa_context, nt_context, record['meta'].CHROM, genomic_context_lower_pos, genomic_context_upper_pos)] = mut_peptide

        # INSERTIONS #
        elif re.search("[0-9]ins", record['aa3_change']):
            print(record)
            print("HERE!")

            ins_rec_aa, bufr, ins_aa3 = record['aa3_change'].partition('ins')
            ins_aa = seq1(ins_aa3)
            print(ins_aa)
            insert_pos_aa = int(ins_rec_aa.split('_')[0][3:]) - 1
            start_aa = seq1(ins_rec_aa.split('_')[0][:3])
            insert_plus_one_pos_aa= int(ins_rec_aa.split('_')[1][3:]) - 1
            start_plus_one_aa = seq1(ins_rec_aa.split('_')[1][:3])

            ins_nt = ''
            start_pos_nt = ''
            stop_pos_nt = ''
            mut_nt = ''

            if re.search('ins', record['nt_change']):
                ins_rec_nt, bufr, ins_nt = record['nt_change'].partition('ins')
                start_pos_nt = int(ins_rec_nt.split('_')[0]) - 1
                stop_pos_nt = int(ins_rec_nt.split('_')[1]) - 1
                mut_nt = ''.join(bufr_nt[:start_pos_nt+1] + [ins_nt] + bufr_nt[stop_pos_nt:])
                print("{}\n{}\n{}\n{}".format(ins_nt, start_pos_nt, stop_pos_nt, mut_nt[start_pos_nt-20:stop_pos_nt+20]))
            elif re.search('dup', record['nt_change']):
                dup_rec_nt, bufr, ins_nt = record['nt_change'].partition('dup')
                if not(re.search('-', dup_rec_nt)):
                    start_pos_nt = int(dup_rec_nt.split('_')[0]) - 1
                    stop_pos_nt = int(dup_rec_nt.split('_')[1]) - 1
                    mut_nt = ''.join(bufr_nt[:start_pos_nt-1] + [ins_nt] + bufr_nt[start_pos_nt:])
                    print("{}\n{}\n{}\n{}".format(ins_nt, start_pos_nt, stop_pos_nt, mut_nt[start_pos_nt-20:stop_pos_nt+20]))
                else:
                    continue



            if bufr_aa[insert_pos_aa] != start_aa or bufr_aa[insert_plus_one_pos_aa] != start_plus_one_aa:
                print("Insertion site amino acid doesn't match (FASTA: {}, snpEff: {}) or insertion site + 1 amino acid doesn't match (FASTA: {}, snpEff {})".format(bufr_aa[insert_pos_aa], start_aa, bufr_aa[insert_plus_one_pos_aa], start_plus_one_aa))
                print("Continuing with next record...")
                continue

            mut_aa = ''.join(bufr_aa[:insert_pos_aa+1] + [ins_aa] + bufr_aa[insert_plus_one_pos_aa:])

            # Should be args.length
            ref_peptide = ''.join(bufr_aa[max(insert_pos_aa-6,0):min(insert_plus_one_pos_aa+7, len(bufr_aa))])
            mut_peptide = ''.join(mut_aa[max(insert_pos_aa-6, 0):min(insert_plus_one_pos_aa+len(ins_aa)+7, len(mut_aa))])


            aa_context = ''.join(mut_aa[max(insert_pos_aa-19, 0):min(insert_plus_one_pos_aa+len(ins_aa)+20, len(mut_aa))])
            #nt_context won't need len(ins_nt) in case of dups.
            nt_context = ''.join(mut_nt[max(start_pos_nt-39, 0):min(stop_pos_nt+len(ins_nt)+40, len(mut_nt))])

            genomic_context_lower_pos = record['meta'].POS - 39
            genomic_context_upper_pos = record['meta'].POS + 40

            print("Reference peptide: {}".format(ref_peptide))
            print("Mutant peptide:    {}".format(mut_peptide))

            if check_mut_pep_in_nt_context(mut_peptide, nt_context):
                md5able_str = "{}".format(':'.join([str(record['meta'].CHROM), str(record['meta'].POS), str(record['transcript']), str(record['meta'].REF), str(record['meta'].ALT)])).encode('utf-8')
                var_md5 = hashlib.md5(md5able_str).hexdigest()[:16]
                mutant_peptides["MD5:{} VARIANT_POS:{}_{} TRANSCRIPT:{} REF:{} ALT:{} INDEL_TYPE:{} PROTEIN_CONTEXT:{} GENOMIC_CONTEXT:{} GENOMIC_CONTEXT_RANGE:{}_{}-{}".format(var_md5, record['meta'].CHROM, record['meta'].POS, record['transcript'], record['meta'].REF, record['meta'].ALT, 'inframe_insertion', aa_context, nt_context, record['meta'].CHROM, genomic_context_lower_pos, genomic_context_upper_pos)] = mut_peptide

        # INSERTIONS/DELETIONS #
        elif re.search("delins", record['aa3_change']):
            del_rec_aa, bufr, ins_aa3 = record['aa3_change'].partition('delins')

            ins_aa = seq1(ins_aa3)
            start_pos_aa = int(del_rec_aa.split('_')[0][3:]) - 1
            start_aa3 = del_rec_aa.split('_')[0][:3]
            start_aa = seq1(start_aa3)
            stop_pos_aa = 0
            stop_aa = 'foo'
            # Dealing with single amino acid deletions
            if re.search('_', del_rec_aa):
                stop_pos_aa = int(del_rec_aa.split('_')[1][3:]) - 1
                stop_aa3 = del_rec_aa.split('_')[1][:3]
                stop_aa = seq1(stop_aa3)
            else:
                stop_pos_aa = start_pos_aa
                stop_aa3 = start_aa3
                stop_aa = start_aa

            print("{} {} {}".format(start_pos_aa, start_aa, start_aa3))
            print("{} {} {}".format(stop_pos_aa, stop_aa, stop_aa3))


            ins_nt = '' # For 'dup' cases...
            mut_nt = ''
            start_pos_nt = 0
            stop_pos_nt = 0
            if re.search('ins', record['nt_change']):
                ins_rec_nt, bufr, ins_nt = record['nt_change'].partition('ins')
                start_pos_nt = int(ins_rec_nt.split('_')[0]) - 1
                stop_pos_nt = int(ins_rec_nt.split('_')[1]) - 1
                mut_nt = ''.join(bufr_nt[:start_pos_nt+1] + [ins_nt] + bufr_nt[stop_pos_nt:])
            elif re.search('dup', record['nt_change']):
                dup_rec_nt, bufr, ins_nt = record['nt_change'].partition('dup')
                start_pos_nt = int(dup_rec_nt.split('_')[0]) - 1
                stop_pos_nt = int(dup_rec_nt.split('_')[1]) - 1
                mut_nt = ''.join(bufr_nt[:start_pos_nt-1] + [ins_nt] + bufr_nt[start_pos_nt:])
#            print("{}\n{}\n{}\n{}".format(mut_nt, start_pos_nt, stop_pos_nt, mut_nt[start_pos_nt-20:stop_pos_nt+20]))

            if bufr_aa[start_pos_aa] != start_aa or bufr_aa[stop_pos_aa] != stop_aa:
                print("Either first amino acid doesn't match (FASTA: {}, snpEff: {}) or last amino acid doesn't match (FASTA: {}, snpEff: {}).".format(bufr_aa[start_pos_aa], start_aa, bufr_aa[stop_pos_aa], stop_aa))
                print("Continuing with next record...")
                continue


            print("{}".format(''.join(bufr_aa[start_pos_aa-10:start_pos_aa+10])))
            del_mut_aa = ''.join(bufr_aa[:start_pos_aa] + bufr_aa[stop_pos_aa+1:])
            print("{}".format(''.join(del_mut_aa[start_pos_aa-10:start_pos_aa+10])))
            delins_mut_aa = ''.join(del_mut_aa[:start_pos_aa] + ins_aa + del_mut_aa[start_pos_aa:])

            aa_context = ''.join(delins_mut_aa[max(start_pos_aa-19, 0):min(stop_pos_aa+len(ins_aa)+20, len(delins_mut_aa))])
            #nt_context won't need len(ins_nt) in case of dups.
            nt_context = ''.join(mut_nt[max(start_pos_nt-39, 0):min(stop_pos_nt+len(ins_nt)+40, len(mut_nt))])

            genomic_context_lower_pos = record['meta'].POS - 39
            genomic_context_upper_pos = record['meta'].POS + 40

            mut_peptide = ''.join(delins_mut_aa[max(start_pos_aa-7, 0):min(start_pos_aa+len(ins_aa)+7, len(delins_mut_aa))])
            ref_peptide = ''.join(bufr_aa[max(start_pos_aa-7,0):min(stop_pos_aa+len(ins_aa)+7, len(bufr_aa))])

            print("Reference peptide: {}".format(ref_peptide))
            print("mutant peptide:    {}".format(mut_peptide))
            print("Protein context:   {}".format(aa_context))
            print("Genomic context:   {}".format(nt_context))

#            if check_mut_pep_in_nt_context(mut_peptide, nt_context):
            md5able_str = "{}".format(':'.join([str(record['meta'].CHROM), str(record['meta'].POS), str(record['transcript']), str(record['meta'].REF), str(record['meta'].ALT)])).encode('utf-8')
            var_md5 = hashlib.md5(md5able_str).hexdigest()[:16]
            mutant_peptides["MD5:{} VARIANT_POS:{}_{} TRANSCRIPT:{} REF:{} ALT:{} INDEL_TYPE:{} PROTEIN_CONTEXT:{} GENOMIC_CONTEXT:{} GENOMIC_CONTEXT_RANGE:{}_{}-{}".format(var_md5, record['meta'].CHROM, record['meta'].POS, record['transcript'], record['meta'].REF, record['meta'].ALT, 'inframe_indel', aa_context, nt_context, record['meta'].CHROM, genomic_context_lower_pos, genomic_context_upper_pos)] = mut_peptide

        # DUPLICATIONS #
        elif re.search("dup", record['aa3_change']):
            dup_rec_aa = record['aa3_change'].strip('dup')
            start_pos_aa = int(dup_rec_aa.split('_')[0][3:]) - 1
            start_aa3 = dup_rec_aa.split('_')[0][:3]
            start_aa = seq1(start_aa3)
            stop_pos_aa = 0
            stop_aa3 = 'foo'
            stop_aa = 'foo'
            # Dealing with single amino acid deletions
            if re.search('_', dup_rec_aa):
                stop_pos_aa = int(dup_rec_aa.split('_')[1][3:]) - 1
                stop_aa3 = dup_rec_aa.split('_')[1][:3]
                stop_aa = seq1(stop_aa3)
            else:
                stop_pos_aa = start_pos_aa
                stop_aa3 = start_aa3
                stop_aa = start_aa

            print("{} {} {}".format(start_pos_aa, start_aa, start_aa3))
            print("{} {} {}".format(stop_pos_aa, stop_aa, stop_aa3))

            if re.search('ins', record['nt_change']):
                ins_rec_nt, bufr, ins_nt = record['nt_change'].partition('ins')
                start_pos_nt = int(ins_rec_nt.split('_')[0]) - 1
                stop_pos_nt = int(ins_rec_nt.split('_')[1]) - 1
                mut_nt = ''.join(bufr_nt[:start_pos_nt+1] + [ins_nt] + bufr_nt[stop_pos_nt:])
            elif re.search('dup', record['nt_change']):
                dup_rec_nt, bufr, ins_nt = record['nt_change'].partition('dup')
                start_pos_nt = int(dup_rec_nt.split('_')[0]) - 1
                stop_pos_nt = int(dup_rec_nt.split('_')[1]) - 1
                mut_nt = ''.join(bufr_nt[:start_pos_nt] + [ins_nt] + bufr_nt[start_pos_nt:])

            if start_pos_aa > len(bufr_aa) or stop_pos_aa > len(bufr_aa):
                print("The mutation occurs outside of the peptide sequence. Check transcript versions.")
                continue

            if bufr_aa[start_pos_aa] != start_aa or bufr_aa[stop_pos_aa] != stop_aa:
                print("Either first amino acid doesn't match (FASTA: {}, snpEff: {}) or last amino acid doesn't match (FASTA: {}, snpEff: {}).".format(bufr_aa[start_pos_aa], start_aa, bufr_aa[stop_pos_aa], stop_aa))
                print("Continuing with next record...")
                continue
            duped_aa = bufr_aa[start_pos_aa:stop_pos_aa+1]
            mut_aa = bufr_aa[:start_pos_aa] + duped_aa + bufr_aa[start_pos_aa:]

            mut_peptide = ''.join(mut_aa[max(stop_pos_aa-6, 0):min(start_pos_aa+len(duped_aa)+7, len(mut_aa))])

            aa_context = ''.join(mut_aa[max(start_pos_aa-19, 0):min(stop_pos_aa+len(duped_aa)+20, len(mut_aa))])
            #nt_context won't need len(ins_nt) in case of dups.
            rt_context = ''.join(bufr_nt[max(start_pos_nt-39, 0):min(stop_pos_nt+len(ins_nt)+40, len(bufr_nt))])
            nt_context = ''.join(mut_nt[max(start_pos_nt-39, 0):min(stop_pos_nt+len(ins_nt)+40, len(mut_nt))])

            genomic_context_lower_pos = record['meta'].POS - 39
            genomic_context_upper_pos = record['meta'].POS + 40

            print("ReferenceGenomic context: {}".format(rt_context))
            print("Mutant Genomic context:   {}".format(nt_context))
            print("Duped seq:       {}".format(''.join(duped_aa)))
            print("Reference protein context: {}".format(''.join(bufr_aa[start_pos_aa-15:stop_pos_aa+15])))
            print("Mutant protein context:    {}".format(''.join(mut_aa[start_pos_aa-15:stop_pos_aa+15])))
            print("Mutant peptide:            {}".format(mut_peptide))

            if check_mut_pep_in_nt_context(mut_peptide, nt_context):
                md5able_str = "{}".format(':'.join([str(record['meta'].CHROM), str(record['meta'].POS), str(record['transcript']), str(record['meta'].REF), str(record['meta'].ALT)])).encode('utf-8')
                var_md5 = hashlib.md5(md5able_str).hexdigest()[:16]
                mutant_peptides["MD5:{} VARIANT_POS:{}_{} TRANSCRIPT:{} REF:{} ALT:{} INDEL_TYPE:{} PROTEIN_CONTEXT:{} GENOMIC_CONTEXT:{} GENOMIC_CONTEXT_RANGE:{}_{}-{}".format(var_md5, record['meta'].CHROM, record['meta'].POS, record['transcript'], record['meta'].REF, record['meta'].ALT, 'inframe_dup', aa_context, nt_context, record['meta'].CHROM, genomic_context_lower_pos, genomic_context_upper_pos)] = mut_peptide

    # FRAMESHIFTS #
    for entry in frameshift_indels.keys():
        record = frameshift_indels[entry]
        tx_no_version = record['transcript'].split('.')[0]
        if record['transcript'] not in tx_to_aa.keys() and tx_no_version not in tx_nv_to_aa.keys():
            continue
        tx_ref_nt = []
        if record['transcript'] in tx_to_cds.keys():
            tx_ref_nt = list(tx_to_cds[record['transcript']])
        elif tx_no_version in tx_nv_to_cds.keys():
            tx_refnt = list(tx_nv_to_cds[tx_no_version])
        bufr_nt = tx_ref_nt[:]

        # DELETIONS #
        if re.search("del", record['nt_change']):
            print("Frameshift deletion {}".format(record))
            del_rec_nt, buffer, del_nt = record['nt_change'].partition('del')
            start_pos_nt = 0
            stop_pos_nt = 0
            start_base = 'A'
            stop_base = 'A'
            if re.search('_', del_rec_nt):
                 start_pos_nt = int(del_rec_nt.split('_')[0]) - 1
                 stop_pos_nt = int(del_rec_nt.split('_')[1]) - 1
                 if del_nt:
                     start_nt = del_nt[0]
                     stop_nt = del_nt[-1]
            else:
                start_pos_nt = int(del_rec_nt) - 1
                stop_pos_nt = start_pos_nt
                if del_nt:
                    start_nt = del_nt[0]
                    stop_nt = start_base
            print("{} {}".format(start_pos_nt, start_nt))
            print("{} {}".format(stop_pos_nt, stop_nt))
            if stop_pos_nt > len(bufr_nt):
                print("The mutation occurs outside of the trancsript sequence. Check transcript versions.")
                continue

            mut_nt = bufr_nt[:start_pos_nt] + bufr_nt[stop_pos_nt+1:]
            print("Reference NT: {}".format(''.join(bufr_nt[start_pos_nt-10:stop_pos_nt+10])))
            print("Mutant    NT: {}".format(''.join(mut_nt[start_pos_nt-10:stop_pos_nt+10])))
            ref_aa = str(Seq(''.join(bufr_nt)).translate(to_stop=True, cds=False))
            mut_aa = str(Seq(''.join(mut_nt)).translate(to_stop=True, cds=False))

            start_pos_aa = int(start_pos_nt/3)
            stop_pos_aa = int(stop_pos_nt/3)

            print("{} {}".format(start_pos_aa, stop_pos_aa))

            ref_peptide = ref_aa[max(start_pos_aa-7, 0):]
            mut_peptide = mut_aa[max(start_pos_aa-7, 0):]
            nt_context = ''.join(mut_nt[start_pos_nt-10:stop_pos_nt+10])
#            nt_context = ''.join(mut_nt[max(start_pos_nt-24, 0):])

            print("Mutant genomic context:    {}".format(nt_context))
            print("Reference peptide: {}".format(ref_peptide))
            print("Mutant peptide:    {}".format(mut_peptide))

            genomic_context_lower_pos = record['meta'].POS - 39
            genomic_context_upper_pos = record['meta'].POS + 40


#            if check_mut_pep_in_nt_context(mut_peptide, mut_nt):
            md5able_str = "{}".format(':'.join([str(record['meta'].CHROM), str(record['meta'].POS), str(record['transcript']), str(record['meta'].REF), str(record['meta'].ALT)])).encode('utf-8')
            var_md5 = hashlib.md5(md5able_str).hexdigest()[:16]
            mutant_peptides["MD5:{} VARIANT_POS:{}_{} TRANSCRIPT:{} REF:{} ALT:{} INDEL_TYPE:{} PROTEIN_CONTEXT:{} GENOMIC_CONTEXT:{} GENOMIC_CONTEXT_RANGE:{}_{}-{}".format(var_md5, record['meta'].CHROM, record['meta'].POS, record['transcript'], record['meta'].REF, record['meta'].ALT, 'frameshift_deletion', 'NA', nt_context, record['meta'].CHROM, genomic_context_lower_pos, genomic_context_upper_pos)] = mut_peptide

        # INSERTIONS #
        if re.search("ins", record['nt_change']):
            print("Frameshift insertion {}".format(record))
            ins_rec_nt, buffer, ins_nt = record['nt_change'].partition('ins')
            start_pos_nt = 0
            stop_pos_nt = 0
            start_nt = ''
            stop_nt = ''
            if re.search('_', ins_rec_nt) and not(re.search('\*', ins_rec_nt)):
                 #Need to deal with insertions that extend beyond the annotated transcript.
                 start_pos_nt = int(ins_rec_nt.split('_')[0]) - 1
                 stop_pos_nt = int(ins_rec_nt.split('_')[1]) - 1
                 if ins_nt:
                     start_nt = ins_nt[0]
                     stop_nt = ins_nt[-1]
            elif not(re.search('\*', ins_rec_nt)):
                start_pos_nt = int(ins_rec_nt) - 1
                stop_pos_nt = start_pos_nt
                if ins_nt:
                    start_nt = ins_seq[0]
                    stop_nt = start_nt
            else:
                pass
            print("{} {}".format(start_pos_nt, start_nt))
            print("{} {}".format(stop_pos_nt, stop_nt))
            if stop_pos_nt > len(bufr_nt):
                print("The mutation occurs outside of the trancsript sequence. Check transcript versions.")
                continue
            mut_nt = bufr_nt[:start_pos_nt + 1] + [ins_nt] + bufr_nt[stop_pos_nt:]
            print("Mut seq: {}".format(''.join(mut_nt)))

            ref_aa = str(Seq(''.join(bufr_nt)).translate(to_stop=True, cds=False))
            mut_aa = str(Seq(''.join(mut_nt)).translate(to_stop=True, cds=False))

            start_pos_aa = int(start_pos_nt/3)
            stop_pos_aa = int(stop_pos_nt/3)

            print("{} {}".format(start_pos_aa, stop_pos_aa))

            print("Reference affected region: {}".format(ref_aa[start_pos_aa-5:stop_pos_aa+5]))
            print("Mutant affected region:    {}".format(mut_aa[start_pos_aa-5:stop_pos_aa+5]))

            mut_peptide = mut_aa[max(start_pos_aa-7,0):]

#            nt_context = ''.join(mut_nt[max(start_pos_nt-24, 0):])
            nt_context = ''.join(mut_nt[start_pos_nt-10:stop_pos_nt+10])


            genomic_context_lower_pos = record['meta'].POS - 39
            genomic_context_upper_pos = record['meta'].POS + 40

#            if check_mut_pep_in_nt_context(mut_peptide, nt_context):
            md5able_str = "{}".format(':'.join([str(record['meta'].CHROM), str(record['meta'].POS), str(record['transcript']), str(record['meta'].REF), str(record['meta'].ALT)])).encode('utf-8')
            var_md5 = hashlib.md5(md5able_str).hexdigest()[:16]
            mutant_peptides["MD5:{} VARIANT_POS:{}_{} TRANSCRIPT:{} REF:{} ALT:{} INDEL_TYPE:{} PROTEIN_CONTEXT:{} GENOMIC_CONTEXT:{} GENOMIC_CONTEXT_RANGE:{}_{}-{}".format(var_md5, record['meta'].CHROM, record['meta'].POS, record['transcript'], record['meta'].REF, record['meta'].ALT, 'frameshift_insertion', 'NA', nt_context, record['meta'].CHROM, genomic_context_lower_pos, genomic_context_upper_pos)] = mut_peptide

        # DUPLICATIONS #
        if re.search("dup", record['nt_change']):
            print("Frameshift duplication {}".format(record))
            dup_rec_nt, buffer, dup_nt = record['nt_change'].partition('dup')
            start_pos_nt = 0
            stop_pos_nt = 0
            start_nt = 'A'
            stop_nt = 'A'
            if not(re.search('-', dup_rec_nt)):
                if re.search('_', dup_rec_nt):
                     start_pos_nt = int(dup_rec_nt.split('_')[0]) - 1
                     start_nt = dup_nt[0]
                     stop_pos_nt = int(dup_rec_nt.split('_')[1]) - 1
                     stop_nt = dup_nt[-1]
                else:
                    start_pos_nt = int(dup_rec_nt) - 1
                    start_nt = dup_nt[0]
                    stop_pos_nt = start_pos_nt
                    stop_nt = start_nt
            else:
                continue

            print("{} {}".format(start_pos_nt, start_nt))
            print("{} {}".format(stop_pos_nt, stop_nt))

            if stop_pos_nt > len(bufr_nt):
                print("The mutation occurs outside of the trancsript sequence. Check transcript versions.")
                continue

            mut_nt = ''.join(bufr_nt[:start_pos_nt] + [dup_nt] +  bufr_nt[stop_pos_nt:])

            ref_aa = str(Seq(''.join(bufr_nt)).translate(to_stop=True, cds=False))
            mut_aa = str(Seq(''.join(mut_nt)).translate(to_stop=True, cds=False))

            start_pos_aa = int(start_pos_nt/3)
            stop_pos_aa = int(stop_pos_nt/3)

            print("{} {}".format(start_pos_aa, stop_pos_aa))

            print("Reference affected region: {}".format(ref_aa[start_pos_aa-5:stop_pos_aa+5]))
            print("Mutant affected region: {}".format(mut_aa[stop_pos_aa-5:stop_pos_aa+5]))

            mut_peptide = mut_aa[max(start_pos_aa-7,0):]
#            nt_context = mut_nt[max(start_pos_nt-24,0):]
            nt_context = ''.join(mut_nt[start_pos_nt-10:stop_pos_nt+10])

            genomic_context_lower_pos = record['meta'].POS - 39
            genomic_context_upper_pos = record['meta'].POS + 40


#            if check_mut_pep_in_nt_context(mut_peptide, nt_context):
            md5able_str = "{}".format(':'.join([str(record['meta'].CHROM), str(record['meta'].POS), str(record['transcript']), str(record['meta'].REF), str(record['meta'].ALT)])).encode('utf-8')
            var_md5 = hashlib.md5(md5able_str).hexdigest()[:16]
            mutant_peptides["MD5:{} VARIANT_POS:{}_{} TRANSCRIPT:{} REF:{} ALT:{} INDEL_TYPE:{} PROTEIN_CONTEXT:{} GENOMIC_CONTEXT:{} GENOMIC_CONTEXT_RANGE:{}_{}-{}".format(var_md5, record['meta'].CHROM, record['meta'].POS, record['transcript'], record['meta'].REF, record['meta'].ALT, 'frameshift_duplication', 'NA', nt_context, record['meta'].CHROM, genomic_context_lower_pos, genomic_context_upper_pos)] = mut_peptide

    with open(args.output, 'w') as ofo:
        for k, v in mutant_peptides.items():
            ofo.write('>{}\n{}\n'.format(k, v))


def check_mut_pep_in_nt_context(mut_peptide, nt_context):
    sanity_check = True
    all_contexts = [nt_context, nt_context[1:], nt_context[2:]]
    all_orfs =[str(Seq(j).translate(cds=False)) for j in all_contexts]
    if not(any([re.search(mut_peptide, orf) for orf in all_orfs])):
        print("Mutant peptide is not within ORFs of NT context. Stopping.")
        sanity_check = False
    return sanity_check

def check_rna_coverage(args):
    """
    Takes a FASTA file with the required header information (see documentation)
    Appends new tag to header regarding RNA coverage in reads. This will
    optionally filter recrods that do not have RNA coverage.

    RNA_COVERAGE: (YES/NO)
    PERC_RNA_READS_W_VAR: (0.0-100.0)
    """
    rna_bam = pysam.AlignmentFile(args.tumor_rna_bam, "rb")

    with open(args.peptide_fasta) as fo:
        for record in SeqIO.parse(fo, 'fasta'):
            print(record)

    tx_to_cds = load_tx_cds(args)

    # This is a fairly shallow check
    # Some quick pseudocode:
    # for record in FASTA:
    #      get tx cds and ensure that txomic context is not present in transcript.
    #      get variant position
    #      get genomic range (var_pos +/- 10 bp)
    #      get txome_context(center +/- 10 bp)
    #      get reads mapping to region from pysam.fetch()
    #      use re.search to determine what proportion of reads contain the txome context ONCE
    #      Get YES/NO and proportion

    new_header_and_seqs = {}

    with open(args.peptide_fasta) as fo:
        for record in SeqIO.parse(fo, 'fasta'):
            #print(record)

            record_meta = {i.split(':')[0]:i.split(':')[1] for i in record.description.split(' ')}
            #print(record_meta)

            genome_coord_chr = record_meta['VARIANT_POS'].split('_')[0]
            genome_coord_lo = int(record_meta['VARIANT_POS'].split('_')[1]) - 10
            genome_coord_hi = int(record_meta['VARIANT_POS'].split('_')[1]) + 10

            total_context = record_meta['GENOMIC_CONTEXT']
            midpoint = int(len(total_context)/2)
            expected_seq = total_context[midpoint-5:midpoint+5]

            covered_reads = []
            for read in rna_bam.fetch(genome_coord_chr, genome_coord_lo, genome_coord_hi):
                covered_reads.append(read.query_sequence)
#            print("Number of reads: {}".format(len(covered_reads)))
            positive_hits = 0
            if len(covered_reads) > 0:
                for read in covered_reads:
#                    print(expected_seq)
#                    print(read)
#                    print(bool(re.search(expected_seq, read)))
                    if bool(re.search(expected_seq, read)):
                        print(expected_seq)
                        print(read)
                        print(re.search(expected_seq, read))
                        positive_hits += 1
                print("Positive hits: {}".format(positive_hits))
                prop_var = positive_hits/(len(covered_reads) + 0.0)
                if prop_var > 0.01:
                    print(positive_hits)
                    print(prop_var)
                    #print("Reads with mutant: {}".format(positive_hits))
                    record_meta['TOTAL_RNA_COVERGE'] = len(covered_reads)
                    record_meta['PROPORTION_RNA_READS_WITH_VARIANT'] = positive_hits/len(covered_reads)
                    new_header = str(record.description) + " {}:{} {}:{:.2f}".format('TOTAL_RNA_COVERAGE', len(covered_reads), 'PROPORTION_RNA_VARIANT_READS', prop_var)
                    new_header_and_seqs[new_header] = record.seq
            else:
                pass
#                new_header = str(record.description) + " {}:{} {}:{}".format('TOTAL_RNA_COVERAGE', '0', 'PROPORTION_RNA_VARIANT_READS', '0.00')



    with open(args.output, 'w') as ofo:
        for k, v in new_header_and_seqs.items():
            ofo.write(">{}\n{}\n".format(k, v))







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
        print("Expression Threshold: {}".format(tx_threshold))
    else:
        tx_threshold = float(args.abundance_threshold)
        print("Expression Threshold: {}".format(tx_threshold))
    expressed_txids = get_expressed_txs(args, tx_threshold)
    print("# of expressed transcripts: {}".format(len(expressed_txids)))
    filtered_records = filter_vcf_by_expression(args, expressed_txids)
    print("# of filtered records: {}".format(len(filtered_records)))
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
    threshold = np.percentile([x for x in counts if float(x) > 0], int(args.percentile))
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
                count = np.log2(float(line[count_column]) + 1)
#                print("{} {}".format(count, threshold))
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
            if transcript in expressed_txids and record not in filtered_records:
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
                line.extend(['Reference_Aff(nM)', 'Reference_Peptide', 'Agretopocity'])
                header = ','.join(line)
            elif len(line) > 14 and line[0] not in ['Pos', 'Protein']:
                # Hardcoding a 1000 nM ceiling for speed purposes.
                if float(line[15]) < 1000:
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
    with open(args.mutant_peptides) as fo:
        for line in fo.readlines():
            if line.startswith('>'):
                line = line.rstrip()
                bufr_dict = {i.split(':')[0]: i.split(':')[1] for i in line.split(' ')[1:]}
                line = line.split(' ')
                checksum = "MD5_{}".format(line[0].lstrip('>').split(':')[1][:-5])
                checksum_to_meta_map[checksum] = bufr_dict

    tx_to_tpm = {}
    tx_to_log2tpm = {}
    tx_to_uqlog2tpm = {}
    with open(args.quants) as qfo:
        tpm_col_idx = ''
        tx_col_idx = ''
        for line_idx, line in enumerate(qfo.readlines()):
            if line_idx == 0:
                tpm_col_idx = line.split('\t').index('TPM')
                tx_col_idx = line.split('\t').index('Name')
            else:
                line = line.split('\t')
                tx_to_tpm[line[tx_col_idx].split('.')[0]] = float(line[tpm_col_idx])

    for k,v in tx_to_tpm.items():
        tx_to_log2tpm[k] = np.log2(v + 1)

    uq = np.percentile(tx_to_log2tpm.values(), 75)

    for k,v in tx_to_log2tpm.items():
        tx_to_uqlog2tpm[k] = v/uq

    tx_to_gene = {}
    with open(args.gtf) as fo:
        for line in fo.readlines():
            line = line.split('\t')
            if len(line) > 3 and line[2] == 'transcript':
                gene_name = str(line).split('gene_name "')[1].split('"')[0]
                tx_id = str(line).split('transcript_id "')[1].split('"')[0]
                tx_id_no_version = str(line).split('transcript_id "')[1].split('"')[0].split('.')[0]
                tx_to_gene[tx_id] = gene_name
                tx_to_gene[tx_id_no_version] = gene_name


    var_to_ccf = {}
    with open(args.cancer_cell_fraction) as fo:
        for line_idx, line in enumerate(fo.readlines()):
            if line_idx != 0:
                line = line.rstrip().replace("b'", '').replace("'", '').split()
                var = line[0]
                ccf = line[3]
                var_to_ccf[var] = ccf

    header = []

    new_lines = []
    with open(args.binding_affinities) as fo:
        for line_idx, line in enumerate(fo.readlines()):
            line = line.rstrip().split(',')
            if line_idx == 0:
                header = line
                try:
                    header.remove('BindLevel')
                except:
                    pass
                header.insert(0, 'AntigenSource')
                header.extend(['GeneName', 'TranscriptIdentifier', 'VariantPosition', 'ReferenceAllele', 'AlternateAllele', 'TPM', 'log2(TPM+1)', 'UQ(log2(TPM+1))', 'CCF', 'ProteinContext', 'NucleotideContext', 'VariantType', 'TotalRnaReadCoverage', 'PercentRnaReadsWithVariant'])
                header = header[:4] + header[9:]
            elif len(line) < 16:
                continue
            else:
                csum = line[10]
                if csum not in checksum_to_meta_map.keys():
                    continue
                tx_id = checksum_to_meta_map[csum]['TRANSCRIPT']
                if float(line[15]) < 1000 and tx_id.split('.')[0] in tx_to_gene.keys() and tx_id.split('.')[0] in tx_to_tpm.keys() and line[10] in checksum_to_meta_map.keys():
                    gene_name = tx_to_gene[tx_id.split('.')[0]]
                    var_pos = checksum_to_meta_map[csum]['VARIANT_POS'].replace('_', ':')
                    ref = checksum_to_meta_map[csum]['REF']
                    alt = checksum_to_meta_map[csum]['ALT'].replace('[','').replace(']','')
                    tpm = tx_to_tpm[tx_id.split('.')[0]]
                    log2tpm = tx_to_log2tpm[tx_id.split('.')[0]]
                    uqlog2tpm = tx_to_uqlog2tpm[tx_id.split('.')[0]]
                    aa_context = checksum_to_meta_map[csum]['PROTEIN_CONTEXT']
                    nt_context = checksum_to_meta_map[csum]['GENOMIC_CONTEXT']
                    var_type = checksum_to_meta_map[csum]['SNV_TYPE']
                    rna_coverage = 'NA'
                    proportion_rna_var = 'NA'
                    rna_coverage = checksum_to_meta_map[csum]['TOTAL_RNA_COVERAGE']
                    proportion_rna_var = checksum_to_meta_map[csum]['PROPORTION_RNA_VARIANT_READS']
                    ccf = 'NA'
                    if var_pos in var_to_ccf.keys():
                        ccf = var_to_ccf[var_pos]
                    if 'SB' in line:
                        line.remove('<=')
                        line.remove('SB')
                    if 'WB' in line:
                       line.remove('<=')
                       line.remove('WB')
                    line.insert(0, 'SNV')
                    line.extend([gene_name, tx_id, var_pos, ref, alt, str(tpm), str(log2tpm), str(uqlog2tpm), ccf, aa_context, nt_context, var_type, rna_coverage, proportion_rna_var])
                    line = line[:4] + line[9:]
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
                line = line.rstrip()
                bufr_dict = {i.split(':')[0]: i.split(':')[1] for i in line.split(' ')[1:]}
                line = line.split(' ')
                checksum = "MD5_{}".format(line[0].lstrip('>').split(':')[1][:-5])
#                checksum = line[0].lstrip('>')[:-1]
                checksum_to_meta_map[checksum] = bufr_dict

    tx_to_tpm = {}
    tx_to_log2tpm = {}
    tx_to_uqlog2tpm = {}
    with open(args.quants) as qfo:
        tpm_col_idx = ''
        tx_col_idx = ''
        for line_idx, line in enumerate(qfo.readlines()):
            if line_idx == 0:
                tpm_col_idx = line.split('\t').index('TPM')
                tx_col_idx = line.split('\t').index('Name')
            else:
                line = line.split('\t')
                tx_to_tpm[line[tx_col_idx].split('.')[0]] = float(line[tpm_col_idx])

    for k,v in tx_to_tpm.items():
        tx_to_log2tpm[k] = np.log2(v + 1)

    uq = np.percentile(tx_to_log2tpm.values(), 75)

    for k,v in tx_to_log2tpm.items():
        tx_to_uqlog2tpm[k] = v/uq

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
            #print(line)
            line = line.rstrip()
            line = line.split()
            if not(line):
                continue
            if line[0] == 'Pos':
                header = line
                header.remove('BindLevel')
                header.extend(['GeneName', 'TranscriptIdentifier', 'VariantPosition', 'ReferenceAllele', 'AlternateAllele', 'TPM', 'log2(TPM+1)', 'UQ(log2(TPM+1))', 'CCF', 'ProteinContext', 'NucleotideContext', 'VariantType', 'TotalRnaReadCoverage', 'PercentRnaReadsWithVariant'])
                header = header[:4] + header[9:]
                header.insert(0, 'AntigenSource')
            elif len(line) < 16:
                continue
            if line[10] not in checksum_to_meta_map.keys():
                continue
            csum = line[10]
            tx_id = checksum_to_meta_map[csum]['TRANSCRIPT']
            if tx_id.split('.')[0] in tx_to_gene.keys() and tx_id.split('.')[0] in tx_to_tpm.keys() and float(line[15]) < 1000:
            #if tx_id.split(',')[0] in tx_to_gene.keys() and tx_id.split('.')[0] in tx_to_tpm.keys():
                gene_name = tx_to_gene[tx_id.split('.')[0]]
                var_pos = checksum_to_meta_map[csum]['VARIANT_POS'].replace('_', ':')
                ref = checksum_to_meta_map[csum]['REF']
                alt = checksum_to_meta_map[csum]['ALT'].replace('[','').replace(']','')
                tpm = tx_to_tpm[tx_id.split('.')[0]]
                log2tpm = tx_to_log2tpm[tx_id.split('.')[0]]
                uqlog2tpm = tx_to_uqlog2tpm[tx_id.split('.')[0]]
                ccf = 'NA'
                aa_context = checksum_to_meta_map[csum]['PROTEIN_CONTEXT']
                nt_context = checksum_to_meta_map[csum]['GENOMIC_CONTEXT']
                var_type = checksum_to_meta_map[csum]['INDEL_TYPE']
                rna_coverage = 'NA'
                rna_proportion_vars = 'NA'
                rna_coverage = checksum_to_meta_map[csum]['TOTAL_RNA_COVERAGE']
                rna_proportion_vars = checksum_to_meta_map[csum]['PROPORTION_RNA_VARIANT_READS']
                if var_pos in var_to_ccf.keys():
                    ccf = var_to_ccf[var_pos]
                if 'SB' in line:
                    line.remove('<=')
                    line.remove('SB')
                if 'WB' in line:
                    line.remove('<=')
                    line.remove('WB')
                line.extend([gene_name, tx_id, var_pos, ref, alt, str(tpm), str(log2tpm), str(uqlog2tpm), ccf, aa_context, nt_context, var_type, rna_coverage, rna_proportion_vars])
                line.insert(0, 'InDel')
                line = line[:4] + line[9:]
                new_lines.append(line)


    with open(args.output, 'w') as ofo:
        ofo.write("{}\n".format('\t'.join(header)))
        for new_line in new_lines:
            ofo.write("{}\n".format('\t'.join(new_line)))

    return new_lines



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
#    tx_threshold = 0
    if not(args.abundance_threshold):
        print("No abundance threshold!")
#        tx_threshold = get_tx_threshold(args, tx_abundances)
        tx_threshold = np.mean(tx_abundances) + (5*np.std(tx_abundances))
        print("Transcript threshold (mean + 5*std): {}".format(tx_threshold))
    else:
        print("Abundance threshold: {}".format(args.abundance_threshold))
        tx_threshold = float(args.abundance_threshold)
    expressed_hervs = get_expressed_txs(args, tx_threshold)
    print("tx_thresshold: {}".format(tx_threshold))
    print("# of expressed transcripts: {}".format(len(expressed_hervs)))
    print("Some expressed transcripts: {}".format(expressed_hervs[:10]))

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
    with open(args.expressed_hervs) as fo:
        for line in fo.readlines():
            line = line.rstrip()
            expressed_hervs.append(line)

    for seq_record in SeqIO.parse(args.herv_ref, "fasta"):
        if seq_record.description in expressed_hervs:
            expressed_hervs_seqs[seq_record.description] = seq_record.seq
            expressed_hervs_seqs["{}_rc".format(seq_record.description)] = seq_record.seq.reverse_complement()

    for i in [0, 1, 2]:
        for id, seq in expressed_hervs_seqs.items():
            #Defaulting to False here. May change to True.
            aa_seq = seq[i:].translate(to_stop=False)
            expressed_hervs_aas["{}_{}".format(id, i)] = aa_seq


    with open(args.output, 'w') as ofo:
        for k, v in expressed_hervs_aas.items():
            if len(v) > 8:
                ofo.write(">{}\n{}\n".format(k, v))


def add_herv_metadata(args):
    """
    """
    header = ''
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

    herv_expression = {}

    output_lines = []

    with open(args.binding_affinities) as fo:
        for line_idx, line in enumerate(fo.readlines()):
            line = line.rstrip()
            line = line.split()
            if len(line) >= 16 and line[0] == 'Pos':
                header = line
                try:
                    header.remove('BindLevel')
                except:
                    pass
                header.extend(['TPM'])
                header.insert(0, 'AntigenSource')
                header = header[:4] + header[9:]
            elif len(line) < 16 or line[0] in ['Protein', 'Pos']:
                continue
            elif float(line[15]) < 1000:
                if 'SB' in line:
                    line.remove('<=')
                    line.remove('SB')
                if 'WB' in line:
                    line.remove('<=')
                    line.remove('WB')
                tpm = tx_to_tpm[line[10]]
                line.extend([tpm])
                line = line[:4] + line[9:]
                line.insert(0, 'ERV')
                output_lines.append(line)

    if header and output_lines:
        with open(args.output, 'w') as ofo:
            ofo.write("{}\n".format('\t'.join(header)))
            for line in output_lines:
                ofo.write("{}\n".format('\t'.join(line)))
    else:
        with open(args.output, 'w') as ofo:
            ofo.write("No viable ERV peptides. Try loosening parameters.")

    return output_lines

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

    # Probably want to prefill relevant_germline_vars with the appropriate self
    # transcripts. Loop through missense and indels and apply to transcripts where appropriate.
    #vcf_reader = vcf.Reader(open(args.vcf), 'r', compressed=True)
    #vcf_reader = vcf.Reader(open(args.germline_vcf), 'r')

    # Extract missense and indels. Safe using missense, conservative indels,
    # and disruptive indels. I would be cautious about using frameshift
    # variants.

    # Best bet here is to use extract_missense_snvs(),
    # extract_conservative_inframe_indels(), extract_disruptive_indels().
    # Probably want to use extrace_frameshift_indels() just to allow for
    # filtering of selfs with missenses.

    print("Extracting missense variants from germline VCF...")
    missense_snvs = extract_missense_snvs(args.germline_vcf)
    print(len(missense_snvs.keys()))
    print("Extracting conservative inframe indels from germline VCF...")
    conserv_inframe_indels = extract_conservative_inframe_indels(args.germline_vcf)
    print(len(conserv_inframe_indels.keys()))
    print("Extracting disruptive inframe indels from germline VCF...")
    disrupt_inframe_indels = extract_disruptive_inframe_indels(args.germline_vcf)
    print(len(disrupt_inframe_indels.keys()))
    print("Extracting frameshift indels from germline VCF...")
    frameshift_indels = extract_frameshift_indels(args.germline_vcf)
    print(len(frameshift_indels.keys()))

    print("Number of expressed self-antigen transcripts: {}".format(len(expressed_self_seqs.keys())))

    frameshift_transcripts = [record['transcript'] for entry, record in frameshift_indels.items()]
    print(len(frameshift_transcripts))

    for self_tx in expressed_self_seqs.keys():
        if self_tx in frameshift_transcripts:
            del(expressed_self_seqs[self_tx])

    print("Number of expressed self-antigen transcripts without frameshifts: {}".format(len(expressed_self_seqs.keys())))

    all_vars = dict(missense_snvs, **conserv_inframe_indels)
    all_vars = dict(all_vars, **disrupt_inframe_indels)

    print(len(all_vars.keys()))

    for entry, record in all_vars.items():
        if record['transcript'] in expressed_self_seqs.keys() or record['transcript'].split('.')[0] in expressed_self_seqs.keys():
            if record['transcript'] not in relevant_germline_vars.keys():
                relevant_germline_vars[record['transcript'].split('.')[0]] = [record]
            else:
                relevant_germline_vars[record['transcript'].split('.')[0]].append(record)

#    print(relevant_germline_vars)


    tx_to_peps = {}
    for expressed_self, expressed_self_seq in expressed_self_seqs.items():
#        print(expressed_self)
        if expressed_self in relevant_germline_vars.keys():
            relevant_vars = relevant_germline_vars[expressed_self]
            print("{}\n{}".format(expressed_self, relevant_vars))
            if not(relevant_vars):
                print("No relevant variants, emitting...")
                tx_to_peps[expressed_self] = expressed_self_seq

            indels = [i for i in relevant_vars if 'aa3_change' in i.keys()]
            if not indels:
                print("No indels detected. Proceeding with missense variants.")
                # TODO: Current code assumes homozygosity of germline variants
                # which is a poor assuming. Check the GT from the record.
                germ_aa = list(str(expressed_self_seq))
                germline_aas = [germ_aa]
                for record in relevant_vars:
                    print(record)
                    if len(germ_aa) != int(record['aa_len']):
                        print("transcript {} shows different lengths between amino acid fasta ({}) and snpEff annotations! ({})".format(tx, record['aa_len'], len(germ_aa)))
                        continue
                    aa_pos = int(record['aa_pos']) - 1
                    germ_aa[aa_pos] = record['alt_aa']
                    print("{}".format(expressed_self_seq[aa_pos-5:aa_pos+5]))
                    print("{}".format(''.join(germ_aa)[aa_pos-5:aa_pos+5]))
                    #print(record['meta'].FORMAT['GT'])
                print("Applied all missense variants, emitting...")
                tx_to_peps[expressed_self] = expressed_self_seq
            else:
                pass
                # This will need to be revamped. The mechanism of how indels
                # are handled right now is a bit complex and needs to be
                # improved before I have faith in this.
#                seq_windows = []
#                window = 21
#                step = 8
#                capture_index = 0
#                for pos_idx in range(len(expressed_self_seq)):
#                    if pos_idx == capture_index:
#                        seq_windows.append([str(expressed_self_seq[pos_idx:(pos_idx + window)]), pos_idx+ 1, pos_idx + window + 1])
#                        capture_index += step
#
#                for seq_window in seq_windows:
#                    print("Seq window: {}".format(seq_window))
#                    seq_window_start = seq_window[1]
#                    seq_window_stop = seq_window[2]
#                    found_vars = 0
#                    for relevant_var in relevant_vars:
#                        print(relevant_var)
#                        if 'aa_pos' in record.keys():
#                            relevant_var_pos = record['aa_pos']
#                        elif re.search("del$", record['aa3_change']):
#                            del_rec_aa = record['aa3_change'].strip('del')
#                            start_pos_aa = int(del_rec_aa.split('_')[0][3:]) - 1
#                            start_aa3 = del_rec_aa.split('_')[0][:3]
#                            start_aa = seq1(start_aa3)
#                            stop_pos_aa = 0
#                            stop_aa3 = 'foo'
#                            stop_aa = 'foo'
#                            # Handling single amino acid deletions
#                            if re.search('_', del_rec_aa):
#                                stop_pos_aa = int(del_rec_aa.split('_')[1][3:]) - 1
#                                stop_aa3 = del_rec_aa.split('_')[1][:3]
#                                stop_aa = seq1(stop_aa3)
#                             else:
#                                 stop_pos_aa = start_pos_aa
#                                 stop_aa3 = start_aa3
#                                 stop_aa = start_aa
#                        elif
#                        print("Relevant var pos: {}".format(relevant_var_pos))
#                        if int(relevant_var_pos) > int(seq_window_start) and int(relevant_var_pos) < int(seq_window_stop):
#                            found_vars = 1
#                            print("Var is within window.")
#                           # print(seq_window)
#                           # print(relevant_var)
#                        relative_var_pos = int(relevant_var_pos) - int(seq_window_start)
#                        print(relative_var_pos)
#                        new_seq = list(seq_window[0])
#                        print(new_seq)
#                        if new_seq[relative_var_pos] != relevant_var[2]:
#                            print("The reference amino acid isn't what's expected.")
#                        new_seq[relative_var_pos] = relevant_var[3]
#                        if relevant_var[-1] in ['0/1', '0|1']:
#                            print("Variant is heterozygote, grabbing both genotypes.")
#                            tx_to_peps["{}_{}_het_alt0".format(expressed_self, seq_window_start)] = seq_window[0]
#                            tx_to_peps["{}_{}_het_alt1".format(expressed_self, seq_window_start)] = ''.join(new_seq)
#                        if relevant_var[-1] in ['1/1', '1|1']:
#                            tx_to_peps["{}_{}_hom_alt".format(expressed_self, seq_window_start)] = ''.join(new_seq)
#                        print("old seq: {}\nnew_seq: {}".format(seq_window[0], ''.join(new_seq)))
#                if found_vars == 0:
#                    print("{}\t{}".format(expressed_self, expressed_self_seq))
#                    tx_to_peps["{}_{}".format(expressed_self, seq_window_start)] = seq_window[0]
#        else:
#            print("No variants")
#            # If there are no applicable variants, simply emit the reference sequence.
#            tx_to_peps[expressed_self] = expressed_self_seq
#
    with open(args.output, 'w') as ofo:
        for expressed_self, expressed_self_seq in sorted(tx_to_peps.items()):
            ofo.write(">{}\n{}\n".format(expressed_self, expressed_self_seq))


def add_self_antigen_metadata(args):
    """
    """
    tx_to_tpm = {}
    tx_to_log2tpm = {}
    tx_to_uqlog2tpm = {}
    with open(args.quants) as qfo:
        tpm_col_idx = ''
        tx_col_idx = ''
        for line_idx, line in enumerate(qfo.readlines()):
            if line_idx == 0:
                tpm_col_idx = line.split('\t').index('TPM')
                tx_col_idx = line.split('\t').index('Name')
            else:
                line = line.split('\t')
                tx_to_tpm[line[tx_col_idx].split('.')[0]] = float(line[tpm_col_idx])

    for k,v in tx_to_tpm.items():
        tx_to_log2tpm[k] = np.log2(v + 1)

    uq = np.percentile(tx_to_log2tpm.values(), 75)

    for k,v in tx_to_log2tpm.items():
        tx_to_uqlog2tpm[k] = v/uq

#    tx_to_tpm = {}
#    with open(args.quants) as qfo:
#        tpm_col_idx = ''
#        tx_col_idx = ''
#        for line_idx, line in enumerate(qfo.readlines()):
#            if line_idx == 0:
#                print(line)
#                tpm_col_idx = line.split('\t').index('TPM')
#                tx_col_idx = line.split('\t').index('Name')
#            else:
#                line = line.split('\t')
#                tx_to_tpm[line[tx_col_idx].split('.')[0]] = line[tpm_col_idx]

    tx_to_gene = {}
    with open(args.gtf) as gtfo:
        for line in gtfo.readlines():
            line = line.split('\t')
            if len(line) > 3 and line[2] == 'transcript':
                gene_name = str(line).split('gene_name "')[1].split('"')[0]
                tx_id = str(line).split('transcript_id "')[1].split('"')[0].partition('.')[0]
                tx_to_gene[tx_id] = gene_name

    output_lines = []

    header = ''

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
                header.extend(['GeneName', 'TranscriptIdentifier', 'TPM', 'log2(TPM+1)', 'UQ(log2(TPM+1))'])
                header.insert(0, 'AntigenSource')
                header = header[:4] + header[9:]
            elif len(line) < 16 or line[0] in ['Protein', 'Pos']:
                continue
            elif float(line[15]) < 1000:
                if 'SB' in line:
                    line.remove('<=')
                    line.remove('SB')
                if 'WB' in line:
                    line.remove('<=')
                    line.remove('WB')
                tx_id = line[10]
                gene_name = tx_to_gene[line[10]]
                # The below line is simply for TCGA-LAML self-antigen filtering
                # and should be removed eventually.
                tpm = tx_to_tpm[line[10]]
                log2tpm = tx_to_log2tpm[tx_id.split('.')[0]]
                uqlog2tpm = tx_to_uqlog2tpm[tx_id.split('.')[0]]
                print(line[10])
                print(tpm)
                line.extend([gene_name, tx_id,str(tpm), str(log2tpm), str(uqlog2tpm)])
                line.insert(0, 'Self-Antigen')
                line = line[:4] + line[9:]
                output_lines.append(line)

    with open(args.output, 'w') as ofo:
        if header:
            ofo.write("{}\n".format('\t'.join(header)))
            for line in output_lines:
                ofo.write("{}\n".format('\t'.join(line)))
        else:
            ofo.write("No viable fusion peptides. Try loosening parameters.")

    return output_lines


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
#    viral_counts = {raw_refs[i]:raw_counts[i] for i in range(len(raw_refs)) if int(raw_counts[i]) > int(args.min_threshold)}
    all_viral_counts = {raw_refs[i]:int(raw_counts[i]) for i in range(len(raw_refs))}
    threshold = ''
    if not args.min_threshold:
        counts = [float(x) for x in all_viral_counts.values()]
        threshold = np.mean(counts) + (3*np.std(counts))
        print("Count threshold (mean + 3*std): {}".format(threshold))
    else:
        threshold= args.min_threshold
    expressed_viruses = {k:v for k, v in all_viral_counts.items() if int(v) > int(threshold)}
    with open(args.output, 'w') as ofo:
        ofo.write("virus_ref\tread_count\n")
        for k, v in expressed_viruses.items():
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

    for seq_record in SeqIO.parse(args.fasta, "fasta"):
        map_info = seq_record.description.split('|')[1].split(' ')[0]
        virus_id, unneeded, protein_id = map_info.partition('_cds_')
        protein_id = re.sub(r'_1$', '', protein_id)
        if virus_id not in virus_to_cds.keys():
            virus_to_cds[virus_id] = [protein_id]
        else:
            virus_to_cds[virus_id].append(protein_id)

        viral_seqs[protein_id] = seq_record.seq.translate()



    with open(args.fasta) as fo:
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
                header.extend(['VirusIdentifier', 'ReadCount'])
                header.insert(0, 'AntigenSource')
                header = header[:4] + header[9:]
            elif len(line) < 16 or line[0] in ['Protein', 'Pos']:
                continue
            elif float(line[15]) < 1000:
                if 'SB' in line:
                    line.remove('<=')
                    line.remove('SB')
                if 'WB' in line:
                    line.remove('<=')
                    line.remove('WB')
                #This should be extended to replace and _[0-9] to .[0-9]
                orig_index = line[10].rfind('_')
                virus_id = cds_to_virus[line[10][:orig_index] + '.' + line[10][orig_index+1:]]
                read_counts = expressed_viruses[virus_id]
                print(read_counts)
                line.extend([virus_id, read_counts])
                line.insert(0, 'Virus')
                line = line[:4] + line[9:]
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
                     if line[header_map['PROT_FUSION_TYPE']] == 'INFRAME':
                         start_prot = (int(line[header_map['CDS_LEFT_RANGE']].split('-')[1]) - 1)/3
                         peptide = line[header_map['FUSION_TRANSL']][start_prot - 8:start_prot + 8]
                         if not(re.search('\*', peptide)):
                             valid_fusions[line[header_map['#FusionName']]] = peptide
                     elif line[header_map['PROT_FUSION_TYPE']] == 'FRAMESHIFT':
                         start_prot = (int(line[header_map['CDS_LEFT_RANGE']].split('-')[1]) - 1)/3
                         init_peptide = line[header_map['FUSION_TRANSL']][start_prot - 8:]
                         print(init_peptide)
                         if re.search('\*', init_peptide):
                             stop_codon = init_peptide.index('*')
                             peptide = init_peptide[:stop_codon]
                             valid_fusions[line[header_map['#FusionName']]] = peptide
                         else:
                             valid_fusions[line[header_map['#FusionName']]] = init_peptide

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
                fusion_metadata[line[0][:15]] = {}
                for elem_idx, elem in enumerate(line[1:]):
                    fusion_metadata[line[0][:15]][col_to_idx[elem_idx + 1]] = elem

    print(fusion_metadata)

    extensions = ['JunctionReadCount', 'SpanningFragCount', 'SpliceType', 'LeftGene',
                  'LeftBreakpoint', 'RightGene', 'RightBreakpoint', 'LargeAnchorSupport',
                  'FFPM', 'PROT_FUSION_TYPE', 'annots']
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
                header_extension = extensions[:]
                header_extension[-1] = 'FusionAnnotations'
                header_extension[-2] = 'FusionType'
                header.extend(header_extension)
                header.insert(0, 'AntigenSource')
                header = header[:4] + header[9:]
            elif len(line) < 16 or line[0] in ['Protein', 'Pos']:
                continue
            elif float(line[15]) < 1000:
                if 'SB' in line:
                    line.remove('<=')
                    line.remove('SB')
                if 'WB' in line:
                    line.remove('<=')
                    line.remove('WB')
                relevant_metadata = fusion_metadata[line[10].replace('_', '.')]
                line.insert(0, 'FusionEvent')
                line_extension = []
                for extension in extensions:
                    line_extension.append(relevant_metadata[extension])
                line.extend(line_extension)
                line = line[:4] + line[9:]
                output_lines.append(line)

    if header and output_lines:
        with open(args.output, 'w') as ofo:
            ofo.write("{}\n".format('\t'.join(header)))
            for line in output_lines:
                ofo.write("{}\n".format('\t'.join(line)))
    else:
        with open(args.output, 'w') as ofo:
            ofo.write("No viable fusion peptides. Try loosening parameters.")

    return output_lines


def add_splice_metadata(args):
    """
    """
    new_header = ['AntigenSource', 'Peptide', 'Chromosome', 'GeneName', 'Strand', 'MHC', 'Aff(nM)', 'BindingProperty', 'NucleotideContext', 'MinimumExpression']
    new_lines = []

    col_idx_map = {}

    mhc_info_map = {}

    with open(args.splice_summary) as fo:
        for line_idx, line in enumerate(fo.readlines()):
            line = line.strip()
            if line_idx == 0:
                header = line.split('\t')
                for col_idx, col in enumerate(header):
                    if not(re.search('HLA', col)):
                        col_idx_map[col] = col_idx
                    else:
                        mhc_info_map[col] = col_idx

                print(col_idx_map)

            else:
                line = line.split('\t')
                base_new_line = ['SpliceVariant', line[col_idx_map['Variant_peptide_sequence']], line[col_idx_map['Chromosome']], line[col_idx_map['Gene']], line[col_idx_map['Strand']], 'ALLELE_PROXY', 'BA_PROXY', 'BP_PROXY', line[col_idx_map['DNA_sequence']], line[col_idx_map['min_expression']]]
                for allele in list(set([i.partition('_')[0] for i in mhc_info_map.keys()])):
                    allele_specific = base_new_line[:]
                    allele_specific[5] = allele
                    allele_specific[6] = line[mhc_info_map["{}_Nm".format(allele)]]
                    allele_specific[7] = line[mhc_info_map["{}_binding_property".format(allele)]]

                if float(allele_specific[6]) < 1000:
                    new_lines.append(allele_specific)
    with open(args.output, 'w') as ofo:
        ofo.write("{}\n".format('\t'.join(new_header)))
        for new_line in new_lines:
            ofo.write("{}\n".format('\t'.join(new_line)))






def add_rna_norms(args):
    """
    """
    manifest_entries = []
    with open(args.manifest) as fo:
        for line in fo.readlines():
            line = line.strip('\n')
            line = line.split('\t')
            trunc_line = [ line[-2], line[20], line[-1], line[19], line[1], line[9]]
#            print(trunc_line)
#            if trunc_line[1].startswith('ar-'):
#                new_line = [line[-2], "{}-CTRL".format(line[20].replace('ar', 'nr')), args.dataset, args.prefix, 'RNA-Seq', 'TRUE']
#                manifest_entries.append('\t'.join(new_line))
            manifest_entries.append('\t'.join(trunc_line))
        new_line = [args.pat_name, "nr-CTRL", args.dataset, args.prefix, 'RNA-Seq', 'TRUE']
        manifest_entries.append('\t'.join(new_line))


    with open(args.output, 'w') as ofo:
        for entry in manifest_entries:
            ofo.write("{}\n".format(entry))
    return manifest_entries


def make_lens_report(args):
    """
    """
    reports = glob(os.path.join(args.metadata_dir, "*metadata.txt"))
    print(reports)

    antigen_sources = ['snv', 'indel', 'viral', 'hervs', 'splice', 'fusion', 'self_antigen']
#    antigen_sources = ['snv', 'indel', 'viral', 'hervs', 'fusion', 'self_antigen']

    #Initiating dataframe with snvs...

    snv_metadata = [i for i in reports if re.search('snv.metadata', i)][0]
    print(snv_metadata)

    report_df = pd.read_csv(snv_metadata, sep='\t')

    # Skipping SNV since it was used to initiate DF.
    for antigen_source in antigen_sources[1:]:
        print(antigen_source)
        report = [i for i in reports if re.search('{}.metadata'.format(antigen_source), i)][0]
        print(report)
        tmp_report = pd.read_csv(report, sep='\t')
        if not(tmp_report.empty):
            report_df = pd.concat([report_df, tmp_report])


    sorted_df = report_df.sort_values(by="Aff(nM)")

    sorted_df.to_csv(args.output, sep='\t', index=False, na_rep='NA')

def make_antigens_barplot(args):
    """
    """
    patients = []
    counts = {}

    reports = glob(os.path.join(args.reports_dir, "*lens_report.txt"))
    print(reports)

    antigen_sources = ['SNV', 'InDel', 'Virus', 'ERV', 'SpliceVariant', 'FusionEvent', 'Self-Antigen']
    for antigen_source in antigen_sources:
        counts[antigen_source] = []

    for patient_report in reports:
        print(patient_report)
        patient_id = patient_report.split('/')[-1].replace('.lens_report.txt', '').partition('-')[2]
        print(patient_id)
        patients.append(patient_id)

        tmp_report = pd.read_csv(patient_report, sep='\t')
        count_df = tmp_report['AntigenSource'].value_counts()
        print(count_df)
        for antigen_source in counts.keys():
            if antigen_source in count_df.index.values:
                counts[antigen_source].append(np.log2(count_df[antigen_source]))
            else:
                counts[antigen_source].append(0)

    print(patients)
    print(counts)


    df = pd.DataFrame(counts, index=patients)
    print(df)

    a_series = (df != 0).any(axis=1)
    df = df.loc[a_series]

    df_sum = df.sum(axis=1)
    print(df_sum)
    df_sum.sort_values(ascending=True, inplace=True)
    df = df.reindex(df_sum.index)

    ax = df.plot.bar(stacked=True, colormap='viridis')
    ax.set_xlabel("Patients")
    ax.set_ylabel("Log-transformed Count of Predicted Peptides\n(<1000 nM binding affinity)")
    ax.set_title("Tumor Antigen Sources among TCGA-LAML patients")
    ax.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off

    fig = ax.get_figure()
    fig.tight_layout()
    fig.savefig(args.output)

def add_tcga_data(args):
    """
    """

    abbrev_to_type = load_tcga_dicts()

#    if args.tumor_type.lower() not in [x.lower() for x in abbrev_to_type.values()]:
#        if args.tumor_type.upper() in abbrev_to_type.keys():
#            args.tumor_type = abbrev_to_type

    tumor_spec_metrics = load_tcga_tx_summ(args.tumor_type, args.tcga_transcript_summary)
    tumor_spec_metrics_nov = {}

    all_tumor_spec_means = [v['mean'] for k,v in tumor_spec_metrics.items()]
    print(all_tumor_spec_means)

    initial_ks = tumor_spec_metrics.keys()
    for initial_k in initial_ks:
        trunc_k = initial_k.split('.')[0]
        tumor_spec_metrics_nov[trunc_k] = tumor_spec_metrics[initial_k]

    new_output_lines = []

    col_idx_map = {}
    with open(args.report) as fo:
       for line_idx, line in enumerate(fo.readlines()):
           line = line.rstrip().split('\t')
           if line_idx == 0:
               for col_idx, col in enumerate(line):
                   col_idx_map[col] = col_idx
               line.extend(['TCGA_TPM_Mean', 'TCGA_TPM_Median', 'TCGA_TPM_Max', 'TCGA_TPM_IQR', 'TCGA_TPM_Percentile'])
               new_output_lines.append(line)
           else:
               if line[col_idx_map['AntigenSource']] in ['SNV', 'InDel']:
                   tx_id = line[col_idx_map['TranscriptIdentifier']].split('.')[0]
                   if tx_id in tumor_spec_metrics_nov.keys():
                       tx_id_mean = float(tumor_spec_metrics_nov[tx_id]['mean'])
                       under_sum = sum([1 for x in all_tumor_spec_means if float(x) < float(tumor_spec_metrics_nov[tx_id]['mean'])])
                       percentile = (sum([1 for x in all_tumor_spec_means if float(x) < float(tumor_spec_metrics_nov[tx_id]['mean'])]) / float(len(all_tumor_spec_means)))
                       line.extend([tumor_spec_metrics_nov[tx_id]['mean'], tumor_spec_metrics_nov[tx_id]['median'], tumor_spec_metrics_nov[tx_id]['max'], tumor_spec_metrics_nov[tx_id]['iqr'], percentile])
               else:
                   line.extend(['NA', 'NA', 'NA', 'NA', 'NA'])
               new_output_lines.append(line)
    with open(args.output, 'w') as ofo:
        for line in new_output_lines:
            ofo.write("{}\n".format('\t'.join([str(x) for x in line])))






def load_tcga_tx_summ(tumor_type, summ_file):
    """
    """
    print(tumor_type)
    col_idx_map = {}
    tumor_spec_metrics = {}
    with open(summ_file) as fo:
       for line_idx, line in enumerate(fo.readlines()):
           line = line.rstrip().split('\t')
           if line_idx == 0:
               for col_idx, col in enumerate(line):
                   col_idx_map[col] = col_idx
           else:
               if line[col_idx_map['tumor_type']] == tumor_type:
                   tumor_spec_metrics[line[col_idx_map['identifier']]] = {}
                   tumor_spec_metrics[line[col_idx_map['identifier']]]['mean'] = line[col_idx_map['mean']]
                   tumor_spec_metrics[line[col_idx_map['identifier']]]['median'] = line[col_idx_map['median']]
                   tumor_spec_metrics[line[col_idx_map['identifier']]]['max'] = line[col_idx_map['max']]
                   tumor_spec_metrics[line[col_idx_map['identifier']]]['iqr'] = line[col_idx_map['iqr']]


    return tumor_spec_metrics









def main():
    args = get_args()
    #print(args)
    if args.command == 'make-snv-peptides':
        make_snv_peptides(args)
    if args.command == 'make-indel-peptides':
        make_indel_peptides(args)
    if args.command == 'filter-expressed-variants':
        expressed_variants(args)
    if args.command == 'check-rna-coverage':
        check_rna_coverage(args)
    if args.command == 'filter-isolated-variants':
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
    if args.command == 'filter-expressed-hervs':
        expressed_hervs(args)
    if args.command == 'make-herv-peptides':
        make_herv_peptides(args)
    if args.command == 'add-herv-metadata':
        add_herv_metadata(args)
    if args.command == 'filter-expressed-self-genes':
        expressed_self_genes(args)
    if args.command == 'make-self-antigen-peptides':
        make_self_antigen_peptides(args)
    if args.command == 'add-self-antigen-metadata':
        add_self_antigen_metadata(args)
    if args.command == 'filter-expressed-viruses':
        filter_virdetect_by_counts(args)
    if args.command == 'make-viral-peptides':
        make_viral_peptides(args)
    if args.command == 'add-viral-metadata':
        add_viral_metadata(args)
    if args.command == 'make-fusion-peptides':
        make_fusion_peptides(args)
    if args.command == 'add-fusion-metadata':
        add_fusion_metadata(args)
    if args.command == 'add-splice-metadata':
        add_splice_metadata(args)
    if args.command == 'consolidate-multiqc-stats':
        consolidate_multiqc_stats(args)
    if args.command == 'add-rna-normals':
        add_rna_norms(args)
    if args.command == 'make-lens-report':
        make_lens_report(args)
    if args.command == 'make-antigens-barplot':
        make_antigens_barplot(args)
    if args.command == 'add-tcga-data':
        add_tcga_data(args)

if __name__=='__main__':
    main()
