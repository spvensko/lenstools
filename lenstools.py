#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio.SeqUtils import seq1
from Bio.SeqUtils import Seq
import csv
from glob import glob
import hashlib
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
import pysam
import pandas as pd
from pprint import pprint
from pysam.libcalignmentfile import IteratorColumnRegion
import re
import scipy
from scipy import stats
import sys
import vcf


def get_args():
    """
    """
    parser = argparse.ArgumentParser(prog="lenstools",
                                     description="lenstools")
    subparsers = parser.add_subparsers(dest='command')


    # Subparser for generating peptides by coding InDels
    parser_mk_indel_peps_cntxt = subparsers.add_parser('make-indel-peptides-context',
                                                       help="Make InDel peptides FASTA.")
    parser_mk_indel_peps_cntxt.add_argument('-vts', '--var-tx-seqs',
                                            help="Directory with var-specific transcript sequences.",
                                            required=True)
    parser_mk_indel_peps_cntxt.add_argument('-st', '--somatic-txs',
                                            help="Expressed transcripts harboring somatic InDels.",
                                            required=True)
    parser_mk_indel_peps_cntxt.add_argument('-sv', '--somatic-vcf',
                                            help="Contains filtered expressed somatic variants of interest.",
                                            required=True)
    parser_mk_indel_peps_cntxt.add_argument('-g', '--gtf',
                                            help="GTF with gene annotations.",
                                            required=True)
    parser_mk_indel_peps_cntxt.add_argument('-l', '--length',
                                            help="Emitted peptide length. (default: 11)",
                                            default=11)
    parser_mk_indel_peps_cntxt.add_argument('-n', '--context-nt-length',
                                            help="Emitted contextual nucleotide sequence length (default: 39)",
                                            default=39)
    parser_mk_indel_peps_cntxt.add_argument('-o', '--output',
                                            help="Output mutant peptides FASTA.",
                                            required=True)
    parser_mk_indel_peps_cntxt.add_argument('--nt-output',
                                            help="Output mutant nucleotides FASTA.",
                                            required=True)
    
    
    # Subparser for generating peptides by missense SNVs
    parser_mk_snv_peps_cntxt = subparsers.add_parser('make-snv-peptides-context',
                                                     help="Make SNV peptides FASTA.")
    parser_mk_snv_peps_cntxt.add_argument('-vts', '--var-tx-seqs',
                                          help="Directory with var-specific transcript sequences.",
                                          required=True)
    parser_mk_snv_peps_cntxt.add_argument('-st', '--somatic-txs',
                                          help="Expressed transcripts harboring somatic SNVs.",
                                          required=True)
    parser_mk_snv_peps_cntxt.add_argument('-sv', '--somatic-vcf',
                                          help="Contains filtered expressed somatic variants of interest.",
                                          required=True)
    parser_mk_snv_peps_cntxt.add_argument('-g', '--gtf',
                                          help="GTF with gene annotations.",
                                          required=True)
    parser_mk_snv_peps_cntxt.add_argument('-l', '--length',
                                          help="Emitted peptide length. (default: 11)",
                                          default=11)
    parser_mk_snv_peps_cntxt.add_argument('-n', '--context-nt-length',
                                          help="Emitted contextual nucleotide sequence length (default: 39)",
                                          default=39)
    parser_mk_snv_peps_cntxt.add_argument('-o', '--mt-output',
                                          help="Output mutant peptides FASTA.",
                                          required=True)
    parser_mk_snv_peps_cntxt.add_argument('--nt-output',
                                          help="Output mutant nucleotides FASTA.",
                                          required=True)
    parser_mk_snv_peps_cntxt.add_argument('-w', '--wt-output',
                                          help="Output wildtype peptides FASTA.",
                                          required=True)
    

    # Subparser for making ERV peptides
    parser_mk_erv_peps = subparsers.add_parser('make-erv-peptides',
                                               help="make ERV peptides FASTA.")
    parser_mk_erv_peps.add_argument('-e', '--expressed-ervs',
                                    help="File containg list of Expressed ERVs",
                                    required=True)
    parser_mk_erv_peps.add_argument('-r', '--patient-ervs-fasta',
                                    help="Patient ERV FASTA (with homozygous germline vars).",
                                    required=True)
    parser_mk_erv_peps.add_argument('-s', '--species',
                                    help="Species (e.g. mm). Default: hs",
                                    default='hs')
    parser_mk_erv_peps.add_argument('-g', '--geve-reference',
                                    help="gEVE reference file",
                                    required=True)
    parser_mk_erv_peps.add_argument('-o', '--output',
                                    help="Output ERV peptides FASTA.",
                                    required=True)
    parser_mk_erv_peps.add_argument('-n', '--nt-output',
                                    help="Output ERV neoleotide FASTA.",
                                    required=True)


    # Subparser for making viral peptides
    parser_mk_viral_peps = subparsers.add_parser('make-viral-peptides',
                                                  help="Make viral peptides FASTA.")
    parser_mk_viral_peps.add_argument('-f', '--patient-viral-fasta',
                                      help="Patient viral fasta (with homozygous germline vars).",
                                      required=True)
    parser_mk_viral_peps.add_argument('-o', '--output',
                                      help="Output viral peptides FASTA.",
                                      required=True)


    # Subparser for making self-antigen peptides
    parser_mk_self_peps = subparsers.add_parser('make-self-antigen-peptides',
                                                help="Make CTA/self-antigen peptides FASTA.")
    parser_mk_self_peps.add_argument('-s', '--selfs-seqs-fasta',
                                     help="CTA/Self-antigens FASTA (with homozygous germline vars).",
                                     required=True)
    parser_mk_self_peps.add_argument('-g', '--gtf',
                                     help="GTF with gene annotations.",
                                     required=True)
    parser_mk_self_peps.add_argument('-e', '--expressed-selfs',
                                     help="File containing list of expressed CTA/self-antigen genes.",
                                     required=True)
    parser_mk_self_peps.add_argument('-o', '--output',
                                     help="Output CTA/self-antigens peptide FASTA.",
                                     required=True)
    parser_mk_self_peps.add_argument('-n', '--nt-output',
                                     help="Output CTA/self-antigens nucleotide FASTA.",
                                     required=True)


    # Subparser for making fusion peptides                                                          
    parser_mk_fusion_peps_cntxt = subparsers.add_parser('make-fusion-peptides-context',                     
                                                        help="Make fusion peptides FASTA.")         
    parser_mk_fusion_peps_cntxt.add_argument('-f', '--fusions',                                     
                                             help="Predicted fusions (STARFusion abridged coding effect format).",         
                                             required=True)                                         
    parser_mk_fusion_peps_cntxt.add_argument('-g', '--gtf',                                         
                                             help="GTF with gene annotations.",         
                                             required=True)                                         
    parser_mk_fusion_peps_cntxt.add_argument('-e', '--exons-fasta',                                    
                                             help="FASTA of fusion transcripts exons (with homozygous germline vars)", 
                                             required=True)                                         
    parser_mk_fusion_peps_cntxt.add_argument('-t', '--fusion-txs',                                    
                                             help="Expressed fusion transcripts (extracted from STARFusion output).",   
                                             required=True)                                         
    parser_mk_fusion_peps_cntxt.add_argument('-o', '--output',                                      
                                             help='Output fusion peptides FASTA.',                  
                                             required=True)                                         
    parser_mk_fusion_peps_cntxt.add_argument('-n', '--nt-output',                                   
                                             help='Output fusion nucleotide FASTA.',                        
                                             required=True)  


    # Subparser for adding SNV metadata
    parser_add_snv_metadata = subparsers.add_parser('add-snv-metadata',
                                                    help="Add SNV peptides metadata.")
    parser_add_snv_metadata.add_argument('-m', '--mutant-peptides',
                                         help="Mutant SNV peptides FASTA.",
                                         required=True)
    parser_add_snv_metadata.add_argument('-q', '--quants',
                                         help="Tumor ranscript abundance file (Salmon format).",
                                         required=True)
    parser_add_snv_metadata.add_argument('-c', '--cancer-cell-fraction',
                                         help="Cancer cell fraction file (PyClone-VI format).",
                                         required=True)
    parser_add_snv_metadata.add_argument('-g', '--gtf',
                                         help="GTF with gene annotations",
                                         required=True)
    parser_add_snv_metadata.add_argument('-b', '--binding-affinities',
                                         help="Binding affinities file (NetMHCpan-4.1b format).",
                                         required=True)
    parser_add_snv_metadata.add_argument('-o', '--output',
                                         help="Output SNV report.",
                                         required=True)


    # Subparser for adding InDel metadata
    parser_add_indel_metadata = subparsers.add_parser('add-indel-metadata',
                                                      help="Add InDel peptides metadata.")
    parser_add_indel_metadata.add_argument('-m', '--mutant-peptides',
                                           help="Mutant InDel peptides FASTA.",
                                           required=True)
    parser_add_indel_metadata.add_argument('-q', '--quants',
                                           help="Tumor ranscript abundance file (Salmon format).",
                                           required=True)
    parser_add_indel_metadata.add_argument('-c', '--cancer-cell-fraction',
                                           help="Cancer cell fraction file (PyClone-VI format).",
                                           required=True)
    parser_add_indel_metadata.add_argument('-g', '--gtf',
                                           help="GTF file with gene annotations.",
                                           required=True)
    parser_add_indel_metadata.add_argument('-b', '--binding-affinities',
                                           help="Binding affinities file (netMHCpan format).",
                                           required=True)
    parser_add_indel_metadata.add_argument('-o', '--output',
                                           help="Output InDel report.",
                                           required=True)


    # Subparser for adding ERV metadata
    parser_add_erv_metadata = subparsers.add_parser('add-erv-metadata',
                                                     help="Add ERV peptides metadata.")
    parser_add_erv_metadata.add_argument('-p', '--peptides',
                                         help="Annotated peptide file.",
                                         required=True)
    parser_add_erv_metadata.add_argument('-b', '--binding-affinities',
                                         help="Binding affinities file (NetMHCpan-4.1b format).",
                                         required=True)
    parser_add_erv_metadata.add_argument('-q', '--quants',
                                         help="Tumor transcript abundance file (Salmon format).",
                                         required=True)
    #This VCF comes from the simple consensus workflow.
    parser_add_erv_metadata.add_argument('-v', '--patient-vcf',
                                         help="Patient ERV VCF",
                                         required=True)
    parser_add_erv_metadata.add_argument('-n', '--nt',
                                         help="Patient-specific ERV FASTA",
                                         required=True)
    parser_add_erv_metadata.add_argument('-d', '--geve-data',
                                         help="External gEVE data.",
                                         required=True)
    parser_add_erv_metadata.add_argument('-t', '--trim-chr-prefix',
                                         help="Trim the chr prefix from gEVE references.",
                                         action='store_true')
    parser_add_erv_metadata.add_argument('-o', '--output',
                                         help="Output ERV report",
                                         required=True)


    # Subparser for getting viral metadata
    parser_add_viral_metadata = subparsers.add_parser('add-viral-metadata',
                                                      help="Add viral peptides metadata.")
    parser_add_viral_metadata.add_argument('-b', '--binding-affinities',
                                           help="Binding affinities data (netMHCpan-4.1b format).",
                                           required=True)
    parser_add_viral_metadata.add_argument('-r', '--viral-cds-ref',
                                           help="Viral CDS reference FASTA (used for VirDetect).",
                                           required=True)
    parser_add_viral_metadata.add_argument('-p', '--viral-pep-ref',
                                           help="Viral peptide reference FASTA (used for VirDetect).",
                                           required=True)
    parser_add_viral_metadata.add_argument('-o', '--output',
                                           help="Output viral report.",
                                           required=True)


    # Subparser for getting self-antigen metadata
    parser_add_self_metadata = subparsers.add_parser('add-self-antigen-metadata',
                                                     help="Add self-antigen peptides metadata.")
    parser_add_self_metadata.add_argument('-q', '--quants',
                                          help="Tumor transcript abundance file (Salmon format).",
                                          required=True)
    parser_add_self_metadata.add_argument('-b', '--binding-affinities',
                                          help="Binding affinitities data (netMHCpan-4.1b format)",
                                          required=True)
    parser_add_self_metadata.add_argument('-f', '--fasta',
                                          help="CTA/Self-antigen Peptide FASTA",
                                          required=True)
    parser_add_self_metadata.add_argument('-g', '--gtf',
                                          help="GTF with gene annotation.s",
                                          required=True)
    parser_add_self_metadata.add_argument('-l', '--gene-list',
                                          help="File with CTA/self-antigen gene list.",
                                          required=True)
    parser_add_self_metadata.add_argument('-o', '--output',
                                          help="Output CTA/self-antigen report.",
                                          required=True)


    # Subparser for adding fusion metadata
    parser_add_fusion_metadata = subparsers.add_parser('add-fusion-metadata',
                                                       help="Add metadata to fusion binding affinity data.")
    parser_add_fusion_metadata.add_argument('-b', '--binding-affinities',
                                            help="Binding affinities data (netMHCpan-4.1b format).",
                                            required=True)
    parser_add_fusion_metadata.add_argument('-f', '--fusions',
                                            help="Predicted fusions (STARFusion format).",
                                            required=True)
    parser_add_fusion_metadata.add_argument('-a', '--fasta',
                                            help="FASTA with fusion peptides.",
                                            required=True)
    parser_add_fusion_metadata.add_argument('-o', '--output',
                                            help="Output file.",
                                            required=True)


    # Subparser for filtering variants for expression
    parser_expressed_variants = subparsers.add_parser('filter-expressed-variants',
                                                 help="Filter variants for expression.")
    parser_expressed_variants.add_argument('-q', '--quants',
                                           help="Transcript abundance file (Salmon format)).",
                                           required=True)
    parser_expressed_variants.add_argument('-m', '--metric',
                                           help="Expression metric (default: TPM).",
                                           default='TPM')
    parser_expressed_variants.add_argument('-r', '--exclude-zeros',
                                           help="Exclude zeros when calcualating expression percentile. (default: False)",
                                           action='store_true')
    parser_expressed_variants.add_argument('-p', '--percentile',
                                           help="Expression percentile for filtering (default: 90).",
                                           default=90)
    parser_expressed_variants.add_argument('-t', '--abundance-threshold',
                                           help="Expression threshold for filtering (default: 0).",
                                           default=0)
    parser_expressed_variants.add_argument('-v', '--vcf',
                                           help="Annotated (snpEff) VCF.",
                                           required=True)
    parser_expressed_variants.add_argument('-o', '--output',
                                           help="Output file.",
                                           required=True)
    parser_expressed_variants.add_argument('-s', '--somatic-txs',
                                           help="Somatic transcripst output file.",
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
    parser_isolated_variants.add_argument('-a', '--allow-silent-and-homozygous',
                                          help="Allow silent and homozygous germline variants near variants. (default: False)",
                                          action="store_true")
    parser_isolated_variants.add_argument('-o', '--output',
                                          help="Output file.",
                                          required=True)


    # Subparser for filtering expressed ERVs through TPM.
    parser_expressed_ervs = subparsers.add_parser('filter-expressed-ervs',
                                                   help="Filter expressed ERVs.")
    parser_expressed_ervs.add_argument('-q', '--quants',
                                       help="Transcript abundance file (Salmon format).",
                                       required=True)
    parser_expressed_ervs.add_argument('-m', '--metric',
                                       help="Column for expression from abundance file. (default: TPM)",
                                       default='TPM')
    parser_expressed_ervs.add_argument('-r', '--exclude-zeros',
                                       help="Exclude zeros for expression percentile.",
                                       action='store_false')
    parser_expressed_ervs.add_argument('-z', '--abundance-threshold',
                                       help="Expression threshold for filtering. (default: 10)",
                                       default=10)
    parser_expressed_ervs.add_argument('-t', '--trim-chr-prefix',
                                       help="Trim the chr prefix from gEVE references.",
                                       action='store_true')
    parser_expressed_ervs.add_argument('-o', '--output',
                                       help="Output file.",
                                       required=True)
   
 
    # Subparser for filtering expressed ERVs through RNA coverage.
    parser_erv_rna_cov = subparsers.add_parser('check-erv-rna-coverage',
                                                help="Filter for ERVs with sufficient average coverage.")
    parser_erv_rna_cov.add_argument('-c', '--coverage-file',
                                    help="Coverage file created by samtools cov.",
                                    required=True)
    parser_erv_rna_cov.add_argument('-e', '--expressed-hervs',
                                    help="File containing list of expressed ERVs.",
                                    default="TPM")
    parser_erv_rna_cov.add_argument('-m', '--mean-depth',
                                    help="Required mean depth for ERV to be classified as expressed (default: 0).",
                                    default=0)
    parser_erv_rna_cov.add_argument('-p', '--coverage',
                                    help="Required cov for ERV to be classified as expressed (default: 0).",
                                    default=0)
    parser_erv_rna_cov.add_argument('-o', '--output',
                                    help="Output file.",
                                    required=True)
   
 
    # Subparser for filtering expressed viruses through RNA coverage.
    parser_virus_rna_cov = subparsers.add_parser('check-virus-rna-cov',
                                                      help="Filter for viruses with sufficient average cov.")
    parser_virus_rna_cov.add_argument('-c', '--cov-file',
                                      help="Coverage file created by samtools cov.",
                                      required=True)
    parser_virus_rna_cov.add_argument('-e', '--expressed-viruses',
                                      help="File containing list of expressed viruses.",
                                      required=True)
    parser_virus_rna_cov.add_argument('-m', '--mean-depth',
                                      help="Required mean depth for virus to be classified as expressed (default: 0).",
                                      default=0)
    parser_virus_rna_cov.add_argument('-p', '--cov',
                                      help="Required cov for virus to be classified as expressed (default: 0).",
                                      default=0)
    parser_virus_rna_cov.add_argument('-o', '--output',
                                      help="Output file.",
                                      required=True)
   
 
    # Subparser for creating expressed ERV bed file.
    parser_expressed_ervs_bed = subparsers.add_parser('get-expressed-ervs-bed',
                                                      help="Create BED file for expressed ERVs.")
    parser_expressed_ervs_bed.add_argument('-e', '--expressed-ervs',
                                           help="File containing list of expressed ERVs.",
                                           required=True)
    parser_expressed_ervs_bed.add_argument('-r', '--geve-reference',
                                           help="gEVE reference file",
                                           required=True)
    parser_expressed_ervs_bed.add_argument('-o', '--output',
                                           help="Output file.",
                                           required=True)
   
 
    # Subparser for creating expressed self-antigen/CTA BED file
    parser_expressed_selfs_bed = subparsers.add_parser('get-expressed-selfs-bed',
                                                       help="Create BED file for expressed self-antigen/CTAs..")
    parser_expressed_selfs_bed.add_argument('-e', '--expressed-selfs',
                                            help="File with expressed CTA/self-antigen gene list",
                                            required=True)
    parser_expressed_selfs_bed.add_argument('-g', '--gff',
                                            help="GFF File.",
                                            required=True)
    parser_expressed_selfs_bed.add_argument('-o', '--output',
                                            help="Output file.",
                                            required=True)
   
 
    # Subparser for creating expressed self-antigen/CTA BED file
    parser_expressed_txs_bed = subparsers.add_parser('get-expressed-transcripts-bed',
                                                       help="Create BED file for expressed transripts.")
    parser_expressed_txs_bed.add_argument('-t', '--expressed-transcripts',
                                          help="File containing list of expressed transcripts",
                                          required=True)
    parser_expressed_txs_bed.add_argument('-g', '--gff',
                                          help="GFF with gene annotations.",
                                          required=True)
    parser_expressed_txs_bed.add_argument('-o', '--output',
                                          help="Output file.",
                                          required=True)
   
 
    # Subparser for filtering expressed viruses
    parser_expressed_viral_bed = subparsers.add_parser('get-expressed-viral-bed',
                                                       help="Filter expressed hERVs.")
    parser_expressed_viral_bed.add_argument('-e', '--expressed-viruses',
                                            help="Expressed viruses file.",
                                            required=True)
    parser_expressed_viral_bed.add_argument('-r', '--viral-cds-ref',
                                            help="Viral coding sequence (CDS) FASTA.",
                                            required=True)
    parser_expressed_viral_bed.add_argument('-o', '--output',
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
   
 
    # Subparser for filtering virdetect outputs for expressed viruses
    parser_filter_viral_cds = subparsers.add_parser('filter-expressed-viral-cds',
                                                   help="Generate list of expressed viral coding sequences")
    parser_filter_viral_cds.add_argument('-q', '--viral-cds-quants',
                                         help="Viral counts file from VirDetect.",
                                         required=True)
    parser_filter_viral_cds.add_argument('-r', '--viral-cds-ref',
                                         help="Viral reference FASTA (used for VirDetect).",
                                         required=True)
    parser_filter_viral_cds.add_argument('-e', '--expressed-viruses',
                                         help="List of expressed viruses.",
                                         required=True)
    parser_filter_viral_cds.add_argument('-o', '--output',
                                         help="Output file.",
                                         required=True)


    # Subparser for filtering self-antigens for expression
    parser_expressed_self_genes = subparsers.add_parser('filter-expressed-self-genes',
                                                        help="Filter expressed CTA/self-antigen genes.")
    parser_expressed_self_genes.add_argument('-q', '--quants',
                                             help="Transcript abundance file (Salmon format).",
                                             required=True)
    parser_expressed_self_genes.add_argument('-m', '--metric',
                                             help="Expression metric. (default: TPM)",
                                             default='TPM')
    parser_expressed_self_genes.add_argument('-r', '--exclude-zeros',
                                             help="Exclude zeros for expression percentile. (default: False)",
                                             action='store_true')
    parser_expressed_self_genes.add_argument('-p', '--percentile',
                                             help="Expression percentile for filtering expression (default: 90).",
                                             default=90)
    parser_expressed_self_genes.add_argument('-t', '--abundance-threshold',
                                             help="Expression abundance threshold for filtering expression (default: 0)",
                                             default=0)
    parser_expressed_self_genes.add_argument('-g', '--gene-list',
                                             help="File containing CTA/Self-antigen genes.",
                                             required=True)
    parser_expressed_self_genes.add_argument('-f', '--gtf',
                                             help="GTF file.",
                                             required=True)
    parser_expressed_self_genes.add_argument('-o', '--output',
                                             help="Output file.",
                                             required=True)


    # Subparser for calculalting agretopicity (mut BA/wt BA)
    parser_calc_agreto = subparsers.add_parser('calculate-agretopicity',
                                               help="Calcuate agreotopicity (mut BA/wt BA).")
    parser_calc_agreto.add_argument('-w', '--wt-fasta',
                                    help="Wildtype peptide FASTA (LENSTools format).",
                                    required=True)
    parser_calc_agreto.add_argument('-m', '--mt-fasta',
                                    help="Mutant peptide FASTA. (LENSTools format)",
                                    required=True)
    parser_calc_agreto.add_argument('-o', '--output',
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
                                      help="Manifest file with original patient data.",
                                      required=True)
    parser_add_rna_norms.add_argument('-p', '--prefix',
                                      help="RNA normal sample Run Name (e.g. nr-CTRL).",
                                      required=True)
    parser_add_rna_norms.add_argument('-d', '--dataset',
                                      help="RNA normal sample dataset.",
                                      required=True)
    parser_add_rna_norms.add_argument('-n', '--pat_name',
                                      help="RNA normal sample Patient Name (e.g. CTRL).",
                                      required=True)
    parser_add_rna_norms.add_argument('-o', '--output',
                                      help="Output manifest file.",
                                      required=True)


    # Subparser for adding splice metadata
    parser_add_splice_meta = subparsers.add_parser('add-splice-metadata',
                                                 help="Add splice metadata.")
    parser_add_splice_meta.add_argument('-s', '--splice-summary',
                                        help="Splice summary file from NeoSplice.",
                                        required=True)
    parser_add_splice_meta.add_argument('-o', '--output',
                                        help="Output file.",
                                        required=True)


    # Subparser for making LENS report
    parser_mk_lens_report = subparsers.add_parser('make-lens-report')
    parser_mk_lens_report.add_argument('--metadata-dir', '-d',
                                       help="Path to metadata directory.",
                                       required=True)
    parser_mk_lens_report.add_argument('--output', '-o',
                                       help="Output file.",
                                       required=True)


    # Subparser for making pan-patient antigen source barplot
    parser_mk_antigen_plot = subparsers.add_parser('make-antigens-barplot')
    parser_mk_antigen_plot.add_argument('--reports-dir', '-d',
                                        help="Path containing patient reports.",
                                        required=True)
    parser_mk_antigen_plot.add_argument('--output', '-o',
                                        help="Output file.",
                                        required=True)


    # Subparser for adding TCGA data to SNV or InDel LENS report
    parser_add_tcga_data = subparsers.add_parser('add-tcga-data')
    parser_add_tcga_data.add_argument('--report', '-r',
                                      help="SNV or InDel LENS Report",
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


    # Subparser for determining read count support for Viral, ERV, and self peptides.
    parser_get_pep_rd_cnt = subparsers.add_parser('get-peptide-read-count')
    parser_get_pep_rd_cnt.add_argument('--netmhcpan', '-n',
                                       help="NetMHCpan file with peptides of interest.",
                                       required=True)
    parser_get_pep_rd_cnt.add_argument('--bam', '-b',
                                       help="RNA-Sequencing BAM.",
                                       required=True)
    parser_get_pep_rd_cnt.add_argument('--cds-fasta', '-c',
                                       help="FASTA with peptide nucleotide coding sequence.",
                                       required=True)
    parser_get_pep_rd_cnt.add_argument('--output', '-o',
                                       help="Output file.",
                                       required=True)
   
 
    # Subparser for determining read count support for SNV-derived peptides.
    parser_get_snv_pep_rd_cnt = subparsers.add_parser('get-snv-peptide-read-count')
    parser_get_snv_pep_rd_cnt.add_argument('--netmhcpan', '-n',
                                           help="NetMHCpan file with peptides of interest.",
                                           required=True)
    parser_get_snv_pep_rd_cnt.add_argument('--bam', '-b',
                                           help="RNA-Sequencing BAM.",
                                           required=True)
    parser_get_snv_pep_rd_cnt.add_argument('--nt-fasta', '-c',
                                           help="FASTA file peptide nucleotide coding sequence.",
                                           required=True)
    parser_get_snv_pep_rd_cnt.add_argument('--gtf', '-g',
                                           help="GTF with gene annotations.",
                                           required=True)
    parser_get_snv_pep_rd_cnt.add_argument('--output', '-o',
                                           help="Output file.",
                                           required=True)
   
 
    # Subparser for determining read count support for InDel peptides.
    parser_get_indel_pep_rd_cnt = subparsers.add_parser('get-indel-peptide-read-count')
    parser_get_indel_pep_rd_cnt.add_argument('--netmhcpan', '-n',
                                             help="NetMHCpan file with peptides of interest.",
                                             required=True)
    parser_get_indel_pep_rd_cnt.add_argument('--bam', '-b',
                                             help="RNA-Sequencing BAM",
                                             required=True)
    parser_get_indel_pep_rd_cnt.add_argument('--nt-fasta', '-c',
                                             help="FASTA file with peptide nucleotide coding sequence",
                                             required=True)
    parser_get_indel_pep_rd_cnt.add_argument('--gtf', '-g',
                                             help="GTF with gene annotations",
                                             required=True)
    parser_get_indel_pep_rd_cnt.add_argument('--output', '-o',
                                             help="Output file.",
                                             required=True)

    
    # Subparser for determining read count support for fusion peptides.
    parser_get_fus_pep_rd_cnt = subparsers.add_parser('get-fusion-peptide-read-count')
    parser_get_fus_pep_rd_cnt.add_argument('--netmhcpan', '-n',
                                           help="NetMHCpan file with peptides of interest.",
                                           required=True)
    parser_get_fus_pep_rd_cnt.add_argument('--fusions', '-f',
                                           help="Fusions (STARFusion format)",
                                           required=True)
    parser_get_fus_pep_rd_cnt.add_argument('--nt-fasta', '-t',
                                           help="FASTA file with peptide nucleotide coding sequences",
                                           required=True)
    parser_get_fus_pep_rd_cnt.add_argument('--fusion-reads', '-r',
                                           help="Junction or discordant reads supporting fusions",
                                           required=True)
    parser_get_fus_pep_rd_cnt.add_argument('--output', '-o',
                                           help="Output file.",
                                           required=True)
   
 
    # Subparser for determining read count support for splice peptides.
    parser_get_splc_pep_rd_cnt = subparsers.add_parser('get-splice-peptide-read-count')
    parser_get_splc_pep_rd_cnt.add_argument('--neosplice-summary', '-n',
                                            help="NeoSplice summary file",
                                            required=True)
    parser_get_splc_pep_rd_cnt.add_argument('--bam', '-b',
                                            help="RNA-Sequencing BAM file",
                                            required=True)
    parser_get_splc_pep_rd_cnt.add_argument('--output', '-o',
                                            help="Output file.",
                                            required=True)
   
 
    # Subparser for filtering peptides that exist in wildtype sample
    parser_fltr_mut_peps = subparsers.add_parser('filter-mutant-peptides')
    parser_fltr_mut_peps.add_argument('--wt-netmhcpan', '-w',
                                      help="Wildtype NetMHCpan",
                                      required=True)
    parser_fltr_mut_peps.add_argument('--mt-netmhcpan', '-m',
                                      help="Mutant NetMHCpan",
                                      required=True)
    parser_fltr_mut_peps.add_argument('--wt-output', '-wo',
                                      help="Wildtype filtered output",
                                      required=False)
    parser_fltr_mut_peps.add_argument('--mt-output', '-mo',
                                      help="Mutant filtered output",
                                      required=False)
    parser_fltr_mut_peps.add_argument('--max-peptide-length', '-l',
                                      help="Maximum peptide length (default: 11)",
                                      default=11)


    return parser.parse_args()


def load_tcga_dicts():
    """
    Convert TCGA tumor type abbreviations to complete name.
    """
    abbrev_to_tumor_type = {
    "LAML": "Acute Myeloid Leukemia",
    "ACC":  "Adrenocortical carcinoma",
    "BLCA": "Bladder Urothelial Carcinoma",
    "LGG":  "Brain Lower Grade Glioma",
    "BRCA": "Breast invasive carcinoma",
    "CESC": "Cervical squamous cell carcinoma and endocervical adenocarcinoma",
    "CHOL": "Cholangiocarcinoma",
    "LCML": "Chronic Myelogenous Leukemia",
    "COAD": "Colon adenocarcinoma",
    "CNTL": "Controls",
    "ESCA": "Esophageal carcinoma",
    "FPPP": "FFPE Pilot Phase II",
    "GBM":  "Glioblastoma multiforme",
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
    "OV":   "Ovarian serous cystadenocarcinoma",
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
    "UCS":  "Uterine Carcinosarcoma",
    "UCEC": "Uterine Corpus Endometrial Carcinoma",
    "UVM":  "Uveal Melanoma"}

    return abbrev_to_tumor_type


def load_vcf(input_vcf):
    """
    Uses VCF module to load VCF file.

    Args:
        input_vcf: Input VCF file.

    Returns:
        vcf_reader (vcf Reader object): loaded VCF.
    """
    vcf_reader = ''
    if input_vcf.endswith('gz'):
        vcf_reader = vcf.Reader(open(input_vcf), 'r', compressed=True)
    else:
        vcf_reader = vcf.Reader(open(input_vcf), 'r')
    return vcf_reader


def extract_missense_snvs(input_vcf):
    """
    Extracts missense SNVs from annotated somatic VCF file.
    Assumes annotated VCF follows ANN standard:
    http://pcingola.github.io/SnpEff/adds/VCFannotationformat_v1.0.pdf

    Args:
        input_vcf: Annotated VCF.

    Returns:
        missense_records (dict): keys are <chr>:<pos> and values are
                                 dictionaries with variant metadata.
    """
    missense_records = {}

    vcf_reader = load_vcf(input_vcf)

    for record in vcf_reader:
        annotations = [x for x in record.INFO['ANN']]
        missense_records["{}:{}".format(record.CHROM, record.POS)] = []
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

                cdna_pos = effects[11].split('/')[0]

                aa3_change = effects[10].lstrip('p.')
                aa3_change = re.split('(\d+)', aa3_change)

                ref_aa3 = aa3_change[0]
                alt_aa3 = aa3_change[2]

                aa_pos = aa3_change[1]

                ref_aa = seq1(ref_aa3)
                alt_aa = seq1(alt_aa3)
                aa_len = effects[13].split('/')[1]

                codon_pos = int((int(effects[13].split('/')[0])*3) - 2)

                var_meta = {'transcript': transcript,
                            'cdna_pos': cdna_pos,
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

                missense_records["{}:{}".format(record.CHROM, record.POS)].append(var_meta)
    
    return missense_records


def extract_conservative_inframe_indels(input_vcf):
    """
    Extracts conservative inframe InDels from annotated somatic VCF file.

    Args:
        input_vcf: Annotated VCF.

    Returns:
        conserv_inframe_indels(dict): keys are <chr>:<pos> and values are
                                      dictionaries with variant metadata.
    """
    conserv_inframe_indels = {}

    vcf_reader = load_vcf(input_vcf)

    for record in vcf_reader:
        annotations = [x for x in record.INFO['ANN']]
        conserv_inframe_indels["{}:{}".format(record.CHROM, record.POS)] = []
        for annotation in annotations:
            effects = annotation.split('|')
            if re.search('^conservative_inframe_[a-z]+$', effects[1]) and effects[13]:
                transcript = effects[6]
                aa3_change = effects[10].lstrip('p.')
                aa_len = effects[13].split('/')[1]
                nt_change = effects[9].lstrip('c.')
                nt_len = effects[12].split('/')[1]
                cdna_pos = effects[11].split('/')[0]
                var_meta = {'transcript': transcript,
                            'aa3_change': aa3_change,                    
                            'nt_change': nt_change,                      
                            'aa_len': aa_len,                            
                            'nt_len': nt_len,                            
                            'cdna_pos': cdna_pos,                        
                            'meta': record}                      
                conserv_inframe_indels["{}:{}".format(record.CHROM, record.POS)].append(var_meta)

    # Removing InDels without annotation metadata.
    ks_to_del = []
    for k,v in conserv_inframe_indels.items():
        if not v:
            ks_to_del.append(k)
    for k in ks_to_del:
        del(conserv_inframe_indels[k])
    return conserv_inframe_indels


def extract_disruptive_inframe_indels(input_vcf):
    """
    Extracts disruptive inframe InDels from annotated somatic VCF file.

    Args:
        input_vcf: Annotated VCF.

    Returns:
        disrupt_inframe_indels(dict): keys are <chr>:<pos> and values are
                                      dictionaries with variant metadata.
    Focusing on deletions now, need to find a good example of a conservative_inframe_insertion
    """
    disrupt_inframe_indels = {}

    vcf_reader = load_vcf(input_vcf)

    for record in vcf_reader:
        annotations = [x for x in record.INFO['ANN']]
        disrupt_inframe_indels["{}:{}".format(record.CHROM, record.POS)] = []
        for annotation in annotations:
            effects = annotation.split('|')
            if re.search('^disruptive_inframe_[a-z]+$', effects[1]) and effects[13]:
                transcript = effects[6]
                aa3_change = effects[10].lstrip('p.')
                aa_len = effects[13].split('/')[1]
                nt_change = effects[9].lstrip('c.')
                nt_len = effects[12].split('/')[1]
                cdna_pos = effects[11].split('/')[0]
                var_meta = {'transcript': transcript,
                            'aa3_change': aa3_change,                    
                            'nt_change': nt_change,                      
                            'aa_len': aa_len,                            
                            'nt_len': nt_len,                            
                            'cdna_pos': cdna_pos,                        
                            'meta': record} 
                disrupt_inframe_indels["{}:{}".format(record.CHROM, record.POS)].append(var_meta)

    # Removing InDels without annotation metadata.
    ks_to_del = []
    for k,v in disrupt_inframe_indels.items():
        if not v:
            ks_to_del.append(k)
    for k in ks_to_del:
        del(disrupt_inframe_indels[k])
    return disrupt_inframe_indels


def extract_frameshift_indels(input_vcf):
    """
    Focusing on deletions now, need to find a good example of a conservative_inframe_insertion
    """
    frameshift_indels = {}
    if input_vcf.endswith('gz'):
        vcf_reader = vcf.Reader(open(input_vcf), 'r', compressed=True)
    else:
        vcf_reader = vcf.Reader(open(input_vcf), 'r')

    for record in vcf_reader:
        annotations = [x for x in record.INFO['ANN']]
        frameshift_indels["{}:{}".format(record.CHROM, record.POS)] = []
        for annotation in annotations:
            effects = annotation.split('|')
            if re.search('^frameshift_variant$', effects[1]) and effects[13] and effects[7] == 'protein_coding' and effects[-1] == '':
                transcript = effects[6]
                aa_len = effects[13].split('/')[1]
                nt_change = effects[9].lstrip('c.')
                nt_len = effects[12].split('/')[1]
                cdna_pos = effects[11].split('/')[0]
                var_meta = {'transcript': transcript,
                            'nt_change': nt_change,                   
                            'aa_len': aa_len,                         
                            'nt_len': nt_len,                         
                            'cdna_pos': cdna_pos,                     
                            'meta': record}
                frameshift_indels["{}:{}".format(record.CHROM, record.POS)].append(var_meta)

    # Removing InDels without annotation metadata.
    ks_to_del = []
    for k,v in frameshift_indels.items():
        if not v:
            ks_to_del.append(k)
    for k in ks_to_del:
        del(frameshift_indels[k])
    return frameshift_indels


def make_snv_peptides_context(args):
    """
    """
    mutant_peptides = {}
    mutant_seqs = {}
    wildtype_peptides = {}

    print("Loading missense SNVs...")
    missense_snvs = extract_missense_snvs(args.somatic_vcf)
    print("Loaded missense SNVs.")

    somatic_txs = []

    print("Loading transcripts harboring somatic variants...")
    with open(args.somatic_txs) as fo:
        for line in fo.readlines():
            somatic_txs.append(line.strip())
    print("Loaded transcripts harboring somatic variants.")

    print("Loading variant transcripts metadata...")
    variant_txs_metadata = {} 
   
    for variant_tx in somatic_txs:
        variant_txs_metadata[variant_tx] = {}
        variant_txs_metadata[variant_tx]['cds'] = []
        variant_txs_metadata[variant_tx]['strand'] = ''

    with open(args.gtf) as fo:
        for line in fo.readlines():
            print(line.split('\t'))
            print(line.split('\t')[8].split('; '))
            print(line.split('\t')[8].split('; ')[2].replace('"', '').replace('transcript_id ', ''))
#            for variant_tx in somatic_txs:
#                if re.search(variant_tx, line):
#                    if re.search('\tCDS\t', line):
            tx_id_idx = ''
            meta_entries = line.split('\t')[8].split('; ')
            for meta_entry_idx, meta_entry in enumerate(meta_entries):
                if re.search('transcript_id ', meta_entry):
                    tx_id_idx = meta_entry_idx
                    break
            variant_tx = line.split('\t')[8].split('; ')[meta_entry_idx].replace('"', '').replace('transcript_id ', '')
            chr = line.split('\t')[0]
            start = line.split('\t')[3]
            stop = line.split('\t')[4]
            strand = line.split('\t')[6] 
            variant_txs_metadata[variant_tx]['strand'] = strand
            coords = "{}:{}-{}".format(chr, start, stop)
            if coords not in variant_txs_metadata[variant_tx]['cds']:
                variant_txs_metadata[variant_tx]['cds'].append("{}:{}-{}".format(chr, start, stop))
    print("Loaded variant transcripts metadata.")

#    expressed_txs_nts = {}
#    expressed_txs_aas = {}

    for entry in missense_snvs.keys():
        for record in missense_snvs[entry]:
#            variant_tx = record['transcript'].partition('.')[0]
            variant_tx = record['transcript']
            if variant_tx not in somatic_txs:
                continue
            norm_nts = []
            norm_aas = []
            tumor_nts = []
            tumor_aas = []
            print("Creating peptides for record: {}".format(record))
            record_coords = "{}_{}".format(record['meta'].CHROM, record['meta'].POS)
            print("Record coordinates: {}".format(record_coords))
            print("Record transcript: {}".format(record['transcript']))

            exon_seqs = {}
            exon_seqs['norm'] = {}
            exon_seqs['tumor'] = {}
            print("Loading normal exonic sequences...")
#            if glob(os.path.join(args.var_tx_seqs, '*{}_{}.normal.fa'.format(record['transcript'].partition('.')[0], record_coords))) and glob(os.path.join(args.var_tx_seqs, '*{}_{}.tumor.fa'.format(record['transcript'].partition('.')[0], record_coords))): 
            if glob(os.path.join(args.var_tx_seqs, '*{}_{}.normal.fa'.format(record['transcript'], record_coords))) and glob(os.path.join(args.var_tx_seqs, '*{}_{}.tumor.fa'.format(record['transcript'], record_coords))): 
#                for seq_record in SeqIO.parse(glob(os.path.join(args.var_tx_seqs, '*{}_{}.normal.fa'.format(record['transcript'].partition('.')[0], record_coords)))[0], "fasta"):
                for seq_record in SeqIO.parse(glob(os.path.join(args.var_tx_seqs, '*{}_{}.normal.fa'.format(record['transcript'], record_coords)))[0], "fasta"):
                    exon_seqs['norm'][seq_record.description] = seq_record.seq
                print("Loaded normal exonic sequences.")
                print("Loading tumor exonic sequences...")
#                for seq_record in SeqIO.parse(glob(os.path.join(args.var_tx_seqs, '*{}_{}.tumor.fa'.format(record['transcript'].partition('.')[0], record_coords)))[0], "fasta"):
                for seq_record in SeqIO.parse(glob(os.path.join(args.var_tx_seqs, '*{}_{}.tumor.fa'.format(record['transcript'], record_coords)))[0], "fasta"):
                    exon_seqs['tumor'][seq_record.description] = seq_record.seq
                print("Loaded tumor exonic sequences.")
            else:
                continue

            variant_tx_seqs = {}

            samp_types = ['norm', 'tumor']
            for i in samp_types:
                variant_tx_seqs[i] = ''
            if variant_txs_metadata[variant_tx]['strand'] == '+':
                print("Transcript is on positive strand.")
                for i in samp_types:
                    for cds in sorted(variant_txs_metadata[variant_tx]['cds']):
                        variant_tx_seqs[i] += exon_seqs[i][cds]
            elif variant_txs_metadata[variant_tx]['strand'] == '-':
                print("Transcript is on negative strand.")
                for i in samp_types:
                    for cds in sorted(variant_txs_metadata[variant_tx]['cds'], reverse=True):
                        variant_tx_seqs[i] += exon_seqs[i][cds].reverse_complement()

            # By default, snv is expected to be in the 30th position.
            # upstream_offset and downstream_offset are to properly localize
            # reading frame.
            snv_pos = 30
            upstream_offset = 0
            downstream_offset = 2
            if (float(record['nt_pos']) - 1) % 3 == 0:
                print("Direct translate")
                upstream_offset += 0
                downstream_offset += 2
            elif (float(record['nt_pos']) - 2) % 3 == 0:
                print("Minus one")
                upstream_offset += -1
                downstream_offset += 1
                snv_pos += 1
            elif (float(record['nt_pos']) - 3) % 3 == 0:
                print("Minus two")
                upstream_offset += -2
                downstream_offset += 0
                snv_pos += 2
            for i in samp_types:
                target_seq = ''
                # Walking through to find position
                steps = 1
                seq = variant_tx_seqs[i]
                print("Sample: {}".format(i))
                snv_pos = 30
                for pos_idx, pos in enumerate(seq):
                    if int(steps) == int(record['nt_pos']):
                        if (int(pos_idx) - 30) < 0:
                            print("Variant is within first 30 bases of transcript. Setting snv_pos to pos_idx...")
                            snv_pos = pos_idx
                        break
                    elif re.search('[A-Z]', pos):
                        steps += 1

                # Add bases upstream or downstream to compensate for neighboring deletions. 
                suff_aa_len = False
                to_add_left = 0
                to_add_right = 0
                while not suff_aa_len:
                    target_seq = list(seq[max(0, pos_idx-30+upstream_offset-to_add_left):min(pos_idx+31 + downstream_offset + to_add_right, len(seq))])
                    print(''.join(target_seq))
                    reduced_target_seq = str(target_seq).replace('X', '')
                                   
                    if len(reduced_target_seq) > 63:
                        print("Sufficient AA length!")
                        suff_aa_len = True
                    else:
                        print("Insufficient AA length!")
                        to_add_left = target_seq[:snv_pos].count('X')
                        to_add_right = target_seq[snv_pos:].count('X')
                        print("To add left: {}".format(to_add_left))
                        print("To add right: {}".format(to_add_right))
         
                if re.search('tumor', i): 
                    if target_seq[snv_pos] == record['alt_nt'] or (record['ref_nt'] in iupac_conversion(target_seq[snv_pos]) and record['alt_nt'] in iupac_conversion(target_seq[snv_pos])):
                        print("Before conversion: {}".format(''.join(target_seq)))
                        list_target_seq = list(target_seq)
                        list_target_seq[snv_pos] = record['alt_nt']
                        converted_seq = ''.join(list_target_seq)
                        print("After  conversion: {}".format(converted_seq))
                        tumor_nts.append(''.join(converted_seq).replace('X', ''))
                    else:
                        print("Cannot add this transcript. The variant position differs between annotations and transcript sequence. Check your references.")
                else:
                    if target_seq[snv_pos] == record['ref_nt']:
                        norm_nts.append(''.join(target_seq))
                    else:
                        print(target_seq[snv_pos])
                        print(target_seq[snv_pos-10:snv_pos+10])
                        print(record['ref_nt'])
    
            for seq in norm_nts:
                norm_aas.append(Seq(''.join(seq).replace('X', '')).translate()) 
            for seq in tumor_nts:
                tumor_aas.append(Seq(''.join(seq).replace('X', '')).translate()) 
                         
            if len(list(set(tumor_aas))) != 1 or len(list(norm_aas)) == 0:
                print("Something is wrong. Check your tumor sequences.")
                print("Something is wrong. Check your transcript versions.")
            else:
                # Making checksum identifier 
                hash_components = [str(record['meta'].CHROM), str(record['meta'].POS), str(record['transcript']), str(record['meta'].REF), str(record['meta'].ALT)]
                var_md5 = hashlib.md5("{}".format(':'.join(hash_components).encode('utf-8'))).hexdigest()[:16]
                #peptide_entry
                mutant_peptides["MD5:{} VARIANT_POS:{}_{} TX_POS:{} TRANSCRIPT:{} REF:{} ALT:{} SNV_TYPE:{} PROTEIN_CONTEXT:NA GENOMIC_CONTEXT:NA".format(var_md5, record['meta'].CHROM, record['meta'].POS, record['nt_pos'], record['transcript'], record['meta'].REF, record['meta'].ALT, 'missense')] = tumor_aas[0]
                #nt_entry
                mutant_seqs["MD5:{} VARIANT_POS:{}_{} TX_POS:{} CDNA_POS:{} TRANSCRIPT:{} REF:{} ALT:{} SNV_TYPE:{} PROTEIN_CONTEXT:NA GENOMIC_CONTEXT:NA".format(var_md5, record['meta'].CHROM, record['meta'].POS, record['nt_pos'], record['cdna_pos'], record['transcript'], record['meta'].REF, record['meta'].ALT, 'missense')] = tumor_nts[0]
                if len(list(set(norm_aas))) > 1:
                    for norm_aa_idx, norm_aa in enumerate(norm_aas): 
                        #wt entry (multiple alleles)
                        wildtype_peptides["MD5{}:{} VARIANT_POS:{}_{} TRANSCRIPT:{} REF:{} ALT:{} SNV_TYPE:{} PROTEIN_CONTEXT:NA GENOMIC_CONTEXT:NA".format(norm_aa_idx, var_md5, record['meta'].CHROM, record['meta'].POS, record['transcript'], record['meta'].REF, record['meta'].ALT, 'missense')] = norm_aa
                else: 
                    #wt entry (single allele)
                    wildtype_peptides["MD5:{} VARIANT_POS:{}_{} TRANSCRIPT:{} REF:{} ALT:{} SNV_TYPE:{} PROTEIN_CONTEXT:NA GENOMIC_CONTEXT:NA".format(var_md5, record['meta'].CHROM, record['meta'].POS, record['transcript'], record['meta'].REF, record['meta'].ALT, 'missense')] = norm_aas[0]


    with open(args.mt_output, 'w') as ofo:
        for k, v in mutant_peptides.items():
            ofo.write('>{}\n{}\n'.format(k, v))
    with open(args.nt_output, 'w') as ofo:
        for k, v in mutant_seqs.items():
            ofo.write('>{}\n{}\n'.format(k, v))
    with open(args.wt_output, 'w') as ofo:
        for k, v in wildtype_peptides.items():
            ofo.write('>{}\n{}\n'.format(k, v))


# Remove
#def make_snv_peptides(args):
#    """
#    Make peptide sequences derived from mutant (somatic) and normal (germline) sequences.
#    """
#    tx_to_aa = load_tx_aas(args)
#    tx_to_cds = load_tx_cds(args)
#
#
#    mutant_peptides = {}
#    reference_peptides = {}
#
#    for entry in missense_snvs.keys():
#        record = missense_snvs[entry]
#        print(record)
#        tx = record['transcript']
#        tx_no_version = tx.split('.')[0]
#        if tx not in tx_to_aa.keys() and tx_no_version not in tx_to_aa.keys():
#            print("{} cannot be found in the transcript-to-amino acid dictionary. Continuing to next record...".format(record['transcript']))
#            continue
#        if tx not in tx_to_cds.keys() and tx_no_version not in tx_to_cds.keys():
#            print("{} cannot be found in the transcript-to-coding sequence dictionary. Continuing to next record...".format(record['transcript']))
#            continue
#
#        try:
#            bufr_aa = list(tx_to_aa[tx])
#        except:
#            bufr_aa = list(tx_to_aa[tx_no_version])
#        try:
#            bufr_nt = list(tx_to_cds[tx])
#        except:
#            bufr_nt = list(tx_to_cds[tx_no_version])
#
#        if len(bufr_aa) != int(record['aa_len']):
#            print("transcript {} shows different lengths between amino acid fasta ({}) and snpEff annotations! ({})".format(tx, record['aa_len'], len(bufr_aa)))
#        if len(bufr_nt) != int(record['nt_len']):
#            print("transcript {} shows different lengths between coding sequence fasta ({}) and snpEff annotations! ({})".format(tx, record['nt_len'], len(bufr_nt)))
#
#        # Deep copy to mut_aa...
#        # This is the step where the annotated germline reference should be
#        # utilized...
#        mut_aa = bufr_aa[:]
#        mut_nt = bufr_nt[:]
#        aa_pos = int(record['aa_pos']) - 1
#        nt_pos = int(record['nt_pos']) - 1
#        # This check should be performed prior to incorporating germline variants.
#
#        if bufr_aa[aa_pos] != record['ref_aa']:
#            print("Reference amino acid from snpEff annotation ({}) doesn't match amino acid from amino acid FASTA ({})! Continuing to next record...".format(record['ref_aa'], bufr_aa[aa_pos]))
#            continue
#        if bufr_nt[nt_pos] != record['ref_nt']:
#            print("Reference nucleotide from snpEff annotation ({}) doesn't match nucleotide from amino acid FASTA ({})! Continuing to next record...".format(record['ref_nt'], bufr_nt[nt_pos]))
#            continue
#
#        # Applying the mutated amino acid
#        mut_aa[aa_pos] = record['alt_aa']
#        mut_nt[nt_pos] = record['alt_nt']
#
#        ref_peptide = ''.join(bufr_aa[max(aa_pos-args.length+1,0):min(aa_pos+args.length, len(bufr_aa))])
#        mut_peptide = ''.join(mut_aa[max(aa_pos-args.length+1, 0):min(aa_pos+args.length, len(mut_aa))])
#
#        print("Reference peptide: {}".format(ref_peptide))
#        print("Mutant peptide:    {}".format(mut_peptide))
#
#        codon_pos = record['codon_pos'] - 1
#        ref_nuc = ''.join(bufr_nt[max(codon_pos-args.context_nt_length, 0):min(codon_pos+args.context_nt_length, len(bufr_nt)) + 3])
#        mut_nuc = ''.join(mut_nt[max(codon_pos-args.context_nt_length, 0):min(codon_pos+args.context_nt_length, len(mut_nt)) + 3])
#
#        print("Reference genomic context: {}".format(ref_nuc))
#        print("Mutant genomic context:    {}".format(mut_nuc))
#
#        header_mut_peptide = ''.join(mut_aa[max(aa_pos-13, 0):min(aa_pos+14, len(mut_aa))])
#
#        #var_md5 is being used to create a unique identifier for the resulting mutant peptide to overcome netMHCpan's length limitations.
#        var_md5 = hashlib.md5("{}".format(':'.join([str(record['meta'].CHROM), str(record['meta'].POS), str(tx), str(record['meta'].REF), str(record['meta'].ALT)])).encode('utf-8')).hexdigest()[:16]
#        mutant_peptides["MD5:{} VARIANT_POS:{}_{} TRANSCRIPT:{} REF:{} ALT:{} SNV_TYPE:{} PROTEIN_CONTEXT:{} GENOMIC_CONTEXT:{}".format(var_md5, record['meta'].CHROM, record['meta'].POS, tx, record['meta'].REF, record['meta'].ALT, 'missense', header_mut_peptide, mut_nuc)] = mut_peptide
#        reference_peptides["MD5:{} VARIANT_POS:{}_{} TRANSCRIPT:{} REF:{} ALT:{}".format(var_md5, record['meta'].CHROM, record['meta'].POS, tx, record['meta'].REF, record['meta'].ALT)] = ref_peptide
#
#    with open(args.mt_output, 'w') as ofo:
#        for k, v in mutant_peptides.items():
#            ofo.write('>{}\n{}\n'.format(k, v))
#    with open(args.wt_output, 'w') as ofo:
#        for k, v in reference_peptides.items():
#            ofo.write('>{}\n{}\n'.format(k, v))


# Remove
#def get_snv_genomic_context(args):
#    """
#    """
#    tx_to_cds = load_tx_cds(args)
#    missense_snvs = extract_missense_snvs(args.vcf)
#    emitted_nucs = {}
#    for record in missense_snvs:
#        if record[0] not in tx_to_cds.keys():
#            continue
#        tx_ref_seq = list(tx_to_cds[record[0]])
##        if len(tx_ref_seq) != int(record[2]):
##            print("transcript {} shows different lengths between peptide fasta ({}) and snpeff! ({})".format(record[0], record[2], len(tx_ref_seq)))
#
#        mut_seq = tx_ref_seq[:]
#        pos = int(record[1]) - 1
##        if mut_seq[pos] != record[3]:
##            print("Reference amino acid doesn't match! Something has gone horribly wrong.")
##            continue
#        # Have to be careful here... we want to start at the first base of the affected codon
#        mut_seq[pos] = record[4]
#        codon_pos = record[5] - 1
#        ref_nuc = ''.join(tx_ref_seq[max(codon_pos-args.length-1, 0):min(codon_pos+args.length, len(tx_ref_seq)) + 3])
#        mut_nuc = ''.join(mut_seq[max(codon_pos-args.length, 0):min(codon_pos+args.length, len(mut_seq)) + 3])
#        print("{}\n{}\n{}".format(record[0],ref_nuc, mut_nuc))
#        print(':'.join([str(record[-1].CHROM), str(record[-1].POS), str(record[0]), str(record[-1].REF), str(record[-1].ALT)]))
#        var_md5 = hashlib.md5("{}".format(':'.join([str(record[-1].CHROM), str(record[-1].POS), str(record[0]), str(record[-1].REF), str(record[-1].ALT)])).encode('utf-8')).hexdigest()[:16]
#        emitted_nucs["{}\t{}:{}\t{}\t{}\t{}".format(var_md5, record[-1].CHROM, record[-1].POS, record[0], record[-1].REF, record[-1].ALT)] = mut_nuc
#
#    with open(args.output, 'w') as ofo:
#        for k, v in emitted_nucs.items():
#            ofo.write('>{}\t{}\n'.format(k, v))

def make_indel_peptides_context(args):
    """
    """
    somatic_txs = []

    print("Loading all transcripts harboring somatic variants...")
    with open(args.somatic_txs) as fo:
        for line in fo.readlines():
            somatic_txs.append(line.strip())
    print("Loaded all transcripts harboring somatic variants.")

    print("Somatic transcripts (harboring somatic vars): {}".format(somatic_txs))

    # Loading actual somatic variants of interest (all InDels)
    # We will need to retrieve the nucleotide sequence records as well.
    print("Getting indels...")
    conserv_inframe_indels = extract_conservative_inframe_indels(args.somatic_vcf)
    print("Retrieved {} conservative inframe indels...".format(len(conserv_inframe_indels)))
    disrupt_inframe_indels = extract_disruptive_inframe_indels(args.somatic_vcf)
    print("Retrieved {} disruptive inframe indels...".format(len(disrupt_inframe_indels)))
    frameshift_indels = extract_frameshift_indels(args.somatic_vcf)
    print("Retrieved {} frameshift indels...".format(len(frameshift_indels)))
    
    inframe_indels = dict(conserv_inframe_indels, **disrupt_inframe_indels)

    mutant_peptides = {}
    mutant_seqs = {}

    # Dictionary for holding exon coding sequences for combining to create
    # transcript sequences.
    variant_txs_metadata = {} 
   
    for expressed_tx in somatic_txs:
        variant_txs_metadata[expressed_tx] = {}
        variant_txs_metadata[expressed_tx]['cds'] = []
        variant_txs_metadata[expressed_tx]['strand'] = ''

    with open(args.gtf) as fo:
        for line in fo.readlines():
            print(line.split('\t'))
            print(line.split('\t')[8].split('; '))
            print(line.split('\t')[8].split('; ')[2].replace('"', '').replace('transcript_id ', ''))
#            for variant_tx in somatic_txs:
#                if re.search(variant_tx, line):
#                    if re.search('\tCDS\t', line):
            meta_entries = line.split('\t')[8].split('; ')
            for meta_entry_idx, meta_entry in enumerate(meta_entries):
                if re.search('transcript_id', meta_entry):
                    tx_id_idx = meta_entry_idx
                    break
            variant_tx = line.split('\t')[8].split('; ')[meta_entry_idx].replace('"', '').replace('transcript_id ', '')
            chr = line.split('\t')[0]
            start = line.split('\t')[3]
            stop = line.split('\t')[4]
            strand = line.split('\t')[6] 
            variant_txs_metadata[variant_tx]['strand'] = strand
            coords = "{}:{}-{}".format(chr, start, stop)
            if coords not in variant_txs_metadata[variant_tx]['cds']:
                variant_txs_metadata[variant_tx]['cds'].append("{}:{}-{}".format(chr, start, stop))
    print("Loaded variant transcripts metadata.")
#    with open(args.gtf) as fo:
#        for line in fo.readlines():
#            for expressed_tx in somatic_txs:
#                if re.search(expressed_tx, line):
#                    if re.search('\tCDS\t', line):
#                        chr = line.split('\t')[0]
#                        start = line.split('\t')[3]
#                        stop = line.split('\t')[4]
#                        strand = line.split('\t')[6] 
#                        variant_txs_metadata[expressed_tx]['strand'] = strand
#                        coords = "{}:{}-{}".format(chr, start, stop)
#                        if coords not in variant_txs_metadata[expressed_tx]['cds']:
#                            variant_txs_metadata[expressed_tx]['cds'].append("{}:{}-{}".format(chr, start, stop))

    expressed_txs_nts = {}

    tumor_aas = []
    for entry in inframe_indels.keys():
        for record in inframe_indels[entry]:
#            variant_tx = record['transcript'].partition('.')[0]
            variant_tx = record['transcript']
            if variant_tx not in somatic_txs:
                continue
            print("NEW INFRAME RECORD: {}".format(record))
            record_coords = "{}_{}".format(record['meta'].CHROM, record['meta'].POS)
            print(record_coords)
            print(variant_tx)
            strands = ['tumor']
            
            print(glob(os.path.join(args.var_tx_seqs, '*{}_{}.tumor.fa'.format(record['transcript'].partition('.')[0], record_coords))))
            exon_seqs = {}
            exon_seqs['tumor'] = {}
            for seq_record in SeqIO.parse(glob(os.path.join(args.var_tx_seqs, '*{}_{}.tumor.fa'.format(record['transcript'].partition('.')[0], record_coords)))[0], "fasta"):
                exon_seqs['tumor'][seq_record.description] = seq_record.seq
            
            variant_tx_seqs = {}

            samp_types = ['tumor']
            for i in samp_types:
                variant_tx_seqs[i] = ''
            if variant_txs_metadata[variant_tx]['strand'] == '+':
                print('positive_strand')
                for i in samp_types:
                    for cds in sorted(variant_txs_metadata[variant_tx]['cds']):
                        variant_tx_seqs[i] += exon_seqs[i][cds]
            elif variant_txs_metadata[variant_tx]['strand'] == '-':
                print('negative_strand')
                for i in samp_types:
                    for cds in sorted(variant_txs_metadata[variant_tx]['cds'], reverse=True):
                        variant_tx_seqs[i] += exon_seqs[i][cds].reverse_complement()

            indel_start = 0
            indel_stop = 0
            indel_len = 0

            print(record)
            
            # DELETIONS #
            if re.search("del", record['nt_change']):
                print("INFRAME DELETION: {}:{}".format(record['meta'].CHROM, record['meta'].POS))
                del_rec_coords, buffer, del_rec_nt = record['nt_change'].partition('del')
                if '_' in del_rec_coords:
                    indel_start = del_rec_coords.split('_')[0]
                    indel_stop = del_rec_coords.split('_')[1]
                    indel_len = int(indel_stop) - int(indel_start) + 2
                else:
                    indel_start = del_rec_coords
                    indel_stop = del_rec_coords
                    indel_len = 1

            # INSERTIONS #
            elif re.search("ins", record['nt_change']):
                print("INFRAME INSERTION: {}:{}".format(record['meta'].CHROM, record['meta'].POS))
                ins_rec_nt, bufr, ins_nt = record['nt_change'].partition('ins')
                indel_start = int(ins_rec_nt.split('_')[0]) - 1
                indel_stop = int(ins_rec_nt.split('_')[1]) - 1
                indel_len = int(indel_stop) - int(indel_start) + 2

#            # INSERTIONS/DELETIONS #
#            elif re.search("delins", record['aa3_change']):
#                print("INFRAME INSERTION/DELETION: {}:{}".format(record['meta'].CHROM, record['meta'].POS))
#              
#                pass
#                del_rec_aa, bufr, ins_aa3 = record['aa3_change'].partition('delins')
#    
#                ins_aa = seq1(ins_aa3)
#                start_pos_aa = int(del_rec_aa.split('_')[0][3:]) - 1
#                start_aa3 = del_rec_aa.split('_')[0][:3]
#                start_aa = seq1(start_aa3)
#                stop_pos_aa = 0
#                stop_aa = 'foo'

            # DUPLICATIONS #
            elif re.search("dup", record['nt_change']):
                dup_rec_nt, bufr, ins_nt = record['nt_change'].partition('dup')
                if '_' in dup_rec_nt:
                    indel_start = dup_rec_nt.split('_')[0]
                    indel_stop = dup_rec_nt.split('_')[1]
                    indel_len = int(indel_stop) - int(indel_start) + 2
                else:
                    indel_start = dup_rec_nt
                    indel_stop = dup_rec_nt
                    indel_len = 1
#                indel_start = int(dup_rec_nt.split('_')[0]) - 1
#                indel_stop = int(dup_rec_nt.split('_')[1]) - 1
#                indel_len = int(indel_stop) - int(indel_start) + 1
  
            print("{} {} {} {}".format(indel_start, indel_stop, indel_len, record['nt_change']))


          
                
            snv_pos = 30
            upstream_offset = 0
            downstream_offset = 2
            if (float(indel_start) - 1) % 3 == 0:
                print("Direct translate")
                upstream_offset += 0
                downstream_offset += 2
            elif (float(indel_start) - 2) % 3 == 0:
                print("Minus one")
                upstream_offset += -1
                downstream_offset += 1
                snv_pos += 1
            elif (float(indel_start) - 3) % 3 == 0:
                print("Minus two")
                upstream_offset += -2
                downstream_offset += 0
                snv_pos += 2
#            variant_tx_seqs = ''
#            if variant_txs_metadata[variant_tx]['strand'] == '+':
#                print('positive_strand')
#                for i in strands:
#                    for cds in sorted(variant_txs_metadata[variant_tx]['cds']):
#                        variant_tx_seqs += exon_seqs[i][cds]
#            elif variant_txs_metadata[variant_tx]['strand'] == '-':
#                print('negative_strand')
#                for i in strands:
#                    for cds in sorted(variant_txs_metadata[variant_tx]['cds'], reverse=True):
#                        variant_tx_seqs += exon_seqs[i][cds].reverse_complement()

            print(variant_tx_seqs)
            target_seq = ''
            # Walking through to find position
            steps = 1
            seq = variant_tx_seqs['tumor']#record['transcript'].partition('.')[0]][i]
            for pos_idx, pos in enumerate(list(seq)):
                if int(steps) == int(indel_start):
                    print(steps)
                    print(indel_start)
                    print(pos_idx)
                    if (int(pos_idx) - 30) < 0:
                        snv_pos = pos_idx
                    target_seq = list(seq[max(0, pos_idx-30+upstream_offset):min(pos_idx+31 + downstream_offset + indel_len, len(seq))])
                    break
                elif re.search('[A-Z]', pos):
                    steps += 1  
            target_seq = ''.join(target_seq)
            target_seq = target_seq.replace('X', '')
            print("Target seq: {}".format(target_seq))
            print("Target AA: {}".format(Seq(target_seq).translate()))
            tumor_aas = []
            tumor_nts = []
            tumor_aas.append("{}".format(Seq(target_seq).translate(to_stop=True)))
            tumor_nts.append("{}".format(Seq(target_seq)))
            if len(list(set(tumor_aas))) > 1:
                print("Something is wrong. We shouldn't have two tumor alleles.")
            else:
                var_md5 = hashlib.md5("{}".format(':'.join([str(record['meta'].CHROM), str(record['meta'].POS), str(record['transcript']), str(record['meta'].REF), str(record['meta'].ALT)])).encode('utf-8')).hexdigest()[:16]
                mutant_peptides["MD5:{} VARIANT_POS:{}_{} TRANSCRIPT:{} REF:{} ALT:{} INDEL_TYPE:{} PROTEIN_CONTEXT:NA GENOMIC_CONTEXT:NA".format(var_md5, record['meta'].CHROM, record['meta'].POS, record['transcript'], record['meta'].REF, record['meta'].ALT, 'inframe')] = tumor_aas[0]
                mutant_seqs["MD5:{} VARIANT_POS:{}_{} CDNA_POS:{} TRANSCRIPT:{} REF:{} ALT:{} INDEL_TYPE:{} PROTEIN_CONTEXT:NA GENOMIC_CONTEXT:NA".format(var_md5, record['meta'].CHROM, record['meta'].POS, record['cdna_pos'], record['transcript'], record['meta'].REF, record['meta'].ALT, 'inframe')] = tumor_nts[0]
                         
    
#                    stop_aa3 = dup_rec_aa.split('_')[1][:3]
#                    stop_aa = seq1(stop_aa3)
#                else:
#                    stop_pos_aa = start_pos_aa
#                    stop_aa3 = start_aa3
#                    stop_aa = start_aa
#    
#                print("{} {} {}".format(start_pos_aa, start_aa, start_aa3))
#                print("{} {} {}".format(stop_pos_aa, stop_aa, stop_aa3))
#    
#                if re.search('ins', record['nt_change']):
#                    ins_rec_nt, bufr, ins_nt = record['nt_change'].partition('ins')
#                    start_pos_nt = int(ins_rec_nt.split('_')[0]) - 1
#                    stop_pos_nt = int(ins_rec_nt.split('_')[1]) - 1
#                    mut_nt = ''.join(bufr_nt[:start_pos_nt+1] + [ins_nt] + bufr_nt[stop_pos_nt:])
#
#    # Do stuff here to put these peptides into the mutant_peptides dict.


    for entry in frameshift_indels.keys():
        for record in frameshift_indels[entry]:
            variant_tx = record['transcript'].partition('.')[0]
            if variant_tx not in somatic_txs:
                continue
            print("NEW FRAMESHIFT RECORD: {}".format(record))
            record_coords = "{}_{}".format(record['meta'].CHROM, record['meta'].POS)
            print(record_coords)
            print(variant_tx)
            strands = ['tumor']
            
            print(glob(os.path.join(args.var_tx_seqs, '*{}_{}.tumor.fa'.format(record['transcript'].partition('.')[0], record_coords))))
            exon_seqs = {}
            exon_seqs['tumor'] = {}
            for seq_record in SeqIO.parse(glob(os.path.join(args.var_tx_seqs, '*{}_{}.tumor.fa'.format(record['transcript'].partition('.')[0], record_coords)))[0], "fasta"):
                exon_seqs['tumor'][seq_record.description] = seq_record.seq

            print(exon_seqs)
            
            variant_tx_seqs = {}

#            strands = ['norm', 'tumor']
            strands = ['tumor']
            for i in strands:
                variant_tx_seqs[i] = ''
            if variant_txs_metadata[variant_tx]['strand'] == '+':
                print('positive_strand')
                for i in strands:
                    for cds in sorted(variant_txs_metadata[variant_tx]['cds']):
                        print(cds)
                        variant_tx_seqs[i] += exon_seqs[i][cds]
            elif variant_txs_metadata[variant_tx]['strand'] == '-':
                print('negative_strand')
                for i in strands:
                    for cds in sorted(variant_txs_metadata[variant_tx]['cds'], reverse=True):
                        print(cds)
                        variant_tx_seqs[i] += exon_seqs[i][cds].reverse_complement()

            print(variant_tx_seqs)

            indel_start = 0
            indel_stop = 0
            indel_len = 0

            # DELETIONS #
#            if re.search("del$", record['aa3_change']):
            if re.search("del", record['nt_change']):
                print("FRAMESHIFT DELETION: {}:{}".format(record['meta'].CHROM, record['meta'].POS))
                del_rec_coords, buffer, del_rec_nt = record['nt_change'].partition('del')

                if '_' in del_rec_coords:
                    indel_start = del_rec_coords.split('_')[0]
                    indel_stop = del_rec_coords.split('_')[1]
                    indel_len = int(indel_stop) - int(indel_start) + 2
                else:
                    indel_start = del_rec_coords
                    indel_stop = del_rec_coords
                    indel_len = 1
            # INSERTIONS #
#            elif re.search("[0-9]ins", record['aa3_change']):
            elif re.search("ins", record['nt_change']):
                print("INFRAME INSERTION: {}:{}".format(record['meta'].CHROM, record['meta'].POS))
                ins_rec_nt, bufr, ins_nt = record['nt_change'].partition('ins')

                indel_start = int(ins_rec_nt.split('_')[0]) - 1
                indel_stop = int(ins_rec_nt.split('_')[1]) - 1
                indel_len = int(indel_stop) - int(indel_start) + 2
#            # INSERTIONS/DELETIONS #
#            elif re.search("delins", record['aa3_change']):
#                print("INFRAME INSERTION/DELETION: {}:{}".format(record['meta'].CHROM, record['meta'].POS))
#              
#                pass
#                del_rec_aa, bufr, ins_aa3 = record['aa3_change'].partition('delins')
#    
#                ins_aa = seq1(ins_aa3)
#                start_pos_aa = int(del_rec_aa.split('_')[0][3:]) - 1
#                start_aa3 = del_rec_aa.split('_')[0][:3]
#                start_aa = seq1(start_aa3)
#                stop_pos_aa = 0
#                stop_aa = 'foo'

            # DUPLICATIONS #
            elif re.search("dup", record['nt_change']):
                dup_rec_coords, bufr, ins_nt = record['nt_change'].partition('dup')
                #indel_start = int(dup_rec_nt.split('_')[0]) - 1
                #indel_stop = int(dup_rec_nt.split('_')[1]) - 1
                #indel_len = int(indel_stop) - int(indel_start) + 1
                if '_' in dup_rec_coords:
                    indel_start = dup_rec_coords.split('_')[0]
                    indel_stop = dup_rec_coords.split('_')[1]
                    indel_len = int(indel_stop) - int(indel_start) + 2
                else:
                    indel_start = dup_rec_coords
                    indel_stop = dup_rec_coords
                    indel_len = 1
  
            print("{} {} {} {}".format(indel_start, indel_stop, indel_len, record['nt_change']))


          
                
            snv_pos = 30
            upstream_offset = 0
            downstream_offset = 2
            if (float(indel_start) - 1) % 3 == 0:
                print("Direct translate")
                upstream_offset += 0
                downstream_offset += 2
            elif (float(indel_start) - 2) % 3 == 0:
                print("Minus one")
                upstream_offset += -1
                downstream_offset += 1
                snv_pos += 1
            elif (float(indel_start) - 3) % 3 == 0:
                print("Minus two")
                upstream_offset += -2
                downstream_offset += 0
                snv_pos += 2
            print(variant_tx_seqs)
            target_seq = ''
            # Walking through to find position
            steps = 1
            seq = variant_tx_seqs['tumor']#record['transcript'].partition('.')[0]][i]
            for pos_idx, pos in enumerate(list(seq)):
                if int(steps) == int(indel_start):
                    print(steps)
                    print(indel_start)
                    print(pos_idx)
                    if (int(pos_idx) - 30) < 0:
                        snv_pos = pos_idx
#                    target_seq = list(seq[max(0, pos_idx-30+upstream_offset):min(pos_idx+31 + downstream_offset + indel_len, len(seq))])
                    target_seq = list(seq[max(0, pos_idx-30+upstream_offset):])
                    break
                elif re.search('[A-Z]', pos):
                    steps += 1  
            target_seq = ''.join(target_seq)
            target_seq = target_seq.replace('X', '')
            print("Target seq: {}".format(target_seq))
            tumor_aas = []
            tumor_nts = []
            tumor_aas.append("{}".format(Seq(target_seq).translate(to_stop=True)))
            tumor_nts.append("{}".format(Seq(target_seq)))
            print("Target AA: {}".format(Seq(target_seq).translate(to_stop=True)))
            if len(list(set(tumor_aas))) > 1:
                print("Something is wrong. We shouldn't have two tumor alleles.")
            else:
                var_md5 = hashlib.md5("{}".format(':'.join([str(record['meta'].CHROM), str(record['meta'].POS), str(record['transcript']), str(record['meta'].REF), str(record['meta'].ALT)])).encode('utf-8')).hexdigest()[:16]
                mutant_peptides["MD5:{} VARIANT_POS:{}_{} TRANSCRIPT:{} REF:{} ALT:{} INDEL_TYPE:{} PROTEIN_CONTEXT:NA GENOMIC_CONTEXT:NA".format(var_md5, record['meta'].CHROM, record['meta'].POS, record['transcript'], record['meta'].REF, record['meta'].ALT, 'frameshift')] = tumor_aas[0]
                mutant_seqs["MD5:{} VARIANT_POS:{}_{} CDNA_POS:{} TRANSCRIPT:{} REF:{} ALT:{} INDEL_TYPE:{} PROTEIN_CONTEXT:NA GENOMIC_CONTEXT:NA".format(var_md5, record['meta'].CHROM, record['meta'].POS, record['cdna_pos'], record['transcript'], record['meta'].REF, record['meta'].ALT, 'inframe')] = tumor_nts[0]
                #mutant_seqs["MD5:{} VARIANT_POS:{}_{} TRANSCRIPT:{} REF:{} ALT:{} INDEL_TYPE:{} PROTEIN_CONTEXT:NA GENOMIC_CONTEXT:NA".format(var_md5, record['meta'].CHROM, record['meta'].POS, record['transcript'], record['meta'].REF, record['meta'].ALT, 'frameshift')] = tumor_nts[0]

#    sys.exit()
#            # DELETIONS #
#            print(record)
#            if re.search("del", record['nt_change']):
#                print("Frameshift deletion {}".format(record))
#                del_rec_nt, buffer, del_nt = record['nt_change'].partition('del')
#                start_pos_nt = 0
#                stop_pos_nt = 0
#                start_base = 'A'
#                stop_base = 'A'
#                if re.search('_', del_rec_nt):
#                     start_pos_nt = int(del_rec_nt.split('_')[0]) - 1
#                     stop_pos_nt = int(del_rec_nt.split('_')[1]) - 1
#                     if del_nt:
#                         start_nt = del_nt[0]
#                         stop_nt = del_nt[-1]
#                else:
#                    start_pos_nt = int(del_rec_nt) - 1
#                    stop_pos_nt = start_pos_nt
#                    if del_nt:
#                        start_nt = del_nt[0]
#                        stop_nt = start_base
#                print("{} {}".format(start_pos_nt, start_nt))
#                print("{} {}".format(stop_pos_nt, stop_nt))
#            # INSERTIONS #
#            if re.search("ins", record['nt_change']):
#                print("Frameshift insertion {}".format(record))
#                ins_rec_nt, buffer, ins_nt = record['nt_change'].partition('ins')
#                start_pos_nt = 0
#                stop_pos_nt = 0
#                start_nt = ''
#                stop_nt = ''
#                if re.search('_', ins_rec_nt) and not(re.search('\*', ins_rec_nt)):
#                     #Need to deal with insertions that extend beyond the annotated transcript.
#                     start_pos_nt = int(ins_rec_nt.split('_')[0]) - 1
#                     stop_pos_nt = int(ins_rec_nt.split('_')[1]) - 1
#                     if ins_nt:
#                         start_nt = ins_nt[0]
#                         stop_nt = ins_nt[-1]
#                elif not(re.search('\*', ins_rec_nt)):
#                    start_pos_nt = int(ins_rec_nt) - 1
#                    stop_pos_nt = start_pos_nt
#                    if ins_nt:
#                        start_nt = ins_seq[0]
#                        stop_nt = start_nt
#                else:
#                    pass
#                print("{} {}".format(start_pos_nt, start_nt))
#                print("{} {}".format(stop_pos_nt, stop_nt))
#            # DUPLICATIONS #
#            if re.search("dup", record['nt_change']):
#                print("Frameshift duplication {}".format(record))
#                dup_rec_nt, buffer, dup_nt = record['nt_change'].partition('dup')
#                start_pos_nt = 0
#                stop_pos_nt = 0
#                start_nt = 'A'
#                stop_nt = 'A'
#                if not(re.search('-', dup_rec_nt)):
#                    if re.search('_', dup_rec_nt):
#                         start_pos_nt = int(dup_rec_nt.split('_')[0]) - 1
#                         start_nt = dup_nt[0]
#                         stop_pos_nt = int(dup_rec_nt.split('_')[1]) - 1
#                         stop_nt = dup_nt[-1]
#                    else:
#                        start_pos_nt = int(dup_rec_nt) - 1
#                        start_nt = dup_nt[0]
#                        stop_pos_nt = start_pos_nt
#                        stop_nt = start_nt
#                else:
#                    continue
#    
#                print("{} {}".format(start_pos_nt, start_nt))
#                print("{} {}".format(stop_pos_nt, stop_nt))
#    
#            upstream_offset = 0
#            if (float(start_pos_nt)-1) % 3 == 0:
#                print("Direct translate")
#                upstream_offset += 0
#            elif (float(start_pos_nt) - 2) % 3 == 0:
#                print("Minus one")
#                upstream_offset += -1
#            elif (float(start_pos_nt) - 3) % 3 == 0:
#                print("Minus two")
#                upstream_offset += -2
#            for i in strands:
#    #Account for indels here within 33 bp.
#                target_seq = ''
#                # Walking through to find position
#                steps = 1
#                seq = expressed_txs_nts[record['transcript'].partition('.')[0]][i]
#                for pos_idx, pos in enumerate(seq):
#                   if int(steps) == int(start_pos_nt):
#                        
#                       target_seq = list(seq[pos_idx-27+upstream_offset:])
#                       break
#                   if pos_idx == 0:
#                       if pos != 'x' and pos in ['A', 'C', 'G', 'T'] or re.search('[acgtryswkmbdhv]', pos) and re.search('[ACGT]', seq[pos_idx+1]):
#                           steps += 1  
#                   else:
#                       if pos != 'x' and pos in ['A', 'C', 'G', 'T'] or re.search('[acgtryswkmbdhv]', pos) and (re.search('[ACGT]', seq[pos_idx-1]) and re.search('[ACGT]', seq[pos_idx+1])):
#                           steps +=1 
#                print(''.join(target_seq)[:60]) 
#                tumor_aa = Seq(''.join(target_seq)).translate(to_stop=True)
#                tumor_aas.append(tumor_aa)
#                
#        

    with open(args.output, 'w') as ofo:
        for k, v in mutant_peptides.items():
            ofo.write('>{}\n{}\n'.format(k, v))
    with open(args.nt_output, 'w') as ofo:
        for k, v in mutant_seqs.items():
            ofo.write('>{}\n{}\n'.format(k, v))


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
#            nt_context = ''.join(mut_nt[start_pos_nt-10:stop_pos_nt+10])
            nt_context = ''.join(mut_nt[max(start_pos_nt-24, 0):])

            print("Mutant genomic context:    {}".format(nt_context))
            print("Reference peptide: {}".format(ref_peptide))
            print("Mutant peptide:    {}".format(mut_peptide))

            genomic_context_lower_pos = record['meta'].POS - 39
            genomic_context_upper_pos = record['meta'].POS + 40


#            if check_mut_pep_in_nt_context(mut_peptide, mut_nt):
            md5able_str = "{}".format(':'.join([str(record['meta'].CHROM), str(record['meta'].POS), str(record['transcript']), str(record['meta'].REF), str(record['meta'].ALT)])).encode('utf-8')
            var_md5 = hashlib.md5(md5able_str).hexdigest()[:16]
            mutant_peptides["MD5:{} VARIANT_POS:{}_{} TRANSCRIPT:{} REF:{} ALT:{} INDEL_TYPE:{} PROTEIN_CONTEXT:{} GENOMIC_CONTEXT:{} GENOMIC_CONTEXT_RANGE:{}_{}-{}".format(var_md5, record['meta'].CHROM, record['meta'].POS, record['transcript'], record['meta'].REF, record['meta'].ALT, 'frameshift_deletion', 'NA', 'NA', record['meta'].CHROM, genomic_context_lower_pos, genomic_context_upper_pos)] = mut_peptide

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



    if new_header_and_seqs:
        with open(args.output, 'w') as ofo:
            for k, v in new_header_and_seqs.items():
                ofo.write(">{}\n{}\n".format(k, v))







def pileup_truncated(bam,contig, start=None, stop=None):
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
    filtered_records, somatic_transcripts = filter_vcf_by_expression(args, expressed_txids)
    print("# of filtered records: {}".format(len(filtered_records)))
    write_expressed_vcf(args, filtered_records)
    with open(args.somatic_txs, 'w') as ofo:
       for transcript in somatic_transcripts:
           ofo.write("{}\n".format(transcript))


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
#                    expressed_txids.append(line[txid_column].split('.')[0])
                    expressed_txids.append(line[txid_column])
    return expressed_txids


def filter_vcf_by_expression(args, expressed_txids):
    """
    """
    somatic_transcripts = []
    filtered_records = []
    vcf_reader = ''
    input_vcf = args.vcf
    if input_vcf.endswith('gz'):
        vcf_reader = vcf.Reader(open(input_vcf), 'r', compressed=True)
    else:
        vcf_reader = vcf.Reader(open(input_vcf), 'r')
    for record in vcf_reader:
        annotations = [x for x in record.INFO['ANN']]
        for annotation in annotations:
            effects = annotation.split('|')
#            transcript = effects[6].partition('.')[0]
            transcript = effects[6]
            #print("Transcript: {}".format(transcript))
            if transcript in expressed_txids and record not in filtered_records:
                filtered_records.append(record)
                somatic_transcripts.append(transcript)
       
    return filtered_records, somatic_transcripts


def write_expressed_vcf(args, filtered_records):
    """
    """
    input_vcf = args.vcf
    vcf_reader = ''
    if input_vcf.endswith('gz'):
        vcf_reader = vcf.Reader(open(input_vcf), 'r', compressed=True)
    else:
        vcf_reader = vcf.Reader(open(input_vcf), 'r')
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


#Remove
#def write_crowded_vcf(args, filtered_records):
#    """
#    refactor
#    """
#    vcf_reader = vcf.Reader(filename=args.somatic_vcf, compressed=True)
#    vcf_writer = vcf.Writer(open(args.crowded_vars_output, 'w'), vcf_reader)
#    for filtered_record in filtered_records:
#        vcf_writer.write_record(filtered_record)


def isolated_variants(args):
    """
    """
    cand_vars = get_candidate_variants(args)
    print("# of candidate variants: {}".format(len(cand_vars)))
    isolated_vars = get_all_isolated_variants(args)
    print("# of isolated variants: {}".format(len(isolated_vars)))
    isolated_cand_vars = []
    crowded_cand_vars = []
    for record in cand_vars:
        if record in isolated_vars:
            isolated_cand_vars.append(record)
        else:
            crowded_cand_vars.append(record)
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
            line = line.rstrip('\n').split()
            line = [x for x in line if x]
            if len(line) > 14 and line[0] not in ['Pos', 'Protein']:
                k = "{}_{}_{}_{}".format(line[1], line[10], line[0], len(line[9]))
                wt_nms[k] = line[15]
                wt_peptides[k] = line[9]

    mts_w_agreto = []

    header = ''

    with open(args.mt_fasta) as mto:
        for line in mto.readlines():
            line = line.rstrip('\n').split()
            line = [x for x in line if x]
            if len(line) > 14 and line[0] == 'Pos':
                line.extend(['wildtype_binding_affinity', 'wildtype_peptide', 'agretopicity'])
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


def filter_mutant_peptides(args):
    """
    """
    threshold = int(args.max_peptide_length)
    print(threshold)
    mt_hits = {}
    wt_hits = {}
    with open(args.mt_netmhcpan) as fo:
        header = ''
        with open(args.mt_output, 'w') as ofo:
            for line_idx, line in enumerate(fo.readlines()):
                line = line.rstrip().split()
                if len(line) > 0 and line[0] == 'Pos' and not header :
                    header = line
                    print(header)
                    header.remove('BindLevel')
                    ofo.write("{}\n".format('\t'.join(header)))
                elif len(line) >= 16 and line[0] not in ['Protein', 'Pos']:
                    if 'SB' in line:
                        line.remove('<=')
                        line.remove('SB')
                    if 'WB' in line:
                        line.remove('<=')
                        line.remove('WB')
                    if int(line[0]) < threshold:
                        if (int(line[0]) + len(line[2])) > threshold:
                            ofo.write("{}\n".format('\t'.join(line)))
                    elif int(line[0]) == threshold: 
                        ofo.write("{}\n".format('\t'.join(line)))
    
    with open(args.wt_netmhcpan) as fo:
        header = ''
        with open(args.wt_output, 'w') as ofo:
            for line_idx, line in enumerate(fo.readlines()):
                line = line.rstrip().split()
                if len(line) > 0 and line[0] == 'Pos' and not header :
                    header = line
                    header.remove('BindLevel')
                    ofo.write("{}\n".format('\t'.join(header)))
                elif len(line) >= 16 and line[0] not in ['Protein', 'Pos']:
                    if 'SB' in line:
                        line.remove('<=')
                        line.remove('SB')
                    if 'WB' in line:
                        line.remove('<=')
                        line.remove('WB')
                    if int(line[0]) < threshold:
                        if (int(line[0]) + len(line[2])) > threshold:
                            ofo.write("{}\n".format('\t'.join(line)))
                    elif int(line[0]) == threshold: 
                        ofo.write("{}\n".format('\t'.join(line)))
    
              

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

    uq = np.percentile(list(tx_to_log2tpm.values()), 75)

    for k,v in tx_to_log2tpm.items():
        tx_to_uqlog2tpm[k] = v/uq

    tx_to_gene = {}
    with open(args.gtf) as fo:
        for line in fo.readlines():
            line = line.split('\t')
            if len(line) > 3 and line[2] == 'transcript' and line[1] != 'geve':
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
    header_extension = ['gene_name', 'transcript_identifier',
                        'variant_position', 'reference_allele',
                        'alternate_allele', 'tpm', 'log2(tpm+1)',
                        'uq(log2(tpm+1))', 'cancer_cell_fraction',
                        'protein_context', 'nucleotide_context',
                        'variant_type']

    new_lines = []
    with open(args.binding_affinities) as fo:
        for line_idx, line in enumerate(fo.readlines()):
            line = line.rstrip().split('\t')
            if line_idx == 0:
                header = line
                try:
                    header.remove('BindLevel')
                except:
                    pass
                print(header)
                header.insert(0, 'antigen_source')
                header.extend(header_extension)
                header = header[:5] + header[11:]
                header[1] = 'pos'                                                                   
                header[2] = 'mhc_allele'                                                            
                header[3] = 'peptide'                                                               
                header[4] = 'peptide_core'                                                          
                header[5] = 'internal_identifier'                                                   
                header[6] = 'score_el'                                                              
                header[7] = 'percent_rank_el'                                                       
                header[8] = 'score_ba'                                                              
                header[9] = 'percent_rank_ba'                                                       
                header[10] = 'binding_affinity'  
                header[11] = 'reference_binding_affinity'  
                header[12] = 'reference_peptide'  
                header[13] = 'agretopicity'  
                header[14] = 'total_rna_reads_full_overlap'  
                header[15] = 'reads_with_peptide'  
                header[16] = 'proportion_rna_reads_with_peptide'  
            elif len(line) < 16:
                continue
            else:
                csum = line[10]
                if csum not in checksum_to_meta_map.keys():
                    continue
                tx_id = checksum_to_meta_map[csum]['TRANSCRIPT']
                if float(line[15]) < 1000 and tx_id.split('.')[0] in tx_to_gene.keys() and tx_id.split('.')[0] in tx_to_tpm.keys() and line[10] in checksum_to_meta_map.keys():
                    line.insert(0, 'SNV')
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
#                    rna_coverage = 'NA'
#                    proportion_rna_var = 'NA'
#                    rna_coverage = checksum_to_meta_map[csum]['TOTAL_RNA_COVERAGE']
#                    proportion_rna_var = checksum_to_meta_map[csum]['PROPORTION_RNA_VARIANT_READS']
                    ccf = 'NA'
                    if var_pos in var_to_ccf.keys():
                        ccf = var_to_ccf[var_pos]
                    if 'SB' in line:
                        line.remove('<=')
                        line.remove('SB')
                    if 'WB' in line:
                       line.remove('<=')
                       line.remove('WB')
                    line.extend([gene_name, tx_id, var_pos, ref, alt, str(tpm), str(log2tpm), str(uqlog2tpm), ccf, aa_context, nt_context, var_type])
                    line = line[:5] + line[11:]
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

    uq = np.percentile(list(tx_to_log2tpm.values()), 75)

    for k,v in tx_to_log2tpm.items():
        tx_to_uqlog2tpm[k] = v/uq

    tx_to_gene = {}
    with open(args.gtf) as gtfo:
        for line in gtfo.readlines():
            line = line.split('\t')
            if len(line) > 3 and line[2] == 'transcript' and line[1] != 'geve':
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
    header_extension = ['gene_name', 'transcript_identifier',
                        'variant_position', 'reference_allele',
                        'alternate_allele', 'tpm', 'log2(tpm+1)',
                        'uq(log2(tpm+1))', 'cancer_cell_fraction',
                        'protein_context', 'nucleotide_context',
                        'variant_type']

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
                header.insert(0, 'antigen_source')
                try:
                    header.remove('BindLevel')
                except:
                    pass
                header.extend(header_extension)
                header = header[:5] + header[11:]
                header[1] = 'pos'                                                                   
                header[2] = 'mhc_allele'                                                            
                header[3] = 'peptide'                                                               
                header[4] = 'peptide_core'                                                          
                header[5] = 'internal_identifier'                                                   
                header[6] = 'score_el'                                                              
                header[7] = 'percent_rank_el'                                                       
                header[8] = 'score_ba'                                                              
                header[9] = 'percent_rank_ba'                                                       
                header[10] = 'binding_affinity'  
            elif len(line) < 16:
                continue
            if line[10] not in checksum_to_meta_map.keys():
                continue
            csum = line[10]
            tx_id = checksum_to_meta_map[csum]['TRANSCRIPT']
            if tx_id.split('.')[0] in tx_to_gene.keys() and tx_id.split('.')[0] in tx_to_tpm.keys() and float(line[15]) < 1000:
            #if tx_id.split(',')[0] in tx_to_gene.keys() and tx_id.split('.')[0] in tx_to_tpm.keys():
                line.insert(0, 'InDel')
                line = line[:5] + line[11:]
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
                #rna_coverage = 'NA'
                #rna_proportion_vars = 'NA'
                #rna_coverage = checksum_to_meta_map[csum]['TOTAL_RNA_COVERAGE']
                #rna_proportion_vars = checksum_to_meta_map[csum]['PROPORTION_RNA_VARIANT_READS']
                if var_pos in var_to_ccf.keys():
                    ccf = var_to_ccf[var_pos]
                if 'SB' in line:
                    line.remove('<=')
                    line.remove('SB')
                if 'WB' in line:
                    line.remove('<=')
                    line.remove('WB')
                line.extend([gene_name, tx_id, var_pos, ref, alt, str(tpm), str(log2tpm), str(uqlog2tpm), ccf, aa_context, nt_context, var_type])
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
#                if len(line[3]) > 0 and len(line[4]) > 0:
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
#        print([re.search('-ad-', i.sample) for i in record.samples])
#        print([i.sample for i in record.samples if re.search('-ad-', i.sample)])
        call = record.genotype([i.sample for i in record.samples if re.search('-ad-', i.sample)][0])
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
            if maj_count != 'NA' and min_count != 'NA':
                print("{}".format(','.join([chr, start, stop, maj_count, min_count])))
                if chr not in copy_no.keys():
                    copy_no[chr] = {}
                if "{}-{}".format(start, stop) not in copy_no[chr].keys():
                    copy_no[chr]["{}-{}".format(start, stop)] = {}
                copy_no[chr]["{}-{}".format(start, stop)]["maj_count"] = maj_count
                copy_no[chr]["{}-{}".format(start, stop)]["min_count"] = min_count

    #This is pointless code, remove it during refactoring.
    easily_parsable = []

    for var in vars.keys():
        print(var)
        chr, pos = var.split('_')
#        if var in depths.keys():
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
                if var in depths.keys():
                    ref_depth = str(depths[var]['ref_depth'])
                    alt_depth = str(depths[var]['alt_depth'])
                    ref_cn = str(copy_no[chr][segment]["maj_count"])
                    alt_cn = str(copy_no[chr][segment]["min_count"])
                    var = var.replace('_', ':')
                    pyclone_inp.append("{}\n".format('\t'.join([var, args.samp_id, ref_depth, alt_depth, ref_cn, alt_cn, '2', '0.001', cellularity])))
    for var in easily_parsable:
        var = "{}:{}".format(var[0], var[1])
        if var not in [i.split('\t')[0] for i in pyclone_inp]:
            pyclone_inp.append("{}\n".format('\t'.join([var, args.samp_id, '0', '0', '2', '0', '2', '0.001', cellularity])))

    with open(args.output, 'a') as ofo:
        for i in pyclone_inp:
            ofo.write(i)


def expressed_hervs(args):
    """
    """
    tpms = {}
    quant_map = {}
    all_quant_meta = {}
    with open(args.quants) as fo:
        for line_idx, line in enumerate(fo.readlines()):
            line = line.rstrip().split('\t')
            if line_idx == 0:
                for col_idx, col in enumerate(line):
                    quant_map[col] = col_idx
            else:
                name = line[quant_map['Name']]
                if re.search('[Mmus|Hsap].*[\+\-]$', name):
                    tpms[name] = line[quant_map['TPM']]
                    all_quant_meta[name] = ','.join(line[1:])

    tx_abundances = [float(x) for x in tpms.values()]
    print(tx_abundances)
 
    if not(args.abundance_threshold):
        print("No abundance threshold!")
        tx_threshold = np.mean(tx_abundances) + (5*np.std(tx_abundances))
        print("Transcript threshold (mean + 5*std): {}".format(tx_threshold))
    else:
        print("Abundance threshold: {}".format(args.abundance_threshold))
        tx_threshold = float(args.abundance_threshold)
    print("tx_thresshold: {}".format(tx_threshold))
#    print("# of expressed transcripts: {}".format(len(expressed_hervs)))
#    print("Some expressed transcripts: {}".format(expressed_hervs[:10]))

    with open(args.output, 'w') as ofo:
        ofo.write('Name,Tumor_CPM,Norm_CPM,log2(Tumor_CPM+1)-log2(Norm_CPM+1),Length,EffectiveLength,TPM,NumReads\n')
        for herv, tpm in tpms.items():
            if float(tpm) > float(tx_threshold):
                if (args.trim_chr_prefix):
                    herv_adj = herv.replace('chr', '')
                ofo.write("{},NA,NA,NA,{}\n".format(herv_adj, all_quant_meta[herv]))

def get_expressed_selfs_bed(args):
    """
    """
    expressed_self_transcripts = []
    with open(args.expressed_selfs) as fo:
        with open(args.output, 'w') as ofo:
            for line in fo.readlines():
                expressed_self_transcripts.append(line.rstrip().split(':')[1])

    with open(args.output, 'w') as ofo:
        with open(args.gff) as fo:
            for line in fo.readlines():
                for expressed_self_transcript in expressed_self_transcripts:
                    if re.search(expressed_self_transcript, line):
#                        for segment in ['CDS', 'five_prime_UTR', 'three_prime_UTR']:
                        if re.search('\tCDS\t', line):
                            chr = line.split('\t')[0]
                            start = line.split('\t')[3]
                            stop = line.split('\t')[4]
                            ofo.write("{}\t{}\t{}\n".format(chr, start, stop))

def get_expressed_transcripts_bed(args):
    """
    """
    expressed_transcripts = []
    with open(args.expressed_transcripts) as fo:
        with open(args.output, 'w') as ofo:
            for line in fo.readlines():
                expressed_transcripts.append(line.rstrip())

    with open(args.output, 'w') as ofo:
        with open(args.gff) as fo:
            for line in fo.readlines():
                line = line.split('\t')
                if len(line) > 1 and line[2] == 'CDS':
                    for expressed_transcript in expressed_transcripts:
                        if re.search(expressed_transcript, line[8]):
                            print(line)
#                        if re.search('\tCDS\t', line):
#                            chr = line.split('\t')[0]
#                            start = line.split('\t')[3]
#                            stop = line.split('\t')[4]
                            chr = line[0]
                            start = line[3]
                            stop = line[4]
                            ofo.write("{}\t{}\t{}\n".format(chr, start, stop))

def get_expressed_ervs_bed(args):
    """
    """
    full_to_post_met = {}

    with open(args.geve_reference) as fo:
        for line_idx, line in enumerate(fo.readlines()):
            if line_idx != 0:
                line = line.split('\t')
                if line[8] != '-':
                    full_to_post_met[line[0]] = line[8].replace('.M', '')



    print(full_to_post_met)
    
    with open(args.expressed_ervs) as fo:
        with open(args.output, 'w') as ofo:
            for line_idx, line in enumerate(fo.readlines()):
                if line_idx != 0:
                    print("line: {}".format(line))
                    if line.split(',')[0] in full_to_post_met.keys():
                        post_met = full_to_post_met[line.split(',')[0]]
                        print("post_met: {}".format(post_met))
                        herv = post_met.split(',')[0]
                        chr = post_met.split('.')[1] 
                        first = post_met.split('.')[2] 
                        second = post_met.split('.')[3] 
#                    herv = line.split(',')[0]
#                    chr = herv.split('.')[1] 
#                    first = herv.split('.')[2] 
#                    second = herv.split('.')[3] 
                        ofo.write("{}\t{}\t{}\n".format(chr, first, second))

def get_expressed_viral_bed(args):
    """
    """

    viral_stops = {}
 
    for seq_record in SeqIO.parse(args.viral_cds_ref, "fasta"):
        viral_stops[seq_record.description.split()[0]] = len(seq_record.seq)

    print(viral_stops)


    with open(args.expressed_viruses) as fo:
        with open(args.output, 'w') as ofo:
            for line_idx, line in enumerate(fo.readlines()):
                ofo.write("{}\t{}\t{}\n".format(line.split()[0], 1, viral_stops[line.split()[0]]))


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


def make_erv_peptides(args):
    """
    """
    post_met_to_full = {}

    with open(args.geve_reference) as fo:
        for line_idx, line in enumerate(fo.readlines()):
            if line_idx != 0:
                line = line.split('\t')
                post_met = line[8].replace('.M', '') 
                print(post_met)
                if len(post_met.split('.')) > 1 :
                    chr = post_met.split('.')[1] 
                    first = post_met.split('.')[2] 
                    second = post_met.split('.')[3] 
                    post_met_to_full["{}:{}-{}".format(chr, first, second)] = line[0]

    expressed_hervs = []
    expressed_hervs_metadata = {}


    expressed_herv_col_map = {}

    with open(args.expressed_ervs) as fo:
        for line_idx, line in enumerate(fo.readlines()):
            line = line.rstrip().split(',')
            if line_idx == 0: 
                for idx, i in enumerate(line): 
                    expressed_herv_col_map[i] = idx
                print(expressed_herv_col_map)
            else:
                herv_id = line[expressed_herv_col_map['Name']]
                expressed_hervs.append(herv_id)
                expressed_hervs_metadata[herv_id] = {}
                expressed_hervs_metadata[herv_id]['tumor_cpm'] = line[expressed_herv_col_map['Tumor_CPM']]
                expressed_hervs_metadata[herv_id]['norm_cpm'] = line[expressed_herv_col_map['Norm_CPM']]
                expressed_hervs_metadata[herv_id]['delta'] = line[expressed_herv_col_map['log2(Tumor_CPM+1)-log2(Norm_CPM+1)']]
                expressed_hervs_metadata[herv_id]['tpm'] = line[expressed_herv_col_map['TPM']]
                expressed_hervs_metadata[herv_id]['numreads'] = line[expressed_herv_col_map['NumReads']]

    print(expressed_hervs_metadata)

    expressed_hervs_seqs = {}

    for seq_record in SeqIO.parse(args.patient_ervs_fasta, "fasta"):
        print(seq_record)
        chr = seq_record.description.split(':')[0]
        start = seq_record.description.split(':')[1].split('-')[0]
        stop = seq_record.description.split(':')[1].split('-')[1]
        if args.species != 'mm':
            new_seq_record_description = 'Hsap38.{}.{}.{}'.format(chr, start, stop)
        else:
            new_seq_record_description = 'Mmus38.{}.{}.{}'.format(chr, start, stop)
        print(new_seq_record_description)
        regexed_herv_id = post_met_to_full["{}:{}-{}".format(chr, start, stop)]
#        regexed_herv_id = [x for x in expressed_hervs_metadata.keys() if re.search(new_seq_record_description, x)][0]
#        regexed_herv_id = [x for x in expressed_hervs_metadata.keys() if re.search(new_seq_record_description, x)][0]
        #regexed_herv_id = [re.findall(new_seq_record_description + '+', x)[0] for x in expressed_hervs_metadata.keys() if re.findall(new_seq_record_description + '+', x)]
        print(regexed_herv_id)
        expressed_hervs_seqs[regexed_herv_id] = seq_record.seq

    print(expressed_hervs_seqs)

    expressed_hervs_aas = {}

    for id, seq in expressed_hervs_seqs.items():
        if id.endswith('-'):
            aa_seq = seq.reverse_complement().translate(to_stop=True)
            expressed_hervs_seqs[id] = seq.reverse_complement()
        else:
            aa_seq = seq.translate(to_stop=True)
        expressed_hervs_aas["{}".format(id)] = aa_seq

    print(expressed_hervs_aas)
    

    with open(args.output, 'w') as ofo:
        for k, v in expressed_hervs_aas.items():
#            header = "MD5:{} NAME:{} TUMOR_CPM:{} NORM_CPM:{} DELTA:{} TPM:{} NUMREADS:{} RNA_COVERAGE_RANGE:{} POTENTIAL_GERMLINE_HETS:{}".format(hashlib.md5(str(k).encode('utf-8')).hexdigest()[:16], k, expressed_hervs_metadata[k]['tumor_cpm'], expressed_hervs_metadata[k]['norm_cpm'], expressed_hervs_metadata[k]['delta'], expressed_hervs_metadata[k]['tpm'], expressed_hervs_metadata[k]['numreads'], orf_rna_coverage[k], ','.join([str(x) for x in orf_het_sites[k]]))
            header = "MD5:{} NAME:{} TUMOR_CPM:{} NORM_CPM:{} DELTA:{} TPM:{} NUMREADS:{}".format(hashlib.md5(str(k).encode('utf-8')).hexdigest()[:16], k, expressed_hervs_metadata[k]['tumor_cpm'], expressed_hervs_metadata[k]['norm_cpm'], expressed_hervs_metadata[k]['delta'], expressed_hervs_metadata[k]['tpm'], expressed_hervs_metadata[k]['numreads'])
            ofo.write(">{}\n{}\n".format(header, v))

    with open(args.nt_output, 'w') as ofo:
        for k, v in expressed_hervs_seqs.items():
            header = "MD5:{} NAME:{}".format(hashlib.md5(str(k).encode('utf-8')).hexdigest()[:16], k)
            ofo.write(">{}\n{}\n".format(header, v))#''.join(orf_patient_seq[k])))


#def make_herv_peptides(args):
#    """
#    """
#    expressed_hervs = []
#    expressed_hervs_metadata = {}
#    expressed_hervs_aas = {}
#
#    expressed_herv_col_map = {}
#
#    with open(args.expressed_hervs) as fo:
#        for line_idx, line in enumerate(fo.readlines()):
#            line = line.rstrip().split(',')
#            if line_idx == 0: 
#                for idx, i in enumerate(line): 
#                    expressed_herv_col_map[i] = idx
#                print(expressed_herv_col_map)
#            else:
#                herv_id = line[expressed_herv_col_map['Name']]
#                expressed_hervs.append(herv_id)
#                expressed_hervs_metadata[herv_id] = {}
#                expressed_hervs_metadata[herv_id]['tumor_cpm'] = line[expressed_herv_col_map['Tumor_CPM']]
#                expressed_hervs_metadata[herv_id]['norm_cpm'] = line[expressed_herv_col_map['Norm_CPM']]
#                expressed_hervs_metadata[herv_id]['delta'] = line[expressed_herv_col_map['log2(Tumor_CPM+1)-log2(Norm_CPM+1)']]
#                expressed_hervs_metadata[herv_id]['tpm'] = line[expressed_herv_col_map['TPM']]
#                expressed_hervs_metadata[herv_id]['numreads'] = line[expressed_herv_col_map['NumReads']]
#
#    expressed_hervs_seqs = {}
#
#
#    for seq_record in SeqIO.parse(args.herv_ref, "fasta"):
#        if seq_record.description in expressed_hervs: 
#            expressed_hervs_seqs[seq_record.description] = seq_record.seq
#
#    alns = pysam.AlignmentFile(args.tumor_bam, 'rb')
#
#    #This is a bit crude right now. It _only_ considers the major allele at
#    #each position for the nucleotide sequence. I will include something in the
#    #header to list potentially heterozygous sites. This can be translated within
#    #the herv_metadata function.
#
#    orf_rna_coverage = {}
#    orf_het_sites = {}
#    orf_patient_seq = {}
#
#
#    for herv in expressed_hervs_metadata.keys():
#        orf_patient_seq[herv] = []
#        orf_het_sites[herv] = []
#         
#        print(herv)
#        print(expressed_hervs_metadata[herv])
#        herv_chr = herv.split('.')[1] 
#        herv_start = int(herv.split('.')[2])
#        herv_stop = int(herv.split('.')[3])
#        print("{}\t{}\t{}".format(herv_chr, herv_start, herv_stop))
#        herv_seq = expressed_hervs_seqs[herv]
#        if herv[-1] == '+':
#            tmp_info = pileup_truncated(alns, herv_chr, herv_start-1, herv_stop)
#            pileup_info = [(x.reference_pos, [i for i in x.get_query_sequences() if i]) for x in tmp_info]
#        elif herv[-1] == '-':
#            tmp_info = pileup_truncated(alns, herv_chr, herv_start-1, herv_stop)
#            pileup_info = reversed([(x.reference_pos, [complement(i) for i in x.get_query_sequences() if i]) for x in tmp_info])
#            
#        else:
#            sys.exit("Cannot determine ERV orientation.")
#
#        covered_pos = [] 
#        for pos_idx, pos in enumerate(pileup_info):
#            
#            query_seqs = [x.lower() for x in pos[1]]
#            if not query_seqs:
#                orf_patient_seq[herv].append(herv_seq[pos_idx])
#            if query_seqs:
#                covered_pos.append(pos[0])
#                counts= {x:query_seqs.count(x) for x in query_seqs}
#                most_freq_finder = lambda x: scipy.stats.mode(x)[0][0]
#                most_freq_base = most_freq_finder(query_seqs)
#                orf_patient_seq[herv].append(most_freq_base)
#                mismatch = ''
#                if herv_seq[pos_idx] != most_freq_base:
#                    mismatch = "Mismatch"
#                print("{}\t{}\t{}\t{}".format(pos[0], herv_seq[pos_idx], most_freq_base, mismatch))
#                most_freq_base_freq = float(query_seqs.count(most_freq_base))/len(query_seqs)
#                if float(most_freq_base_freq) < 0.80:
#                    orf_het_sites[herv].append(pos[0])
#                    print("WARNING: Possible germline mutation!")
#                    print("{}\t{}\t{}\t{}\t{}".format(pos[0], herv_seq[pos_idx], most_freq_base, query_seqs.count(most_freq_base)/len(query_seqs), counts))
#
#
#        orf_rna_coverage[herv] = "{},{}".format(covered_pos[0]+1, covered_pos[-1]+1)
#        print(herv_seq)
#        print(''.join(orf_patient_seq[herv]))
#        print(orf_het_sites[herv])
#        print(orf_rna_coverage[herv])
#        
#
#    for id, seq in expressed_hervs_seqs.items():
#        aa_seq = Seq(''.join(orf_patient_seq[id])).translate(to_stop=True)
#        expressed_hervs_aas["{}".format(id)] = aa_seq
#
#    with open(args.output, 'w') as ofo:
#        for k, v in expressed_hervs_aas.items():
#            header = "MD5:{} NAME:{} TUMOR_CPM:{} NORM_CPM:{} DELTA:{} TPM:{} NUMREADS:{} RNA_COVERAGE_RANGE:{} POTENTIAL_GERMLINE_HETS:{}".format(hashlib.md5(str(k).encode('utf-8')).hexdigest()[:16], k, expressed_hervs_metadata[k]['tumor_cpm'], expressed_hervs_metadata[k]['norm_cpm'], expressed_hervs_metadata[k]['delta'], expressed_hervs_metadata[k]['tpm'], expressed_hervs_metadata[k]['numreads'], orf_rna_coverage[k], ','.join([str(x) for x in orf_het_sites[k]]))
#            ofo.write(">{}\n{}\n".format(header, v))
#
#    with open(args.nt_output, 'w') as ofo:
#        for k, v in expressed_hervs_aas.items():
#            header = "MD5:{}".format(hashlib.md5(str(k).encode('utf-8')).hexdigest()[:16])
#            ofo.write(">{}\n{}\n".format(header, ''.join(orf_patient_seq[k])))

def complement(i):
    comp = {'a': 't',
            't': 'a',
            'g': 'c',
            'c': 'g'}
    return comp[i.lower()]


def add_erv_metadata(args):
    """
    """

    het_sites = {}

    vcf_reader = ''
    vcf_reader = vcf.Reader(compressed=True, filename=args.patient_vcf)

    for record in vcf_reader:
        hets = record.get_hets()
        for het_call in hets:
            if het_call.site.CHROM in het_sites.keys():
                het_sites[het_call.site.CHROM].append(het_call.site.POS)
            else:
                het_sites[het_call.site.CHROM] = [het_call.site.POS]

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
                if line[tx_col_idx].startswith('Hsap'):
                    tx_to_tpm[line[tx_col_idx]] = float(line[tpm_col_idx])                
                if line[tx_col_idx].startswith('Mmus'):
                    tx_to_tpm[line[tx_col_idx].replace('chr', '')] = float(line[tpm_col_idx])                
                else:                                             
                    tx_to_tpm[line[tx_col_idx].split('.')[0]] = float(line[tpm_col_idx])                
                                                                                                    
    for k,v in tx_to_tpm.items():                                                                   
        tx_to_log2tpm[k] = np.log2(v + 1)                                                           
                                                                                                    
    uq = np.percentile(list(tx_to_log2tpm.values()), 75)                                            
                                                                                                    
    for k,v in tx_to_log2tpm.items():                                                               
        tx_to_uqlog2tpm[k] = v/uq   


    checksum_to_meta_map = {}
    with open(args.peptides) as fo:
        for line in fo.readlines():
            if line.startswith('>'):
                line = line.rstrip()
                bufr_dict = {i.split(':')[0]: i.split(':')[1] for i in line.split(' ')[1:]}
                line = line.split(' ')
                checksum = "MD5_{}".format(line[0].lstrip('>').split(':')[1][:-5])
                checksum_to_meta_map[checksum] = bufr_dict


    hsap_external_data = {}
    hsap_idx_to_col = {}
    with open(args.geve_data) as fo:
        for line_idx, line in enumerate(fo.readlines()):
            if line_idx == 0:
                line = line.rstrip().split('\t')
                #print(line)
                for col_idx, col in enumerate(line):
                    hsap_idx_to_col[col_idx] = col
                #print(hsap_idx_to_col)
            else:
                line = line.rstrip().split('\t')
                #print(line)
                id = line[0]
                hsap_external_data[id] = {}
                for val_idx, val in enumerate(line):
                    #print("{}\t{}".format(val_idx, val))
                    if val_idx > 0:
                        #if re.search('chr', id):
                        #    id = id.replace('chr', '')
                        hsap_external_data[id][hsap_idx_to_col[val_idx]] = val
                print(hsap_external_data[id])

    ### Determining gEVE neighbors
    hsap_neighbors = {}
    hsap_raw = []
    for expressed_erv in checksum_to_meta_map.keys():
        hsap_raw.append(checksum_to_meta_map[expressed_erv]['NAME'])

    hsap_deplexed = {}

    hsap_chrs = [x.split('.')[1] for x in hsap_raw]

    for chr in hsap_chrs:
        hsap_deplexed[chr] = {}
        hsap_deplexed[chr]['+'] = []
        hsap_deplexed[chr]['-'] = []

    ###Populating the hsap_deplexed dict
    for hsap in hsap_raw:
        chr = hsap.split('.')[1]
        start = hsap.split('.')[2]
        stop = hsap.split('.')[3]
        strand= hsap.split('.')[4]

        hsap_deplexed[chr][strand].append([start,stop])

    for hsap in hsap_raw:
        chr = hsap.split('.')[1]
        start = hsap.split('.')[2]
        stop = hsap.split('.')[3]
        strand= hsap.split('.')[4]
        midpoint = (int(start) + int(stop))/2

        candidate_neighbors = hsap_deplexed[chr][strand]

        if candidate_neighbors:
            midpoint_diffs = []
            for candidate in candidate_neighbors:
                if start != candidate[0] and stop != candidate[1]:
                    cand_midpoint = (int(candidate[0]) + int(candidate[1]))/2
                    mid_to_cand_mid_dist = abs(midpoint - cand_midpoint)
                    if mid_to_cand_mid_dist < 10000:
                        midpoint_diffs.append(mid_to_cand_mid_dist)
                if midpoint_diffs:
                    desired_midpoint = min(midpoint_diffs)
                 
                    for candidate in candidate_neighbors:
                        cand_midpoint = (int(candidate[0]) + int(candidate[1]))/2
                        if cand_midpoint == desired_midpoint:
                            hsap_neighbors[hsap] = {}
                            hsap_neighbors[hsap]['neighbor'] = "Hsap38.{}.{}.{}.{}".format(chr, candidate[0], candidate[1], strand)
                            hsap_neighbors[hsap]['distance'] = abs(midpoint - cand_midpoint)


    pat_nts = {}
    pat_aas = {}

    pat_seqs = {}
    for seq_record in SeqIO.parse(args.peptides, "fasta"):
        print(seq_record.id)
        orf = seq_record.id.replace(':', '_')[:15]
        pat_aas[orf] = seq_record.seq

    for seq_record in SeqIO.parse(args.nt, "fasta"):
        geve = seq_record.id.replace(':', '_')[:15]
        pat_nts[geve] = seq_record.seq
    


    header = ''

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
                header.insert(0, 'antigen_source')
                header = header[:5] + header[11:]
                header[1] = 'pos'
                header[2] = 'mhc_allele'
                header[3] = 'peptide'
                header[4] = 'peptide_core'
                header[5] = 'internal_identifier'
                header[6] = 'score_el'
                header[7] = 'percent_rank_el'
                header[8] = 'score_ba'
                header[9] = 'percent_rank_ba'
                header[10] = 'binding_affinity'
                header.extend(['geve_orf'])
                header.extend(['tumor_cpm'])
                header.extend(['norm_cpm'])
                header.extend(['log2(tumor_cpm+1)-log2(norm_cpm+1)'])
                header.extend(['tpm'])
                header.extend(['log2(tpm+1)'])
                header.extend(['uq(log2(tpm+1))'])
                header.extend(['geve_viral_blast'])
                header.extend(['geve_retrotector'])
                header.extend(['geve_expressed_neighbor_orf'])
                header.extend(['geve_expressed_neighbor_orf_distance'])
                header.extend(['geve_peptide_downstream_of_start_codon'])
                header.extend(['geve_peptide_potential_germ_het'])
                header.extend(['geve_orf_protein_length'])
                header.extend(['geve_met_orf_protein_length'])
                header.extend(['nucleotide_context'])
                header.extend(['protein_context'])

            elif len(line) < 16 or line[0] in ['Protein', 'Pos']:
                continue
            elif float(line[15]) > 0:
                print(line) 
                line.insert(0, 'ERV')
                line = line[:5] + line[11:]
                if 'SB' in line:
                    line.remove('<=')
                    line.remove('SB')
                if 'WB' in line:
                    line.remove('<=')
                    line.remove('WB')
                tumor_cpm = checksum_to_meta_map[line[5]]['TUMOR_CPM']
                norm_cpm = checksum_to_meta_map[line[5]]['NORM_CPM']
                delta = checksum_to_meta_map[line[5]]['DELTA']
                tpm = checksum_to_meta_map[line[5]]['TPM']
                name = checksum_to_meta_map[line[5]]['NAME']
                log2tpm = tx_to_log2tpm[name]
                uqlog2tpm = tx_to_uqlog2tpm[name]
                #print(hsap_external_data[name])
                viral_blast = hsap_external_data[name]['Viral_BLAST']
                retrotector = hsap_external_data[name]['RetroTector']
                line.extend([name])
                line.extend([tumor_cpm])
                line.extend([norm_cpm])
                line.extend([delta])
                line.extend([tpm])
                line.extend([log2tpm])
                line.extend([uqlog2tpm])
                if viral_blast != '-':
                    line.extend([viral_blast])
                else:
                    line.extend(['NA'])
                if retrotector != '.':
                    line.extend([retrotector])
                else:
                    line.extend(['NA'])
                # Neighbor data
                if name not in hsap_neighbors.keys():
                    line.extend(['NA', 'NA'])
                else:
                    line.extend([hsap_neighbors[name]['neighbor'], hsap_neighbors[name]['distance']])

                # Downstream from start codon
                print(name)
                print(hsap_external_data[name])
                print(hsap_external_data[name]['MetID'])
                if hsap_external_data[name]['strand'] == '+':
                    if hsap_external_data[name]['MetID'] != '-':
                        met_start = hsap_external_data[name]['MetID'].split('.')[2]
                        aa_start = hsap_external_data[name]['start']
                        start_codon_idx = int(met_start) - int(aa_start)
                        if int(start_codon_idx) > int(line[1]):
                            line.extend(['False'])
                        else:
                            line.extend(['True'])
                    else:
                        line.extend(['N/A'])

                    peptide_start = int(hsap_external_data[name]['start']) + int(start_codon_idx)
                    peptide_stop = peptide_start + len(line[3])
                   
                    has_het = 0
                    chr = name.split('.')[1]
                    if chr in het_sites.keys():
                        for het in het_sites[chr]:
                            if het in range(peptide_start, peptide_stop):
                                has_het = 1

                    if has_het:
                        line.extend(['True'])
                    else:
                        line.extend(['False'])
                    
          
                elif hsap_external_data[name]['strand'] == '-':
                    if hsap_external_data[name]['MetID'] != '-':
                        met_start = hsap_external_data[name]['MetID'].split('.')[3]
                        aa_start = hsap_external_data[name]['end']
                        start_codon_idx = int(aa_start) - int(met_start)
                        if int(start_codon_idx) > int(line[1]):
                            line.extend(['False'])
                        else:
                            line.extend(['True'])
                    else:
                        line.extend(['N/A'])
                    
                    peptide_start = int(hsap_external_data[name]['end']) - int(start_codon_idx)
                    peptide_stop = peptide_start - len(line[3])
                   

#                    rna_coverage_start = int(checksum_to_meta_map[line[5]]['POTENTIAL_GERMLINE_HETS'].strip('(').strip(')').split(',')[1])
#                    rna_coverage_stop = int(checksum_to_meta_map[line[5]]['POTENTIAL_GERMLINE_HETS'].strip('(').strip(')').split(',')[0])

                    has_het = 0
                    chr = name.split('.')[1]
                    if chr in het_sites.keys():
                        for het in het_sites[chr]:
                            if het in range(peptide_stop, peptide_start):
                                has_het = 1

                    if has_het:
                        line.extend(['True'])
                    else:
                        line.extend(['False'])

                line.extend([hsap_external_data[name]['AA_length']])
                line.extend([hsap_external_data[name]['Met_AA_length']])

                aa_start_idx = str(pat_aas[line[5]]).index(line[3])
                aa_stop_idx = aa_start_idx + len(line[3])

                adjusted_start_idx = aa_start_idx - 8
                adjusted_stop_idx = aa_stop_idx + 8
                nt_start_idx = adjusted_start_idx * 3
                nt_stop_idx = adjusted_stop_idx * 3

                protein_context = str(pat_aas[line[5]][max(0, adjusted_start_idx):min(adjusted_stop_idx, len(pat_aas[line[5]]))])
                nuc_context = str(pat_nts[line[5]][max(0, nt_start_idx):min(nt_stop_idx, len(pat_nts[line[5]]))])
                line.extend([nuc_context])
                line.extend([protein_context])

                print(line)
                output_lines.append(line)

    print(output_lines)

    if header and output_lines:
        with open(args.output, 'w') as ofo:
            ofo.write("{}\n".format('\t'.join(header)))
            for line in output_lines:
                str_line = [str(x) for x in line]
                ofo.write("{}\n".format('\t'.join(str_line)))
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
            gene_list.append(line.rstrip().split('\t')[0])
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
            if len(line) > 3 and line[2] == 'transcript' and 'gene_name' in str(line) and (re.search('transcript_type "protein_coding"', '\t'.join(line)) or re.search('transcript_biotype "protein_coding"', '\t'.join(line))):
                gene_name = str(line).split('gene_name "')[1].split('"')[0]
                tx_id = str(line).split('transcript_id "')[1].split('"')[0].split('.')[0]
                print("{}\t{}".format(gene_name, tx_id))
                tx_to_gene[tx_id] = gene_name



    for tx, gene in tx_to_gene.items():
        if gene in gene_list and tx in tx in expressed_txids:
            print(gene, tx)
            expressed_selfs.append("{}:{}".format(gene, tx))

    if expressed_selfs:
        with open(args.output, 'w') as ofo:
            for expressed_self in expressed_selfs:
                ofo.write("{}\n".format(expressed_self))


def make_self_antigen_peptides(args):
    """
    """
    print(args)
    expressed_selfs_exon_seqs = {}
    for seq_record in SeqIO.parse(args.selfs_seqs_fasta, "fasta"):
        expressed_selfs_exon_seqs[seq_record.description] = seq_record.seq

    print(expressed_selfs_exon_seqs.keys())

    expressed_selfs = []
    with open(args.expressed_selfs) as fo:
        for line in fo.readlines():
            expressed_selfs.append(line.split(':')[1].rstrip('\n'))
  
    expressed_selfs_tx_metadata = {} 

    for expressed_self in expressed_selfs:
        print(expressed_self)
        expressed_selfs_tx_metadata[expressed_self] = {}
        expressed_selfs_tx_metadata[expressed_self]['cds'] = []
        expressed_selfs_tx_metadata[expressed_self]['UTR'] = []
        expressed_selfs_tx_metadata[expressed_self]['strand'] = ''


    with open(args.gtf) as fo:
        for line in fo.readlines():
            for expressed_self in expressed_selfs:
                if re.search(expressed_self, line):
                    if re.search('\tCDS\t', line):
                        chr = line.split('\t')[0]
                        start = line.split('\t')[3]
                        stop = line.split('\t')[4]
                        strand = line.split('\t')[6] 
                        expressed_selfs_tx_metadata[expressed_self]['cds'].append("{}:{}-{}".format(chr, start, stop))
                        expressed_selfs_tx_metadata[expressed_self]['strand'] = strand
                    #Assumes GENCODE GTF here...
                    elif re.search('\tUTR|five_prime_utr\t', line):
                        print(line)
                        chr = line.split('\t')[0]
                        start = line.split('\t')[3]
                        stop = line.split('\t')[4]
                        strand = line.split('\t')[6] 
                        expressed_selfs_tx_metadata[expressed_self]['UTR'].append("{}:{}-{}".format(chr, start, stop))

    expressed_selfs_nts = {}
    expressed_selfs_utr_buffers = {}

    for expressed_self in expressed_selfs_tx_metadata.keys():
        expressed_self_tx_seq = ''
        print(expressed_self)
        if expressed_selfs_tx_metadata[expressed_self]['strand'] == '+':
            print('positive_strand')
            print(sorted(expressed_selfs_tx_metadata[expressed_self]['cds'])[0])
            utr_buffer = 0
            if 'UTR' in expressed_selfs_tx_metadata[expressed_self].keys():
                for utr in expressed_selfs_tx_metadata[expressed_self]['UTR']:
                    print(utr)
                    print(utr.split(':')[1].split('-')[1])
                    print(sorted(expressed_selfs_tx_metadata[expressed_self]['cds'])[0].split(':')[1].split('-')[0])
                    if int(utr.split(':')[1].split('-')[1]) < int(sorted(expressed_selfs_tx_metadata[expressed_self]['cds'])[0].split(':')[1].split('-')[0]):
                        print("Adding!")
                        utr_buffer += (int(utr.split(':')[1].split('-')[1])  + 1 - int(utr.split(':')[1].split('-')[0]))
            for cds in sorted(expressed_selfs_tx_metadata[expressed_self]['cds']):
                print(cds)
                expressed_self_tx_seq += expressed_selfs_exon_seqs[cds]
        elif expressed_selfs_tx_metadata[expressed_self]['strand'] == '-':
            print('negative_strand')
            print(expressed_self)
            utr_buffer = 0
            if 'UTR' in expressed_selfs_tx_metadata[expressed_self].keys():
                for utr in expressed_selfs_tx_metadata[expressed_self]['UTR']: 
                    print("CDS: {}".format(sorted(expressed_selfs_tx_metadata[expressed_self]['cds'], reverse=True)[0]))
                    print("UTR: {}".format(utr))
                    if int(utr.split(':')[1].split('-')[0]) > int(sorted(expressed_selfs_tx_metadata[expressed_self]['cds'], reverse=True)[0].split(':')[1].split('-')[1]):
                        utr_buffer += (int(utr.split(':')[1].split('-')[1]) + 1  - int(utr.split(':')[1].split('-')[0]))
                        print("Adding!")
            for cds in sorted(expressed_selfs_tx_metadata[expressed_self]['cds']):
                print(cds)
                expressed_self_tx_seq += expressed_selfs_exon_seqs[cds]
            expressed_self_tx_seq = str(expressed_self_tx_seq.reverse_complement())
#        if expressed_selfs_tx_metadata[expressed_self]['strand'] == '-':
        print(expressed_self_tx_seq)
        expressed_selfs_nts[expressed_self] = expressed_self_tx_seq
        expressed_selfs_utr_buffers[expressed_self] = utr_buffer
                
   
     
        
#    with open(args.self_seqs) as eso:
#        for line in eso.readlines():
#            expressed_selfs.append(line.rstrip())

#    tx_to_aa = load_tx_aas(args)
#    tx_to_aa_trunc = {k.partition('.')[0]:v for k, v in tx_to_aa.items()}


    # Need to incorporate the germline variants here.

#    for expressed_self in expressed_selfs:
#        tx = expressed_self.partition(':')[2]
#        if tx in tx_to_aa_trunc.keys():
#            aa_seq = tx_to_aa_trunc[tx]
#            expressed_self_seqs[tx] = aa_seq

#    relevant_germline_vars = {}

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

#    print("Extracting missense variants from germline VCF...")
#    missense_snvs = extract_missense_snvs(args.germline_vcf)
#    print(len(missense_snvs.keys()))
#    print("Extracting conservative inframe indels from germline VCF...")
#    conserv_inframe_indels = extract_conservative_inframe_indels(args.germline_vcf)
#    print(len(conserv_inframe_indels.keys()))
#    print("Extracting disruptive inframe indels from germline VCF...")
#    disrupt_inframe_indels = extract_disruptive_inframe_indels(args.germline_vcf)
#    print(len(disrupt_inframe_indels.keys()))
#    print("Extracting frameshift indels from germline VCF...")
#    frameshift_indels = extract_frameshift_indels(args.germline_vcf)
#    print(len(frameshift_indels.keys()))

#    print("Number of expressed self-antigen transcripts: {}".format(len(expressed_self_seqs.keys())))

#    frameshift_transcripts = [record['transcript'] for entry, record in frameshift_indels.items()]
#    print(len(frameshift_transcripts))

#    for self_tx in expressed_self_seqs.keys():
#        if self_tx in frameshift_transcripts:
#            del(expressed_self_seqs[self_tx])

#    print("Number of expressed self-antigen transcripts without frameshifts: {}".format(len(expressed_self_seqs.keys())))

#    all_vars = dict(missense_snvs, **conserv_inframe_indels)
#    all_vars = dict(all_vars, **disrupt_inframe_indels)

#    print(len(all_vars.keys()))

#    for entry, record in all_vars.items():
#        if record['transcript'] in expressed_self_seqs.keys() or record['transcript'].split('.')[0] in expressed_self_seqs.keys():
#            if record['transcript'] not in relevant_germline_vars.keys():
#                relevant_germline_vars[record['transcript'].split('.')[0]] = [record]
#            else:
#                relevant_germline_vars[record['transcript'].split('.')[0]].append(record)

#    print(relevant_germline_vars)


    selfs_peps = {}
    for self, self_seq in expressed_selfs_nts.items():
        seq = Seq(str(self_seq))
        pep_seq = seq.translate(to_stop=True)
        if pep_seq.startswith('M'):
            selfs_peps[self] = pep_seq
#        if expressed_self in relevant_germline_vars.keys():
#            relevant_vars = relevant_germline_vars[expressed_self]
#            print("{}\n{}".format(expressed_self, relevant_vars))
#            if not(relevant_vars):
#                print("No relevant variants, emitting...")
#                tx_to_peps[expressed_self] = expressed_self_seq
#
#            indels = [i for i in relevant_vars if 'aa3_change' in i.keys()]
#            if not indels:
#                print("No indels detected. Proceeding with missense variants.")
#                # TODO: Current code assumes homozygosity of germline variants
#                # which is a poor assuming. Check the GT from the record.
#                germ_aa = list(str(expressed_self_seq))
#                germline_aas = [germ_aa]
#                for record in relevant_vars:
#                    print(record)
#                    if len(germ_aa) != int(record['aa_len']):
#                        print("transcript {} shows different lengths between amino acid fasta ({}) and snpEff annotations! ({})".format(tx, record['aa_len'], len(germ_aa)))
#                        continue
#                    aa_pos = int(record['aa_pos']) - 1
#                    germ_aa[aa_pos] = record['alt_aa']
#                    print("{}".format(expressed_self_seq[aa_pos-5:aa_pos+5]))
#                    print("{}".format(''.join(germ_aa)[aa_pos-5:aa_pos+5]))
#                    #print(record['meta'].FORMAT['GT'])
#                print("Applied all missense variants, emitting...")
#                tx_to_peps[expressed_self] = expressed_self_seq
#            else:
#                pass
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
        for expressed_self, expressed_self_seq in sorted(selfs_peps.items()):
            ofo.write(">MD5:{} NAME:{}\n{}\n".format(hashlib.md5("{}".format(expressed_self).encode('utf-8')).hexdigest()[:16], expressed_self, expressed_self_seq))
    
    with open(args.nt_output, 'w') as ofo:
        for self, self_seq in expressed_selfs_nts.items():
            seq = Seq(str(self_seq))
            ofo.write(">MD5:{} NAME:{} UTR_BUFFER:{}\n{}\n".format(hashlib.md5("{}".format(self).encode('utf-8')).hexdigest()[:16], self, expressed_selfs_utr_buffers[self], seq))


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

    uq = np.percentile(list(tx_to_log2tpm.values()), 75)

    for k,v in tx_to_log2tpm.items():
        tx_to_uqlog2tpm[k] = v/uq
    

    md5_to_tx = {}

    for seq_record in SeqIO.parse(args.fasta, "fasta"):
        id = seq_record.id
        print("ID: {}".format(id))
        print("Description: {}".format(seq_record.description))
        md5 = seq_record.description.split(' ')[0].replace(':', '_')[:15]
        tx = seq_record.description.split(' ')[1].replace('NAME:', '')
        print(md5)
        print(tx)
        md5_to_tx[md5] = tx
      

    

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
            if len(line) > 3 and line[2] == 'transcript' and line[1] != 'geve':
                gene_name = str(line).split('gene_name "')[1].split('"')[0]
                tx_id = str(line).split('transcript_id "')[1].split('"')[0].partition('.')[0]
                tx_to_gene[tx_id] = gene_name

    additional_notes = {}
    with open(args.gene_list) as glo:
        for line in glo.readlines():
            line = line.rstrip().split('\t')
            if len(line) > 1:
                additional_notes[line[0]] = line[1]
            else:
                additional_notes[line[0]] = 'NA'

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
                header.extend(['gene_name', 'transcript_identifier', 'tpm', 'log2(tpm+1)', 'uq(log2(tpm+1))', 'additional_notes'])
                header.insert(0, 'antigen_source')
                header = header[:5] + header[11:]
                header[1] = 'pos'
                header[2] = 'mhc_allele'
                header[3] = 'peptide'
                header[4] = 'peptide_core'
                header[5] = 'internal_identifier'
                header[6] = 'score_el'
                header[7] = 'percent_rank_el'
                header[8] = 'score_ba'
                header[9] = 'percent_rank_ba'
                header[10] = 'binding_affinity'
            elif len(line) < 16 or line[0] in ['Protein', 'Pos']:
                continue
            elif float(line[15]) < 1000:
                if 'SB' in line:
                    line.remove('<=')
                    line.remove('SB')
                if 'WB' in line:
                    line.remove('<=')
                    line.remove('WB')
                line.insert(0, 'Self-Antigen')
                line = line[:5] + line[11:]
                tx_id = line[5]
                if re.search('MD5', line[5]):
                    tx_id = md5_to_tx[line[5]]
                gene_name = tx_to_gene[tx_id]
                # The below line is simply for TCGA-LAML self-antigen filtering
                # and should be removed eventually.
                tpm = tx_to_tpm[tx_id]
                log2tpm = tx_to_log2tpm[tx_id.split('.')[0]]
                uqlog2tpm = tx_to_uqlog2tpm[tx_id.split('.')[0]]
                line.extend([gene_name, tx_id,str(tpm), str(log2tpm), str(uqlog2tpm), additional_notes[gene_name]])
                output_lines.append(line)

    with open(args.output, 'w') as ofo:
        if header:
            ofo.write("{}\n".format('\t'.join(header)))
            for line in output_lines:
                ofo.write("{}\n".format('\t'.join(line)))
        else:
            ofo.write("No viable fusion peptides. Try loosening parameters.")

    return output_lines


def filter_viral_cds(args):
    """
    """
    raw_counts = []                                                                                 
    raw_refs = []                                                                                   
    with open(args.viral_cds_quants) as vco:                                                            
        raw_counts = vco.readlines()[0].rstrip().split()[1:]                                        
    with open(args.viral_cds_ref) as vro:                                                               
        for line in vro.readlines():                                                                
            if line.startswith('>'):                                                                
                line = line.rstrip('\n').lstrip('>')                                                
                raw_refs.append(line.split(' ')[0])                                                 
#    viral_counts = {raw_refs[i]:raw_counts[i] for i in range(len(raw_refs)) if int(raw_counts[i]) > int(args.min_threshold)}
    all_viral_cds_counts = {raw_refs[i]:int(raw_counts[i]) for i in range(len(raw_refs)) if int(raw_counts[i]) > 0} 
    print(all_viral_cds_counts)

    expressed_viruses = []
    with open(args.expressed_viruses) as fo:
        for line_idx, line in enumerate(fo.readlines()):
            if line_idx > 0:
                virus_id = line.split(' ')[0].split('|')[3]
                expressed_viruses.append(virus_id)

    print(expressed_viruses)


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
                raw_refs.append(line.split(' ')[0])
#    viral_counts = {raw_refs[i]:raw_counts[i] for i in range(len(raw_refs)) if int(raw_counts[i]) > int(args.min_threshold)}
    all_viral_counts = {raw_refs[i]:int(raw_counts[i]) for i in range(len(raw_refs))}
    pprint(all_viral_counts)
    threshold = ''
    if not args.min_threshold:
        counts = [float(x) for x in all_viral_counts.values()]
        threshold = np.mean(counts) + (20*np.std(counts))
        print("Count threshold (mean + 20*std): {}".format(threshold))
    else:
        threshold= args.min_threshold
    expressed_viruses = {k:v for k, v in all_viral_counts.items() if int(v) > int(threshold)}
    if expressed_viruses:
        with open(args.output, 'w') as ofo:
            ofo.write("virus_ref\tread_count\n")
            for k, v in expressed_viruses.items():
                ofo.write("{}\t{}\n".format(k, v))


def make_viral_peptides(args):
    """
    """
#    alns = pysam.AlignmentFile(args.viral_bam, 'rb')
#
#    expressed_viruses_quant = {}
#    virus_to_contig = {}
#    viral_peptide_seqs = {}
#    expressed_viral_proteins = {}
#
#    with open(args.expressed_viruses) as fo:
#        for line_idx, line in enumerate(fo.readlines()):
#            if line_idx != 0:
#                line = line.rstrip().split('\t')
#                virus_to_contig[line[0].split('|')[3]] = line[0]
#                expressed_viruses_quant[line[0].split('|')[3]] = line[1]
#
    patient_nts = {}
#    virus_to_cds = {}
#    cds_to_virus = {}
#
    for seq_record in SeqIO.parse(args.fasta, "fasta"):
#        map_info = seq_record.description.split('|')[1].split(' ')[0]
#        virus_id, unneeded, protein_id = map_info.partition('_cds_')
#        protein_id = re.sub(r'_1$', '', protein_id)
#
        patient_nts[seq_record.description.split()[0].split('_cds_')[1][:-2]] = seq_record.seq

#        if virus_id not in virus_to_cds.keys():
#            virus_to_cds[virus_id] = [protein_id]
#        else:
#            virus_to_cds[virus_id].append(protein_id)
#
#        cds_to_virus[protein_id] = virus_id
#
#    cds_coords = {}
#
#    with open(args.viral_gff) as fo:
#        for line in fo.readlines():
#            line = line.rstrip().split('\t')
#            meta = line[-1]
#            protein_id = {x.split('=')[0]: x.split('=')[1] for x in meta.split(';')}['protein_id']
#            cds_coords[protein_id] = (line[3], line[4])
#
#
#
##    with open(args.fasta) as fo:
##        for line in fo.readlines():
##            if line.startswith('>'):
##                line = line.split('|')[1].split(' ')[0]
##                virus_id, unneeded, protein_id = line.partition('_cds_')
##                protein_id = re.sub(r'_1$', '', protein_id)
##                if virus_id not in virus_to_cds.keys():
##                    virus_to_cds[virus_id] = [protein_id]
##                else:
##                    virus_to_cds[virus_id].append(protein_id)
#
##    for seq_record in SeqIO.parse(args.viral_pep_ref, "fasta"):
##        peptide_id = seq_record.description.split(' ')[0]
##        viral_peptide_seqs[peptide_id] = seq_record.seq
#
#    expressed_viral_proteins = {}
#
#    for expressed_virus in expressed_viruses_quant.keys():
#        print(expressed_virus)
#        if expressed_virus in virus_to_cds.keys():
#             for protein in virus_to_cds[expressed_virus]:
#                expressed_viral_proteins[protein] = protein_nts[protein]
#        else:
#            print("Cannot find virus {} in cds data.".format(expressed_virus))
#
#    print(expressed_viral_proteins)
#
#    # evp = expressed viral protein
#    evp_rna_coverage = {}                                                                          
#    evp_het_sites = {}                                                                             
#    evp_patient_seq = {}    
#
#
#    for evp, evp_nts in expressed_viral_proteins.items():
#        evp_contig = virus_to_contig[cds_to_virus[evp]]
#        print(evp_contig)
#        evp_start = int(cds_coords[evp][0])
#        evp_stop = int(cds_coords[evp][1])
#        print(evp_start, evp_stop)
#        evp_per_pos_coverage_nt_level = alns.count_coverage(evp_contig, evp_start, evp_stop)
#        evp_per_pos_coverage_total = []
#        for base in evp_per_pos_coverage_nt_level[0]:
#            evp_per_pos_coverage_total.append(0)
#        for base in evp_per_pos_coverage_nt_level:
#            for pos_idx, count in enumerate(base):
#                evp_per_pos_coverage_total[pos_idx] = evp_per_pos_coverage_total[pos_idx] + count
#        print(evp_per_pos_coverage_total)
#        evp_coverage_avg = float(sum(evp_per_pos_coverage_total))/len(evp_per_pos_coverage_total)
#        print(evp_coverage_avg)
#
#        # 25x coverage is an abritrary threshold for now.
#        if evp_coverage_avg >= 25:
#            evp_patient_seq[evp] = []                                                                  
#            evp_het_sites[evp] = []  
#            evp_seq = evp_nts
#            for i in alns.pileup(contig=vpro_contig):
#                print(i)
#            tmp_info = pileup_truncated(alns, evp_contig, evp_start, evp_stop)                    
#            pileup_info = [(x.reference_pos, [i for i in x.get_query_sequences() if i]) for x in tmp_info]
#                                                                                                    
#            covered_pos = []                                                                            
#            for pos_idx, pos in enumerate(pileup_info):                                                 
#                                                                                                    
#                query_seqs = [x.lower() for x in pos[1]]                                                
#                print(query_seqs)
#                if not query_seqs:                                                                      
#                    evp_patient_seq[evp].append(evp_nts[pos_idx])                                     
#                if query_seqs:                                                                          
#                    covered_pos.append(pos[0])                                                          
#                    counts= {x:query_seqs.count(x) for x in query_seqs}                                 
#                    most_freq_finder = lambda x: scipy.stats.mode(x)[0][0]                              
#                    most_freq_base = most_freq_finder(query_seqs)                                       
#                    evp_patient_seq[evp].append(most_freq_base)                                        
#                    mismatch = ''                                                                       
#                    if vpro_seq[pos_idx] != most_freq_base:                                             
#                        mismatch = "Mismatch"                                                           
#                    print("{}\t{}\t{}\t{}".format(pos[0], herv_seq[pos_idx], most_freq_base, mismatch)) 
#                    most_freq_base_freq = float(query_seqs.count(most_freq_base))/len(query_seqs)       
#                    if float(most_freq_base_freq) < 0.80:                                               
#                        evp_het_sites[evp].append(pos[0])                                              
#                        print("WARNING: Possible germline mutation!")                                   
#                        print("{}\t{}\t{}\t{}\t{}".format(pos[0], herv_seq[pos_idx], most_freq_base, query_seqs.count(most_freq_base)/len(query_seqs), counts))
#
#            evp_rna_coverage[evp] = "{},{}".format(covered_pos[0]+1, covered_pos[-1]+1)                
#                             

    patient_aas = {}                                                                      
                                                                                                    
    for id, seq in patient_nts.items():                                                    
        aa_seq = seq.translate(to_stop=True)                          
        patient_aas[id] = aa_seq    


    with open(args.output, 'w') as ofo:
        for id, seq in patient_aas.items():
            ofo.write(">{}\n{}\n".format(id.partition('_1:1')[0], seq))
    
#    with open(args.nt_output, 'w') as ofo:
#        for id, seq in patient_nts:
#            ofo.write(">{}\n{}\n".format(id, seq))


def add_viral_metadata(args):
    """
    """
#    expressed_viruses = {}
    cds_to_virus = {}

    with open(args.viral_cds_ref) as fo:
        for line in fo.readlines():
            if line.startswith('>'):
                line = line.split('|')[1].split(' ')[0]
                virus_id, trash, protein_id = line.partition('_cds_')
                protein_id = re.sub(r'_1$', '', protein_id)
                cds_to_virus[protein_id] = virus_id

#    with open(args.viral_quants) as fo:
#        for line_idx, line in enumerate(fo.readlines()):
#            line = line.rstrip().split('\t')
#            expressed_viruses[line[0].split('_cds_')[1][:-2]] = line[1]

#    print(expressed_viruses)

    output_lines = []
    header = []


    viral_protein_names = {}
    viral_names = {}
    with open(args.viral_pep_ref) as fo:
        for line in fo.readlines():
            if line.startswith('>'):
                line = line.rstrip()
                part = line.partition(' ')
                viral_protein_id = part[0].lstrip('>')

                id_and_name = part[2].partition('[')
                viral_protein_name = id_and_name[0]
                virus_name = id_and_name[2]
                viral_protein_names[viral_protein_id] = viral_protein_name.rstrip(' ')
                viral_names[viral_protein_id] = virus_name.rstrip(']')

    print(viral_names)

    with open(args.binding_affinities) as fo:
        for line_idx, line in enumerate(fo.readlines()):
            line = line.rstrip().split()
            if len(line) >= 16 and line[0] == 'Pos':
                header = line
                header.insert(0, 'antigen_source')
                header = header[:5] + header[11:]
                try:
                    header.remove('BindLevel')
                except:
                    pass
#                header.extend(['virus_protein_identifier', 'virus_protein_name', 'virus_identifier', 'virus_name', 'rna_read_count'])
                header.extend(['virus_protein_identifier', 'virus_protein_name', 'virus_identifier', 'virus_name'])
                header[1] = 'pos'
                header[2] = 'mhc_allele'
                header[3] = 'peptide'
                header[4] = 'peptide_core'
                header[5] = 'internal_identifier'
                header[6] = 'score_el'
                header[7] = 'percent_rank_el'
                header[8] = 'score_ba'
                header[9] = 'percent_rank_ba'
                header[10] = 'binding_affinity'
            elif len(line) < 16 or line[0] in ['Protein', 'Pos']:
                continue
            elif float(line[15]) < 1000:
                if 'SB' in line:
                    line.remove('<=')
                    line.remove('SB')
                if 'WB' in line:
                    line.remove('<=')
                    line.remove('WB')
                line.insert(0, 'Virus')
                line = line[:5] + line[11:]
                #This should be extended to replace and _[0-9] to .[0-9]
                orig_index = line[5].rfind('_')
                virus_cds = "NA"
                virus_id = "NA"
                viral_name = "NA"
                viral_protein_name = "NA" 
                read_counts = "NA"
                virus_cds = line[5][:orig_index] + '.' + line[5][orig_index+1:]
                print(virus_cds)
                virus_id = cds_to_virus[virus_cds]
#                read_counts = expressed_viruses[virus_cds]
                viral_name = viral_names[virus_cds]
                viral_protein_name = viral_protein_names[virus_cds]
#                line.extend([virus_cds, viral_protein_name, virus_id, viral_name, read_counts])
                line.extend([virus_cds, viral_protein_name, virus_id, viral_name])
                output_lines.append(line)

    with open(args.output, 'w') as ofo:
        ofo.write("{}\n".format('\t'.join(header)))
        for line in output_lines:
            ofo.write("{}\n".format('\t'.join(line)))

def make_fusion_peptides_context(args):
    """
    """

    fusion_txs = []
    with open(args.fusion_txs) as fo:
        for line in fo.readlines():
            fusion_txs.append(line.rstrip())

    print("Fusion transcripts: {}".format(fusion_txs))

    exon_seqs = {}
#    for seq_record in SeqIO.parse(glob(os.path.join(args.var_tx_seqs, '*{}_{}.normal.fa'.format(record['transcript'].partition('.')[0], record_coords)))[0], "fasta"):
    for seq_record in SeqIO.parse(args.exons_fasta, "fasta"):#glob(os.path.join(args.var_tx_seqs, '*{}_{}.normal.fa'.format(record['transcript'].partition('.')[0], record_coords)))[0], "fasta"):
        exon_seqs[seq_record.description] = seq_record.seq
   
    fusion_txs_metadata = {} 
    for fusion_tx in fusion_txs:
        fusion_txs_metadata[fusion_tx] = {}
        fusion_txs_metadata[fusion_tx]['cds'] = []
        fusion_txs_metadata[fusion_tx]['strand'] = ''
    with open(args.gtf) as fo:
        for line in fo.readlines():
            for fusion_tx in fusion_txs:
                if re.search(fusion_tx, line):
                    if re.search('\tCDS\t', line):
                        print(fusion_tx)
                        chr = line.split('\t')[0]
                        start = line.split('\t')[3]
                        stop = line.split('\t')[4]
                        strand = line.split('\t')[6] 
                        coords = "{}:{}-{}".format(chr, start, stop)
                        print(coords)
                        print(strand)
                        fusion_txs_metadata[fusion_tx]['strand'] = strand
                        print("Actual dict: {}".format(fusion_txs_metadata[fusion_tx]['strand']))
                        if coords not in fusion_txs_metadata[fusion_tx]['cds']:
                            print("Adding CDS")
                            fusion_txs_metadata[fusion_tx]['cds'].append("{}:{}-{}".format(chr, start, stop))
                            print("Actual dict: {}".format(fusion_txs_metadata[fusion_tx]['cds']))   

    pprint(fusion_txs_metadata)   
 
    fusion_tx_seqs = {}
    for fusion_tx in fusion_txs:
        fusion_tx_seqs[fusion_tx] = ''
        print(fusion_tx)
        if fusion_txs_metadata[fusion_tx]['strand'] == '+':
            print('positive_strand')
            for cds in sorted(fusion_txs_metadata[fusion_tx]['cds']):
                fusion_tx_seqs[fusion_tx] += exon_seqs[cds]
        elif fusion_txs_metadata[fusion_tx]['strand'] == '-':
            print('negative_strand')
            for cds in sorted(fusion_txs_metadata[fusion_tx]['cds'], reverse=True):
                fusion_tx_seqs[fusion_tx] += exon_seqs[cds].reverse_complement()

    pprint(fusion_tx_seqs)


    print("Loaded variant transcripts metadata.")
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
                    fused_tx = ''
                    left_tx = line[header_map['CDS_LEFT_ID']]
                    right_tx = line[header_map['CDS_RIGHT_ID']]
                    if left_tx in fusion_tx_seqs and right_tx in fusion_tx_seqs and fusion_tx_seqs[left_tx] and fusion_tx_seqs[right_tx]:
                        left_tx_end_pos = int(line[header_map['CDS_LEFT_RANGE']].split('-')[1])
                        right_tx_start_pos = int(line[header_map['CDS_RIGHT_RANGE']].split('-')[0]) - 1
                        print("Both transcripts observed.")
                        print("{}\t{}".format(left_tx, left_tx_end_pos))
                        print("{}\t{}".format(right_tx, right_tx_start_pos))
                        print(line[header_map['#FusionName']])
                        try:
                            fused_tx = fusion_tx_seqs[left_tx][left_tx_end_pos].lower() + fusion_tx_seqs[right_tx][right_tx_start_pos:].upper()
                            print(fused_tx)
                        except:
                            print("Indexing issues with fused_tx.")
                    start_prot = 0
                    nuc_sub = 48
                    print(line[header_map['#FusionName']])
                    print(line[header_map['CDS_LEFT_RANGE']])
                    print(line[header_map['CDS_LEFT_RANGE']].split('-')[1])
                    valid_fusions[line[header_map['#FusionName']]] = {}
                    if float(line[header_map['CDS_LEFT_RANGE']].split('-')[1]) % 3 == 0:
                        print("Direct translate")
                        start_prot = int(int(int(line[header_map['CDS_LEFT_RANGE']].split('-')[1]))/3)
                    elif (float(line[header_map['CDS_LEFT_RANGE']].split('-')[1]) - 1) % 3 == 0:
                        print("Minus one")
                        start_prot = int(int(int(line[header_map['CDS_LEFT_RANGE']].split('-')[1]) - 1)/3)
                        nuc_sub = nuc_sub + 1
                    elif (float(line[header_map['CDS_LEFT_RANGE']].split('-')[1]) - 2) % 3 == 0:
                        print("Minus two")
                        start_prot = int(int(int(line[header_map['CDS_LEFT_RANGE']].split('-')[1]) - 2)/3)
                        nuc_sub = nuc_sub + 2
                    else:
                        print("No luck")
                    print(start_prot)
                    full_nuc_seq = ''
                    full_pep_seq = ''
                    if fused_tx:
                        full_nuc_seq = fused_tx
                        full_pep_seq = str(full_nuc_seq.translate())
                    else:
                        full_nuc_seq = line[header_map['FUSION_CDS']]
                        full_pep_seq = line[header_map['FUSION_TRANSL']]
                    if line[header_map['PROT_FUSION_TYPE']] == 'INFRAME':
                        print('inframe')
#                        peptide = line[header_map['FUSION_TRANSL']][start_prot - 8:start_prot + 8]
#                        peptide_context = line[header_map['FUSION_TRANSL']][start_prot - 16:start_prot + 16]
#                        nuc_context = line[header_map['FUSION_CDS']][int(line[header_map['CDS_LEFT_RANGE']].split('-')[1]) - nuc_sub:int(line[header_map['CDS_LEFT_RANGE']].split('-')[1]) + nuc_sub]
#                        nuc_seq = line[header_map['FUSION_CDS']][int(line[header_map['CDS_LEFT_RANGE']].split('-')[1]) - (nuc_sub - 24):int(line[header_map['CDS_LEFT_RANGE']].split('-')[1]) + (nuc_sub - 24)]
                        peptide = full_pep_seq[start_prot - 8:start_prot + 8]
                        peptide_context = full_pep_seq[start_prot - 16:start_prot + 16]
                        nuc_context = full_nuc_seq[int(line[header_map['CDS_LEFT_RANGE']].split('-')[1]) - nuc_sub:int(line[header_map['CDS_LEFT_RANGE']].split('-')[1]) + nuc_sub]
                        nuc_seq = full_nuc_seq[int(line[header_map['CDS_LEFT_RANGE']].split('-')[1]) - (nuc_sub - 24):int(line[header_map['CDS_LEFT_RANGE']].split('-')[1]) + (nuc_sub - 24)]
                        print("peptide: {}".format(peptide))
                        print("peptide_context: {}".format(peptide_context))
                        print("nuc_context: {}".format(nuc_context))
                        if not(re.search('\*', peptide)):
                            valid_fusions[line[header_map['#FusionName']]]['peptide'] = peptide
                            valid_fusions[line[header_map['#FusionName']]]['peptide_context'] = peptide_context
                            valid_fusions[line[header_map['#FusionName']]]['nuc_context'] = nuc_context
                            valid_fusions[line[header_map['#FusionName']]]['nuc_seq'] = nuc_seq
                        else:
                            valid_fusions[line[header_map['#FusionName']]]['peptide'] = peptide
                            valid_fusions[line[header_map['#FusionName']]]['peptide_context'] = peptide_context
                            valid_fusions[line[header_map['#FusionName']]]['nuc_seq'] = nuc_seq
                    elif line[header_map['PROT_FUSION_TYPE']] == 'FRAMESHIFT':
                        print('frameshift')
#                        init_peptide = line[header_map['FUSION_TRANSL']][start_prot - 8:]
#                        init_peptide_context = line[header_map['FUSION_TRANSL']][start_prot - 16:]
#                        init_nuc_context = line[header_map['FUSION_CDS']][int(line[header_map['CDS_LEFT_RANGE']].split('-')[1]) - nuc_sub:]
#                        init_nuc_seq = line[header_map['FUSION_CDS']][int(line[header_map['CDS_LEFT_RANGE']].split('-')[1]) - 24:]
                        init_peptide = full_pep_seq[start_prot - 8:]
                        init_peptide_context = full_pep_seq[start_prot - 16:]
                        init_nuc_context = full_nuc_seq[int(line[header_map['CDS_LEFT_RANGE']].split('-')[1]) - nuc_sub:]
                        init_nuc_seq = full_nuc_seq[int(line[header_map['CDS_LEFT_RANGE']].split('-')[1]) - 24:]
                        print(init_peptide)
                        if re.search('\*', init_peptide):
                            stop_codon = init_peptide.index('*')
                            peptide = init_peptide[:stop_codon]
                            peptide_context = init_peptide_context[:stop_codon+8]
                            valid_fusions[line[header_map['#FusionName']]]['peptide'] = peptide
                            valid_fusions[line[header_map['#FusionName']]]['peptide_context'] = peptide_context
                            valid_fusions[line[header_map['#FusionName']]]['nuc_context'] = init_nuc_context
                            valid_fusions[line[header_map['#FusionName']]]['nuc_seq'] = init_nuc_seq
                        else:
                            valid_fusions[line[header_map['#FusionName']]]['peptide'] = init_peptide
                            valid_fusions[line[header_map['#FusionName']]]['peptide_context'] = init_peptide_context
                            valid_fusions[line[header_map['#FusionName']]]['nuc_context'] = init_nuc_context
                            valid_fusions[line[header_map['#FusionName']]]['nuc_seq'] = init_nuc_seq

    with open(args.output, 'w') as ofo:
        for valid_fusion_id in valid_fusions.keys():
            print(valid_fusion_id)
            meta = valid_fusions[valid_fusion_id]
            print(meta)
            # Context sequences are too long.
            ofo.write(">{} NUCLEOTIDE_CONTEXT:{} PEPTIDE_CONTEXT:{}\n{}\n".format(valid_fusion_id, 'N/A', 'N/A', meta['peptide']))
    
    with open(args.nt_output, 'w') as ofo:
        for valid_fusion_id in valid_fusions.keys():
            print(valid_fusion_id)
            meta = valid_fusions[valid_fusion_id]
            print(meta)
            # Context sequences are too long.
            ofo.write(">{} NUCLEOTIDE_CONTEXT:{} PEPTIDE_CONTEXT:{}\n{}\n".format(valid_fusion_id, 'N/A', 'N/A', meta['nuc_seq']))

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
                     start_prot = 0
                     nuc_sub = 48
                     print(line[header_map['#FusionName']])
                     print(line[header_map['CDS_LEFT_RANGE']])
                     print(line[header_map['CDS_LEFT_RANGE']].split('-')[1])
                     valid_fusions[line[header_map['#FusionName']]] = {}
                     if float(line[header_map['CDS_LEFT_RANGE']].split('-')[1]) % 3 == 0:
                         print("Direct translate")
                         start_prot = int(int(int(line[header_map['CDS_LEFT_RANGE']].split('-')[1]))/3)
                     elif (float(line[header_map['CDS_LEFT_RANGE']].split('-')[1]) - 1) % 3 == 0:
                         print("Minus one")
                         start_prot = int(int(int(line[header_map['CDS_LEFT_RANGE']].split('-')[1]) - 1)/3)
                         nuc_sub = nuc_sub + 1
                     elif (float(line[header_map['CDS_LEFT_RANGE']].split('-')[1]) - 2) % 3 == 0:
                         print("Minus two")
                         start_prot = int(int(int(line[header_map['CDS_LEFT_RANGE']].split('-')[1]) - 2)/3)
                         nuc_sub = nuc_sub + 2
                     else:
                         print("No luck")
                     print(start_prot)
                     if line[header_map['PROT_FUSION_TYPE']] == 'INFRAME':
                         print('inframe')
                         peptide = line[header_map['FUSION_TRANSL']][start_prot - 8:start_prot + 8]
                         peptide_context = line[header_map['FUSION_TRANSL']][start_prot - 16:start_prot + 16]
                         nuc_context = line[header_map['FUSION_CDS']][int(line[header_map['CDS_LEFT_RANGE']].split('-')[1]) - nuc_sub:int(line[header_map['CDS_LEFT_RANGE']].split('-')[1]) + nuc_sub]
                         nuc_seq = line[header_map['FUSION_CDS']][int(line[header_map['CDS_LEFT_RANGE']].split('-')[1]) - (nuc_sub - 24):int(line[header_map['CDS_LEFT_RANGE']].split('-')[1]) + (nuc_sub - 24)]
                         print("peptide: {}".format(peptide))
                         print("peptide_context: {}".format(peptide_context))
                         print("nuc_context: {}".format(nuc_context))
                         if not(re.search('\*', peptide)):
                             valid_fusions[line[header_map['#FusionName']]]['peptide'] = peptide
                             valid_fusions[line[header_map['#FusionName']]]['peptide_context'] = peptide_context
                             valid_fusions[line[header_map['#FusionName']]]['nuc_context'] = nuc_context
                             valid_fusions[line[header_map['#FusionName']]]['nuc_seq'] = nuc_seq
                         else:
                             valid_fusions[line[header_map['#FusionName']]]['peptide'] = peptide
                             valid_fusions[line[header_map['#FusionName']]]['peptide_context'] = peptide_context
                             valid_fusions[line[header_map['#FusionName']]]['nuc_seq'] = nuc_seq
                     elif line[header_map['PROT_FUSION_TYPE']] == 'FRAMESHIFT':
                         print('frameshift')
                         init_peptide = line[header_map['FUSION_TRANSL']][start_prot - 8:]
                         init_peptide_context = line[header_map['FUSION_TRANSL']][start_prot - 16:]
                         init_nuc_context = line[header_map['FUSION_CDS']][int(line[header_map['CDS_LEFT_RANGE']].split('-')[1]) - nuc_sub:]
                         init_nuc_seq = line[header_map['FUSION_CDS']][int(line[header_map['CDS_LEFT_RANGE']].split('-')[1]) - 24:]
                         print(init_peptide)
                         if re.search('\*', init_peptide):
                             stop_codon = init_peptide.index('*')
                             peptide = init_peptide[:stop_codon]
                             peptide_context = init_peptide_context[:stop_codon+8]
                             valid_fusions[line[header_map['#FusionName']]]['peptide'] = peptide
                             valid_fusions[line[header_map['#FusionName']]]['peptide_context'] = peptide_context
                             valid_fusions[line[header_map['#FusionName']]]['nuc_context'] = init_nuc_context
                             valid_fusions[line[header_map['#FusionName']]]['nuc_seq'] = init_nuc_seq
                         else:
                             valid_fusions[line[header_map['#FusionName']]]['peptide'] = init_peptide
                             valid_fusions[line[header_map['#FusionName']]]['peptide_context'] = init_peptide_context
                             valid_fusions[line[header_map['#FusionName']]]['nuc_context'] = init_nuc_context
                             valid_fusions[line[header_map['#FusionName']]]['nuc_seq'] = init_nuc_seq

     with open(args.output, 'w') as ofo:
         for valid_fusion_id in valid_fusions.keys():
             print(valid_fusion_id)
             meta = valid_fusions[valid_fusion_id]
             print(meta)
             # Context sequences are too long.
             ofo.write(">{} NUCLEOTIDE_CONTEXT:{} PEPTIDE_CONTEXT:{}\n{}\n".format(valid_fusion_id, 'N/A', 'N/A', meta['peptide']))
     
     with open(args.nt_output, 'w') as ofo:
         for valid_fusion_id in valid_fusions.keys():
             print(valid_fusion_id)
             meta = valid_fusions[valid_fusion_id]
             print(meta)
             # Context sequences are too long.
             ofo.write(">{} NUCLEOTIDE_CONTEXT:{} PEPTIDE_CONTEXT:{}\n{}\n".format(valid_fusion_id, 'N/A', 'N/A', meta['nuc_seq']))


def add_fusion_metadata(args):
    """
    """

    fusion_metadata = {}
#    with open(args.fusions) as fo:
#        col_to_idx = {}
#        for line_idx, line in enumerate(fo.readlines()):
#            line = line.rstrip().split()
#            if line_idx ==  0:
#                for col_idx, col in enumerate(line):
#                    col_to_idx[col_idx] = col
#            else:
#                fusion_metadata[line[0][:15]] = {}
#                for elem_idx, elem in enumerate(line[1:]):
#                    fusion_metadata[line[0][:15]][col_to_idx[elem_idx + 1]] = elem

#    with open(args.fusions) as fo:
#        for line in fo.readlines():
#            if line.startswith('>'):
#                line = line.lstrip('>').rstrip().split(' ')
#                nuc_context = line[1].split(':')[1] 
#                protein_context = line[2].split(':')[1] 
#                fusion_metadata[line[0][:15]]['nuc_context'] = nuc_context
#                fusion_metadata[line[0][:15]]['protein_context'] = protein_context


    metadata = ['JunctionReadCount', 'SpanningFragCount', 'SpliceType', 'LeftGene',
                'LeftBreakpoint', 'LeftBreakDinuc', 'LeftBreakEntropy', 'RightGene', 
                 'RightBreakpoint', 'RightBreakDinuc', 'RightBreakEntropy',
                 'LargeAnchorSupport', 'FFPM', 'PROT_FUSION_TYPE', 'annots']

    header_extension = ['fusion_junction_read_count',
                        'fusion_spanning_frag_count', 'fusion_splice_type',
                        'fusion_left_gene', 'fusion_left_breakpoint',
                        'fusion_left_break_dinuc', 'fusion_left_break_entropy',
                        'fusion_right_gene', 'fusion_right_breakpoint',
                        'fusion_right_break_dinuc',
                        'fusion_right_break_entropy',
                        'fusion_large_anchor_supprot', 'ffpm', 'fusion_type',
                        'fusion_annotations', 'nucleotide_context',
                        'protein_context']


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
#                header_extension = extensions[:]
                header.extend(header_extension)
                header.insert(0, 'antigen_source')
#                print(header)
                header = header[:5] + header[11:]
#                print(header)
#                sys.exit(0)
                header[1] = 'pos'
                header[2] = 'mhc_allele'
                header[3] = 'peptide'
                header[4] = 'peptide_core'
                header[5] = 'internal_identifier'
                header[6] = 'score_el'
                header[7] = 'percent_rank_el'
                header[8] = 'score_ba'
                header[9] = 'percent_rank_ba'
                header[10] = 'binding_affinity'
            elif len(line) < 16 or line[0] in ['Protein', 'Pos']:
                continue
            elif float(line[15]) < 1000:
                if 'SB' in line:
                    line.remove('<=')
                    line.remove('SB')
                if 'WB' in line:
                    line.remove('<=')
                    line.remove('WB')
                line.insert(0, 'FusionEvent')
                line = line[:5] + line[11:]
                print(line[5])
                print(line[3])
#                if '_' in line[5] and line[5].index('_') > line[5].index('-'):
#                    print(fusion_metadata[line[5].rsplit('_')[0].replace('_', '.')]['protein_context'])
#                else:
#                    print(fusion_metadata[line[5].replace('_', '.')]['protein_context'])

# Need to fix context code...
#                if re.search(line[3], fusion_metadata[line[5].rsplit('_')[0].replace('_', '.')]['protein_context']):
                
#                    idx_start = fusion_metadata[line[5].rsplit('_')[0].replace('_', '.')]['protein_context'].index(line[3])
#                    idx_stop = idx_start + len(line[3])
#                    adjusted_protein_idx_start = idx_start - 8
#                    adjusted_protein_idx_stop = idx_stop + 8
#                    nuc_idx_start = adjusted_protein_idx_start * 3
#                    nuc_idx_stop = adjusted_protein_idx_stop * 3
  
#                    print(adjusted_protein_idx_start)
#                    print(adjusted_protein_idx_stop)
#                    print(nuc_idx_start)
#                    print(nuc_idx_stop)

#                print(line[5])
#                print(line[5].rsplit('_')[0])
#                protein_context = fusion_metadata[line[5].rsplit('_')[0].replace('_', '.')]['protein_context'][max(0, adjusted_protein_idx_start):min(adjusted_protein_idx_stop, len(fusion_metadata[line[5].rsplit('_')[0].replace('_', '.')]['protein_context']))]
#                nuc_context = fusion_metadata[line[5].rsplit('_')[0].replace('_', '.')]['nuc_context'][max(0, nuc_idx_start):min(nuc_idx_stop, len(fusion_metadata[line[5].rsplit('_')[0].replace('_', '.')]['nuc_context']))]
                protein_context = "N/A"
                nuc_context = "N/A"

#                print("Peptide: {}".format(line[3]))
#                print("Full protein context: {}".format(fusion_metadata[line[5].rsplit('_')[0].replace('_', '.')]['protein_context']))
#                print("Protein context: {}".format(protein_context))
#                print("Nucleotide context: {}".format(nuc_context))

#                relevant_metadata = {}
#                if '_' in line[5] and line[5].index('_') > line[5].index('-'):
#                    relevant_metadata = fusion_metadata[line[5].rsplit('_')[0].replace('_', '.')]
#                else:
#                    relevant_metadata = fusion_metadata[line[5].replace('_', '.')]
                line_extension = []
#                for metadatum in metadata:
#                    line_extension.append(relevant_metadata[metadatum])
                line_extension.append(nuc_context)
                line_extension.append(protein_context)
                line.extend(line_extension)
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
    new_header = ['antigen_source', 'peptide', 'chromosome', 'gene_name', 'strand', 'mhc_allele', 'binding_affinity', 'nucleotide_context', 'neosplice_min_expression', 'reads_with_peptide']
    new_lines = []

    col_idx_map = {}

    mhc_info_map = {}

    with open(args.splice_summary) as fo:
        for line_idx, line in enumerate(fo.readlines()):
            line = line.strip()
            if line_idx == 0:
                header = line.split('\t')
                for col_idx, col in enumerate(header):
                    if col.endswith('_y'):
                        pass
                    elif col.endswith('_x'):
                        col = col.rstrip('_x')
                    if not(re.search('HLA', col)):
                        col_idx_map[col] = col_idx
                    else:
                        mhc_info_map[col] = col_idx

                print(col_idx_map)

            else:
                line = line.split('\t')
                base_new_line = ['SpliceVariant', line[col_idx_map['Variant_peptide_sequence']], line[col_idx_map['Chromosome']], line[col_idx_map['Gene']], line[col_idx_map['Strand']], 'ALLELE_PROXY', 'BA_PROXY', line[col_idx_map['DNA_sequence']], line[col_idx_map['min_expression']], line[col_idx_map['reads_with_peptide']]]
                for allele in list(set([i.partition('_')[0] for i in mhc_info_map.keys()])):
                    print(allele)
                    allele_specific = base_new_line[:]
                    allele_specific[5] = allele
                    allele_specific[6] = line[mhc_info_map["{}_Nm".format(allele)]]
                    print(allele_specific)
#                    allele_specific[7] = line[mhc_info_map["{}_binding_property".format(allele)]]

                    if float(allele_specific[6]) < 1000:
                        new_lines.append(allele_specific)
    with open(args.output, 'w') as ofo:
        ofo.write("{}\n".format('\t'.join(new_header)))
        for new_line in new_lines:
            ofo.write("{}\n".format('\t'.join(new_line)))






def add_rna_norms(args):
    """
    """
    hmap = {}
    manifest_entries = []
    with open(args.manifest) as fo:
        for line_idx, line in enumerate(fo.readlines()):
            if line_idx == 0:
                line = line.strip('\n').split('\t')
                for col_idx, col in enumerate(line):
                    hmap[col] = col_idx
                print(hmap)
            else:
                line = line.strip('\n')
                line = line.split('\t')
                trunc_line = [line[hmap['Patient_Name']], line[hmap['Run_Name']], line[hmap['Dataset']], line[hmap['File_Prefix']], line[hmap['Sequencing_Method']], line[hmap['Normal']]]
                manifest_entries.append('\t'.join(trunc_line))
        new_line = [args.pat_name, "nr-{}".format(args.pat_name), args.dataset, args.prefix, 'RNA-Seq', 'TRUE']
        manifest_entries.append('\t'.join(new_line))


    with open(args.output, 'w') as ofo:
        ofo.write('{}\n'.format('\t'.join(['Patient_Name', 'Run_Name', 'Dataset', 'File_Prefix', 'Sequencing_Method', 'Normal'])))
        for entry in manifest_entries:
            ofo.write("{}\n".format(entry))
    return manifest_entries


def make_lens_report(args):
    """
    """
    reports = glob(os.path.join(args.metadata_dir, "*metadata*"))
    print(reports)

    antigen_sources = ['self_antigen', 'indel', 'hervs', 'splice', 'snv', 'fusion', 'viral']

    #Initiating dataframe with snvs...

#    init_metadata = [i for i in reports if not re.search('fusion|indel', i) and not(pd.read_csv(i, sep='\t').empty)][0]
#    init_metadata = [i for i in reports if not pd.read_csv(i, sep='\t').empty][0]
 
#    report_df = pd.read_csv(init_metadata, sep='\t')

#    for report in [i for i in reports if i not in [init_metadata]]:
    report_df = pd.DataFrame()
    for report in reports:
        tmp_report = ''
        print(report)
        try:
            tmp_report = pd.read_csv(report, sep='\t')
        except:
            continue
        if re.search('self_antigen', report):
            print("SELF_ANTIGEN HERE!!!!")
            tmp_report = tmp_report[tmp_report["binding_affinity"] <= 50.0]
        if not(tmp_report.empty) and report_df.empty:
            report_df = tmp_report
        elif not(tmp_report.empty):
            report_df = pd.concat([report_df, tmp_report])


    #sorted_df = report_df.sort_values(by="Aff(nM)")
    #print(sorted_df)

    #Rearranging column headers.
    #if 'SpliceType' in sorted_df.columns:
    #    sorted_df = sorted_df[['MHC', 'Peptide', 'Icore', 'Aff(nM)', 'Score_BA', '%Rank_BA', 'Score_EL', '%Rank_EL', 'Pos', 'AntigenSource', 'VariantType', 'CCF', 'Agretopocity', 'Reference_Aff(nM)', 'Reference_Peptide', 'ReferenceAllele', 'AlternateAllele', 'Chromosome', 'VariantPosition', 'Strand', 'GeneName', 'MinimumExpression', 'TPM', 'ReadCount', 'PercentRnaReadsWithVariant', 'TotalRnaReadCoverage', 'TranscriptIdentifier', 'VirusIdentifier', 'BindingProperty', 'NucleotideContext', 'ProteinContext', 'JunctionReadCount', 'SpanningFragCount', 'SpliceType', 'LeftGene', 'RightGene', 'LargeAnchorSupport', 'FFPM', 'Identity']]
    #else:
    #    sorted_df = sorted_df[['MHC', 'Peptide', 'Icore', 'Aff(nM)', 'Score_BA', '%Rank_BA', 'Score_EL', '%Rank_EL', 'Pos', 'AntigenSource', 'VariantType', 'CCF', 'Agretopocity', 'Reference_Aff(nM)', 'Reference_Peptide', 'ReferenceAllele', 'AlternateAllele', 'Chromosome', 'VariantPosition', 'Strand', 'GeneName', 'MinimumExpression', 'TPM', 'ReadCount', 'PercentRnaReadsWithVariant', 'TotalRnaReadCoverage', 'TranscriptIdentifier', 'VirusIdentifier', 'BindingProperty', 'NucleotideContext', 'ProteinContext', 'Identity']]

    #sorted_df.to_csv(args.output, sep='\t', index=False, na_rep='NA')

    print(report_df.dtypes)
    report_df.astype({'binding_affinity': 'float64'})
    print(report_df.dtypes)

    filtered_df = report_df[report_df["binding_affinity"] <= 50.0]
    filtered_df.to_csv(args.output, sep='\t', index=False, na_rep='NA')


#####k
#    reports = glob(os.path.join(args.metadata_dir, "*metadata.txt"))
#    print(reports)
#
#    antigen_sources = ['snv', 'indel', 'viral', 'hervs', 'splice', 'fusion', 'self_antigen']
##    antigen_sources = ['snv', 'indel', 'viral', 'hervs', 'fusion', 'self_antigen']
#
#    #Initiating dataframe with snvs...
#
#    snv_metadata = [i for i in reports if re.search('snv.metadata', i)][0]
#    print(snv_metadata)
#
#    report_df = pd.read_csv(snv_metadata, sep='\t')
#
#    # Skipping SNV since it was used to initiate DF.
#    for antigen_source in antigen_sources[1:]:
#        print(antigen_source)
#        report = [i for i in reports if re.search('{}.metadata'.format(antigen_source), i)][0]
#        print(report)
#        tmp_report = pd.read_csv(report, sep='\t')
#        if not(tmp_report.empty):
#            report_df = pd.concat([report_df, tmp_report])


#    sorted_df = report_df.sort_values(by="Aff(nM)")

#    sorted_df.to_csv(args.output, sep='\t', index=False, na_rep='NA')

def make_antigens_barplot(args):
    """
    """
    patients = []
    counts = {}

    reports = glob(os.path.join(args.reports_dir, "*lens_report*.txt"))
    print(reports)

    antigen_sources = ['SNV', 'InDel', 'Virus', 'ERV', 'SpliceVariant', 'FusionEvent', 'Self-Antigen']
#    antigen_sources = ['SNV', 'InDel', 'Virus', 'ERV', 'SpliceVariant', 'FusionEvent']
    for antigen_source in antigen_sources:
        counts[antigen_source] = []

    for patient_report in reports:
        print(patient_report)
        patient_id = patient_report.split('/')[-1].replace('.lens_report.txt', '').partition('-')[2]
        print(patient_id)
        patients.append(patient_id)

        tmp_report = pd.read_csv(patient_report, sep='\t')
        count_df = tmp_report['antigen_source'].value_counts()
        print(count_df)
        for antigen_source in counts.keys():
            if antigen_source in count_df.index.values:
#                counts[antigen_source].append(np.log2(count_df[antigen_source]))
                counts[antigen_source].append(count_df[antigen_source])
            else:
                counts[antigen_source].append(0)
    counts['CTA/Self-Antigen'] = counts['Self-Antigen']
    del counts['Self-Antigen']

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

#    ax = df.plot.bar(stacked=True, colormap='viridis', log=True)
    ax = df.plot.bar(stacked=True, colormap='viridis')
    ax.set_xlabel("Replicates")
    ax.set_ylabel("Count of Predicted Peptides\n(<50 nM binding affinity)")
    ax.set_title("Tumor Antigen Sources among BBN963 cell line replicates")
    ax.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
    ax.legend(loc=2, prop={'size': 5.5}) 

    fig = ax.get_figure()
    fig.tight_layout()
    fig.savefig(args.output, format='pdf')

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
               line.extend(['tcga_tpm_mean', 'tcga_tpm_median', 'tcga_tpm_max', 'tcga_tpm_iqr', 'tcga_tpm_percentile'])
               new_output_lines.append(line)
           else:
               if line[col_idx_map['antigen_source']] in ['SNV', 'InDel']:
                   tx_id = line[col_idx_map['transcript_identifier']].split('.')[0]
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


def check_erv_rna_coverage(args):
    """
    """
    rna_covered_hervs = []
    expressed_hervs = {}
    out_header = ''
    with open(args.expressed_hervs) as fo:
        for line_idx, line in enumerate(fo.readlines()):
            if line_idx == 0:
                out_header = line
            else:
                line = line.rstrip().split(',')
                expressed_hervs[line[0]] = ','.join(line[1:])
    print(expressed_hervs)


    col_map = {}
    with open(args.coverage_file) as fo:
        for line_idx, line in enumerate(fo.readlines()):
            if line_idx == 0:
                line = line.strip().split('\t')
                for col_idx, col in enumerate(line):
                    col_map[col] = col_idx
            else:
                line = line.strip().split('\t')
                samp = line[col_map['#rname']]
                if samp in expressed_hervs.keys():
                    print(line)
                    mean_depth = line[col_map['meandepth']]
                    coverage = line[col_map['coverage']]
                    print("okay!")
                    print("{}\t{}\t{}\t{}".format(mean_depth, args.mean_depth, coverage, args.coverage))
                    if float(coverage) > float(args.coverage) and float(mean_depth) > float(args.mean_depth):
#                        rna_covered_hervs.append("{}\t{}".format(samp, line[col_map['numreads']]))
                        rna_covered_hervs.append("{},{}".format(samp, expressed_hervs[samp]))

    if rna_covered_hervs:
        with open(args.output, 'w') as fo:
            fo.write(out_header)
            for i in rna_covered_hervs:
                fo.write("{}\n".format(i))

def check_virus_rna_coverage(args):
    """
    """
    rna_covered_viruses = []
    expressed_viruses = []
    with open(args.expressed_viruses) as fo:
        for line_idx, line in enumerate(fo.readlines()):
            expressed_viruses.append(line.rstrip().split('\t')[0])
    print(expressed_viruses)


    col_map = {}
    with open(args.coverage_file) as fo:
        for line_idx, line in enumerate(fo.readlines()):
            if line_idx == 0:
                line = line.strip().split('\t')
                for col_idx, col in enumerate(line):
                    col_map[col] = col_idx
            else:
                line = line.strip().split('\t')
                samp = line[col_map['#rname']]
                if samp in expressed_viruses:
                    print(line)
                    mean_depth = line[col_map['meandepth']]
                    coverage = line[col_map['coverage']]
                    if float(coverage) > float(args.coverage) and float(mean_depth) > float(args.mean_depth):
                        rna_covered_viruses.append("{}\t{}".format(samp, line[col_map['numreads']]))

    if rna_covered_viruses:
        with open(args.output, 'w') as fo:
            for i in rna_covered_viruses:
                fo.write("{}\n".format(i))


def get_snv_peptide_read_count(args):
    """
    This accepts netmhcpan inputs, extracts the peptide and postiion within the
    contig (e.g. viral, self, or ERV coding sequence) and determines how many reads
    fall onto that region and what proportion of reads may code for that peptide.
    """
    rna_bam = pysam.AlignmentFile(args.bam, "rb")

    contig_seqs = {}

    md5_to_contig_map = {}

    md5_meta = {}

    tx_meta = {}

    
    for seq_record in SeqIO.parse(args.nt_fasta, "fasta"):
        id = seq_record.id
        print("ID: {}".format(id))
        print("Description: {}".format(seq_record.description))
        if re.search('MD5', id):
            print("Found MD5")
            id = seq_record.description
            split_id = id.split(' ')
            md5_meta[split_id[0].replace(':', '_')[:15]] = {}
            print("split id: {}".format(split_id))
            md5_meta[split_id[0].replace(':', '_')[:15]]['cdna_pos'] = split_id[3].replace('CDNA_POS:', '')
            md5_meta[split_id[0].replace(':', '_')[:15]]['tx'] = split_id[4].replace('TRANSCRIPT:', '')
            md5_meta[split_id[0].replace(':', '_')[:15]]['seq'] = seq_record.seq
#            id = split_id[0].replace(':', '_')
        elif re.search('ENST|ENSMUST', id):
            print("Found ENST|ENSMUST")
            id = seq_record.description
            split_id = id.split()
            print(split_id)
            tx_to_utr_buffer[split_id[0]] = split_id[2].replace('UTR_BUFFER:', '')
            contig_seqs[split_id[0]] = seq_record.seq
        else:
            id = id.partition(':')[0]
            contig_seqs[id] = seq_record.seq


    print(md5_meta)
    printable_lines = []

    txs = [md5_meta[x]['tx'] for x in md5_meta.keys()]

    print(txs)


    with open(args.gtf) as fo:
        for line in fo.readlines():
            for tx in txs:
                if re.search(tx, line):
                    if re.search('\tCDS\t', line):
                        strand = line.split('\t')[6] 
                        tx_meta[tx] = strand
                elif re.search(tx.split('.')[0], line):
                    if re.search('\tCDS\t', line):
                        strand = line.split('\t')[6] 
                        tx_meta[tx] = strand
    print("Loaded transcripts metadata.")
   
    header = [] 
    with open(args.netmhcpan) as fo:
        for line_idx, line in enumerate(fo.readlines()):
            line = line.rstrip().split(',')
            if len(line) > 16 and line[0] == 'Pos':
                header = line
            elif len(line) > 16 and line[0] not in ['Protein', 'Pos']:
                if 'SB' in line:
                    line.remove('<=')
                    line.remove('SB')
                if 'WB' in line:
                    line.remove('<=')
                    line.remove('WB')
                if float(line[15]) < 1000:
                    print("{}\t{}\t{}".format(line[0], line[2], line[10]))
                    lo_coord = int(line[0]) * 3
                    hi_coord = lo_coord + (int(len(line[2]) * 3))
                    print("{}\t{}".format(lo_coord, hi_coord))
                    matched_id = ''
                    if re.search('MD5_', line[10]):
                        print(line[10])
                        print(md5_to_contig_map)
                        matched_id = md5_meta[line[10]]['tx']
                    print("{}".format(matched_id))
                    peptide_nt_seq = str(md5_meta[line[10]]['seq'][lo_coord-3:hi_coord-3])
                    print(peptide_nt_seq)
                    print(lo_coord)
                    print(hi_coord)
                                                                                                    
                    lo_coord_gnm = int(md5_meta[line[10]]['cdna_pos']) - 30 + lo_coord - 4
                    hi_coord_gnm = int(md5_meta[line[10]]['cdna_pos']) - 30 + hi_coord - 4

                    print("{}\t{}".format(lo_coord_gnm, hi_coord_gnm))
                    if lo_coord_gnm < 0:
                        lo_coord_gnm = 0

                    covered_reads = []
                    contig_id = matched_id
                    contig_id_no_suffix = contig_id.split('.')[0]
                    bam_contig_id = contig_id
                    if contig_id not in list(rna_bam.references) and contig_id_no_suffix not in list(rna_bam.references):
                        print("Contig not discovered in references.")
                        continue
                    if contig_id_no_suffix in list(rna_bam.references):
                        #Test further. Appears mouse reference doesn't use tx version.
                        bam_contig_id = contig_id_no_suffix
                    for read in rna_bam.fetch(bam_contig_id, lo_coord_gnm, hi_coord_gnm):
                        if read.get_overlap(lo_coord_gnm, hi_coord_gnm) >= len(peptide_nt_seq):
                            if bool(re.search(peptide_nt_seq, read.query_sequence)):
                                covered_reads.append(read)
                    if len(covered_reads) > 0:
                        for read in covered_reads:
                            print("{}\t{}".format(read.get_overlap(lo_coord_gnm, hi_coord_gnm), len(peptide_nt_seq)))
                            print("Read: {}".format(read))
                            print("Query Seq: {}".format(read.query_sequence))
                            print("Ref positions: {}".format(read.get_reference_positions()))
                            print("Variant Seq: {}".format(peptide_nt_seq))
                            if bool(re.search(peptide_nt_seq, read.query_sequence)):
                                print("Hit!")
#                    print("Number of full overlap reads: {}".format(len(covered_reads)))
                    positive_hits = 0
                    print(tx_meta.keys()[:10])
                    tx_strand = ''
                    if contig_id in tx_meta.keys(): 
                        tx_strand = tx_meta[contig_id]
                    else:
                        tx_strand = tx_meta[contig_id.split('.')[0]]
                    if tx_strand == '-':
                        peptide_nt_seq = str(Seq(peptide_nt_seq).reverse_complement())
                    if len(covered_reads) > 0:
                        for read in covered_reads:
                            if bool(re.search(peptide_nt_seq, read.query_sequence)):
                                positive_hits += 1
                        print("Positive hits: {}".format(positive_hits))
                        prop_var = positive_hits/(len(covered_reads) + 0.0)
                        if prop_var > 0:
                            line.extend([str(len(covered_reads)), str(positive_hits), str(prop_var)])
#                            printable_lines.append(line.extend([positive_hits, prop_var]))
                            printable_lines.append(line)
    print(printable_lines)

    header.extend(['total_reads_full_overlap', 'reads_with_peptide', 'proportion_full_overlap_reads_with_peptide'])
    try:
        header.remove('BindLevel')
    except:
        pass
    with open(args.output, 'w') as ofo:
        print(header)
        ofo.write('{}\n'.format('\t'.join(header)))
        for line in printable_lines:
            ofo.write('{}\n'.format('\t'.join(line)))

def get_indel_peptide_read_count(args):
    """
    This accepts netmhcpan inputs, extracts the peptide and postiion within the
    contig (e.g. viral, self, or ERV coding sequence) and determines how many reads
    fall onto that region and what proportion of reads may code for that peptide.
    """
    rna_bam = pysam.AlignmentFile(args.bam, "rb")

    contig_seqs = {}

    md5_to_contig_map = {}

    md5_meta = {}

    tx_meta = {}

    
    for seq_record in SeqIO.parse(args.nt_fasta, "fasta"):
        id = seq_record.id
        print("ID: {}".format(id))
        print("Description: {}".format(seq_record.description))
        if re.search('MD5', id):
            print("Found MD5")
            id = seq_record.description
            split_id = id.split(' ')
            md5_meta[split_id[0].replace(':', '_')[:15]] = {}
            print("split id: {}".format(split_id))
            md5_meta[split_id[0].replace(':', '_')[:15]]['cdna_pos'] = split_id[2].replace('CDNA_POS:', '')
            md5_meta[split_id[0].replace(':', '_')[:15]]['tx'] = split_id[3].replace('TRANSCRIPT:', '')
            md5_meta[split_id[0].replace(':', '_')[:15]]['seq'] = seq_record.seq
#            id = split_id[0].replace(':', '_')
        elif re.search('ENST|ENSMUST', id):
            print("Found ENST|ENSMUST")
            id = seq_record.description
            split_id = id.split()
            print(split_id)
            tx_to_utr_buffer[split_id[0]] = split_id[1].replace('UTR_BUFFER:', '')
            contig_seqs[split_id[0]] = seq_record.seq
        else:
            id = id.partition(':')[0]
            contig_seqs[id] = seq_record.seq


    print(md5_meta)
    printable_lines = []

    txs = [md5_meta[x]['tx'] for x in md5_meta.keys()]

    print(txs)


    with open(args.gtf) as fo:
        for line in fo.readlines():
            for tx in txs:
                if re.search(tx, line):
                    if re.search('\tCDS\t', line):
                        strand = line.split('\t')[6] 
                        tx_meta[tx] = strand
                elif re.search(tx.split('.')[0], line):
                    if re.search('\tCDS\t', line):
                        strand = line.split('\t')[6] 
                        tx_meta[tx] = strand
    print("Loaded transcripts metadata.")
#    with open(args.gtf) as fo:
#        for line in fo.readlines():
#            for tx in txs:
#                if re.search(tx, line):
#                    if re.search('\tCDS\t', line):
#                        strand = line.split('\t')[6] 
#                        tx_meta[tx] = strand
#    print("Loaded transcripts metadata.")
   
    header = [] 
    with open(args.netmhcpan) as fo:
        for line_idx, line in enumerate(fo.readlines()):
            line = line.rstrip().split()
            #print(line)
            if len(line) > 16 and line[0] == 'Pos':
                header = line
            elif len(line) > 16 and line[0] not in ['Protein', 'Pos']:
                if 'SB' in line:
                    line.remove('<=')
                    line.remove('SB')
                if 'WB' in line:
                    line.remove('<=')
                    line.remove('WB')
                if float(line[15]) < 1000:
                    print("{}\t{}\t{}".format(line[0], line[2], line[10]))
                    lo_coord = int(line[0]) * 3
                    hi_coord = lo_coord + (int(len(line[2]) * 3))
                    print("{}\t{}".format(lo_coord, hi_coord))
                    matched_id = ''
                    if re.search('MD5_', line[10]):
                        print(line[10])
                        print(md5_to_contig_map)
                        matched_id = md5_meta[line[10]]['tx']
                    print("{}".format(matched_id))
                    peptide_nt_seq = str(md5_meta[line[10]]['seq'][lo_coord-3:hi_coord-3])
                    print(peptide_nt_seq)
                    print(lo_coord)
                    print(hi_coord)
                                                                                                    
                    lo_coord_gnm = int(md5_meta[line[10]]['cdna_pos']) - 30 + lo_coord - 4
                    hi_coord_gnm = int(md5_meta[line[10]]['cdna_pos']) - 30 + hi_coord - 4

                    print("{}\t{}".format(lo_coord_gnm, hi_coord_gnm))

                    covered_reads = []
#                    contig_id = matched_id
#                    if contig_id not in list(rna_bam.references):
#                        print("Contig not discovered in references.")
#                        continue
#                    for read in rna_bam.fetch(contig_id, lo_coord_gnm, hi_coord_gnm):
                    contig_id = matched_id
                    contig_id_no_suffix = contig_id.split('.')[0]
                    bam_contig_id = contig_id
                    if contig_id not in list(rna_bam.references) and contig_id_no_suffix not in list(rna_bam.references):
                        print("Contig not discovered in references.")
                        continue
                    if contig_id_no_suffix in list(rna_bam.references):
                        #Test further. Appears mouse reference doesn't use tx version.
                        bam_contig_id = contig_id_no_suffix
                    for read in rna_bam.fetch(bam_contig_id, lo_coord_gnm, hi_coord_gnm):
                        if read.get_overlap(lo_coord_gnm, hi_coord_gnm) >= len(peptide_nt_seq):
                            covered_reads.append(read)
                    if len(covered_reads) > 0:
                        for read in covered_reads:
                            print("{}\t{}".format(read.get_overlap(lo_coord_gnm, hi_coord_gnm), len(peptide_nt_seq)))
                            print("Read: {}".format(read))
                            print("Query Seq: {}".format(read.query_sequence))
                            print("Ref positions: {}".format(read.get_reference_positions()))
                            print("Variant Seq: {}".format(peptide_nt_seq))
#                    print("Number of full overlap reads: {}".format(len(covered_reads)))
                    positive_hits = 0
                    tx_strand = ''
                    if contig_id in tx_meta.keys(): 
                        tx_strand = tx_meta[contig_id]
                    else:
                        tx_strand = tx_meta[contig_id.split('.')[0]]
                    #tx_strand = tx_meta[contig_id]
                    if tx_strand == '-':
                        peptide_nt_seq = str(Seq(peptide_nt_seq).reverse_complement())
                    if len(covered_reads) > 0:
                        for read in covered_reads:
                            if bool(re.search(peptide_nt_seq, read.query_sequence)):
                                positive_hits += 1
                        print("Positive hits: {}".format(positive_hits))
                        prop_var = positive_hits/(len(covered_reads) + 0.0)
                        if prop_var > 0:
                            line.extend([str(len(covered_reads)), str(positive_hits), str(prop_var)])
#                            printable_lines.append(line.extend([positive_hits, prop_var]))
                            printable_lines.append(line)
    print(printable_lines)

    header.extend(['total_reads_full_overlap', 'reads_with_peptide', 'proportion_full_overlap_reads_with_peptide'])
    try:
        header.remove('BindLevel')
    except:
        pass
    with open(args.output, 'w') as ofo:
        print(header)
        ofo.write('{}\n'.format('\t'.join(header)))
        for line in printable_lines:
            ofo.write('{}\n'.format('\t'.join(line)))



def get_peptide_read_count(args):
    """
    This accepts netmhcpan inputs, extracts the peptide and postiion within the
    contig (e.g. viral, self, or ERV coding sequence) and determines how many reads
    fall onto that region and what proportion of reads may code for that peptide.
    """
    rna_bam = pysam.AlignmentFile(args.bam, "rb")

    contig_seqs = {}

    md5_to_contig_map = {}

    tx_to_utr_buffer = {}
    
    for seq_record in SeqIO.parse(args.cds_fasta, "fasta"):
        id = seq_record.id
        print(seq_record.description)
        if re.search('MD5', id):
            id = seq_record.description
            split_id = id.split(' ')
            print(split_id)
            md5_to_contig_map[split_id[0].replace(':', '_')[:15]] = split_id[1].replace('NAME:', '')
            if re.search('UTR_BUFFER', id):
                tx_to_utr_buffer[split_id[1].replace('NAME:','')] = split_id[2].replace('UTR_BUFFER:', '')
            id = split_id[0].replace(':', '_')
#            contig_seqs[id[:15]] = seq_record.seq
            contig_seqs[split_id[1].replace('NAME:', '')] = seq_record.seq
        elif re.search('ENST|ENSMUST', id):
            id = seq_record.description
            split_id = id.split()
            print(split_id)
            tx_to_utr_buffer[split_id[0]] = split_id[1].replace('UTR_BUFFER:', '')
            contig_seqs[split_id[0]] = seq_record.seq
        else:
            id = id.partition(':')[0]
            contig_seqs[id] = seq_record.seq


    print(md5_to_contig_map)
    printable_lines = []
   
    header = [] 
    with open(args.netmhcpan) as fo:
        for line_idx, line in enumerate(fo.readlines()):
            line = line.rstrip().split()
            if len(line) > 16 and line[0] == 'Pos':
                header = line
            elif len(line) > 16 and line[0] not in ['Protein', 'Pos']:
                if 'SB' in line:
                    line.remove('<=')
                    line.remove('SB')
                if 'WB' in line:
                    line.remove('<=')
                    line.remove('WB')
                if float(line[15]) < 1000:
                    print("{}\t{}\t{}".format(line[0], line[2], line[10]))
                    lo_coord = int(line[0]) * 3
                    hi_coord = lo_coord + (int(len(line[2]) * 3))
                    print("{}\t{}".format(lo_coord, hi_coord))
                    print(contig_seqs.keys())
                    matched_id = ''
                    if re.search('MD5_', line[10]):
                        print(line[10])
                        print(md5_to_contig_map)
                        matched_id = md5_to_contig_map[line[10]]
                    else:
                        print(line[10])
                        print(contig_seqs.keys())
                        matched_id = [x for x in contig_seqs.keys() if re.search(line[10], x.replace('.', '_'))][0]
#                    matched_id = [x for x in contig_seqs.keys() if re.search(line[10], x.replace('.', '_'))][0]
                    if re.search('MD5_', matched_id):
                        matched_id = matched_id[:15]
                    print("{}".format(matched_id))
                    peptide_nt_seq = str(contig_seqs[matched_id][lo_coord-3:hi_coord-3])
                    print(peptide_nt_seq)
                    print(lo_coord)
                    print(hi_coord)
                    covered_reads = []
                    contig_id = matched_id
                    if re.search('MD5_', matched_id):
                        contig_id = md5_to_contig_map[line[10]] 
                    if re.match('ENST|ENSMUST', matched_id):
                        contig_id = [x for x in rna_bam.references if re.search(matched_id, x)][0]
                        lo_coord = lo_coord + int(tx_to_utr_buffer[matched_id])
                        hi_coord = hi_coord + int(tx_to_utr_buffer[matched_id])
                    for read in rna_bam.fetch(contig_id, lo_coord-3, hi_coord-3):
                        print("{}\t{}".format(read.get_overlap(lo_coord-3, hi_coord-3), len(peptide_nt_seq)))
                        if read.get_overlap(lo_coord-3, hi_coord-3) == len(peptide_nt_seq):
                            print(read)
                            print(read.query_sequence)
                            print(read.get_reference_positions())
                            print(peptide_nt_seq)
                            covered_reads.append(read.query_sequence)
                            if bool(re.search(peptide_nt_seq, read.query_sequence)):
                                print("hit!")
#                                print(read.query_sequence)
#                                print(peptide_nt_seq)
                    print("Number of full overlap reads: {}".format(len(covered_reads)))
                    positive_hits = 0
                    if len(covered_reads) > 0:
                        for read in covered_reads:
#                            print(peptide_nt_seq)
#                            print(read)
#                            print(re.search(peptide_nt_seq, read))
                            if bool(re.search(peptide_nt_seq, read)):
                                positive_hits += 1
                        print("Positive hits: {}".format(positive_hits))
                        prop_var = positive_hits/(len(covered_reads) + 0.0)
                        if prop_var > 0:
                            line.extend([str(len(covered_reads)), str(positive_hits), str(prop_var)])
#                            printable_lines.append(line.extend([positive_hits, prop_var]))
                            printable_lines.append(line)
    print(printable_lines)

    header.extend(['total_reads_full_overlap', 'reads_with_peptide', 'proportion_full_overlap_reads_with_peptide'])
    header.remove('BindLevel')
    with open(args.output, 'w') as ofo:
        print(header)
        ofo.write('{}\n'.format('\t'.join(header)))
        for line in printable_lines:
            ofo.write('{}\n'.format('\t'.join(line)))
                        

def get_fusion_peptide_read_count(args):
    """
    This accepts netmhcpan inputs, extracts the peptide and postiion within the
    contig (e.g. viral, self, or ERV coding sequence) and determines how many reads
    fall onto that region and what proportion of reads may code for that peptide.
    """

    contig_seqs = {}

    for seq_record in SeqIO.parse(args.nt_fasta, "fasta"):
        id = '{}'.format(seq_record.description[:15].replace(' ', '_').replace('.', '_')) #The _NUCLEOTIDE addition is required to match to the NETMHCpan outputs.
        contig_seqs[id] = seq_record.seq

    print("contig_seqs keys")
    print(contig_seqs.keys())


    reads_and_seqs = {}
    for seq_record in SeqIO.parse(args.fusion_reads, "fastq"):
        print("NEW RECORD")
        new_id = ''
        if seq_record.id.endswith('/1'):
            new_id = seq_record.id.strip('/1')
        elif seq_record.id.endswith('/2'):
            new_id = seq_record.id.strip('/2')
        print(seq_record.id)
        print(new_id)
        print("{}\t{}".format(seq_record.id, new_id))
        print(seq_record.seq)
        if new_id not in reads_and_seqs.keys():
            reads_and_seqs[new_id] = [str(seq_record.seq)]
        else:
            reads_and_seqs[new_id].append(str(seq_record.seq))

    fusion_reads = {}
    with open(args.fusions) as fo:
        for line_idx, line in enumerate(fo.readlines()):
            if line_idx > 0:
                line = line.split()
                id = '{}_NUCLEOTIDE'.format(line[0][:15].replace(' ', '_').replace('.', '_'))[:15]
                reads = ['{}'.format(x.partition('@')[2]) for x in line[8].split(',')]
                reads.extend(['{}'.format(x.partition('@')[2]) for x in line[8].split(',')])
                reads.extend(['{}'.format(x.partition('@')[2]) for x in line[9].split(',')])
                reads.extend(['{}'.format(x.partition('@')[2]) for x in line[9].split(',')])
                fusion_reads[id] = reads
    print("fusion_reads keys")
    print(fusion_reads.keys())

    print(fusion_reads)
    print(reads_and_seqs)

    printable_lines = []
   
    header = [] 
    with open(args.netmhcpan) as fo:
        for line_idx, line in enumerate(fo.readlines()):
            line = line.rstrip().split()
            if len(line) > 16 and line[0] == 'Pos':
                header = line
            elif len(line) > 16 and line[0] not in ['Protein', 'Pos']:
                if 'SB' in line:
                    line.remove('<=')
                    line.remove('SB')
                if 'WB' in line:
                    line.remove('<=')
                    line.remove('WB')
                if float(line[15]) < 1000:
                    print("{}\t{}\t{}".format(line[0], line[2], line[10]))
                    lo_coord = int(line[0]) * 3
                    hi_coord = lo_coord + (int(len(line[2]) * 3))
                    print("{}\t{}".format(lo_coord, hi_coord))
#                    print(contig_seqs.keys())
                    matched_id = line[10].replace('.', '_')
#                    if re.search('MD5_', line[10]):
#                        print(line[10])
#                        print(md5_to_contig_map)
#                        matched_id = md5_to_contig_map[line[10]]
#                    else:
#                        print(line[10])
#                        print(contig_seqs.keys())
#                        matched_id = [x for x in contig_seqs.keys() if re.search(line[10], x.replace('.', '_'))][0]
##                    matched_id = [x for x in contig_seqs.keys() if re.search(line[10], x.replace('.', '_'))][0]
#                    if re.search('MD5_', matched_id):
#                        matched_id = matched_id[:15]
                    print("{}".format(contig_seqs[matched_id]))
                    peptide_nt_seq = str(contig_seqs[matched_id][lo_coord-3:hi_coord-3]).upper()
                    print(peptide_nt_seq)
                    print(lo_coord)
                    print(hi_coord)
                    covered_reads = []
                    reads_with_peptide_nt = 0
                    potential_reads = fusion_reads[matched_id]
                    for read in potential_reads:
                        reads = list(set(reads_and_seqs[read]))
                        print(reads)
                        for subread in reads:
                           if re.search(peptide_nt_seq, subread) or re.search(peptide_nt_seq, str(Seq(subread).reverse_complement())):
                                reads_with_peptide_nt += 1
                                print(reads_with_peptide_nt)
                    print("Potential reads: {}".format(len(potential_reads)))
                    print("Reads with peptide: {}".format(reads_with_peptide_nt))
                    if reads_with_peptide_nt > 0:
                        line.extend(['NA', str(reads_with_peptide_nt), 'NA'])
                        printable_lines.append(line)
    print(printable_lines)

    header.extend(['total_reads_full_overlap', 'reads_with_peptide', 'proportion_full_overlap_reads_with_peptide'])
    if 'BindLevel' in header:
        header.remove('BindLevel')
    with open(args.output, 'w') as ofo:
        print(header)
        ofo.write('{}\n'.format('\t'.join(header)))
        for line in printable_lines:
            ofo.write('{}\n'.format('\t'.join(line)))


def get_splice_peptide_read_count(args):
    """
    This accepts netmhcpan inputs, extracts the peptide and postiion within the
    contig (e.g. viral, self, or ERV coding sequence) and determines how many reads
    fall onto that region and what proportion of reads may code for that peptide.
    """
    rna_bam = pysam.AlignmentFile(args.bam, "rb")

    header = ''

    printable_lines = []

    with open(args.neosplice_summary) as fo:
        for line_idx, line in enumerate(fo.readlines()):
            line = line.rstrip().split('\t')
            if line_idx == 0:
                header = line
            else:
                dna_seq = line[1]
                print(dna_seq)
                chr = line[2]
                strand = line[7]
                coord_lo = int(line[3].split('), (')[0].split(', ')[0].replace('(', ''))
                coord_hi = int(line[3].split('), (')[-1].split(',')[1].replace(' ',''))
                print(coord_lo)
                print(coord_hi)
                if strand == '+':
                    covered_reads = []
                    positive_hits = 0
                    for read in rna_bam.fetch(chr, coord_lo, coord_hi):
#                        print("{}\t{}".format(read.get_overlap(coord_lo, coord_hi), len(dna_seq)))
                        # greater than or equal to instead of equal to due to
                        # looking at a larger range than actual DNA sequence
                        if read.get_overlap(coord_lo, coord_hi) >= len(dna_seq):
                            covered_reads.append(read.query_sequence)
                    if len(covered_reads) > 0:
                        for read in covered_reads:
                            if bool(re.search(dna_seq, read)):
                                print("Peptide NT: {}".format(dna_seq))
                                print("Read: {}".format(read))
                                print(re.search(dna_seq, read))
                                positive_hits += 1
                        print("Positive hits: {}".format(positive_hits))
                        prop_var = positive_hits/(len(covered_reads) + 0.0)
                        if prop_var > 0.01:
                            print(positive_hits)
                            print(prop_var)
                elif strand == '-':
                    covered_reads = []
                    positive_hits = 0
                    coord_lo, coord_hi = coord_hi, coord_lo
                    print("{}\t{}".format(coord_lo, coord_hi))
                    for read in rna_bam.fetch(chr, coord_lo, coord_hi):
#                        print("{}\t{}".format(read.get_overlap(coord_lo, coord_hi), len(dna_seq)))
                        # greater than or equal to instead of equal to due to
                        # looking at a larger range than actual DNA sequence
                        if read.get_overlap(coord_lo, coord_hi) >= len(dna_seq):
                            covered_reads.append(read.query_sequence)
                    if len(covered_reads) > 0:
                        for read in covered_reads:
                            #print(read)
                            #print(Seq(read).reverse_complement())
                            if bool(re.search(dna_seq, str(Seq(read).reverse_complement()))):
                                #print(dna_seq)
                                #print(read)
                                #print(re.search(dna_seq, str(Seq(read).reverse_compement())))
                                print("Peptide NT: {}".format(dna_seq))
                                print("Read: {}".format(read))
                                print(re.search(dna_seq, read))
                                positive_hits += 1
                        print("Positive hits: {}".format(positive_hits))
                        prop_var = positive_hits/(len(covered_reads) + 0.0)
                        if prop_var > 0.01:
                            print(positive_hits)
                            print(prop_var)
                if positive_hits > 0:
                    line.extend(['NA', str(positive_hits), 'NA'])
                    printable_lines.append(line)
    print(printable_lines)

    header.extend(['total_reads_full_overlap', 'reads_with_peptide', 'proportion_full_overlap_reads_with_peptide'])
    with open(args.output, 'w') as ofo:
        print(header)
        ofo.write('{}\n'.format('\t'.join(header)))
        for line in printable_lines:
            ofo.write('{}\n'.format('\t'.join(line)))


def iupac_conversion(base):
    conversion = {'R': ['A', 'G'],
                  'Y': ['C', 'T'],
                  'S': ['G', 'C'],
                  'W': ['A', 'T'],
                  'K': ['G', 'T'],
                  'M': ['A', 'C']}
    if base in conversion.keys():
        return conversion[base]
    else:
        return []


                
                

def main():
    args = get_args()
    #print(args)
    if args.command == 'make-snv-peptides':
        make_snv_peptides(args)
    if args.command == 'make-snv-peptides-context':
        make_snv_peptides_context(args)
    if args.command == 'make-indel-peptides':
        make_indel_peptides(args)
    if args.command == 'make-indel-peptides-context':
        make_indel_peptides_context(args)
    if args.command == 'make-fusion-peptides-context':
        make_fusion_peptides_context(args)
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
    if args.command == 'filter-expressed-ervs':
        expressed_hervs(args)
    if args.command == 'get-expressed-ervs-bed':
        get_expressed_ervs_bed(args)
    if args.command == 'get-expressed-viral-bed':
        get_expressed_viral_bed(args)
    if args.command == 'get-expressed-selfs-bed':
        get_expressed_selfs_bed(args)
    if args.command == 'make-erv-peptides':
        make_erv_peptides(args)
    if args.command == 'add-erv-metadata':
        add_erv_metadata(args)
    if args.command == 'filter-expressed-self-genes':
        expressed_self_genes(args)
    if args.command == 'make-self-antigen-peptides':
        make_self_antigen_peptides(args)
    if args.command == 'add-self-antigen-metadata':
        add_self_antigen_metadata(args)
    if args.command == 'filter-expressed-viruses':
        filter_virdetect_by_counts(args)
    if args.command == 'filter-expressed-viral-cds':
        filter_viral_cds(args)
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
    if args.command == 'check-erv-rna-coverage':
        check_erv_rna_coverage(args)
    if args.command == 'check-virus-rna-coverage':
        check_virus_rna_coverage(args)
    if args.command == 'get-peptide-read-count':
        get_peptide_read_count(args)
    if args.command == 'get-fusion-peptide-read-count':
        get_fusion_peptide_read_count(args)
    if args.command == 'get-splice-peptide-read-count':
        get_splice_peptide_read_count(args)
    if args.command == 'get-snv-peptide-read-count':
        get_snv_peptide_read_count(args)
    if args.command == 'get-indel-peptide-read-count':
        get_indel_peptide_read_count(args)
    if args.command == 'filter-mutant-peptides':
        filter_mutant_peptides(args)
    if args.command == 'get-expressed-transcripts-bed':
        get_expressed_transcripts_bed(args)


if __name__=='__main__':
    main()
