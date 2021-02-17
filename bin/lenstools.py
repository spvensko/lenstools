#!/usr/bin/env python

import argparse
import vcf
from Bio import SeqIO
from Bio.SeqUtils import seq1
import re
import sys
from pprint import pprint
import numpy as np
from scipy import stats
import csv
import pysam
import hashlib

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
                snapshot = pos_total.split('/')[0]
                tlen = pos_total.split('/')[1] 
                if pos != snapshot:
                    sys.exit(1)
                filtered_records.append([transcript, pos, tlen, pre1, post1, chr_pos])
                
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
                ref_peptide = ''.join(tx_ref_seq[max(pos-8,0):min(pos+8, len(mut_seq))])
                mut_peptide = ''.join(mut_seq[max(pos-8, 0):min(pos+8, len(mut_seq))])
                print("{}\n{}\n{}".format(record[0],ref_peptide, mut_peptide))
                var_md5sum = hashlib.md5("{}".format(':'.join([str(record[-1].CHROM), str(record[-1].POS), str(record[0]), str(record[-1].REF), str(record[-1].ALT)]))).hexdigest()[:16]
                emitted_peptides["{} {}:{} {} {} {}".format(var_md5sum, record[-1].CHROM, record[-1].POS, record[0], record[-1].REF, record[-1].ALT)] = mut_peptide
                wildtype_peptides["{} {}:{} {} {} {}".format(var_md5sum, record[-1].CHROM, record[-1].POS, record[0], record[-1].REF, record[-1].ALT)] = ref_peptide

    with open(args.mt_output, 'w') as ofo:
        for k, v in emitted_peptides.items():
            ofo.write('>{}\n{}\n'.format(k, v))
    with open(args.wt_output, 'w') as ofo:
        for k, v in wildtype_peptides.items():
            ofo.write('>{}\n{}\n'.format(k, v))

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
    tx_threhsold = 0
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
            if np.absolute(upstream_record.POS - record.POS) > args.proximity:
                upstream_isolated = True
        else:
            upstream_isolated = True

        if downstream_record.CHROM == record.CHROM:
            if np.absolute(downstream_record.POS - record.POS) > args.proximity:
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
    with open(args.wt_fasta) as wto:
        for line in wto.readlines():
            line = line.rstrip('\n').split(' ')
            line = [x for x in line if x]
            if len(line) > 14 and line[0] not in ['Pos', 'Protein']:
#                print(line)
                k = "{}_{}_{}_{}".format(line[1], line[10], line[0], len(line[9]))
#                print(k)
                wt_nms[k] = line[15]
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
                    result.extend([wt_nms[m], str(float(line[15])/float(wt_nms[m]))])
                else:
                    result.extend(["NA", "NA"])
                mts_w_agreto.append(','.join(result))
            elif len(line) > 14 and line[0] == 'Pos':
                print(header)
                line.extend(['Wildtype Aff(nM)', 'Agretopocity'])
                header = ','.join(line)

    with open(args.output, 'w') as outf:
        outf.write("{}\n".format(header))
        for line in mts_w_agreto:
            outf.write("{}\n".format(line))



def create_lens_report(args):
    """
    """
    pass


def add_metadata(args):
    """
    """
    checksum_to_meta_map = {}
    with open(args.mt_aa) as mto:
        for line in mto.readlines():
            if line.startswith('>'):
                line = line.split(' ')
                checksum_to_meta_map[line[0]] = line[1]

    tx_to_tpm = {}
    with open(args.quant_file) as qfo:
        tpm_col_idx = ''
        tx_col_idx = ''
        for line_idx, line in qfo.readlines():
            if line_idx == 0:
                tpm_col_idx = line.split('\t').index('TPM')
                tx_col_idx = line.split('\t').index('Name')
            else:     
                line = line.split(' ')
                tx_to_tpm[line[tx_col_idx]] = line[tpm_col_idx]

    new_lines = []
    with open(args.mut_nmp) as mno:
        for line in mno.readlines():
            line = line.split('\t')
            print(line)

    with open(args.output) as ofo:
        pass


def make_pyclone_vi_inputs(args):                                                                   
    """                                                                                             
    This requires an isec vcf, a proper vcf (with depth info), and sequenza results info.           
    """                                                                                             
    # This dictionary will be populated with normal and tumor depth information                     
    # from the proper VCFs. It'll initially have keys populated by the isec VCF.                    
    # vars[var]['var_depth'] =, vars[var]['ref_depth'] =, vars[var]['cn'] =                         
    vars = {}                                                                                       
                                                                                                    
                                                                                                    
    with open(args.candidate_vcf) as cvo:                                                           
        for line in cvo.readlines():                                                                
            if line.startswith('#'):                                                                
                pass                                                                                
            else:                                                                                   
                line = line.rstrip().split('\t')                                                    
                if len(line[3]) == 1 and len(line[4]) == 1:                                         
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
                pyclone_inp.append("{}\n".format('\t'.join([var, args.samp_id, ref_depth, alt_depth, ref_cn, alt_cn, '2', '0.001'])))

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

if __name__=='__main__':
    main()
