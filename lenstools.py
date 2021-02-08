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
    parser_expressed_variants.add_argument('-v', '--vcf',
                                      help="Annotated VCF).",
                                      required=True)
    parser_expressed_variants.add_argument('-o', '--output',
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
#            if re.search("\|conserative_inframe_insertion\|", possible_variant) or re.search("\|conserative_inframe_deletion\|", possible_variant):
            split_possible = possible_variant.split('|')
            if re.search('^conservative_inframe_[a-z]+$', split_possible[1]) and split_possible[13] and split_possible[-1] == '':
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
    print(tx_to_aa['ENST00000369529.1'])
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
                emitted_peptides["{} {} {} {} {}".format(record[0], record[-1].CHROM, record[-1].POS, record[-1].REF, record[-1].ALT)] = mut_peptide
                wildtype_peptides["{} {} {} {} {}".format(record[0], record[-1].CHROM, record[-1].POS, record[-1].REF, record[-1].ALT)] = ref_peptide

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
    print(tx_to_aa['ENST00000369529.1'])
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
                    emitted_peptides["{} {} {} {} {}".format(record[0], record[-1].CHROM, record[-1].POS, record[-1].REF, record[-1].ALT)] = mut_peptide

    with open(args.output, 'w') as ofo:
        for k, v in emitted_peptides.items():
            ofo.write('>{}\n{}\n'.format(k, v))

def expressed_variants(args):
    """
    """
    tx_abundances = load_tx_abundances(args)
    tx_threshold = get_tx_threshold(args, tx_abundances)
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
                if count > threshold:
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
                


def main():
    args = get_args()
    #print(args)
    if args.command == 'make-peptides':
        make_peptides(args)    
    if args.command == 'make-indel-peptides':
        make_indel_peptides(args)    
    if args.command == 'expressed-variants':
        expressed_variants(args)   
    if args.command == 'isolated-variants':
        isolated_variants(args) 
    if args.command == 'calculate-agretopicity':
        calculate_agretopicity(args)

if __name__=='__main__':
    main()
