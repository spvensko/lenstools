#!/usr/bin/env python

import argparse
import vcf
from Bio import SeqIO
from Bio.SeqUtils import seq1
import re
import sys
from pprint import pprint

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
    parser_make_peptides.add_argument('-o', '--output',
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
    vcf_reader = vcf.Reader(open(args.vcf), 'r', compressed=True)
    for record in vcf_reader:
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

def make_peptides(args):
    """
    """
    tx_to_aa = load_tx_aas(args)
    print(tx_to_aa['ENST00000369529.1'])
    filtered_records = extract_potential_somatic_vars(args)
    emitted_peptides = {}
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

    with open(args.output, 'w') as ofo:
        for k, v in emitted_peptides.items():
            ofo.write('>{}\n{}\n'.format(k, v))
   


def main():
    args = get_args()
    #print(args)
    if args.command == 'make-peptides':
        make_peptides(args)    

if __name__=='__main__':
    main()
