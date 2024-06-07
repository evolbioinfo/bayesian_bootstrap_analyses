#!/usr/bin/env python3
# coding: utf-8

import sys, getopt
import csv
import logging
import gzip
import os
import lzma
import io
import subprocess
from Bio import SeqIO


def main(argv):
    cpus=1
    try:
        opts, args = getopt.getopt(argv,"hi:n:o:m:t:")
    except getopt.GetoptError:
        print('selectByEPIISL.py -i <sequences.fa> -n <epi isl ids> -m <metadata> -o <outfile> -t <cpus>')
        sys.exit(1)
    for opt, arg in opts:
        if opt == '-h':
            print('selectByEPIISL.py -i <sequences.fa> -n <epi isl ids> -m <metadata> -o <outfile> -t <cpus>')
            sys.exit(0)
        elif opt in ("-i"):
            inputseq = arg
        elif opt in ("-m"):
            inputmeta = arg
        elif opt in ("-o"):
            outfile = arg
        elif opt in ("-t"):
            cpus = arg
        elif opt in ("-n"):
            epiisl = arg

    logging.info(f'Parsing epiisl file {epiisl}')
    epilist=dict()
    with open(epiisl) as file:
        tsv_file = csv.reader(file, delimiter="\t")
        for line in tsv_file:        
            epilist[line[0]]=1
            
    logging.info(f'Parsing metadata file {inputmeta}')
    # We parse the metadata file for getting the gisaid name of the epiisl ids
    proc = subprocess.Popen(["sh","-c", f"tar -x -O -I 'xz -T{cpus}' -f {inputmeta} metadata.tsv"], stdout=subprocess.PIPE)
    tsv_file = csv.reader(io.TextIOWrapper(proc.stdout), delimiter="\t")
    gisaidlist=dict()
    next(tsv_file,None)
    for line in tsv_file:
        id=line[0]
        epiid=line[2]
        if epiid in epilist:
            gisaidlist[id]=epiid

    logging.info(f'Parsing sequence file {inputseq}')
    # We then parse the sequence file to print each sequence in the proper output file

    with open(outfile,'wt') as outf:

        proc = subprocess.Popen(["sh","-c", f"tar -x -O -I 'xz -T{cpus}' -f {inputseq} sequences.fasta"], stdout=subprocess.PIPE)
        for record in SeqIO.parse(io.TextIOWrapper(proc.stdout), "fasta"):
            id=record.description.split("|")[0]
            if id in gisaidlist:
                outf.write(f'>{id}|{gisaidlist[id]}\n')
                outf.write(f'{record.seq}\n')


if __name__ == "__main__":
    logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))
    main(sys.argv[1:])
