#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Gets the sequence from a fasta file (warning : the file must contain only one sequence)
def readFasta(filename):
    filin=open(filename, "r")
    lnum=0
    nb_ref=0
    seq=""
    for line in filin:
        if lnum>0:
            seq+=line.rstrip("\n")
        lnum+=1
    filin.close()
    return seq
