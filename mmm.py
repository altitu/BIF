#!/usr/bin/python2
# -*- coding: utf-8 -*-

# Partie 1 : extraire les k-mers du read et les match avec le génome

import bwt
import tools_karkkainen_sanders as tks
import seedandextend as sae
from filemanager import openFasta
import argparse

def main(refFilename, readsFilename, k, dmax, psa, pr):

#	s = "GTATGATCAGAA$"
#	sa = tks.simple_kark_sort(s)
#	b = bwt.GET_BWT(s, sa)
#	print b
#	respos = bwt.is_Q_in_S(b, bwt.GET_N(b), sa, "AT")
#	print respos

	#s = "ACCCCGTACCCCGTACCCC$"
	genome = openFasta(refFilename, 0) #ok
	reads = openFasta(readsFilename, 0) #ok
#	genome = [['commentaire random'],['GCATGCTTTTGCCGAT']]
#	reads = [['commentaire random'],['ATGC','TTGC']]
	print "resultat des reads:"
	print reads[1]
	s = genome[1][0] + "$"
	print "et pour le genome:"
	print s
	sa = tks.simple_kark_sort(s)
	b = bwt.getBWT(s, sa)
	sa = bwt.subsampleArray(sa, psa)
	print "bwt du génome:" + str(b)
	n = bwt.getN(b)
	ranks = bwt.buildRankArray(b)
	ranks = bwt.subsampleArray(ranks, pr)
#	respos = bwt.is_Q_in_S(b, bwt.GET_N(b), sa, "ACC") #savoir s'il est présent dans la séquence
#	respos = bwt.is_Q_in_S(b, bwt.GET_N(b), sa, "ATGATCAG") #savoir s'il est présent dans la séquence
#	print  "position de match parfait obtenus: " + str(respos)

#	#data = sae.FastaFile()
#	test0 = sae.extends(respos, 4, "ACC", "CCGTACCCC", s, 0, 1, -10)
#	test0 = sae.extends(respos, 0, "ATGATCA", "ATGATCAG", s, 0, 1, -10)
#	print test0

	sae.distributeReads(reads, k, dmax, s, b, sa, psa, n, ranks, pr)
	#dmax = len(read) - max acceptable

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description=
		"Finds the position and differences of reads inside a reference genome")

	parser.add_argument('ref', help="Filename to reference genome")
	parser.add_argument('reads', help="Filename to reads")
	parser.add_argument('k', help="K-mer size", type=int)
	parser.add_argument('dmax', help="Maximum number of allowed substitutions in a match", type=int)
	parser.add_argument('--psa', help="Step between elements in the suffix array", type=int, default=1)
	parser.add_argument('--pr', help="Step between elements in the rank array", type=int, default=1)

	args = parser.parse_args()

	if args.k < 1:
		print "Invalid parameter k"
		exit(-1)

	if args.psa < 1:
		print "Invalid parameter psa"
		exit(-1)

	if args.pr < 1:
		print "Invalid parameter pr"
		exit(-1)

	main(args.ref, args.reads, args.k, args.dmax, args.psa, args.pr)
