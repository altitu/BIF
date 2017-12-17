#!/usr/bin/python2
# -*- coding: utf-8 -*-

# Partie 1 : extraire les k-mers du read et les match avec le génome

import bwt
import tools_karkkainen_sanders as tks
import seedandextend as sae
from filemanager import openFasta
import argparse

def main(refFilename, readsFilename, k, dmax, verb, outputconsole, debug):

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
	print "bwt du génome:" + str(b)
	n = bwt.getN(b)
	ranks = bwt.buildRankArray(b)
#	respos = bwt.is_Q_in_S(b, bwt.GET_N(b), sa, "ACC") #savoir s'il est présent dans la séquence
#	respos = bwt.is_Q_in_S(b, bwt.GET_N(b), sa, "ATGATCAG") #savoir s'il est présent dans la séquence
#	print  "position de match parfait obtenus: " + str(respos)

#	#data = sae.FastaFile()
#	test0 = sae.extends(respos, 4, "ACC", "CCGTACCCC", s, 0, 1, -10)
#	test0 = sae.extends(respos, 0, "ATGATCA", "ATGATCAG", s, 0, 1, -10)
#	print test0

	sae.distributeReads(reads, k, dmax, s, b, sa, n, ranks)
	#dmax = len(read) - max acceptable

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description=
		"Finds the position and differences of reads inside a reference genome")

	parser.add_argument('-rf','--ref', help="Filename to reference genome")
	parser.add_argument('-rds','--reads', help="Filename to reads")
	parser.add_argument('-o','--output', help="Filename to output the result", default="./output.txt")
	parser.add_argument('-k','--kparam', help="K-mer size, default value = 20", default = 20, type=int)
	parser.add_argument('-dmax','--dmaxparam', help="Maximum number of allowed substitutions in a match, default value = 5", default=5, type=int)
	parser.add_argument('-v','--verbosity', action="count", default=0, help="Increase verbosity level")
	parser.add_argument('-oc','--outputconsole', default=True, help="Print also the output to the screen", type=bool)
	parser.add_argument('-db','--debug', default=False, help="help tracking k-mers extension and positions outputs of reads", type=bool) 

	args = parser.parse_args()

#	if args.k < 1:
#		print "Invalid parameter k"
#		exit(-1)
	error = 0
	if args.reads == None or args.ref == None:
		print "error: option ref and reads are mandatory"
		error += 1
	if args.kparam < 1:
		print "error: values <1 for -k won't produce any k-mer"
		error += 1
	if args.dmaxparam < 0:
		print "error: values <0 for -d won't match in any case"
		error += 1

	if error > 0:
		print "\n{} error(s) in total".format(error)
		exit(-1)

	main(args.ref, args.reads, args.kparam, args.dmaxparam, args.verbosity, args.outputconsole, args.debug)
