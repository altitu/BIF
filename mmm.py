#!/usr/bin/python2
# -*- coding: utf-8 -*-

import bwt
import tools_karkkainen_sanders as tks
import seedandextend as sae
from filemanager import openFasta
import argparse


def main(refFilename, readsFilename, k, dmax, psa, pr, verb, outputconsole, debug):

	genome = openFasta(refFilename, 0) 
	reads = openFasta(readsFilename, 0)
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

	sae.distributeReads(reads, k, dmax, s, b, sa, psa, n, ranks, pr)
	#dmax = len(read) - max acceptable

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description=
		"Finds the position and differences of reads inside a reference genome")

	parser.add_argument('rf', help="Filename to reference genome")
	parser.add_argument('rds', help="Filename to reads")
	parser.add_argument('-o','--output', help="Filename to output the result", default="./output.txt")
	parser.add_argument('-k','--kparam', help="K-mer size, default value = 20", default = 20, type=int)
	parser.add_argument('-dmax','--dmaxparam', help="Maximum number of allowed substitutions in a match, default value = 5", default=5, type=int)
	parser.add_argument('--psa', help="Amount of subsampling in the suffix array", type=int, default=1)
	parser.add_argument('--pr', help="Amount of subsampling in the rank array", type=int, default=1)
	parser.add_argument('-v','--verbosity', action="count", default=0, help="Increase verbosity level")
	parser.add_argument('-oc','--outputconsole', default=True, help="Print also the output to the screen", type=bool)
	parser.add_argument('-db','--debug', default=False, help="help tracking k-mers extension and positions outputs of reads", type=bool) 

	args = parser.parse_args()

	error = 0
	if args.kparam < 1:
		print "error: values <1 for -k won't produce any k-mer"
		error += 1
	if args.dmaxparam < 0:
		print "error: values <0 for -d won't match in any case"
		error += 1

	if args.psa < 1:
		print "error: invalid parameter psa"
		error += 1

	if args.pr < 1:
		print "errpr: invalid parameter pr"
		error += 1

	if error > 0:
		print "\n{} error(s) in total".format(error)
		exit(-1)

	main(args.rf, args.rds, args.kparam, args.dmaxparam, args.psa, args.pr, args.verbosity, args.outputconsole, args.debug)
