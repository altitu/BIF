#!/usr/bin/python2
# -*- coding: utf-8 -*-

import bwt
import tools_karkkainen_sanders as tks
import seedandextend as sae
from filemanager import openFasta
import argparse
import time


def main(refFilename, readsFilename, k, dmax, psa, pr, verb, outputconsole, debug, output, bench):
	start = time.time()
	genome = openFasta(refFilename, verb) 
	reads = openFasta(readsFilename, verb)
	if verb >=2 or debug:
		print "read results of reads sequence:"
		print reads[1]
	s = genome[1][0] + "$"
	if verb >= 2 or debug:
		print "read results of genome sequence:"
		print s
	sa = tks.simple_kark_sort(s)
	b = bwt.getBWT(s, sa)
	sa = bwt.subsampleArray(sa, psa)
	if verb >= 1:
		print "genome BWT:" + str(b)
	n = bwt.getN(b)
	ranks = bwt.buildRankArray(b)
	ranks = bwt.subsampleArray(ranks, pr)
	lfmap = bwt.getLFMapping(b)
	
	mid = time.time()

	if bench:
		print "BWT : {} seconds(s)".format(mid-start)

	ffile = open(output, 'w')
	ffile.write(sae.distributeReads(reads, k, dmax, s, b, sa, psa, n, ranks, pr, lfmap, verb, debug, outputconsole))
	ffile.close()
	end = time.time()
	if bench:
		print "benchmark: computations finished in {} second(s)".format(end - start)

	#dmax = len(read) - max acceptable

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description=
		"Finds the position and differences of reads inside a reference genome",epilog=
		"Usage:\n"+
		"\tdefault values for reference and reads:\n"+
		"\t\t./mmm.py ./test1/reference1.fasta ./test1/reads.fasta\n\n"+
		"\tk = 5 and dmax = 20 for reference and reads:\n"+
		"\t\t./mmm.py ./test1/reference1.fasta ./test1/reads.fasta -dmax 20 -k 5\n\n"+
		"\tverbosity, but the output.txt will not be print on screen:\n"+
		"\t\t./mmm.py ./test1/reference1.fasta ./test1/reads.fasta -v -oc\n\n"+
		"\tvery verbose, debug mode on:\n"+
		"\t\t./mmm.py ./test1/reference1.fasta ./test1/reads.fasta -vv -db\n\n"+
		"\trank and suffix arrays both subsampled by a factor 2 and display computation time:\n"+
		"\t\t./mmm.py ./test1/reference1.fasta ./test1/reads.fasta --psa=2 --pr=2 -b\n\n"+
		"\tk = 5, rank subsampled by 2, no output on screen and output in log.txt:\n"+
		"\t\t./mmm.py ./test1/reference1.fasta ./test1/reads.fasta -k 5 --pr 2 -oc -o ./log.txt\n",
		formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('reference', help="Filename to reference genome")
	parser.add_argument('reads', help="Filename to reads")
	parser.add_argument('-o','--output', help="Filename to output the result", default="./output.txt")
	parser.add_argument('-k','--kparam', help="K-mer size, default value = 20", default = 20, type=int)
	parser.add_argument('-dmax','--dmaxparam', help="Maximum number of allowed substitutions in a match, default value = 5", default=5, type=int)
	parser.add_argument('--psa', help="Amount of subsampling in the suffix array", type=int, default=1)
	parser.add_argument('--pr', help="Amount of subsampling in the rank array", type=int, default=1)
	parser.add_argument('-v','--verbosity', action="count", default=0, help="Increase verbosity level")
	parser.add_argument('-oc','--outputconsole', const=False, nargs='?', default=True, help="Print also the output to the screen (True by default)", type=bool)
	parser.add_argument('-db','--debug', const=True, nargs='?', default=False, help="help tracking k-mers extension and positions outputs of reads (False by default)", type=bool)
	parser.add_argument('-b','--benchmark', const=True, nargs='?', default=False, help="print total computation time from the first opening to the final closing", type=bool)

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
		print "error: invalid parameter pr"
		error += 1

	if error > 0:
		print "\n{} error(s) in total".format(error)
		exit(-1)

	main(args.reference, args.reads, args.kparam, args.dmaxparam, args.psa, args.pr, args.verbosity, args.outputconsole, args.debug, args.output, args.benchmark)
