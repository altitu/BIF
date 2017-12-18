#!/usr/bin/python2

import sys
import bwt
import tools_karkkainen_sanders as tks

if __name__ == "__main__":
	psa = int(sys.argv[1])
	pr = int(sys.argv[2])
	s = "GTATGATCAGAA$"
	sa = tks.simple_kark_sort(s)
	b = bwt.getBWT(s, sa)
	sa = bwt.subsampleArray(sa, psa)
	n =  bwt.getN(b)
	ranks = bwt.buildRankArray(b, pr)

	print "BWT =", b
	print "SA =", sa
	print "N =", n
	print "Raw ranks = ", ranks
	print "Ranks ="
	for c in ['$', 'A', 'C', 'G', 'T']:
		st = c + " "
		for i in range(0, len(b)):
			st += str(bwt.rank(ranks, pr, b, c, i)) + " "
		print st


	expectations = [
		("AT", [2,5]),
		("AGAA", [8]),
		("TAT", [1]),
		("GTAC", []),
		("GTATGATCAGAA", [0])
	]
	for e in expectations:
		res = bwt.findSeqInBWT(b, n, ranks, pr, sa, psa, e[0])
		if (len(e[1]) != len(res)):
			print "Different number of results for " + e[0] + " : " + str(res)
		for v in res:
			if v not in e[1]: print "Result not in expectation for " + e[0]
		for v in e[1]:
			if v not in res: print "Expectation not in result for " + e[0]


	source = bwt.getSourceFromBWT(b, ranks, pr, n) + "$"
	if source != s:
		print "Source not matching : " + source


