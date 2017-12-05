#!/usr/bin/python2
# -*- coding: utf-8 -*-

# Partie 1 : extraire les k-mers du read et les match avec le génome

import bwt
import tools_karkkainen_sanders as tks
import prettyprint as pp

def main():

#	s = "GTATGATCAGAA$"
#	sa = tks.simple_kark_sort(s)
#	b = bwt.GET_BWT(s, sa)
#	print b
#	respos = bwt.is_Q_in_S(b, bwt.GET_N(b), sa, "AT")
#	print respos

	s = "ACCCCTACCCCTACCCCG$"
	sa = tks.simple_kark_sort(s)
	b = bwt.GET_BWT(s, sa)
	print b
	respos = bwt.is_Q_in_S(b, bwt.GET_N(b), sa, "ACC")
	print respos

	#data = pp.FastaFile()
	test0 = pp.extends(respos, 1, "ACC", "TACCCC", "ACCCCTACCCCTACCCCG", 0, 1, -10)
	print test0


if __name__ == "__main__":
	main()
