#!/usr/bin/python2
# -*- coding: utf-8 -*-

# Partie 1 : extraire les k-mers du read et les match avec le génome

import bwt
import tools_karkkainen_sanders as tks
import seedandextend as sae

def main():

#	s = "GTATGATCAGAA$"
#	sa = tks.simple_kark_sort(s)
#	b = bwt.GET_BWT(s, sa)
#	print b
#	respos = bwt.is_Q_in_S(b, bwt.GET_N(b), sa, "AT")
#	print respos

	s = "ACCCCGTACCCCGTACCCC$"
	sa = tks.simple_kark_sort(s)
	b = bwt.GET_BWT(s, sa)
	print b
	respos = bwt.is_Q_in_S(b, bwt.GET_N(b), sa, "ACC")
	print respos

	#data = sae.FastaFile()
	test0 = sae.extends(respos, 4, "ACC", "CCGTACCCC", "ACCCCGTACCCCGTACCCC", 0, 1, -10)
	print test0


if __name__ == "__main__":
	main()
