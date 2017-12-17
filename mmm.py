#!/usr/bin/python2
# -*- coding: utf-8 -*-

# Partie 1 : extraire les k-mers du read et les match avec le génome

import bwt
import tools_karkkainen_sanders as tks
import seedandextend as sae
from filemanager import openFasta

#parser les args ici

def main():

#	s = "GTATGATCAGAA$"
#	sa = tks.simple_kark_sort(s)
#	b = bwt.GET_BWT(s, sa)
#	print b
#	respos = bwt.is_Q_in_S(b, bwt.GET_N(b), sa, "AT")
#	print respos

	#s = "ACCCCGTACCCCGTACCCC$"
#	genome = openFasta("./test1/reference1.fasta", 0) #ok
#	reads = openFasta("./test1/reads.fasta", 0) #ok
	genome = [['commentaire random'],['GCATGCTTTTGCCGAT']]
	reads = [['commentaire random'],['ATGC','TTGC']]
	print "resultat des reads:"
	print reads[1]
	s = genome[1][0] + "$"
	print "et pour le genome:"
	print s
	sa = tks.simple_kark_sort(s)
	b = bwt.GET_BWT(s, sa)
	print "bwt du génome:" + str(b)
#	respos = bwt.is_Q_in_S(b, bwt.GET_N(b), sa, "ACC") #savoir s'il est présent dans la séquence
#	respos = bwt.is_Q_in_S(b, bwt.GET_N(b), sa, "ATGATCAG") #savoir s'il est présent dans la séquence
#	print  "position de match parfait obtenus: " + str(respos)

#	#data = sae.FastaFile()
#	test0 = sae.extends(respos, 4, "ACC", "CCGTACCCC", s, 0, 1, -10)
#	test0 = sae.extends(respos, 0, "ATGATCA", "ATGATCAG", s, 0, 1, -10)
#	print test0

	sae.distributeReads(reads, 2, 1, s, b, sa)
	#dmax = len(read) - max acceptable

if __name__ == "__main__":
	main()
