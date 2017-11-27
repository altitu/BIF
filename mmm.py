#!/usr/bin/python
# -*- coding: utf-8 -*-

# Partie 1 : extraire les k-mers du read et les match avec le génome

import bwt
import tools_karkkainen_sanders as tks

def main():
	s = "GTATGATCAGAA$"
	sa = tks.simple_kark_sort(s)
	b = bwt.GET_BWT(s, sa)
	print b
	print bwt.is_Q_in_S(b, bwt.GET_N(b), sa, "AT")


if __name__ == "__main__":
	main()