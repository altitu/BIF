#!/usr/bin/python2
# -*- coding: utf-8 -*-

import tools_karkkainen_sanders as tks
import time

def subsampleArray(array, step):
	return [array[i] for i in range(0, len(array), step)]

# Returns the Burrows-Wheeler transform and suffix array from string
def getBWTAndSA(s, psa):
	sa = tks.simple_kark_sort(s)
	bwt = [0]*len(sa)
	for i in range(len(sa)):
		bwt[i] = s[(sa[i]-1)%len(s)]
	newsa = subsampleArray(sa, psa)
	return (bwt, newsa)

# Returns index associated with character for storing
def getIndex(a):
	if a == '$': return 0
	elif a == 'A': return 1
	elif a == 'C': return 2
	elif a == 'G': return 3
	elif a == 'T': return 4

# Returns the count of each character in the bwt
def getN(bwt):
	n = [0]*5
	for i in range(len(bwt)):
		n[getIndex(bwt[i])] += 1
	return n

# Returns the index in suffix array of character `a` of rank `k`
def lf(n, a, k):
	o = 0
	for i in range(0, getIndex(a)):
		o += n[i]
	return o + k

# Builds the rank array from `bwt`
def buildRankArray(bwt, pr):
	ranks = [[],[],[],[],[]]
	count = [0]*5
	for i in range(len(bwt)):
		ind = getIndex(bwt[i])
		count[ind] += 1
		for j in range(len(ranks)):
			ranks[j].append(count[j])
	for i in range(len(ranks)):
		ranks[i] = subsampleArray(ranks[i], pr)
	return ranks

def rank(ranks, pr, bwt, c, i):
	l = ranks[getIndex(c)]
	if i<0: return 0
	charCount = 0
	prevI = i/pr
	# Get the previous sample and count the number of same characters
	for j in range(i, prevI*pr, -1):
		if bwt[j] == c: charCount += 1
	return l[prevI]+charCount

# Rebuilds source string from `bwt`
# The number of occurences `n` is needed
def getSourceFromBWT(bwt, ranks, pr, n):
	lastChar = '$'
	lastRank = 0
	word = ""
	while True:
		# find position of last character
		pos = lf(n, lastChar, lastRank)
		prevChar = bwt[pos]
		prevRank = rank(ranks, pr, bwt, prevChar, pos-1)
		if (prevChar == '$'): break
		# insert character in front
		word = prevChar + word
		lastChar = prevChar
		lastRank = prevRank
	return word

"""
Returns the position in the reference of the pattern `q`
The number of occurences `n`, the rank array `ranks`, the
suffix array `sa` are needed 
"""
def findSeqInBWT(bwt, n, ranks, pr, sa, psa, q):
	i = len(q)-2
	c = q[i+1]
	# Range in which we search
	r = [lf(n, c, 0), lf(n, c, n[getIndex(c)]-1)]
	while i >= 0:
		prevChar = q[i]
		newranks = [rank(ranks, pr, bwt, prevChar, r[0]-1), rank(ranks, pr, bwt, prevChar, r[1])]
		if newranks[0] == newranks[1]: return []
		r = [lf(n, prevChar, newranks[0]), lf(n, prevChar, newranks[1]-1)]
		i -= 1

	# Find positions in original text
	positions = []
	for ind in range(r[0], r[1]+1):
		# Backtrack the BWT until we find a sample in the suffix array
		i = ind
		backtrack = 0
		while (i%psa) != 0:
			c = bwt[i]
			r = rank(ranks, pr, bwt, c, i-1)
			i = lf(n, c, r)
			backtrack += 1
		# Add numbers of backtracks and wrap
		result = (sa[i/psa] + backtrack)%len(bwt)
		positions.append(result)

	return positions