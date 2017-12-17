#!/usr/bin/python2
# -*- coding: utf-8 -*-

import tools_karkkainen_sanders as tks

# Returns the Burrows-Wheeler transform from string s and suffix array sa
def getBWT(s, sa):
	bwt = [0]*len(sa)
	for i in range(len(sa)):
		bwt[i] = s[(sa[i]-1)%len(s)]
	return bwt

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
def lf(a, k, n):
	o = 0
	for i in range(0, getIndex(a)):
		o += n[i]
	return o + k - 1

# Builds the rank array from `bwt`
def buildRankArray(bwt):
	prec = [0] * 5
	rang = [0] * len(bwt)
	for i in range(len(bwt)):
		index = getIndex(bwt[i])
		prec[index]+=1
		rang[i] = prec[index]
	return rang

# Rebuilds source string from `bwt`
# The number of occurences `n` is needed
def getSourceFromBWT(bwt, ranks, n):
	lastChar = '$'
	lastOcc = 1
	word = ""
	while True:
		# find position of last character
		pos = lf(lastChar, lastOcc, n)
		prevChar = bwt[pos]
		prevOcc = ranks[pos]
		if (prevChar == '$'): break
		# insert character in front
		word = prevChar + word
		lastChar = prevChar
		lastOcc = prevOcc
	return word


"""
Returns the position in the reference of the pattern `q`
The number of occurences `n`, the rank array `ranks`, the
suffix array `sa` are needed 
"""
def findSeqInBWT(bwt, n, ranks, sa, q):
	i = len(q)-2
	c = q[i+1]
	# Range in which we search
	r = [lf(c, 1, n), lf(c, n[getIndex(c)], n)]
	while i >= 0:
		# Find the new range
		r1 = [-1, -1]
		# Find first bound
		for j in range(r[0], r[1]+1):
			if bwt[j] == q[i] :
				r1[0] = lf(q[i], ranks[j], n)
				break
		# Find last bound
		for j in range(r[1], r[0]-1, -1):
			if bwt[j] == q[i] :
				r1[1] = lf(q[i], ranks[j], n)
				break
		# If new range has no elements, no elements found
		if r1[0] == -1: return []
		i -= 1
		r = r1

	positions = [sa[ind] for ind in range(r[0], r[1]+1)]

	return positions