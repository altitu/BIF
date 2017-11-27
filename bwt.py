#!/usr/bin/python

import tools_karkkainen_sanders as tks
import random
import time
import fasta

def LCP_first(s, sa):
	lcp = [0] * len(sa)
	for i in range(1, len(sa)):
		while True:
			if s[sa[i]+lcp[i]] == s[sa[i-1]+lcp[i]]:
				lcp[i] += 1
			else: break
	return lcp

#Q3
def INV(sa):
	inv = [0]*len(sa)
	for i in range(len(sa)):
		inv[sa[i]] = i
	return inv

#Q4
def LCPfromL(s, sa, line, L):
	if (line == 0): return 0
	while s[sa[line]+L] == s[sa[line-1]+L] :
		L += 1
	return L

def LCP_linear(s, sa, inv):
	lcp = [0]*len(sa)
	prec = 0
	for it in range(len(sa)):
		line = inv[it]
		lcp[line] = LCPfromL(s, sa, line, max(prec-1, 0))
		prec = lcp[line]
	return lcp

def LARGEST_REPEAT_one(s, sa, lcp):
	maximum = 0
	pos = 0
	for i in range(len(lcp)):
		if lcp[i] > maximum:
			maximum = lcp[i]
			pos = i
	print("largest = " + str(s[sa[i]:sa[i]+maximum]) + " at " + str(sa[i]))


def LARGEST_REPEAT_all(s, sa, lcp):
	maxlcp = max(lcp)
	positions = [i for i, x in enumerate(lcp) if x == maxlcp]
	print positions

def GET_BWT(s, sa):
	bwt = [0]*len(sa)
	for i in range(len(sa)):
		bwt[i] = s[(sa[i]-1)%len(s)]
	return bwt

def GET_N(bwt):
	n = {'$' : 0, 'A' : 0, 'C' : 0, 'G' : 0, 'T' : 0}
	for i in range(len(bwt)):
		n[bwt[i]] += 1
	return n

def GET_OFFSETS(n):
	o = {}
	o['$'] = 0
	o['A'] = n['$']
	o['C'] = n['A'] + o['A']
	o['G'] = n['C'] + o['C']
	o['T'] = n['G'] + o['G']
	return o

def LF(alpha, k, n):
	o = GET_OFFSETS(n)
	return o[alpha] + k - 1

def build_rank(bwt):
	prec = [0] * 5
	rang = [0] * len(bwt)
	for i in range(len(bwt)):
		if (bwt[i]=='$'):
			prec[0]+=1
			rang[i] = prec[0]
		elif bwt[i] == 'A':
			prec[1]+=1
			rang[i] = prec[1]
		elif bwt[i] == 'C':
			prec[2]+=1
			rang[i] = prec[2]
		elif bwt[i] == 'G':
			prec[3]+=1
			rang[i] = prec[3]
		elif bwt[i] == 'T':
			prec[4]+=1
			rang[i] = prec[4]
	return rang

def rank(l, ranks):
	return ranks[l]

# Reconstruit la séquence à partir de la bwt et les rangs n
def BWT2SEQ(bwt, n):
	ranks = build_rank(bwt)
	lastChar = '$'
	lastOcc = 1
	word = ""
	while True:
		pos = LF(lastChar, lastOcc, n)
		prevChar = bwt[pos]
		prevOcc = rank(pos, ranks)
		if (prevChar == '$'): break
		word = prevChar + word
		lastChar = prevChar
		lastOcc = prevOcc
	return word


# Cherche un pattern Q dans la bwt avec les rangs n
def is_Q_in_S(bwt, n ,q):
	ranks = build_rank(bwt)
	i = len(q)-2
	c = q[i+1]
	r = [LF(c, 1, n), LF(c, n[c], n)]
	while i >= 0:
		r1 = [-1, -1]
		# Find first
		for j in range(r[0], r[1]+1):
			if bwt[j] == q[i] :
				r1[0] = LF(q[i], rank(j, ranks), n)
				break
		# Find last
		for j in range(r[1], r[0]-1, -1):
			if bwt[j] == q[i] :
				r1[1] = LF(q[i], rank(j, ranks), n)
				break
		if r1[0] == -1: return False
		i -= 1
		r = r1
	return True
