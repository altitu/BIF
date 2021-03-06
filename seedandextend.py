#!/usr/bin/python
# -*- coding: utf-8 -*-

import bwt
import time

def flattenUniq(n, isDebug):
	#flattening
	result = []
	res = []

	if isDebug:
		print "starting n:"
		print n

	for i in n:
		for j in i:
			if j != []:
				result.append(j)
	if isDebug:
		print "flattened n:"
		print result

	#uniquify could probably be optimised
	#nlog(n) is not reachable here
	#we do not have a partially ordered set "<"
	#the result is not a well-founded relation
	for i in range(0, len(result)):
		if result[i] != []:
			for j in range(1, len(result)):
					if result[j] != [] and j != i:

						if isDebug:
							print "i= "+str(result[i])+" et j= "+str(result[j])
							print "p1= "+str(i)+" p2= "+str(j)
							print "state: "+str(result)

						if result[i] == result[j]:
							result[j] = []
	#Because we can't use set()
	if isDebug:	
		print "n only:"
		print result

	for i in result:
		if i != []:
			res.append(i)
	if isDebug:
		print "the result:"
		print res

	return res
				

def revComp(string):
	output = ""
	for i in range(len(string)-1, -1, -1):
		if (string[i] == 'A'):
			output += 'T'
		elif (string[i] == 'G'):
			output += 'C'
		elif (string[i] == 'C'):
			output += 'G'
		elif (string[i] == 'T'):
			output += 'A'
		else:
			print "error: maybe $ was left in the string, or it was uncapitalized"
			output += "?"
	
	return output

#pretty print of 2 strings
def compprettyprint(a, b):
	output = ""
	d = 0
	if len(a) != len(b):
		print "ERROR: length of string a != length of string b!"
	for i in range(0, len(a)):
		if a[i] == b[i]:
			output += "|"
		else:
			output += ":"
			d += 1
	return (output, d)

#return the string to be outputed in a specified file (see mmm.py) for the user
def distributeReads(reads, k, dmax, genome, b, sa, psa, narray, ranks, pr, verbosity, debug, oc): #the genome sequence should be provided with an "$" ending
	output = ""
	norm_reads = reads[1] #norm_reads becomes a read array
	result = []
	result_comp = []
	numread = 1

	nb_reads = len(norm_reads)
	firsttime = time.time()

	for read in norm_reads:
		
		if (1.0*numread/nb_reads)*100%10 == 0:
			nowtime = time.time()
			pourcentage = (1.0*numread/nb_reads)*100
			print "{}% finished, time remaining: {} second(s)".format((1.0*numread/nb_reads*100),
					(nowtime-firsttime)* ((100 -pourcentage)/pourcentage))
		
		result = []
		result_comp = []
		read_comp = revComp(read)
		srand = 1
		r = cutread(read, k) #on a les kmers
		r_comp = cutread(read_comp, k)
		
		if verbosity >=2:
			print "the read "+read+" returns the following k-mers:"
			print str(r)+"\n"
			print "the reverse complement (revcomp) of read "+read_comp+" returns the following k-mers:"
			print str(r_comp)+"\n"
		
		i = 0
		for kmere in r:
			start = time.time()
			respos = bwt.findSeqInBWT(b, narray, ranks, pr, sa, psa, kmere)
		##if len(respos) > 0:
			if verbosity >=1:
				print "perfect match position obtained for "+str(kmere)+":"
				print str(respos)+"\n"
			#We compute again the abolutes positions of the kmers with respect to the genome sequence
			rext = extends(respos, i, kmere, read, genome, 0, 1, dmax, verbosity, debug)
			for n in range(0, len(rext)):
				if rext[n] != []:
					rext[n] = respos[n] - i
		##if len(rext) > 0:
			result.append(rext)
			i += 1

		l = 0
		for kmere in r_comp:
			respos_comp = bwt.findSeqInBWT(b, narray, ranks, pr, sa, psa, kmere)
		##if len(respos) > 0:
			if verbosity >=1:
				print "perfect match position obtained for "+str(kmere)+":"
				print str(respos_comp)+"\n"
			#We compute again the abolutes positions of the kmers with respect to the genome sequence
			rext_comp = extends(respos_comp, l, kmere, read_comp, genome, 0, 1, dmax, verbosity, debug)
			for n in range(0, len(rext_comp)):
				if rext_comp[n] != []:
					rext_comp[n] = respos_comp[n] - l
		##if len(rext) > 0:
			result_comp.append(rext_comp)
			l += 1

		norm_flat = flattenUniq(result, debug)
		if verbosity >=2:
			print "list of alignments for this read:"
			print str(norm_flat) + "\n"
		comp_flat = flattenUniq(result_comp, debug)
		if verbosity >=2:
			print "list of alignments for this read revcomp:"
			print str(comp_flat) + "\n"

		if (norm_flat != [] or comp_flat != []):
			output += ">read"+str(numread)+"\n"
		
		numalign = 0
		for i in range(0, len(norm_flat)):
			rescpp = compprettyprint(read,genome[norm_flat[i]:norm_flat[i]+len(read)])
			output += "  >>alignment "+str(numalign)+"\n"
			output += "  #pos="+str(norm_flat[i])+"\n"
			output += "  #strand=+1"+"\n"
			output += "  #d="+str(rescpp[1])+"\n"
			output += "  "+str(read)+"\n"
			output += "  "+rescpp[0]+"\n"
			output += "  "+	genome[norm_flat[i]:norm_flat[i]+len(read)]+"\n"
			numalign +=1

		for l in range(0, len(comp_flat)):
			rescpp_comp = compprettyprint(read_comp,genome[comp_flat[l]:comp_flat[l]+len(read_comp)])
			output += "  >>alignment "+str(numalign)+"\n"
			output += "  #pos="+str(comp_flat[l])+"\n"
			output += "  #strand=-1"+"\n"
			output += "  #d="+str(rescpp_comp[1])+"\n"
			output += "  "+str(read_comp)+"\n"
			output += "  "+rescpp_comp[0]+"\n"
			output += "  "+	genome[comp_flat[l]:comp_flat[l]+len(read_comp)]+"\n"
			numalign +=1

		

		numread += 1
	if oc:
		print output
	return output

#cut reads into k-mers
def cutread(read, k):
	result = []
	lenread = len(read)
	if (k > lenread or k < 1):
		#exception
		print "error: k value should be between of one and a read ({}), both included".format(lenread)
		exit(-1)
	for i in range(0, lenread-k + 1):
		result.append(read[i:k+i])
	return result

#extends output the number of difference between the match found and the read sequence, or [] if that number is > to dmax
#it does not provide any position informations
#such informations are computed in distributeReads()
def extends(tabpos_kmere_sur_genome, poskr, kmere, read, genome, smatch, smismatch, seuil, verbosity, debug):
	tabresult = []
	result = 0
	lenread = len(read)
	lengen = len(genome)
	lenkmere = len(kmere)
	if (tabpos_kmere_sur_genome != []):
		for poskg in tabpos_kmere_sur_genome:
			i = 1
			result = 0
			while ((poskg - i >= 0) and (poskr - i >= 0)):
				if (genome[poskg - i] == read[poskr - i]) :
					result += smatch
				else:
					result += smismatch
				if debug:
					print "aread " + read[poskr - i]
					print "ageno " + genome[poskg - i]
					print "pkr-i " + str(poskr - i)
					print "pkg-i " + str(poskg - i)
				i += 1
			i -= 1

			if ((poskg - (i + 1)) < 0):
				if debug:
					print "brefore kmere"
				result += (poskr - i) * smismatch
			if debug:
				print "preresult: "+str(result)
			i = lenkmere
			result += lenkmere * smatch
			if debug:
				print "kmere " + kmere
			while ((poskg + i < lengen) and (poskr + i < lenread)):
				if (genome[poskg + i] == read[poskr + i]) :
					result += smatch
				else:
					result += smismatch
				if debug:
					print "bread " + read[poskr + i]
					print "bgeno " + genome[poskg + i]
					print "pkr-i " + str(poskr + i)
					print "pkg-i " + str(poskg + i)
				i += 1
			i -= 1

			if (poskg + i + 1 >= lengen):
				if debug:
					print "after kmere"
				result += (lenread - (poskr + i + 1)) * smismatch
			if (result <= seuil) :
				tabresult.append(result)
			else:
				tabresult.append([]) #is not a int

	return tabresult

