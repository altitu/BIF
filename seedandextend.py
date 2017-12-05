#!/usr/bin/python
# -*- coding: utf-8 -*-


#do we treat directly with reverse complement as mentionned in the bonus section of the subject ?

#dmax correspond à 1, -1, dmax

'''poskr = pos_kmere_sur_read, score_match ex:0, score_mismatch ex:1, seuil par le dessus de renvoit ex: -1 exclut'''
def extends(respos, poskr, kmere, read, genome, smatch, smismatch, seuil):
	test = respos[0]
	tabpos_kmere_sur_genome = respos[1]
	tabresult = []
	result = 0
	lenread = len(read)
	lengen = len(genome)
	lenkmere = len(kmere)
	if (test == True):
		for poskg in tabpos_kmere_sur_genome:
			i = 1
			result = 0
			while ((poskg - i >= 0) and (poskr - i >= 0)):
				if (genome[poskg - i] == read[poskr - i]) :
					result += smatch
				else:
					result += smismatch
				print "aread " + read[poskr - i]
				print "ageno " + genome[poskg - i]
				print "pkr-i " + str(poskr - i)
				print "pkg-i " + str(poskg - i)
				i += 1
			i -= 1
			#si il reste alors smismatch * (lenkmere - i)
			#result += max((poskg - i), (poskr - i)) * smismatch
			print "blabla " + str(poskg - (i + 1))
			if ((poskg - (i + 1)) < 0):
				print "ici                                                ici"
				result += (poskr - i) * smismatch
			print result
			i = lenkmere
			result += lenkmere * smatch
			print "kmere " + kmere
			while ((poskg + i < lengen) and (poskr + i < lenread)):
				if (genome[poskg + i] == read[poskr + i]) :
					result += smatch
				else:
					result += smismatch
				print "bread " + read[poskr + i]
				print "bgeno " + genome[poskg + i]
				print "pkr-i " + str(poskr + i)
				print "pkg-i " + str(poskg + i)
				i += 1
			i -= 1
			#result += (lenread - (poskr + i)) * smismatch
			if (poskg + i + 1 >= lengen):
				print "ici"
				result += (lenread - (poskr + i + 1)) * smismatch
			if (result > seuil) :
				tabresult.append(result)
			print " "
	return tabresult

