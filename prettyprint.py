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
				result += lenread - (poskr + i + 1)
			if (result > seuil) :
				tabresult.append(result)
			print " "
	return tabresult
'''
def distrib_kmere(read, ktaille, ?):
	#dire que si le kmere est > taille read warning
	
	AAAA
	AA
	 ^pos1
	GGAAAATT
	  ^pos 1
	GAA erreur de poskr = 1 - i + 1 car i déjà upper

	pourGGAAAATT
		^ pos ?
	on a AAT donc lenkread = 4 et poskr = 1 et lenkmere = 2 => lenkread - (lenkmere + poskr) = 4-3 = 1
	- i + 1

def 
'''			
		
'''
class FastaFile:

	fsave = []
	extracheck = False
#[(int id_read, string src, [(int alignement, int pos, bool strand, int d, string comp)])]

	

	def getNbAlign(id_read):
		return len(self[id_read][2])

	def getNbRead():
		return len(self)

	def addAlign(id_read, pos, strand, d, comp):
		self[id_read][2].append((getNBAlign(id_read)-1, pos, strand, d, comp))
	
	def addRead(id_read, src):
		if (if_read-1 > len(ffile)):
			self.append((id_read, src, 0, []))
		else:
			print("warning: kmer n°"+id_read+" already exists in db")

	def addRead(id_read, isRC, pos, strand, d, src, comp):
		if (if_read-1 > getNbRead()): 
			self.append((id_read, src, [(0, pos, strand, d, comp)]))
		else:
			addAlign(id_read, pos, strand, d, comp)

	def getIdRead(name):
		try:
			ret = self.index(name)
			if(extracheck == True and self.count(name) > 1): #short-circuiting because O(n*|nb kmer|)
				print("warning: db unicity failure for kmer: "+name)
		except:
			return -1
		return ret
'''
'''
	def dumpAlign(id_read):
		#for i in range (0, getNbAlign(id_read)-1):
			
		
	def dumpRead():
		#for i in range(0, 


#test: len read = 16 and len kmer = 5
read0 = "AACGTGCTAGTAGCTG"
resB0 = "|||:|:|||||||:||"
genB0 = "AACCTTCTAGTAGGTG"
#        ^pos 0 genome

read1 = "CTAGCAGCGTAGATGC"
resB1 = ":|||:|::||||||||"
genB1 = "GTAGGAAAGTAGATGC"
#        ^pos 17 genome

#kmer  on read n°0
kmer0 = "AACGT"
res0  = "|||:|"
gen0  = "AACCT"
d0    = 1
pos0  = 0

kmer1 = "CGTGC"
res1  = "|:|:|"
gen1  = "CCTTC"
d1    = 2
pos1  = 2

kmer2 = "AACGT"
res2  = "|||:|"
gen2  = "AACCT"
d2    = 1
pos2  = 0

#kmer on read n°1
kmer3 = "CTAGC"
res3  = ":|||:"
gen3  = "GTAGG"
d3    = 1
pos3  = 0

kmer4 = "CAGCG"
res4  = ":|::|"
gen4  = "GAAAG"
d4    = 2
pos4  = 2

kmer5 = "GATGC"
res5  = "|||||"
gen5  = "GATGC"
d5    = 1
pos5  = 0

test0 = FastaFile()
test0.extracheck = True

'''
'''
on a une bwt
on a une collone f qu'on a pas vraiment

on a la colonne i,qui numérote de 0 à ++

on a aussi le suffixe array

un vrai exemple:
l'exemple de GTATGATCAGAA$
(de l'exemple du TD3)
alors:
i |BWT|F
0 |A|$
1 |A|a
2 |G|a
3 |C|a
4 |G|a
5 |T|a
6 |T|c
7 |A|
8 |T
9 |$
10|G
11|A
12|A

pour la query on fait un walk left (on part de la gauche)
exemple: Q = AGA

on regarde le entre les A, on prend du A du début au A de la fin de la colonne F
on est entre A1 position i = 0 et A = 5 (sur le bloque de lignes entre i=1 et i=6)
On va sur la BWT entre les même position 1 et 6, et depuis le A2 (bwt) et T2 (bwt) et on remonte depuis T2 vers le 1er G (on trouve G2) et on descend depuis a2(bwt) et on trouve G1
on revient sur F, on reprend G1, reprend G2
et on regarde sur les même lignes de cet interval le A de la query "AGA" et on trouve A3 (qui existe aussi sur F quand on revient dessus après)

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pour reconstruire la séquence, on part de la position 0 ($=0) on prend $
on regarde sur la ligne correspondante sur BWT on trouve A1
on utilise LF(alpha,k) pour trouver la position du A de rang 1 dans la colone F
on trouve la position i = 2 sur F
on regarde sur la même ligne la lettre corresppnodante et son rang
on trouve du coup A2 sur BWT (pour l'instant ça va car la fin de la demande est [...]AA$)
LF(A,2) = (environ ~3) 3
on continue et on trouve sur BWT à la position 8 (du coup là où est sur F le A2) G1



dans l'archive on s'attend à trouver un fichier python
un manuel developpeur 
comment est orga le code (plus qu'un fichier python possible) on s'attend à un MMM.py, par exemple le fm index peut être dans un fichier fmindex à côté etc decoupage libre)
c'est une doc pour quelqu'un qui reprendrait le code: quel fonction est appelé à quel moment et est dans quel fichier etc

dans le manuel utilisateur = doc utilisateur = "quel sont les lignes de commandes classsiques, les âram d'entrée, le format en entrée et le format en sortie (mais pas obligé de copier collé) et guider l'utilisateur dans les parametres de notre mappeur: deux fichier en entrer + les reads + la refences
ne pas hésiter à utiliser des package de type getup qui facilitent la gestion des options ( ce sont les tirets + lettre minuscules)
ex: MMM.py -r reads.fasta -g genome.fasta
cela permet l'appel des options dans différents ordres
-k -d
si les utilisateurs ne les mettent pas ça génère des options par défaut
getup ou argparse etc
si pas le temps préciser les lignes de commandes précises
un truc sympa = .o sortie.txt pour préciser le nom de sortie
(des petits trucs)

ex: ... -k 20 -d 2
on va expliquer comment l'utilisateur quel valeur il faut mettre à k et à d selon le résultat qu'on veut:
on s'attent donc à ce que nous nous fassions des tests = jeu de données selon les jeux de tests disposés
ex: si k = 19 on aligne bien des read sinon si k = 20 on y arrive pas:
permet de montrer les petits effets de bords qu'on avait pas vu.

dans test2: il y a des données plus importantes: génome de references qui fait 2 M de paires de bases
jeu avec 1000 read et jeu avec 100 000 read.

On attend à ce que nous mesurions le jeu de test1 et test2: combien de temps à obtenir le mapping selon la variation de k sur le temps de mapping
et sur le résultat
plus k grand plus temps de mapping court mais plus rater alignement!

dans le fichier de test2 on a donné le nombre d'allignement qu'on est sencé obtenir (calculé par prof avec k = tout petit) pour savoir le nombre d'alignement possible
si nous => k = 30 alors on aura peut être que 80% des alignements

pour différentes valeurs de dmax on va pouvoir faire un tableauoù on a la valeur de k, le recall (le pourcentage), et le temps d'exec en secondes

pour dmax = 4:
K      recall(%)       time(s)
4      100       
10     1480/1531=?%

dire à l'utilisateur:
"si on veut un dmax de temps alors mieux utiliser k = à ..."
"si on veut tout alignement alors utiliser un k au moins égal à ..."
valeur importante pour nous: le nombre d'alignement qu'on est sencé trouver dans le fihier deread au total

on aattend aussi dans le manuel dev de detailler les opti faite par nous pour ameliorer les perf, et s'il y a eu des choix d'implementations particulier

il n'y a pas de nombre de pages max mais:
pour le manuel dev: max 3 pages 
utilisateur = 6 pages max de max (pour les tableaux haha)

'''
