#!/usr/bin/python
# -*- coding: utf-8 -*-


def openFasta(name, verbosity):

	if (verbosity > 0):
		#on parse le nom et on affiche le type supposé de fichier
		pos = -1
		posbar = -1
		ln = len(name)
		for l in range(ln-1, -1, -1):
			if (name[l] == '.' and (pos == -1)):
				#print "name: "+name[l]
				pos = l+1
				#print "pos: "+str(pos)
			if (name[l] == '/' and (posbar == -1)):
				#print "name: "+name[l]
				posbar = l+1
				#print "posbar: "+str(posbar)
		#nom sans slash
		nameWObar = name[posbar:ln]

		print("reading "+nameWObar) 
		#print name[pos:ln]
		if (name[pos: ln] == "fasta" or name[pos: ln] == "fas" or name[pos: ln] == "fa"):
			print("(generic FASTA file)")
		elif (name[pos: ln] == "fna"):
			print("(fasta nucleic acid)")
		elif (name[pos: ln] == "ffn"):
			print("(fasta functional nucleotide)")
		elif (name[pos: ln] == "faa"):
			print("(fasta amino acid)")
		elif (name[pos: ln] == "frn"):
			print("(fasta RNA non-coding)")
		else:
			print ("(assuming generic FASTA file)")	
	#we attempt to read the file
	try:
		ffile = open(name, 'r')
	except:
		print "ERROR: file cannot be open"
		#traiter l'expetion
	
	#on regarde après l'espace et on indique la provenance supposé de la base de donnée
	#on parse également les reads du fasta
	#isR = False
	#isN = False
	flushed = True
	modeIdent = False
	identifiers=[]
	identifier=""
	reads=[]
	read=""
	firstLine = False
	for i in ffile.read():
		if(i == '\r'):
			
			#isR = True
			modeIdent = False
		elif(i == '\n'):
			#isR = False
			modeIdent = False
			#if (read != ""):
			#	reads.append(read)
			#	read = ""
			#	flushed = True
			if (identifier != ""):
				identifiers.append(identifier)
				identifier = ""
		elif(i == '>'):
			modeIdent = True
			if (read != ""): #les reads peuvent être splitté en plusieurs lignes contigües
				reads.append(read)
				read = ""
				flushed = True
		elif(i == '|'):
			identifiers.append(identifier) #même si vide
			identifier = ""
		else:
			if (modeIdent == True): #on est dans l'ident
				identifier += i
			else: #on est dans un read
				flushed = False
				read += i
				#print read
	if (flushed == False): #si on ne termine pas par un '\n'
		reads.append(read)
		flushed = True
	if (verbosity > 1):
		#on parse pour estimer la provenance de la BDD
		alreadySaid = False
		#db = identifiers[0].split('|')
		db = identifiers
		#print("db:: "+str(db))
		if (len(db) == 1):
			if (alreadySaid == False):
				print "warning: cannot find NCBI sequence identifier"
				print "\tswitching to second parsing mode"
				alreadySaid = True
			db = identifiers[0].split(' ')
			if len(db) > 1:
				print ("identifier: "+db[0])
				if ((db[1][0] == '[') and (db[1][len(db[1])-1] == ']')):
					line = db[1].split('=')
					try:
						print(line[0][1:len(line[0])] +":\n\t" + line[1][0:len(line[1]-1)])
						ending = ""
						for i in range(2, len(line)):
							ending += line[i] + " "
						print(ending)
					except:
						pass
				else:
					ending = ""
					for i in range(1, len(db)):
							ending += db[i] + " "
					print("\t"+ending+"\n")
			else:
				print "\tsecond parsing mode failed, sequence:"
				print str(db)+"\n"
		else:
			recognized = True
			strings = []
			if(db[0] == "bbs" and len(db) > 1):
				strings.append("GenInfo Backbone Id")
				strings.append("Number: "+str(db[1]))
				strings.append("comment:")
				for i in range(2, len(db)):
					if (db[i] != ""):
						strings.append("\t"+db[i])
			elif(db[0] == "dbj" and len(db) > 2):
				strings.append("DDBJ, DNA Database of Japan")
				strings.append("accession: "+str(db[1]))
				strings.append("locus: "+str(db[2]))
				strings.append("comment:")
				for i in range(3, len(db)):
					if (db[i] != ""):
						strings.append("\t"+db[i])
			elif(db[0] == "emb" and len(db) > 2):
				strings.append("EMBL Data Library")
				strings.append("accession: "+str(db[1]))
				strings.append("locus: "+str(db[2]))
				strings.append("comment:")
				for i in range(3, len(db)):
					if (db[i] != ""):
						strings.append("\t"+db[i])
			elif(db[0] == "gb" and len(db) > 2):
				strings.append("GenBank")
				strings.append("accession: "+str(db[1]))
				strings.append("locus: "+str(db[2]))
				strings.append("comment:")
				for i in range(3, len(db)):
					if (db[i] != ""):
						strings.append("\t"+db[i])
			elif(db[0] == "gnl" and len(db) > 2):
				strings.append("General database identifier")
				strings.append("database: "+str(db[1]))
				strings.append("identifier: "+str(db[2]))
				strings.append("comment:")
				for i in range(3, len(db)):
					if (db[i] != ""):
						strings.append("\t"+db[i])
			elif(db[0] == "lcl" and len(db) > 1):
				strings.append("Local Sequence identifier")
				strings.append("identifier: "+str(db[1]))
				strings.append("comment:")
				for i in range(2, len(db)):
					if (db[i] != ""):
						strings.append("\t"+db[i])
			elif(db[0] == "pat" and len(db) > 2):
				strings.append("Patents")
				strings.append("country: "+str(db[1]))
				strings.append("number: "+str(db[2]))
				strings.append("comment:")
				for i in range(3, len(db)):
					if (db[i] != ""):
						strings.append("\t"+db[i])
			elif(db[0] == "pdb" and len(db) > 2):
				strings.append("Brookhaven Protein Data Bank")
				strings.append("entry: "+str(db[1]))
				strings.append("chain: "+str(db[2]))
				strings.append("comment:")
				for i in range(3, len(db)):
					if (db[i] != ""):
						strings.append("\t"+db[i])
			elif(db[0] == "pir" and len(db) > 2):
				strings.append("NBRF PIR")
				strings.append("entry: "+str(db[2]))
				strings.append("comment:")
				for i in range(3, len(db)):
					if (db[i] != ""):
						strings.append("\t"+db[i])
			elif(db[0] == "prf" and len(db) > 2):
				strings.append("Protein Research Foundation")
				strings.append("name: "+str(db[2]))
				strings.append("comment:")
				for i in range(3, len(db)):
					if (db[i] != ""):
						strings.append("\t"+db[i])
			elif(db[0] == "ref" and len(db) > 2):
				strings.append("NCBI Reference Sequence")
				strings.append("accession: "+str(db[1]))
				strings.append("locus: "+str(db[2]))
				strings.append("comment:")
				for i in range(3, len(db)):
					if (db[i] != ""):
						strings.append("\t"+db[i])
			elif(db[0] == "sp" and len(db) > 2):
				strings.append("SWISS-PROT")
				strings.append("accession: "+str(db[1]))
				strings.append("entry name: "+str(db[2]))
				strings.append("comment:")
				for i in range(3, len(db)):
					if (db[i] != ""):
						strings.append("\t"+db[i])
			else:
				print "unrecognize db provenance"
				recognized = False

			if (recognized == True):
				for i in strings:
					print(i)

	return (identifiers, reads)
			
			
			
	
	#on retournera un tableau de read dans 0 et de commentaire dans 1: (read, comment)*

def main():
	genome = openFasta("./test_perso/test_dbPRF.fasta", 0)
	print genome[1][0]


if (__name__ == '__main__'):
	main()
