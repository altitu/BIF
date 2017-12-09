#!/usr/bin/python
# -*- coding: utf-8 -*-


def openFasta(name):
	#on parse le nom et on affiche le type supposé de fichier	
	pos = -1
	posbar = -1
	ln = len(name)
	for l in range(ln-1, -1):
		if (name[l] == '.'):
			pos = l+1
		if (name[l] == '\\'):
			posbar = l+1
	#nom sans slash
	nameWObar = name[posbar:ln-1]

	print("reading "+nameWObar) 
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
	if (flushed == False):
		reads.append(read)
		flushed = True

	return (identifiers, reads)
			
			
			
	
	#on retournera un tableau de read dans 0 et de commentaire dans 1: (read, comment)*

def main():
	print openFasta("./test_perso/test_readMultiLignes.fasta")

if (__name__ == '__main__'):
	main()
