attention, ce manuel dev est limité à 3 pages

première partie: BW sur reference.fasta

se servir de karkkainen_sanders pour le tableau des suffixes
# voir tp
-> trier suffixe
http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.428.6589&rep=rep1&type=pdf
https://en.wikipedia.org/wiki/List_of_sequence_alignment_software#Short-Read_Sequence_Alignment

pour eviter de tester le même read plusieurs fois on peut =>
utiliser arbre suffixe sur read, comme ça test plusieurs read sur un même alignement
(ne pas oublié les reverses complements => penser à rajouter une colonne dans ce tableau pour le 
strand -> le sens de lecture de la seq et RC (strand = 1 -> normal; -1 -> RC)




19/20 kmer
4/5 dmax
(jeu de test)

d-> pcID (pourcentage d'identité = inverse du pourcentage de mismatch)
si d = 3 alors pcID = 97
