#Parses the KEG htext file format

# from argparse import ArgumentParser
from collections import defaultdict
from collections import namedtuple
import re
import sys
import STD_Functions as STD

"""Hierarchy structure:
	A series of dictionaries created; A dictionary of B sub dictionaries; B dictionary of C sub dictionaries; C dictionary of D keg names;
	Script to parse out the ko00001.keg file for the Kegg ontology numbers to count tags by higher numbers

	To use: Input raw counts file and then the 
"""


if __name__=="__main__":

	if len(sys.argv)<4:
		print "Usage: python kegg_parserer.py ko00001.keg KAAS_output_file.tab TPM_Counts.tab"
		sys.exit(1)

	k_args=sys.argv[1]
	kegg_genes=sys.argv[2]	
	print kegg_genes
	CountFile=sys.argv[3]
	
# """	
#  PARSE THE ko00001.keg file to get dictionary of dictionaries for each of the levels. 
# """
	A_2_B_Dict={}
	B_2_C_Dict={}
	C_2_D_Dict={}
	for line in open(k_args):
		line=re.sub('<[^<]+?>', '', line)
		if line.startswith("A"):
			A=line.strip("A").strip()
		elif line.startswith("B"):
			line=line.strip().strip("B").strip()
			if line!='':
				b=line
				if A in A_2_B_Dict:
					A_2_B_Dict[A]+=[b]
				else:
					A_2_B_Dict[A]=[b]
		elif line.startswith("C"):
			line=line.strip().strip("C").strip()
			if line!='':
				c=line
				if b in B_2_C_Dict:
					B_2_C_Dict[b]+=[c]
				else:
					B_2_C_Dict[b]=[c]
		elif line.startswith("D"):
			line=line.strip().strip("D").strip()
			if line!='':
				d=line.strip().split()[0]
				if c in C_2_D_Dict:
					C_2_D_Dict[c]+=[d]
				else:
					C_2_D_Dict[c]=[d]

# """	
#  PARSE THE KAAS file (alignment to the KEGG directory of transcripts)
# """
	KASS_Gene_Dict={}
	KASS_Gene_Dict['Unknown']=[]
	for line in open(kegg_genes):
		line=line.strip().split()
		if len(line)==1:
			KASS_Gene_Dict['Unknown']+=line
		else:
			KASS=line[1]
			gene=line[0]
			if KASS in KASS_Gene_Dict:
				KASS_Gene_Dict[KASS]+=[gene]
			else: 
				KASS_Gene_Dict[KASS]=[gene]
# """	
#  IMPORT Count Files
# """
	Counts_raw=STD.importDict(CountFile)
	#Calculate tpm
	Counts_tpm=STD.TPMCalc(Counts_raw)
# """	
#  Sum by KEGG Ontology
# """
# 	for key in B_2_C_Dict:
# 		for k in C_2_D_Dict[key]:
# 			for j in KASS_Gene_Dict:
# 	

	Kass_Count={}			
	for j in KASS_Gene_Dict:
		for gID in KASS_Gene_Dict[j]:
			c=Counts_tpm[gID]
			if j in Kass_Count:
				d=Kass_Count[j]
				cd = [c[i] + d[i] for i in range(len(d))]
				Kass_Count[j]=cd
			else:
				Kass_Count[j]=c
	
	C_2_D_Count={}
	for j in C_2_D_Dict:
		for gID in C_2_D_Dict[j]:
			if gID in Kass_Count:
				c=Kass_Count[gID]
			if j in C_2_D_Count:
				d=C_2_D_Count[j]
				cd = [c[i] + d[i] for i in range(len(d))]
				C_2_D_Count[j]=cd
			else:
				C_2_D_Count[j]=c

	C_2_D_Count['Unknown']=Kass_Count['Unknown']

	B_2_C_Count={}
	for j in B_2_C_Dict:
		for gID in B_2_C_Dict[j]:
			if gID in C_2_D_Dict:
				c=C_2_D_Count[gID]
			if j in B_2_C_Count:
				d=B_2_C_Count[j]
				cd = [c[i] + d[i] for i in range(len(d))]
				B_2_C_Count[j]=cd
			else:
				B_2_C_Count[j]=c

	B_2_C_Count['Unknown']=C_2_D_Count['Unknown']

	for x in KASS_Count:
		print x
		print KASS_Count[x]