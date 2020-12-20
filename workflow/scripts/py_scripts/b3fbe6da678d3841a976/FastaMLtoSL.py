#!/usr/bin/env python 
# Created by: Lee Bergstrand 
# Descript: Converts multiline FASTAs to single line FASTAs
#  
# Usage: FastaMLtoSL.py <sequences.faa> 
# Example: FastaMLtoSL.py mySeqs.faa 
#----------------------------------------------------------------------------------------
#===========================================================================================================
#Imports:
	
import sys
import re

class FastaMLtoSL():

	def __init__(self, inFile):
		self.inFile = inFile			
		# Stores file one for input checking.
	
	def returnSLFasta(self) :
# do whatever the script does
# #==========================================================================
		print(">> Opening FASTA file...")
		# Reads sequence file list and stores it as a string object. Safely closes file:
		try:	
			with open(self.inFile,"r") as newFile:
				sequences = newFile.read()
				sequences = re.split("^>", sequences, flags=re.MULTILINE) # Only splits string at the start of a line.
				del sequences[0] # The first fasta in the file is split into an empty empty element and and the first fasta
								# Del removes this empty element.
				newFile.close()
		except IOError:
			print("Failed to open " + self.inFile)
			exit(1)

		print(">> Converting FASTA file from multiline to single line and writing to file.")
		# Conversts multiline fasta to single line. Writes new fasta to file.
		try:	
			with open(self.inFile + ".out","w") as newFasta:
				for fasta in sequences:
					try:
						header, sequence = fasta.split("\n", 1) # Split each fasta into header and sequence.
					except ValueError:
						print(fasta)
					header = ">" + header + "\n" # Replace ">" lost in ">" split, Replace "\n" lost in split directly above.
					sequence = sequence.replace("\n","") + "\n" # Replace newlines in sequence, remember to add one to the end.
					newFasta.write(header + sequence)
				newFasta.close()
		except IOError:
			print("Failed to open " + self.inFile)

		print(">> Done!")
