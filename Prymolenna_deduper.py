#!/usr/bin/env python3
import argparse
import re

def get_args(): # Usage for user input
	parser = argparse.ArgumentParser(description='Deduplicate aligned data based on the position of the aligned reads')
	parser.add_argument("-f", "--file", help = "This needs to be a sorted SAM alignment file. Please provide an ordered SAM file.", required = True)
	parser.add_argument("-p", "--paired", choices=['True', 'False'], help = "Provide this flag if the reads will be paired. WARNING: This program does not tollerate paired end reads.", required=False)
	parser.add_argument("-u", "--umi", help = "This needs to be a FASTQ file of reverse barcode reads. Please provide a list of UMIs", required=False)
	parser.add_argument("-o", "--output", help = "This is the prefix your output files will contain", required = False)
	return parser.parse_args()
args = get_args()

SAM_file = args.file # read in user input and save as a file record
UMI_file = args.umi # read in umi file from user input
paired_check = args.paired # checks for presence of a paired flag
prefix = args.output

if paired_check == 'True':
	print("This program does not tollerate paired end reads. \n\nExiting.")
	quit()

def check_bitflag(flag):
	'''Check for the presense of a bitflag 16, if present means read is on the (-) strand '''
	if((flag & 16) != 16):
		confirmation = False # Checks if the current read is on (+) strand
	elif((flag & 16 ) == 16):
		confirmation = True # Checks if the current read is (-) strand
	return(confirmation)

def check_UMI(qname, list):
	'''Checks if the UMI contained in the QNAME matches a list of known UMIs'''
	sequence = re.search("([A-Z]+$)", qname) # search the qname for the UMI located at the end
	UMI = sequence.group(0) # extract just the sequence from the search object
	if UMI in list: #if the umi matches a known one
		known = True # return the true
	else: # if the umi is not known
		known = False # return false
	return(UMI, known) # output the sequence of the UMI and T/F whether the umi is a known umi

def parse_CIGAR_stranded(cigar_string, start_position):
	'''Parses the CIGAR string for beginning soft clipping and adjusts the starting position'''
	sequence = re.search("(^[0-9]+S)", cigar_string) # search the cigar string for the beginign soft clipping
	if sequence == None:
		clip = 0
		adjusted_pos = start_position-clip
	else:
		clip = int((sequence.group(0))[:-1])
		adjusted_pos = start_position-clip
	return(adjusted_pos)

def parse_CIGAR_unstranded(cigar_string, start_position):
	'''Parses the CIGAR string for beginning soft clipping and adjusts the starting position'''
	end_S = re.findall("([0-9]+)S$", cigar_string)
	D = re.findall("([0-9]+)D", cigar_string)
	N = re.findall("([0-9]+)N", cigar_string)
	M = re.findall("([0-9]+)M", cigar_string)
	addition = 0
	for i in end_S:
		addition+=int(i)
	for i in D:
		addition+=int(i)
	for i in N:
		addition+=int(i)
	for i in M:
		addition+=int(i)
	adjusted_pos = start_position+addition
	return(adjusted_pos)

UMI_handle = open(UMI_file, "r") # open UMI file for reading
SAM_handle = open(SAM_file, "r") # open SAM file for reading
outfile = open(prefix+"_Deduped.sam", "w")
dupe_outfile = open(prefix+"_Duplicates.sam", "w")
trash_outfile = open(prefix+"_UnknownUMI.sam", "w")

UMI_list = [] # initialize empty list to store known UMIs in for comparison later
for line in UMI_handle: # parse through UMI file
	line = line.strip() # strip each line of white space to only contain the UMI sequence
	UMI_list.append(line) # add UMI sequence to list

store_dict = {} # initialize dict to store reads for comparison (key = (strandedness, UMI, adjusted left most starting position), value = sequence)
chr_count = 1 # initialize counter tp keep track of chromsome number
for line in SAM_handle:
	if "@" == line[0]:
		outfile.write(line) # write out beginning lines into the new file
	elif "@" != line[0]:
		line = line.strip() # strip line of whitespace
		line = line.split() # QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL, optional fields
		CHR = line[2] # save chromosome as variable
		if CHR == chr_count: # continue on if the chromosome of the line matched the counter
			QNAME = line[0] # query name contains UMI
			UMI, known = check_UMI(QNAME, UMI_list) # pull out UMI
			if known == False: # trash the line if the UMI is not recognized
				trash_outfile.write('	'.join(line)+"\n") # write to outfile
				continue # exit the most current loop
			else: # otherwise, processs the line
				FLAG = int(line[1]) # contains strandedness confirmation
				flag = check_bitflag(FLAG) # check for strand association
				#print (flag, FLAG)
				if flag == False: # read is on (+) strand
					POS = int(line[3]) # left most position
					CIGAR = str(line[5]) # store the cigar string in a variable
					adj_pos = parse_CIGAR_stranded(CIGAR, POS) # use function to get adjusted left most position by parsing cigar string
					if (CHR, UMI, flag, adj_pos) not in store_dict: # if the unique chr, umi, flag, and adjusted position combo is not in the dictionary
						store_dict[(CHR, UMI, flag, adj_pos)] = line # associate the file line with that unique combo
					elif (CHR, UMI, flag, adj_pos) in store_dict: # if that combo is in the dictionary
						dupe_outfile.write('	'.join(line)+"\n") # write out the encountered line to the dupe file because we already saw one like it
				else: # read is on (-) strand
					POS = int(line[3]) # left most position
					CIGAR = str(line[5]) # store the cigar string in a variable
					adj_pos = parse_CIGAR_unstranded(CIGAR, POS) # use function to get adjusted left most position by parsing cigar string
					if (CHR, UMI, flag, adj_pos) not in store_dict: # if the unique chr, umi, flag, and adjusted position combo is not in the dictionary
						store_dict[(CHR, UMI, flag, adj_pos)] = line # associate the file line with that unique combo
					elif (CHR, UMI, flag, adj_pos) in store_dict: # if that combo is in the dictionary already
						dupe_outfile.write('	'.join(line)+"\n") # write out the encountered line to the dupe file because we already saw one like it
		else: # if the line being parsed has a chromosome name that is different than the one currently being referenced
			chr_count=CHR # make the chromosome contained on the currently processed read the new chromosome
			for key in store_dict: # parse the currently loaded dictionary with the previous chromosome
				outfile.write("{}\n".format('\t'.join(store_dict[key]))) # write out all the keys (lines) to the deduped file
			store_dict = {} # clear the dictionary of the old chr
			# process the current read, don't just pass it because the loop is currently parsing a read
			QNAME = line[0] # query name contains UMI
			UMI, known = check_UMI(QNAME, UMI_list) # pull out umi
			if known == False: # trash the read if the umi is unrecognizable
				trash_outfile.write('	'.join(line)+"\n") # write to trash outfile
				continue # then exit the most immediate loop
			else: # is the UMI is known...
				FLAG = int(line[1]) # contains strandedness confirmation
				flag = check_bitflag(FLAG) # grrab the confirmation off the bitflag of strandedness
				if flag == False: # read is on (+) strand
					POS = int(line[3]) # left most position
					CIGAR = str(line[5])
					adj_pos = parse_CIGAR_stranded(CIGAR, POS)
					store_dict[(CHR, UMI, flag, adj_pos)] = line
				else: # read is on (-) strand
					POS = int(line[3]) # left most position
					CIGAR = str(line[5])
					adj_pos = parse_CIGAR_unstranded(CIGAR, POS)
					store_dict[(CHR, UMI, flag, adj_pos)] = line

for key in store_dict: # eptry the dictioanry one last time beacause there are no new chromosomes to pass
	outfile.write("{}\n".format('\t'.join(store_dict[key]))) # write this out to the deduped file

UMI_handle.close() # close all files currently open for reading or writing
SAM_handle.close()
outfile.close()
dupe_outfile.close()
trash_outfile.close()
