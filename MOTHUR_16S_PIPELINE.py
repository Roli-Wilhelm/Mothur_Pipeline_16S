#!/usr/bin/python
import sys, os, re, getopt, glob, numpy as np, random
import timeit
import itertools

now = timeit.default_timer()

Usage = """
Usage:  ./MOTHUR_16S_PIPELINE.py -o BACTERIA_0.05 -n BACTERIA_0.05 -i IIKFCBR01.fasta -f IIKFCBR01.oligos -p 8

or

Usage:  ./MOTHUR_16S_PIPELINE.py -o BACTERIA_0.05 -n BACTERIA_0.05 -m IIKFCBR01.fasta -f IIKFCBR01.oligos -p 8


REQUIRED ARGUMENTS:
		-o	output directory
	
                -n      the name given to the project

		-p	the number of processors available to use

		EITHER:
                -i      the input file (accepted: .fasta, .sff, .shhh.fasta. MUST SPECIFY if the latter two)
			(if providing .fasta, one must provide a .qual file)

		&&

		-f	the .oligos file

		OR JUST:
		-m	specify directory if you would like to combine multiple .sff files from a specific directory.
			(the directory must also contain all the associated .oligos files for each run)

DEPENDENCY:
		-mothur v.1.32.1 (written using, probably valide for earlier and later versions)
		-Green Genes "Mothur Formatted" Database
		-Silva Bacteria Database
		

OPTIONAL:
		-S <Y>	Input is .sff file
			(Must provide "Y")

		-R <Y>	SELECT IF you've already run .shhh flows and do not want to re-run (obviously b/c it takes forever), choose this option.
			(Required: *.shhh.names, *.shhh.groups, *.shhh.fasta)

		-O      If you get an error regarding your "flow.files", there is a chance you must specify a different option for flow file to use (ex. B). Do so here.

                -b <Y>  Enter LOGGING MODE - saves every command and other important information to "command.log" 
                        (Must provide "Y")

UTILITY:
This script will performs the general steps in the mothur pipeline for fungal ITS sequences. This script will require tailoring for different datasets, but provides a simple archetype as well as a quick implementation of whatever you decide. The script does do "summary.seqs" at various points along the way, but the user will have to manually examine the outputs to see how the cleaning went. There is one point where manual intervention is required and the script will wait until one does so."


NOTES: 
1) PROCESSORS ARE SET TO 6. MANUALLY CHANGE IF DESIRED.
2) Minflows is set to 180 (i.e. the minimum length of sequence will be 180bp)
3) Maxflows is set to 450 (i.e. the maximum length of sequence will be 450bp)
4) If you are merging multiple runs AND if there are overlapping barcodes, the script CAN handle this, HOWEVER
it makes the assumption that the first 9 nucleotides of your primer are not degenerate (i.e. the same for all sequences).
This is b/c the barcodes will be made unique for all of your reads and to ensure short barcode length sequences
within your sequences are not erroneously changed, the script searches for barcode + 9 nucleotides of primer.

"""

if len(sys.argv)<3:
        print Usage
        sys.exit()


# Store input and output file names
NAME=''
INPUT=''
OLIGOS=''
RERUN=''
SFF=''
ORDER=''
PROCESSORS=''
OUTPUT=''
MIX=''
Debug=''
NO_SHHH=''

# Read command line args
myopts, args = getopt.getopt(sys.argv[1:],"n:i:o:S:b:R:f:O:p:m:")

###############################
# o == option
# a == argument passed to the o
###############################
for o, a in myopts:
    if o == '-o':
        OUTPUT= a
    if o == '-n':
        NAME= a
    if o == '-i':
        INPUT= a
    if o == '-f':
        OLIGOS= a
    if o == '-S':
        SFF= a
    if o == '-R':
        RERUN= a
    if o == '-O':
        ORDER= a
    if o == '-p':
        PROCESSORS= a
    if o == '-m':
        MIX= a
    if o == '-b':
        Debug= a

####
# Function to Find Duplicates in List
####
def list_duplicates(seq):
	seen = set()
	# adds all elements it doesn't know yet to seen and all other to seen_twice
	seen_add = seen.add

	# turn the set into a list (as requested)
	seen_twice = set( x for x in seq if x in seen or seen_add(x) )

	return list(seen_twice)

def id_generator(size, chars):
        return ''.join(random.choice(chars) for _ in range(size))

if Debug:
        print "You are in Debugging Mode"
        error = open("command.log", "w")

if len(OUTPUT)>0:
        if os.path.exists('./' + OUTPUT):
                print "\nOutput Folder Exists - Caution: Files May Be Re-Written"
        else:
                os.mkdir(OUTPUT)
if MIX:
	INPUT = []
	if Debug:
       		print "You chose to merge multiple files for input.\n"

	if SFF:
		for FILE in glob.glob("./"+MIX+"/*.sff"):
			INPUT.append(re.sub(".sff", "", FILE))

			if Debug:
			        error.write("The base filename of your input file is:"+re.sub(".sff", "", FILE)+"\n")
	elif RERUN:
		for FILE in glob.glob("./"+MIX+"/*.shhh.fasta"):
			INPUT.append(re.sub(".shhh.fasta", "", FILE))

			if Debug:
			        error.write("The base filename of your input file is:"+re.sub(".shhh.fasta", "", FILE)+"\n")
	else:
		for FILE in glob.glob("./"+MIX+"/*.fasta"):
			INPUT.append(re.sub(".fasta", "", FILE))

			if Debug:
			        error.write("The base filename of your input file is:"+re.sub(".fasta", "", FILE)+"\n")
	if Debug:
       		error.write(str(INPUT)+"\n")
	
else:
	if SFF:
		INPUT = re.sub(".sff", "", INPUT)
		print INPUT

		if Debug:
        		print "The basename of your input file is:\n"
		        error.write(INPUT+"\n")
	elif RERUN:
		INPUT = re.sub(".shhh.fasta", "", INPUT)
		print INPUT

		if Debug:
        		print "The basename of your input file is:\n"
		        error.write(INPUT+"\n")
	else:
		INPUT = re.sub(".fasta", "", INPUT)
		print INPUT

		if Debug:
        		print "The basename of your input file is:\n"
		        error.write(INPUT+"\n")

	OLIGOS = re.sub(".oligos", "", OLIGOS)
	print OLIGOS

#Process and merge output from .sff files
if MIX and SFF:
	for FILE in INPUT:
		if Debug:
			error.write("You are now processing: "+FILE+"\n\n")
		## The chemistry CAN be different for various sequencing runs, MOTHUR provides for this based on the "ORDER" of flows are read. This is sometimes necessary to specify.
        	if ORDER:
			os.system(' '.join([
				"for n in",
				FILE+".sff; do mothur \"# sffinfo(sff=$n, flow=T)\"; done"
			]))
			
			if Debug:
				error.write(' '.join([
                                	"for n in",
                                	FILE+".sff; do mothur \"# sffinfo(sff=$n, flow=T)\"; done"
                        	])+"\n")


			os.system(' '.join([
				"for n in",
				FILE+".flow; do mothur \"# trim.flows(flow=$n,",
				"oligos="+FILE+".oligos,",
				"pdiffs=2, bdiffs=1, minflows=200, maxflows=450,",
	                        "order="+ORDER+",",
				"processors="+PROCESSORS+")\"; done"
			]))

			if Debug:
				error.write(' '.join([
                        	        "for n in",
                                	FILE+".flow; do mothur \"# trim.flows(flow=$n,",
	                                "oligos="+FILE+".oligos,"
        	                        "pdiffs=2, bdiffs=1, minflows=200, maxflows=450,",
                	                "order="+ORDER+",",
                        	        "processors="+PROCESSORS+")\"; done"
	                        ])+"\n")

        	        os.system(' '.join([
                	        "for n in",
                        	FILE+".flow.files; do mothur \"# shhh.flows(file=$n,",
	                        "order="+ORDER+",",
				"processors="+PROCESSORS+")\"; done"
                	]))

			if Debug:
				error.write(' '.join([
                	                "for n in",
        	                        FILE+".flow.files; do mothur \"# shhh.flows(file=$n,",
                        	        "order="+ORDER+",",
	                                "processors="+PROCESSORS+")\"; done"
        	                ])+"\n")

	        else:
        	        os.system(' '.join([
                	        "for n in",
                        	FILE+".sff; do mothur \"# sffinfo(sff=$n, flow=T)\"; done"
	                ]))

			if Debug:
				error.write(' '.join([
	                                "for n in",
        	                        FILE+".sff; do mothur \"# sffinfo(sff=$n, flow=T)\"; done"
                	        ])+"\n")

        	        os.system(' '.join([
                	        "for n in",
                        	FILE+".flow; do mothur \"# trim.flows(flow=$n,",
	                        "oligos="+FILE+".oligos,",
        	                "pdiffs=2, bdiffs=1, minflows=200, maxflows=450,",
				"processors="+PROCESSORS+")\"; done"
	                ]))

			if Debug:
				error.write(' '.join([
                                	"for n in",
	                                FILE+".flow; do mothur \"# trim.flows(flow=$n,",
        	                        "oligos="+FILE+".oligos,"
                	                "pdiffs=2, bdiffs=1, minflows=200, maxflows=450,",
                        	        "processors="+PROCESSORS+")\"; done"
                        	])+"\n")

        	        os.system(' '.join([
                	        "for n in",
                        	FILE+".flow.files; do mothur \"# shhh.flows(file=$n,",
				"processors="+PROCESSORS+")\"; done"
        	        ]))

			if Debug:
				error.write(' '.join([
	                                "for n in",
        	                        FILE+".flow.files; do mothur \"# shhh.flows(file=$n,",
                	                "processors="+PROCESSORS+")\"; done"
                        	])+"\n")

if MIX and not SFF and not RERUN:
	NO_SHHH = "TRUE"
	print "NO_SHHH set to TRUE"

################################################################################
# Deal with the fact that the same barcode sequences may be used in multiple runs
################################################################################

#Strategy is to get a list of duplicates, and run through them by changing the oligo's file and then the .fasta sequence with randomly generated sequences
if MIX:
	#Concatenate OLIGOS Files
	SAMPLE_LOCATION_DICT = {}	

	#Do concatenation
	COMBO_NAME = NAME+"_Combined_Libraries"

	if Debug:
		error.write("Your files hae been concatenated with the base name: "+COMBO_NAME+"\n")

	count = 0
	with open(MIX+'/'+COMBO_NAME+'.oligos', 'w') as outfile:
		for fname in INPUT:
			if not re.search("Combined_Libraries", fname):
				with open(MIX+"/"+fname+".oligos") as infile:
					if count != 0:
						next(infile)
					
					for line in infile:
						line = line.strip("\r\n")

						#Write to new concatenated oligos file
						outfile.write(line+"\n")
						line = line.split()

						#Store which pyrotag sample came from which oligos file
						if np.size(line)>2:
							fname = re.sub("\\./","",fname)
							SAMPLE_LOCATION_DICT[line[2]] = fname 
							
					count = count + 1

	## Check for Duplicates
	DUPLICATE_DICT = {}
	DUPLICATE_LIST = []
	OLIGOS = MIX+'/'+COMBO_NAME

	with open(OLIGOS+".oligos") as infile:
        	#Grab the first 9 nucleotides of the primer for future use
	        PRIMER = infile.readline()
		PRIMER = PRIMER.strip()
        	SPLIT = PRIMER.split()
	        PRIMER_8 = SPLIT[1][0:8]

		#Read in all barcodes & sample IDs
		for line in infile:
			line = line.strip()
			line = line.split()
	
			#Make list of all barcodes
			DUPLICATE_LIST.append(line[1])
	
			#Make dictionary of sample name and barcode
			if not DUPLICATE_DICT.has_key(line[1]):
				DUPLICATE_DICT[line[1]] = [line[2]]
			else:
				DUPLICATE_DICT[line[1]].append(line[2])

	#Get list of duplicates
	DUPLICATE_BARCODES = list_duplicates(DUPLICATE_LIST)

	#Change sequence barcodes
	mock_oligos = open("TEMP.oligos", "w")
	mock_oligos.write(PRIMER+"\n")
	
	with open(OLIGOS+".oligos") as infile:
		next(infile)

		for line in infile:
			DUPLICATE_SEQ = line.strip("\r\n")
			DUPLICATE_SEQ = DUPLICATE_SEQ.split()
			DUPLICATE_SEQ = DUPLICATE_SEQ[1]

			#See if sequence is duplicated
			if DUPLICATE_SEQ in DUPLICATE_BARCODES:

				if np.size(DUPLICATE_DICT[DUPLICATE_SEQ]) > 1:
					DUPLICATE_ID = DUPLICATE_DICT[DUPLICATE_SEQ][0]
	
					print "Now Substituting Barcodes Found in Sample: "+DUPLICATE_ID+" due to overlap with another sample.\n"
						
					#Only substitute barcode if there is > 1 instance
					if re.search(DUPLICATE_ID, line):
						BARCODE_LENGTH = len(DUPLICATE_SEQ)
						NEW_BARCODE = id_generator(BARCODE_LENGTH, "TCGA")
	
						line = re.sub(DUPLICATE_SEQ, NEW_BARCODE, line)
						mock_oligos.write(line)
	
						#Remove element just in case there are more than two duplications (this would mean the duplicate list
						#would contain multiples of hte same DUPLICATE_SEQ and you'll cycle through until that list is exhausted
						del DUPLICATE_DICT[DUPLICATE_SEQ][0]
						
						if Debug:
							error.write("The Length of NO_SHHH is: "+str(len(NO_SHHH))+"\n")
	
						## Find Correct NAME file
						CORRECT_NAME = SAMPLE_LOCATION_DICT[DUPLICATE_ID]
						
						if len(NO_SHHH) > 1:

							#Replace all instances in the fasta file with sed
							os.system(' '.join([
								"sed",
								"-i",
								"\'s/"+DUPLICATE_SEQ+PRIMER_8+"/"+NEW_BARCODE+PRIMER_8+"/g\'",
								MIX+'/'+CORRECT_NAME+".fasta"
							]))
	
							if Debug:
								error.write(' '.join([
									"sed",
									"-i",
									"\'s/"+DUPLICATE_SEQ+PRIMER_8+"/"+NEW_BARCODE+PRIMER_8+"/g\'",
									MIX+'/'+CORRECT_NAME+".fasta"+"\n"
								]))
	
						else:
							os.system(' '.join([
								"sed",
								"-i",
								"\'s/"+DUPLICATE_SEQ+PRIMER_8+"/"+NEW_BARCODE+PRIMER_8+"/g\'",
								MIX+'/'+CORRECT_NAME+".shhh.fasta"
							]))
	
							if Debug:
								error.write(' '.join([
									"sed",
									"-i",
									"\'s/"+DUPLICATE_SEQ+PRIMER_8+"/"+NEW_BARCODE+PRIMER_8+"/g\'",
									MIX+'/'+CORRECT_NAME+".shhh.fasta"+"\n"
								]))
	
						error.write("Finished Processing: "+DUPLICATE_ID+"\n")
					else:
						mock_oligos.write(line)
				else:
					mock_oligos.write(line)
	
			#Write oligos that have nothing to do with duplicates
			else:
				mock_oligos.write(line)
	
		mock_oligos.close()
	
		#Re-name old oligos file and new oligos file
		os.system(' '.join([
			"mv",
			MIX+'/'+COMBO_NAME+".oligos",
			MIX+'/'+COMBO_NAME+".original.oligos"
		]))
	
		os.system(' '.join([
			"mv",
			"TEMP.oligos",
			MIX+'/'+COMBO_NAME+".oligos"
		]))

## Concatenate multiple files into one
if MIX and SFF or MIX and RERUN:

	#FASTA
	with open(MIX+'/'+COMBO_NAME+'.shhh.fasta', 'w') as outfile:
	    for fname in INPUT:
		if not re.search("Combined_Libraries", fname):
	        	with open(fname+".shhh.fasta") as infile:
		            for line in infile:
				line = line.strip("\r\n")
        	        	outfile.write(line+"\n")

	#NAMES
	with open(MIX+'/'+COMBO_NAME+'.shhh.names', 'w') as outfile:
	    for fname in INPUT:
		if not re.search("Combined_Libraries", fname):
	        	with open(fname+".shhh.names") as infile:
		            for line in infile:
				line = line.strip("\r\n")
        	        	outfile.write(line+"\n")

	OLIGOS = MIX+'/'+COMBO_NAME
	INPUT = MIX+'/'+COMBO_NAME
	
	if Debug:
		error.write(str(OLIGOS))

# Do concatenation of other files
elif MIX and not SFF and not RERUN:
	#Do concatenation
	COMBO_NAME = NAME+"_Combined_Libraries"

	#FASTA
	with open(MIX+'/'+COMBO_NAME+'.fasta', 'w') as outfile:
	    for fname in INPUT:
        	with open(fname+".fasta") as infile:
	            for line in infile:
			line = line.strip("\r\n")
        	        outfile.write(line+"\n")

	#NAMES
	with open(MIX+'/'+COMBO_NAME+'.names', 'w') as outfile:
	    for fname in INPUT:
        	with open(fname+".names") as infile:
	            for line in infile:
			line = line.strip("\r\n")
        	        outfile.write(line+"\n")

	OLIGOS = MIX+'/'+COMBO_NAME
	INPUT = MIX+'/'+COMBO_NAME

	if Debug:
		error.write(str(OLIGOS))

if not MIX and SFF:
	## The chemistry CAN be different for various sequencing runs, MOTHUR provides for this based on the "ORDER" of flows are read. This is sometimes necessary to specify.
        if ORDER:
		os.system(' '.join([
			"for n in",
			INPUT+".sff; do mothur \"# sffinfo(sff=$n, flow=T)\"; done"
		]))

		if Debug:
			error.write(' '.join([
                   		"for n in",
	                        INPUT+".sff; do mothur \"# sffinfo(sff=$n, flow=T)\"; done\n"
        	        ]))

		os.system(' '.join([
			"for n in",
			INPUT+".flow; do mothur \"# trim.flows(flow=$n,",
			"oligos="+OLIGOS+".oligos,"
			"pdiffs=2, bdiffs=1, minflows=200, maxflows=450,",
                        "order="+ORDER+",",
			"processors="+PROCESSORS+")\"; done"
		]))

		if Debug:
			error.write(' '.join([
	                        "for n in",
        	                INPUT+".flow; do mothur \"# trim.flows(flow=$n,",
                	        "oligos="+OLIGOS+".oligos,"
                        	"pdiffs=2, bdiffs=1, minflows=200, maxflows=450,",
	                        "order="+ORDER+",",
        	                "processors="+PROCESSORS+")\"; done\n"
                	]))

                os.system(' '.join([
                        "for n in",
                        INPUT+".flow.files; do mothur \"# shhh.flows(file=$n,",
                        "order="+ORDER+",",
			"processors="+PROCESSORS+")\"; done"
                ]))

		if Debug:
			error.write(' '.join([
	                        "for n in",
        	                INPUT+".flow.files; do mothur \"# shhh.flows(file=$n,",
                	        "order="+ORDER+",",
                        	"processors="+PROCESSORS+")\"; done\n"
	                ]))

        else:
                os.system(' '.join([
                        "for n in",
                        INPUT+".sff; do mothur \"# sffinfo(sff=$n, flow=T)\"; done"
                ]))

		if Debug:
			error.write(' '.join([
                   	     "for n in",
	                        INPUT+".sff; do mothur \"# sffinfo(sff=$n, flow=T)\"; done\n"
        	        ]))

                os.system(' '.join([
                        "for n in",
                        INPUT+".flow; do mothur \"# trim.flows(flow=$n,",
                        "oligos="+OLIGOS+".oligos,"
                        "pdiffs=2, bdiffs=1, minflows=200, maxflows=450,",
			"processors="+PROCESSORS+")\"; done"
                ]))

		if Debug:
			error.write(' '.join([
        	                "for n in",
	                        INPUT+".flow; do mothur \"# trim.flows(flow=$n,",
                        	"oligos="+OLIGOS+".oligos,"
                	        "pdiffs=2, bdiffs=1, minflows=200, maxflows=450,",
        	                "processors="+PROCESSORS+")\"; done\n"
	                ]))

                os.system(' '.join([
                        "for n in",
                        INPUT+".flow.files; do mothur \"# shhh.flows(file=$n,",
			"processors="+PROCESSORS+")\"; done"
                ]))

		if Debug:
			error.write(' '.join([
	                        "for n in",
        	                INPUT+".flow.files; do mothur \"# shhh.flows(file=$n,",
                	        "processors="+PROCESSORS+")\"; done\n"
	                ]))

if SFF or RERUN:
	os.system(' '.join([
		"for n in",
		INPUT+".shhh.fasta; do mothur \"# summary.seqs(fasta=$n,",
		"name="+INPUT+".shhh.names,",
		"processors="+PROCESSORS+")\"; done"
	]))

	if Debug:
		error.write(' '.join([
	                "for n in",
        	        INPUT+".shhh.fasta; do mothur \"# summary.seqs(fasta=$n,",
                	"name="+INPUT+".shhh.names,",
	                "processors="+PROCESSORS+")\"; done\n"
        	]))

	os.system(' '.join([
		"for n in",
		INPUT+".shhh.fasta; do mothur \"# trim.seqs(fasta=$n,",
		"oligos="+OLIGOS+".oligos, name="+INPUT+".shhh.names,",
		"maxambig=0, maxhomop=8, bdiffs=1, pdiffs=2, minlength=200,",
		"processors="+PROCESSORS+")\"; done"
	]))

	if Debug:
		error.write(' '.join([
	                "for n in",
        	        INPUT+".shhh.fasta; do mothur \"# trim.seqs(fasta=$n,",
                	"oligos="+OLIGOS+".oligos, name="+INPUT+".shhh.names,",
	                "maxambig=0, maxhomop=8, bdiffs=1, pdiffs=2, minlength=200,",
        	        "processors="+PROCESSORS+")\"; done\n"
	        ]))

	os.system(' '.join([
		"for n in",
		INPUT+".shhh.trim.fasta; do mothur \"# unique.seqs(fasta=$n,",
		"name="+INPUT+".shhh.names)\"; done"
	]))

	if Debug:
		error.write(' '.join([
	                "for n in",
        	        INPUT+".shhh.trim.fasta; do mothur \"# unique.seqs(fasta=$n,",
                	"name="+INPUT+".shhh.names)\"; done\n"
	]))

	os.system(' '.join([
		"for n in",
		INPUT+".shhh.trim.unique.fasta; do mothur \"# align.seqs(fasta=$n,",
		"reference=~/Phylogenetic_Gene_Databases/silva.bacteria/silva.bacteria.fasta, flip=T,",
		"processors="+PROCESSORS+")\"; done"
	]))

	if Debug:
		error.write(' '.join([
	                "for n in",
        	        INPUT+".shhh.trim.unique.fasta; do mothur \"# align.seqs(fasta=$n,",
                	"reference=~/Phylogenetic_Gene_Databases/silva.bacteria/silva.bacteria.fasta, flip=T,",
	                "processors="+PROCESSORS+")\"; done\n"
        	]))

	os.system(' '.join([
		"for n in",
		INPUT+".shhh.trim.unique.align; do mothur \"# summary.seqs(fasta=$n,",
		"name="+INPUT+".shhh.trim.names,",
		"processors="+PROCESSORS+")\"; done"
	]))

	if Debug:
		error.write(' '.join([
	                "for n in",
        	        INPUT+".shhh.trim.unique.align; do mothur \"# summary.seqs(fasta=$n,",
                	"name="+INPUT+".shhh.trim.names,",
	                "processors="+PROCESSORS+")\"; done\n"
        	]))

	#Get user feedback
	print "Based on the summary above, please choose a starting position and minimum length.\n"
	print "ENTER the starting position."
	start = raw_input()
	print "ENTER the minimum length."
	minlength = raw_input()

	if Debug:
		error.write("You selected to trim your sequences from starting point: "+str(start)+" with a minimum length of: "+minlength+".\n")

	os.system(' '.join([
		"for n in",
		INPUT+".shhh.trim.unique.align; do mothur \"# screen.seqs(fasta=$n,",
		"name="+INPUT+".shhh.trim.names,",
		"group="+INPUT+".shhh.groups,",
		"start="+start+",",
		"minlength="+minlength+",",
		"processors="+PROCESSORS+")\"; done"
	]))

	if Debug:
		error.write(' '.join([
	                "for n in",
        	        INPUT+".shhh.trim.unique.align; do mothur \"# screen.seqs(fasta=$n,",
                	"name="+INPUT+".shhh.trim.names,",
	                "group="+INPUT+".shhh.groups,",
        	        "start="+start+",",
                	"minlength="+minlength+",",
	                "processors="+PROCESSORS+")\"; done\n"
        	]))

	os.system(' '.join([
		"for n in",
		INPUT+".shhh.trim.unique.good.align; do mothur \"# filter.seqs(fasta=$n,",
		"vertical=T,",
		"processors="+PROCESSORS+")\"; done"
	]))

	if Debug:
		error.write(' '.join([
	                "for n in",
        	        INPUT+".shhh.trim.unique.good.align; do mothur \"# filter.seqs(fasta=$n,",
                	"vertical=T,",
	                "processors="+PROCESSORS+")\"; done\n"
        	]))

	os.system(' '.join([
		"for n in",
		INPUT+".shhh.trim.unique.good.filter.fasta; do mothur \"# unique.seqs(fasta=$n,",
		"name="+INPUT+".shhh.trim.good.names)\"; done"
	]))

	if Debug:
		error.write(' '.join([
	                "for n in",
        	        INPUT+".shhh.trim.unique.good.filter.fasta; do mothur \"# unique.seqs(fasta=$n,",
                	"name="+INPUT+".shhh.trim.good.names)\"; done\n"
	        ]))

	os.system(' '.join([
		"for n in",
		INPUT+".shhh.trim.unique.good.filter.unique.fasta; do mothur \"# pre.cluster(fasta=$n,",
		"name="+INPUT+".shhh.trim.unique.good.filter.names,",
		"group="+INPUT+".shhh.good.groups, diffs=2,",
		"processors="+PROCESSORS+")\"; done"
	]))

	if Debug:
		error.write(' '.join([
	                "for n in",
        	        INPUT+".shhh.trim.unique.good.filter.unique.fasta; do mothur \"# pre.cluster(fasta=$n,",
                	"name="+INPUT+".shhh.trim.unique.good.filter.names,",
	                "group="+INPUT+".shhh.good.groups, diffs=2,",
        	        "processors="+PROCESSORS+")\"; done\n"
	        ]))


	os.system(' '.join([
		"for n in",
		INPUT+".shhh.trim.unique.good.filter.unique.precluster.fasta; do mothur \"# chimera.uchime(fasta=$n,",
		"name="+INPUT+".shhh.trim.unique.good.filter.unique.precluster.names,",
		"group="+INPUT+".shhh.good.groups,",
		"processors="+PROCESSORS+")\"; done"
	]))

	if Debug:
		error.write(' '.join([
	                "for n in",
        	        INPUT+".shhh.trim.unique.good.filter.unique.precluster.fasta; do mothur \"# chimera.uchime(fasta=$n,",
                	"name="+INPUT+".shhh.trim.unique.good.filter.unique.precluster.names,",
	                "group="+INPUT+".shhh.good.groups,",
        	        "processors="+PROCESSORS+")\"; done\n"
	        ]))

	os.system(' '.join([
		"for n in",
		INPUT+".shhh.trim.unique.good.filter.unique.precluster.uchime.accnos; do mothur \"# remove.seqs(accnos=$n,",
		"fasta="+INPUT+".shhh.trim.unique.good.filter.unique.precluster.fasta,",
		"name="+INPUT+".shhh.trim.unique.good.filter.unique.precluster.names,",
		"group="+INPUT+".shhh.good.groups)\"; done"
	]))

	if Debug:
		error.write(' '.join([
	                "for n in",
        	        INPUT+".shhh.trim.unique.good.filter.unique.precluster.uchime.accnos; do mothur \"# remove.seqs(accnos=$n,",
                	"fasta="+INPUT+".shhh.trim.unique.good.filter.unique.precluster.fasta,",
	                "name="+INPUT+".shhh.trim.unique.good.filter.unique.precluster.names,",
        	        "group="+INPUT+".shhh.good.groups)\"; done\n"
	        ]))

	os.system(' '.join([
		"for n in",
		INPUT+".shhh.trim.unique.good.filter.unique.precluster.pick.fasta; do mothur \"# summary.seqs(fasta=$n,",
		"name="+INPUT+".shhh.trim.unique.good.filter.unique.precluster.pick.names,",
		"processors="+PROCESSORS+")\"; done"
	]))

	if Debug:
		error.write(' '.join([
	                "for n in",
        	        INPUT+".shhh.trim.unique.good.filter.unique.precluster.pick.fasta; do mothur \"# summary.seqs(fasta=$n,",
                	"name="+INPUT+".shhh.trim.unique.good.filter.unique.precluster.pick.names,",
	                "processors="+PROCESSORS+")\"; done\n"
        	]))

	os.system(' '.join([
		"for n in",
		INPUT+".shhh.trim.unique.good.filter.unique.precluster.pick.fasta; do mothur \"# classify.seqs(fasta=$n,",
		"name="+INPUT+".shhh.trim.unique.good.filter.unique.precluster.pick.names,",
		"group="+INPUT+".shhh.good.pick.groups,",
		"template=~/Phylogenetic_Gene_Databases/silva.bacteria/silva.bacteria.fasta,",
		"taxonomy=~/Phylogenetic_Gene_Databases/silva.bacteria/silva.bacteria.silva.tax, cutoff=50,",
		"processors="+PROCESSORS+")\"; done"
	]))

	if Debug:
		error.write(' '.join([
	                "for n in",
        	        INPUT+".shhh.trim.unique.good.filter.unique.precluster.pick.fasta; do mothur \"# classify.seqs(fasta=$n,",
                	"name="+INPUT+".shhh.trim.unique.good.filter.unique.precluster.pick.names,",
	                "group="+INPUT+".shhh.good.pick.groups,",
        	        "template=~/Phylogenetic_Gene_Databases/silva.bacteria/silva.bacteria.fasta,",
                	"taxonomy=~/Phylogenetic_Gene_Databases/silva.bacteria/silva.bacteria.silva.tax, cutoff=50,",
	                "processors="+PROCESSORS+")\"; done\n"
        	]))

	os.system(' '.join([
		"for n in",
		INPUT+".shhh.trim.unique.good.filter.unique.precluster.pick.fasta; do mothur \"# remove.lineage(fasta=$n,",
		"name="+INPUT+".shhh.trim.unique.good.filter.unique.precluster.pick.names,",
		"group="+INPUT+".shhh.good.pick.groups,",
		"taxonomy="+INPUT+".shhh.trim.unique.good.filter.unique.precluster.pick.silva.wang.taxonomy,",
		"taxon=Mitochondria-Cyanobacteria_Chloroplast-Archaea-Eukarya-unknown)\"; done"
	]))

	if Debug:
		error.write(' '.join([
	                "for n in",
                	INPUT+".shhh.trim.unique.good.filter.unique.precluster.pick.fasta; do mothur \"# remove.lineage(fasta=$n,",
        	        "name="+INPUT+".shhh.trim.unique.good.filter.unique.precluster.pick.names,",
	                "group="+INPUT+".shhh.good.pick.groups,",
        	        "taxonomy="+INPUT+".shhh.trim.unique.good.filter.unique.precluster.pick.silva.wang.taxonomy,",
                	"taxon=Mitochondria-Cyanobacteria_Chloroplast-Archaea-Eukarya-unknown)\"; done\n"
	        ]))


	os.system(' '.join([
		"for n in",
		INPUT+".shhh.trim.unique.good.filter.unique.precluster.pick.pick.fasta; do mothur \"# classify.seqs(fasta=$n,",
		"name="+INPUT+".shhh.trim.unique.good.filter.unique.precluster.pick.pick.names,",
		"group="+INPUT+".shhh.good.pick.pick.groups,",
		"template=~/Phylogenetic_Gene_Databases/GreenGenes_MothurFormatted/gg_13_5_99.fasta,",
		"taxonomy=~/Phylogenetic_Gene_Databases/GreenGenes_MothurFormatted/gg_13_5_99.gg.tax, cutoff=80,",
		"processors="+PROCESSORS+")\"; done"
	]))

	if Debug:
		error.write(' '.join([
	                "for n in",
        	        INPUT+".shhh.trim.unique.good.filter.unique.precluster.pick.pick.fasta; do mothur \"# classify.seqs(fasta=$n,",
                	"name="+INPUT+".shhh.trim.unique.good.filter.unique.precluster.pick.pick.names,",
	                "group="+INPUT+".shhh.good.pick.pick.groups,",
        	        "template=~/Phylogenetic_Gene_Databases/GreenGenes_MothurFormatted/gg_13_5_99.fasta,",
                	"taxonomy=~/Phylogenetic_Gene_Databases/GreenGenes_MothurFormatted/gg_13_5_99.gg.tax, cutoff=80,",
	                "processors="+PROCESSORS+")\"; done\n"
        	]))

	##Rename Files to More Readable Names
	os.system(' '.join([
		"cp",
		INPUT+".shhh.trim.unique.good.filter.unique.precluster.pick.pick.names",
		NAME+"_final.names"
	]))

	if Debug:
		error.write(' '.join([
	                "cp",
        	        INPUT+".shhh.trim.unique.good.filter.unique.precluster.pick.pick.names",
                	NAME+"_final.names\n"
	        ]))

	os.system(' '.join([
		"cp",
		INPUT+".shhh.trim.unique.good.filter.unique.precluster.pick.pick.fasta",
		NAME+"_final.fasta"
	]))

	if Debug:
		error.write(' '.join([
           	     "cp",
	                INPUT+".shhh.trim.unique.good.filter.unique.precluster.pick.pick.fasta",
        	        NAME+"_final.fasta\n"
	        ]))

	os.system(' '.join([
		"cp",
		INPUT+".shhh.good.pick.pick.groups",
		NAME+"_final.groups"
	]))

	if Debug:
		error.write(' '.join([
           	     "cp",
	                INPUT+".shhh.good.pick.pick.groups",
        	        NAME+"_final.groups\n"
	        ]))

	os.system(' '.join([
		"cp",
		INPUT+".shhh.trim.unique.good.filter.unique.precluster.pick.silva.wang.pick.taxonomy",
		NAME+"_final.silva.taxonomy"
	]))

	if Debug:
		error.write(' '.join([
	                "cp",
        	        INPUT+".shhh.trim.unique.good.filter.unique.precluster.pick.silva.wang.pick.taxonomy",
                	NAME+"_final.silva.taxonomy\n"
	        ]))

	os.system(' '.join([
		"cp",
		INPUT+".shhh.trim.unique.good.filter.unique.precluster.pick.pick.gg.wang.taxonomy",
		NAME+"_final.gg.taxonomy"
	]))

	if Debug:
		error.write(' '.join([
	                "cp",
        	        INPUT+".shhh.trim.unique.good.filter.unique.precluster.pick.pick.gg.wang.taxonomy",
                	NAME+"_final.gg.taxonomy\n"
	        ]))

	#Make List File
	os.system(' '.join([
		"for n in",
		NAME+"_final.fasta; do mothur \"# dist.seqs(fasta=$n,",
		"cutoff=0.3,",
		"processors="+PROCESSORS+")\"; done"
	]))

	if Debug:
		error.write(' '.join([
	                "for n in",
        	        NAME+"_final.fasta; do mothur \"# dist.seqs(fasta=$n,",
                	"cutoff=0.3,",
	                "processors="+PROCESSORS+")\"; done\n"
        	]))

	os.system(' '.join([
		"for n in",
		NAME+"_final.fasta; do mothur \"# cluster.split(fasta=$n,",
		"name="+NAME+"_final.names,",
		"taxonomy="+NAME+"_final.silva.taxonomy,",
		"splitmethod=classify, taxlevel=3,",
		"processors="+PROCESSORS+")\"; done"
	]))

	if Debug:
		error.write(' '.join([
           	     "for n in",
	                NAME+"_final.fasta; do mothur \"# cluster.split(fasta=$n,",
	      	   	"name="+NAME+"_final.names,",
	                "taxonomy="+NAME+"_final.silva.taxonomy,",
        	        "splitmethod=classify, taxlevel=3,",
	                "processors="+PROCESSORS+")\"; done\n"
        	]))

	#Move all Final Files to Output directory
	os.system(' '.join([
		"mv",
		NAME+"_final.*",
		"./"+OUTPUT+"/"
	]))

	if Debug:
		error.write(' '.join([
	                "mv",
        	        NAME+"_final.*",
                	"./"+OUTPUT+"/\n"
	        ]))

	## Cat all logfiles in order of creation and move
        os.system(' '.join([
        	"cat",
                "$(ls -t mothur.*)",
                ">",
                "./"+OUTPUT+"/"+NAME+".mothur.logfiles"
	]))

	if Debug:
		error.write(' '.join([
	                "cat",
        	        "$(ls -t mothur.*)",
                	">",
	                "./"+OUTPUT+"/"+NAME+".mothur.logfiles\n"
        	]))

	## Provide a version of the GG taxonomy file acceptable for importing into R
	## I do not use Silva for this b/c it does not have consistent delimiting according to taxonomic levels 
        os.system(' '.join([
        	"cp",
                "./"+OUTPUT+"/"+NAME+"_final.gg.taxonomy",
                "./"+OUTPUT+"/"+NAME+"_final.gg.R.taxonomy"
	]))

	if Debug:
		error.write(' '.join([
	                "cp",
        	        "./"+OUTPUT+"/"+NAME+"_final.gg.taxonomy",
                	"./"+OUTPUT+"/"+NAME+"_final.gg.R.taxonomy\n"
	        ]))


        os.system(' '.join([
        	"sed -i 's/\t/;/g'",
                "./"+OUTPUT+"/"+NAME+"_final.gg.R.taxonomy"
	]))

	if Debug:
		error.write(' '.join([
           	     "sed -i 's/\t/;/g'",
                	"./"+OUTPUT+"/"+NAME+"_final.gg.R.taxonomy\n"
	        ]))

        os.system(' '.join([
        	"sed -i 's/;$//g'",
                "./"+OUTPUT+"/"+NAME+"_final.gg.R.taxonomy"
	]))

	if Debug:
		error.write(' '.join([
	                "sed -i 's/;$//g'",
        	        "./"+OUTPUT+"/"+NAME+"_final.gg.R.taxonomy\n"
	        ]))

        os.system(' '.join([
        	"rm",
                "./mothur*.logfile"
	]))

	if Debug:
		error.write(' '.join([
	        	"rm",
        	        "./mothur*.logfile\n"
	        ]))

else:
	os.system(' '.join([
		"for n in",
		INPUT+".fasta; do mothur \"# summary.seqs(fasta=$n,",
		"name="+INPUT+".names,",
		"processors="+PROCESSORS+")\"; done"
	]))

	if Debug:
		error.write(' '.join([
	                "for n in",
        	        INPUT+".fasta; do mothur \"# summary.seqs(fasta=$n,",
                	"name="+INPUT+".names,",
	                "processors="+PROCESSORS+")\"; done\n"
        	]))

	os.system(' '.join([
		"for n in",
		INPUT+".fasta; do mothur \"# trim.seqs(fasta=$n,",
		"oligos="+OLIGOS+".oligos,",
		"qfile="+INPUT+".qual,",
		"maxambig=0, maxhomop=8, bdiffs=1, pdiffs=2, qwindowaverage=35, qwindowsize=50, minlength=200,",
		"processors="+PROCESSORS+")\"; done"
	]))

	if Debug:
		error.write(' '.join([
			"for n in",
			INPUT+".fasta; do mothur \"# trim.seqs(fasta=$n,",
			"oligos="+OLIGOS+".oligos,",
			"qfile="+INPUT+".qual,",
			"maxambig=0, maxhomop=8, bdiffs=1, pdiffs=2, qwindowaverage=35, qwindowsize=50, minlength=200,",
			"processors="+PROCESSORS+")\"; done\n"
		]))

	os.system(' '.join([
		"for n in",
		INPUT+".trim.fasta; do mothur \"# unique.seqs(fasta=$n,",
		"name="+INPUT+".names)\"; done"
	]))

	if Debug:
		error.write(' '.join([
			"for n in",
			INPUT+".trim.fasta; do mothur \"# unique.seqs(fasta=$n,",
			"name="+INPUT+".names)\"; done\n"
		]))

	os.system(' '.join([
		"for n in",
		INPUT+".trim.unique.fasta; do mothur \"# align.seqs(fasta=$n,",
		"reference=~/Phylogenetic_Gene_Databases/silva.bacteria/silva.bacteria.fasta, flip=T,",
		"processors="+PROCESSORS+")\"; done"
	]))

	if Debug:
		error.write(' '.join([
			"for n in",
			INPUT+".trim.unique.fasta; do mothur \"# align.seqs(fasta=$n,",
			"reference=~/Phylogenetic_Gene_Databases/silva.bacteria/silva.bacteria.fasta, flip=T,",
			"processors="+PROCESSORS+")\"; done\n"
		]))

	os.system(' '.join([
		"for n in",
		INPUT+".trim.unique.align; do mothur \"# summary.seqs(fasta=$n,",
		"name="+INPUT+".trim.names,",
		"processors="+PROCESSORS+")\"; done"
	]))

	if Debug:
		error.write(' '.join([
			"for n in",
			INPUT+".trim.unique.align; do mothur \"# summary.seqs(fasta=$n,",
			"name="+INPUT+".trim.names,",
			"processors="+PROCESSORS+")\"; done\n"
		]))

	#Get user feedback
	print "Based on the summary above, please choose a starting position and minimum length.\n"
	print "ENTER the starting position."
	start = raw_input()
	print "ENTER the minimum length."
	minlength = raw_input()

	if Debug:
		error.write("You selected to trim your sequences from starting point: "+str(start)+" with a minimum length of: "+minlength+".\n")

	os.system(' '.join([
		"for n in",
		INPUT+".trim.unique.align; do mothur \"# screen.seqs(fasta=$n,",
		"name="+INPUT+".trim.names,",
		"group="+INPUT+".groups,",
		"start="+start+",",
		"minlength="+minlength+",",
		"processors="+PROCESSORS+")\"; done"
	]))

	if Debug:
		error.write(' '.join([
			"for n in",
			INPUT+".trim.unique.align; do mothur \"# screen.seqs(fasta=$n,",
			"name="+INPUT+".trim.names,",
			"group="+INPUT+".groups,",
			"start="+start+",",
			"minlength="+minlength+",",
			"processors="+PROCESSORS+")\"; done\n"
		]))

	os.system(' '.join([
		"for n in",
		INPUT+".trim.unique.good.align; do mothur \"# filter.seqs(fasta=$n,",
		"vertical=T,",
		"processors="+PROCESSORS+")\"; done"
	]))

	if Debug:
		error.write(' '.join([
			"for n in",
			INPUT+".trim.unique.good.align; do mothur \"# filter.seqs(fasta=$n,",
			"vertical=T,",
			"processors="+PROCESSORS+")\"; done\n"
		]))

	os.system(' '.join([
		"for n in",
		INPUT+".trim.unique.good.filter.fasta; do mothur \"# unique.seqs(fasta=$n,",
		"name="+INPUT+".trim.good.names)\"; done"
	]))

	if Debug:
		error.write(' '.join([
			"for n in",
			INPUT+".trim.unique.good.filter.fasta; do mothur \"# unique.seqs(fasta=$n,",
			"name="+INPUT+".trim.good.names)\"; done\n"
		]))

	os.system(' '.join([
		"for n in",
		INPUT+".trim.unique.good.filter.unique.fasta; do mothur \"# pre.cluster(fasta=$n,",
		"name="+INPUT+".trim.unique.good.filter.names,",
		"group="+INPUT+".good.groups, diffs=2,",
		"processors="+PROCESSORS+")\"; done"
	]))

	if Debug:
		error.write(' '.join([
			"for n in",
			INPUT+".trim.unique.good.filter.unique.fasta; do mothur \"# pre.cluster(fasta=$n,",
			"name="+INPUT+".trim.unique.good.filter.names,",
			"group="+INPUT+".good.groups, diffs=2,",
			"processors="+PROCESSORS+")\"; done\n"
		]))

	os.system(' '.join([
		"for n in",
		INPUT+".trim.unique.good.filter.unique.precluster.fasta; do mothur \"# chimera.uchime(fasta=$n,",
		"name="+INPUT+".trim.unique.good.filter.unique.precluster.names,",
		"group="+INPUT+".good.groups,",
		"processors="+PROCESSORS+")\"; done"
	]))

	if Debug:
		error.write(' '.join([
			"for n in",
			INPUT+".trim.unique.good.filter.unique.precluster.fasta; do mothur \"# chimera.uchime(fasta=$n,",
			"name="+INPUT+".trim.unique.good.filter.unique.precluster.names,",
			"group="+INPUT+".good.groups,",
			"processors="+PROCESSORS+")\"; done\n"
		]))

	os.system(' '.join([
		"for n in",
		INPUT+".trim.unique.good.filter.unique.precluster.uchime.accnos; do mothur \"# remove.seqs(accnos=$n,",
		"fasta="+INPUT+".trim.unique.good.filter.unique.precluster.fasta,",
		"name="+INPUT+".trim.unique.good.filter.unique.precluster.names,",
		"group="+INPUT+".good.groups)\"; done"
	]))

	if Debug:
		error.write(' '.join([
			"for n in",
			INPUT+".trim.unique.good.filter.unique.precluster.uchime.accnos; do mothur \"# remove.seqs(accnos=$n,",
			"fasta="+INPUT+".trim.unique.good.filter.unique.precluster.fasta,",
			"name="+INPUT+".trim.unique.good.filter.unique.precluster.names,",
			"group="+INPUT+".good.groups)\"; done\n"
		]))

	os.system(' '.join([
		"for n in",
		INPUT+".trim.unique.good.filter.unique.precluster.pick.fasta; do mothur \"# summary.seqs(fasta=$n,",
		"name="+INPUT+".trim.unique.good.filter.unique.precluster.pick.names,",
		"processors="+PROCESSORS+")\"; done"
	]))

	if Debug:
		error.write(' '.join([
			"for n in",
			INPUT+".trim.unique.good.filter.unique.precluster.pick.fasta; do mothur \"# summary.seqs(fasta=$n,",
			"name="+INPUT+".trim.unique.good.filter.unique.precluster.pick.names,",
			"processors="+PROCESSORS+")\"; done\n"
		]))

	os.system(' '.join([
		"for n in",
		INPUT+".trim.unique.good.filter.unique.precluster.pick.fasta; do mothur \"# classify.seqs(fasta=$n,",
		"name="+INPUT+".trim.unique.good.filter.unique.precluster.pick.names,",
		"group="+INPUT+".good.pick.groups,",
		"template=~/Phylogenetic_Gene_Databases/silva.bacteria/silva.bacteria.fasta,",
		"taxonomy=~/Phylogenetic_Gene_Databases/silva.bacteria/silva.bacteria.silva.tax, cutoff=50,",
		"processors="+PROCESSORS+")\"; done"
	]))

	if Debug:
		error.write(' '.join([
			"for n in",
			INPUT+".trim.unique.good.filter.unique.precluster.pick.fasta; do mothur \"# classify.seqs(fasta=$n,",
			"name="+INPUT+".trim.unique.good.filter.unique.precluster.pick.names,",
			"group="+INPUT+".good.pick.groups,",
			"template=~/Phylogenetic_Gene_Databases/silva.bacteria/silva.bacteria.fasta,",
			"taxonomy=~/Phylogenetic_Gene_Databases/silva.bacteria/silva.bacteria.silva.tax, cutoff=50,",
			"processors="+PROCESSORS+")\"; done\n"
		]))
		
	os.system(' '.join([
		"for n in",
		INPUT+".trim.unique.good.filter.unique.precluster.pick.fasta; do mothur \"# remove.lineage(fasta=$n,",
		"name="+INPUT+".trim.unique.good.filter.unique.precluster.pick.names,",
		"group="+INPUT+".good.pick.groups,",
		"taxonomy="+INPUT+".trim.unique.good.filter.unique.precluster.pick.silva.wang.taxonomy,",
		"taxon=Mitochondria-Cyanobacteria_Chloroplast-Archaea-Eukarya-unknown)\"; done"
	]))

	if Debug:
		error.write(' '.join([
			"for n in",
			INPUT+".trim.unique.good.filter.unique.precluster.pick.fasta; do mothur \"# remove.lineage(fasta=$n,",
			"name="+INPUT+".trim.unique.good.filter.unique.precluster.pick.names,",
			"group="+INPUT+".good.pick.groups,",
			"taxonomy="+INPUT+".trim.unique.good.filter.unique.precluster.pick.silva.wang.taxonomy,",
			"taxon=Mitochondria-Cyanobacteria_Chloroplast-Archaea-Eukarya-unknown)\"; done\n"
		]))

	os.system(' '.join([
		"for n in",
		INPUT+".trim.unique.good.filter.unique.precluster.pick.pick.fasta; do mothur \"# classify.seqs(fasta=$n,",
		"name="+INPUT+".trim.unique.good.filter.unique.precluster.pick.pick.names,",
		"group="+INPUT+".good.pick.pick.groups,",
		"template=~/Phylogenetic_Gene_Databases/GreenGenes_MothurFormatted/gg_13_5_99.fasta,",
		"taxonomy=~/Phylogenetic_Gene_Databases/GreenGenes_MothurFormatted/gg_13_5_99.gg.tax, cutoff=80,",
		"processors="+PROCESSORS+")\"; done"
	]))

	if Debug:
		error.write(' '.join([
			"for n in",
			INPUT+".trim.unique.good.filter.unique.precluster.pick.pick.fasta; do mothur \"# classify.seqs(fasta=$n,",
			"name="+INPUT+".trim.unique.good.filter.unique.precluster.pick.pick.names,",
			"group="+INPUT+".good.pick.pick.groups,",
			"template=~/Phylogenetic_Gene_Databases/GreenGenes_MothurFormatted/gg_13_5_99.fasta,",
			"taxonomy=~/Phylogenetic_Gene_Databases/GreenGenes_MothurFormatted/gg_13_5_99.gg.tax, cutoff=80,",
			"processors="+PROCESSORS+")\"; done\n"
		]))

	##Rename Files to More Readable Names
	os.system(' '.join([
		"cp",
		INPUT+".trim.unique.good.filter.unique.precluster.pick.pick.names",
		NAME+"_final.names"
	]))

	if Debug:
		error.write(' '.join([
			"cp",
			INPUT+".trim.unique.good.filter.unique.precluster.pick.pick.names",
			NAME+"_final.names\n"
		]))
	
	os.system(' '.join([
		"cp",
		INPUT+".trim.unique.good.filter.unique.precluster.pick.pick.fasta",
		NAME+"_final.fasta"
	]))

	if Debug:
		error.write(' '.join([
			"cp",
			INPUT+".trim.unique.good.filter.unique.precluster.pick.pick.fasta",
			NAME+"_final.fasta\n"
		]))

	os.system(' '.join([
		"cp",
		INPUT+".good.pick.pick.groups",
		NAME+"_final.groups"
	]))

	if Debug:
		error.write(' '.join([
			"cp",
			INPUT+".good.pick.pick.groups",
			NAME+"_final.groups\n"
		]))

	os.system(' '.join([
		"cp",
		INPUT+".trim.unique.good.filter.unique.precluster.pick.silva.wang.pick.taxonomy",
		NAME+"_final.silva.taxonomy"
	]))

	if Debug:
		error.write(' '.join([
			"cp",
			INPUT+".trim.unique.good.filter.unique.precluster.pick.silva.wang.pick.taxonomy",
			NAME+"_final.silva.taxonomy\n"
		]))

	os.system(' '.join([
		"cp",
		INPUT+".trim.unique.good.filter.unique.precluster.pick.pick.gg.wang.taxonomy",
		NAME+"_final.gg.taxonomy"
	]))

	if Debug:
		error.write(' '.join([
			"cp",
			INPUT+".trim.unique.good.filter.unique.precluster.pick.pick.gg.wang.taxonomy",
			NAME+"_final.gg.taxonomy\n"
		]))

	#Make List File
	os.system(' '.join([
		"for n in",
		NAME+"_final.fasta; do mothur \"# dist.seqs(fasta=$n,",
		"cutoff=0.3,",
		"processors="+PROCESSORS+")\"; done"
	]))

	if Debug:
			error.write(' '.join([
			"for n in",
			NAME+"_final.fasta; do mothur \"# dist.seqs(fasta=$n,",
			"cutoff=0.3,",
			"processors="+PROCESSORS+")\"; done\n"
		]))

	os.system(' '.join([
		"for n in",
		NAME+"_final.fasta; do mothur \"# cluster.split(fasta=$n,",
		"name="+NAME+"_final.names,",
		"taxonomy="+NAME+"_final.silva.taxonomy,",
		"splitmethod=classify, taxlevel=3,",
		"processors="+PROCESSORS+")\"; done"
	]))

	if Debug:
		error.write(' '.join([
			"for n in",
			NAME+"_final.fasta; do mothur \"# cluster.split(fasta=$n,",
			"name="+NAME+"_final.names,",
			"taxonomy="+NAME+"_final.silva.taxonomy,",
			"splitmethod=classify, taxlevel=3,",
			"processors="+PROCESSORS+")\"; done\n"
		]))

	#Move all Final Files to Output directory
	os.system(' '.join([
		"mv",
		NAME+"_final.*",
		"./"+OUTPUT+"/"
	]))

	if Debug:
		error.write(' '.join([
			"mv",
			NAME+"_final.*",
			"./"+OUTPUT+"/\n"
		]))

	## Cat all logfiles in order of creation and move
        os.system(' '.join([
        	"cat",
                "$(ls -t mothur.*)",
                ">",
                "./"+OUTPUT+"/"+NAME+".mothur.logfiles"
	]))

	if Debug:
		error.write(' '.join([
			"cat",
			"$(ls -t mothur.*)",
			">",
			"./"+OUTPUT+"/"+NAME+".mothur.logfiles\n"
		]))

	## Provide a version of the GG taxonomy file acceptable for importing into R
	## I do not use Silva for this b/c it does not have consistent delimiting according to taxonomic levels 
        os.system(' '.join([
        	"cp",
                "./"+OUTPUT+"/"+NAME+"_final.gg.taxonomy",
                "./"+OUTPUT+"/"+NAME+"_final.gg.R.taxonomy"
	]))

	if Debug:
		error.write(' '.join([
			"cp",
			"./"+OUTPUT+"/"+NAME+"_final.gg.taxonomy",
			"./"+OUTPUT+"/"+NAME+"_final.gg.R.taxonomy\n"
		]))

        os.system(' '.join([
        	"sed -i 's/\t/;/g'",
                "./"+OUTPUT+"/"+NAME+"_final.gg.R.taxonomy"
	]))

	if Debug:
		error.write(' '.join([
			"sed -i 's/\t/;/g'",
			"./"+OUTPUT+"/"+NAME+"_final.gg.R.taxonomy\n"
		]))

        os.system(' '.join([
        	"sed -i 's/;$//g'",
                "./"+OUTPUT+"/"+NAME+"_final.gg.R.taxonomy"
	]))

	if Debug:
		error.write(' '.join([
			"sed -i 's/;$//g'",
			"./"+OUTPUT+"/"+NAME+"_final.gg.R.taxonomy\n"
		]))

        os.system(' '.join([
        	"rm",
                "./mothur*.logfile"
	]))

	if Debug:
		error.write(' '.join([
	        	"rm",
        	        "./mothur*.logfile\n"
	        ]))

end = timeit.default_timer()

print end - now

