#!/usr/bin/python
import sys, os, re, getopt, glob, numpy as np
import timeit
import itertools

now = timeit.default_timer()

Usage = """
Usage:  ./MOTHUR_16S_PIPELINE.py -o BACTERIA_0.05 -n BACTERIA_0.05 -i IIKFCBR01.fasta -f IIKFCBR01.oligos -p 8

or

Usage:  ./MOTHUR_16S_PIPELINE.py -o BACTERIA_0.05 -n BACTERIA_0.05 -m IIKFCBR01.fasta -p 8


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
			- Provide the FULL PATH (i.e. don't use "." to specify current directory)
			- the directory must also contain all the associated .oligos files for each run)

DEPENDENCY:
		-mothur v.1.32.1 (written using, probably valide for earlier and later versions)
		-Green Genes "Mothur Formatted" Database
		-Silva Bacteria Database
		

OPTIONAL:
		-S <Y>	Input is .sff file (or WAS derived from a .sff file)
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

Usage:  ./MOTHUR_16S_PIPELINE.py -o BACTERIA_0.05 -n BACTERIA_0.05 -i IIKFCBR01.fasta -f IIKFCBR01.oligos -p 8

or

Usage:  ./MOTHUR_16S_PIPELINE.py -o BACTERIA_0.05 -n BACTERIA_0.05 -m IIKFCBR01.fasta -p 8
"""

if len(sys.argv)<3:
        print Usage
        sys.exit()


# Store input and output file names
NAME=''
INPUT=[]
OLIGOS=[]
RERUN=''
SFF=''
ORDER=''
PROCESSORS=''
OUTPUT=''
MIX=''
Debug=''

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

if len(OUTPUT)>0:
        if os.path.exists('./' + OUTPUT):
                print "\nOutput Folder Exists - Caution: Files May Be Re-Written"
        else:
                os.mkdir(OUTPUT)

## Print Debug Info (not really used much in this script)
if Debug:
        print "You are in Debugging Mode"
        
        if len(OUTPUT)>0:
        	error = open(OUTPUT+"/command.log", "w")
        else:
        	error = open("command.log", "w")
        	
if MIX:
	INPUT = []
	if Debug:
       		print "You chose to merge multiple files for input.\n"

	if SFF and not RERUN:
		for FILE in glob.glob(MIX+"/*.sff"):
			FILE = re.sub(MIX, "", FILE)
                        FILE = re.sub("./", "", FILE)
			INPUT.append(re.sub(".sff", "", FILE))

		if Debug:
		        error.write("The base filename of your input file is:"+re.sub(".sff", "", FILE)+"\n")
	elif SFF and RERUN:
		for FILE in glob.glob(MIX+"/*.shhh.fasta"):
                        FILE = re.sub(MIX, "", FILE)
                        FILE = re.sub("./", "", FILE)
			INPUT.append(re.sub(".shhh.fasta", "", FILE))

		if Debug:
		        error.write("The base filename of your input file is:"+re.sub(".shhh.fasta", "", FILE)+"\n")
	else:
		for FILE in glob.glob(MIX+"/*.fasta"):
                        FILE = re.sub(MIX, "", FILE)
                        FILE = re.sub("./", "", FILE)
			INPUT.append(re.sub(".fasta", "", FILE))

		if Debug:
		        error.write("The base filename of your input file is:"+re.sub(".fasta", "", FILE)+"\n")
	if Debug:
       		error.write(str(INPUT)+"\n")
	
else:
	if SFF and not RERUN:
		INPUT = re.sub(".sff", "", INPUT)
		print INPUT

		if Debug:
        		print "The basename of your input file is:\n"
		        error.write(INPUT+"\n")
	elif SFF and RERUN:
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
if SFF and not RERUN:
	for FILE in INPUT:
		if Debug:
			error.write("You are now processing: "+FILE+"\n\n")
		## The chemistry CAN be different for various sequencing runs, MOTHUR provides for this based on the "ORDER" of flows are read. This is sometimes necessary to specify.
	        if MIX:
        		OLIGOS = FILE
        		
        	if ORDER:
	
			os.system(' '.join([
				"for n in",
				FILE+".sff; do mothur \"# sffinfo(sff=$n, flow=T)\"; done"
			]))
			
			if Debug:
				error.write(' '.join([
                                	"for n in",
                                	FILE+".sff; do mothur \"# sffinfo(sff=$n, flow=T)\"; done\n"
                        	]))

			os.system(' '.join([
				"for n in",
				FILE+".flow; do mothur \"# trim.flows(flow=$n,",
				"oligos="+OLIGOS+".oligos,",
				"pdiffs=2, bdiffs=1, minflows=200, maxflows=450,",
	                        "order="+ORDER+",",
				"processors="+PROCESSORS+")\"; done"
			]))

			if Debug:
				error.write(' '.join([
                        	        "for n in",
                                	FILE+".flow; do mothur \"# trim.flows(flow=$n,",
	                                "oligos="+OLIGOS+".oligos,"
        	                        "pdiffs=2, bdiffs=1, minflows=200, maxflows=450,",
                	                "order="+ORDER+",",
                        	        "processors="+PROCESSORS+")\"; done\n"
	                        ]))

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
	                                "processors="+PROCESSORS+")\"; done\n"
        	                ]))

	        else:
        	        os.system(' '.join([
                	        "for n in",
                        	FILE+".sff; do mothur \"# sffinfo(sff=$n, flow=T)\"; done"
	                ]))

			if Debug:
				error.write(' '.join([
	                                "for n in",
        	                        FILE+".sff; do mothur \"# sffinfo(sff=$n, flow=T)\"; done\n"
                	        ]))

        	        os.system(' '.join([
                	        "for n in",
                        	FILE+".flow; do mothur \"# trim.flows(flow=$n,",
	                        "oligos="+OLIGOS+".oligos,",
        	                "pdiffs=2, bdiffs=1, minflows=200, maxflows=450,",
				"processors="+PROCESSORS+")\"; done"
	                ]))

			if Debug:
				error.write(' '.join([
                                	"for n in",
	                                FILE+".flow; do mothur \"# trim.flows(flow=$n,",
        	                        "oligos="+OLIGOS+".oligos,"
                	                "pdiffs=2, bdiffs=1, minflows=200, maxflows=450,",
                        	        "processors="+PROCESSORS+")\"; done\n"
                        	]))

        	        os.system(' '.join([
                	        "for n in",
                        	FILE+".flow.files; do mothur \"# shhh.flows(file=$n,",
				"processors="+PROCESSORS+")\"; done"
        	        ]))

			if Debug:
				error.write(' '.join([
	                                "for n in",
        	                        FILE+".flow.files; do mothur \"# shhh.flows(file=$n,",
                	                "processors="+PROCESSORS+")\"; done\n"
                        	]))

if SFF:
	for FILE in INPUT:
		if MIX:
			OLIGOS = FILE

		os.system(' '.join([
			"for n in",
			FILE+".shhh.fasta; do mothur \"# summary.seqs(fasta=$n,",
			"name="+FILE+".shhh.names,",
			"processors="+PROCESSORS+")\"; done"
		]))
	
		if Debug:
			error.write(' '.join([
				"for n in",
				FILE+".shhh.fasta; do mothur \"# summary.seqs(fasta=$n,",
				"name="+FILE+".shhh.names,",
				"processors="+PROCESSORS+")\"; done\n"
			]))
	
		os.system(' '.join([
			"for n in",
			FILE+".shhh.fasta; do mothur \"# trim.seqs(fasta=$n,",
			"oligos="+OLIGOS+".oligos, name="+FILE+".shhh.names,",
			"maxambig=0, maxhomop=8, bdiffs=1, pdiffs=2, minlength=200,",
			"processors="+PROCESSORS+")\"; done"
		]))
	
		if Debug:
			error.write(' '.join([
				"for n in",
				FILE+".shhh.fasta; do mothur \"# trim.seqs(fasta=$n,",
				"oligos="+OLIGOS+".oligos, name="+FILE+".shhh.names,",
				"maxambig=0, maxhomop=8, bdiffs=1, pdiffs=2, minlength=200,",
				"processors="+PROCESSORS+")\"; done\n"
			]))

	if MIX:
		## Concatenate files
		for FILE in INPUT:
			## Fasta
			os.system(' '.join([ 
				"cat",
				FILE+".shhh.trim.fasta",
				">>",
				NAME+"_Combined_Libraries.shhh.trim.fasta"
			]))

			if Debug:
				error.write(' '.join([
					"cat",
					FILE+".shhh.trim.fasta",
					">>",
					NAME+"_Combined_Libraries.shhh.trim.fasta\n"
				]))
				
			## Names
			os.system(' '.join([ 
				"cat",
				FILE+".shhh.trim.names",
				">>",
				NAME+"_Combined_Libraries.shhh.trim.names"
			]))				

			if Debug:
				error.write(' '.join([
					"cat",
					FILE+".shhh.trim.names",
					">>",
					NAME+"_Combined_Libraries.shhh.trim.names\n"
				]))
				
			## Group	
			os.system(' '.join([ 
				"cat",
				FILE+".shhh.groups",
				">>",
				NAME+"_Combined_Libraries.shhh.groups"
			]))

			if Debug:
				error.write(' '.join([
					"cat",
					FILE+".shhh.groups",
					">>",
					NAME+"_Combined_Libraries.shhh.groups\n"
				]))
				
		INPUT = NAME+"_Combined_Libraries"	

	os.system(' '.join([
		"for n in",
		INPUT+".shhh.trim.fasta; do mothur \"# unique.seqs(fasta=$n,",
		"name="+INPUT+".shhh.trim.names)\"; done"
	]))

	if Debug:
		error.write(' '.join([
	                "for n in",
        	        INPUT+".shhh.trim.fasta; do mothur \"# unique.seqs(fasta=$n,",
                	"name="+INPUT+".shhh.trim.names)\"; done\n"
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
		"name="+INPUT+".shhh.trim.unique.names,",
		"processors="+PROCESSORS+")\"; done"
	]))

	if Debug:
		error.write(' '.join([
	                "for n in",
        	        INPUT+".shhh.trim.unique.align; do mothur \"# summary.seqs(fasta=$n,",
                	"name="+INPUT+".shhh.trim.unique.names,",
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
		"name="+INPUT+".shhh.trim.unique.names,",
		"group="+INPUT+".shhh.groups,",
		"start="+start+",",
		"minlength="+minlength+",",
		"processors="+PROCESSORS+")\"; done"
	]))

	if Debug:
		error.write(' '.join([
	                "for n in",
        	        INPUT+".shhh.trim.unique.align; do mothur \"# screen.seqs(fasta=$n,",
                	"name="+INPUT+".shhh.trim.unique.names,",
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
		"name="+INPUT+".shhh.trim.unique.good.names)\"; done"
	]))

	if Debug:
		error.write(' '.join([
	                "for n in",
        	        INPUT+".shhh.trim.unique.good.filter.fasta; do mothur \"# unique.seqs(fasta=$n,",
                	"name="+INPUT+".shhh.trim.unique.good.names)\"; done\n"
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
	for FILE in INPUT:
		if MIX:
			OLIGOS = FILE

		os.system(' '.join([
			"for n in",
			FILE+".fasta; do mothur \"# summary.seqs(fasta=$n,",
			"name="+FILE+".names,",
			"processors="+PROCESSORS+")\"; done"
		]))
	
		if Debug:
			error.write(' '.join([
				"for n in",
				FILE+".fasta; do mothur \"# summary.seqs(fasta=$n,",
				"name="+FILE+".names,",
				"processors="+PROCESSORS+")\"; done\n"
			]))
	
		os.system(' '.join([
			"for n in",
			FILE+".fasta; do mothur \"# trim.seqs(fasta=$n,",
			"oligos="+OLIGOS+".oligos,",
			"qfile="+FILE+".qual,",
			"maxambig=0, maxhomop=8, bdiffs=1, pdiffs=2, qwindowaverage=35, qwindowsize=50, minlength=200,",
			"processors="+PROCESSORS+")\"; done"
		]))
	
		if Debug:
			error.write(' '.join([
				"for n in",
				FILE+".fasta; do mothur \"# trim.seqs(fasta=$n,",
				"oligos="+OLIGOS+".oligos,",
				"qfile="+FILE+".qual,",
				"maxambig=0, maxhomop=8, bdiffs=1, pdiffs=2, qwindowaverage=35, qwindowsize=50, minlength=200,",
				"processors="+PROCESSORS+")\"; done\n"
			]))
	if MIX:
		## Concatenate files
		for FILE in INPUT:
			## Fasta
			os.system(' '.join([ 
				"cat",
				FILE+".trim.fasta",
				">>",
				NAME+"_Combined_Libraries.trim.fasta"
			]))

			if Debug:
				error.write(' '.join([
					"cat",
					FILE+".trim.fasta",
					">>",
					NAME+"_Combined_Libraries.trim.fasta\n"
				]))
				
			## Names
			os.system(' '.join([ 
				"cat",
				FILE+".trim.names",
				">>",
				NAME+"_Combined_Libraries.trim.names"
			]))				

			if Debug:
				error.write(' '.join([
					"cat",
					FILE+".trim.names",
					">>",
					NAME+"_Combined_Libraries.trim.names\n"
				]))
				
			## Group	
			os.system(' '.join([ 
				"cat",
				FILE+".groups",
				">>",
				NAME+"_Combined_Libraries.groups"
			]))

			if Debug:
				error.write(' '.join([
					"cat",
					FILE+".groups",
					">>",
					NAME+"_Combined_Libraries.groups\n"
				]))
				
		INPUT = NAME+"_Combined_Libraries"
		
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

