#!/usr/bin/python
import sys, os, re, getopt, glob, numpy as np
import timeit
import itertools

now = timeit.default_timer()

Usage = """
Usage:  ./MOTHUR_16S_PIPELINE.py -o BACTERIA_0.05 -n BACTERIA_0.05 -i IIKFCBR01.fasta -f IIKFCBR01.oligos -p 8

REQUIRED ARGUMENTS:
		-o	output directory
	
                -n      the name given to the project

                -i      the input file (accepted: .fasta, .sff, .shhh.fasta. MUST SPECIFY if the latter two)
			(if providing .fasta, one must provide a .qual file)

		-f	the .oligos file

		-p	the number of processors available to use

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

                -b <V>  Enter DEBUG MODE - saves debugging information to "error.log" 
                        (Must provide "V")

UTILITY:
This script will performs the general steps in the mothur pipeline for fungal ITS sequences. This script will require tailoring for different datasets, but provides a simple archetype as well as a quick implementation of whatever you decide. The script does do "summary.seqs" at various points along the way, but the user will have to manually examine the outputs to see how the cleaning went. There is one point where manual intervention is required and the script will wait until one does so."


NOTES: 
1) PROCESSORS ARE SET TO 6. MANUALLY CHANGE IF DESIRED.
2) Minflows is set to 180 (i.e. the minimum length of sequence will be 180bp)
3) Maxflows is set to 450 (i.e. the maximum length of sequence will be 450bp)

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
Debug=''

# Read command line args
myopts, args = getopt.getopt(sys.argv[1:],"n:i:o:S:b:R:f:O:p:")

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
    if o == '-b':
        Debug= a

if Debug:
        print "You are in Debugging Mode"
        error = open("error.log", "w")

if len(OUTPUT)>0:
        if os.path.exists('./' + OUTPUT):
                print "\nOutput Folder Exists - Caution: Files May Be Re-Written"
        else:
                os.mkdir(OUTPUT)
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


if SFF:
	## The chemistry CAN be different for various sequencing runs, MOTHUR provides for this based on the "ORDER" of flows are read. This is sometimes necessary to specify.
        if ORDER:
		os.system(' '.join([
			"for n in",
			INPUT+".sff; do mothur \"# sffinfo(sff=$n, flow=T)\"; done"
		]))

		os.system(' '.join([
			"for n in",
			INPUT+".flow; do mothur \"# trim.flows(flow=$n,",
			"oligos="+OLIGOS+".oligos,"
			"pdiffs=2, bdiffs=1, minflows=200, maxflows=450,",
                        "order="+ORDER+",",
			"processors="+PROCESSORS+")\"; done"
		]))

                os.system(' '.join([
                        "for n in",
                        INPUT+".flow.files; do mothur \"# shhh.flows(file=$n,",
                        "order="+ORDER+",",
			"processors="+PROCESSORS+")\"; done"
                ]))

        else:
                os.system(' '.join([
                        "for n in",
                        INPUT+".sff; do mothur \"# sffinfo(sff=$n, flow=T)\"; done"
                ]))

                os.system(' '.join([
                        "for n in",
                        INPUT+".flow; do mothur \"# trim.flows(flow=$n,",
                        "oligos="+OLIGOS+".oligos,"
                        "pdiffs=2, bdiffs=1, minflows=200, maxflows=450,",
			"processors="+PROCESSORS+")\"; done"
                ]))

                os.system(' '.join([
                        "for n in",
                        INPUT+".flow.files; do mothur \"# shhh.flows(file=$n,",
			"processors="+PROCESSORS+")\"; done"
                ]))

if SFF or RERUN:
	os.system(' '.join([
		"for n in",
		INPUT+".shhh.fasta; do mothur \"# summary.seqs(fasta=$n,",
		"name="+INPUT+".shhh.names,",
		"processors="+PROCESSORS+")\"; done"
	]))

	os.system(' '.join([
		"for n in",
		INPUT+".shhh.fasta; do mothur \"# trim.seqs(fasta=$n,",
		"oligos="+OLIGOS+".oligos, name="+INPUT+".shhh.names,",
		"maxambig=0, maxhomop=8, bdiffs=1, pdiffs=2, minlength=200,",
		"processors="+PROCESSORS+")\"; done"
	]))

	os.system(' '.join([
		"for n in",
		INPUT+".shhh.trim.fasta; do mothur \"# unique.seqs(fasta=$n,",
		"name="+INPUT+".shhh.names)\"; done"
	]))

	os.system(' '.join([
		"for n in",
		INPUT+".shhh.trim.unique.fasta; do mothur \"# align.seqs(fasta=$n,",
		"reference=~/Phylogenetic_Gene_Databases/silva.bacteria/silva.bacteria.fasta, flip=T,",
		"processors="+PROCESSORS+")\"; done"
	]))

	os.system(' '.join([
		"for n in",
		INPUT+".shhh.trim.unique.align; do mothur \"# summary.seqs(fasta=$n,",
		"name="+INPUT+".shhh.trim.names,",
		"processors="+PROCESSORS+")\"; done"
	]))

	#Get user feedback
	print "Based on the summary above, please choose a starting position and minimum length.\n"
	print "ENTER the starting position."
	start = raw_input()
	print "ENTER the minimum length."
	minlength = raw_input()

	os.system(' '.join([
		"for n in",
		INPUT+".shhh.trim.unique.align; do mothur \"# screen.seqs(fasta=$n,",
		"name="+INPUT+".shhh.trim.names,",
		"group="+INPUT+".shhh.groups,",
		"start="+start+",",
		"minlength="+minlength+",",
		"processors="+PROCESSORS+")\"; done"
	]))

	os.system(' '.join([
		"for n in",
		INPUT+".shhh.trim.unique.good.align; do mothur \"# filter.seqs(fasta=$n,",
		"vertical=T,",
		"processors="+PROCESSORS+")\"; done"
	]))

	os.system(' '.join([
		"for n in",
		INPUT+".shhh.trim.unique.good.filter.fasta; do mothur \"# unique.seqs(fasta=$n,",
		"name="+INPUT+".shhh.trim.good.names)\"; done"
	]))

	os.system(' '.join([
		"for n in",
		INPUT+".shhh.trim.unique.good.filter.unique.fasta; do mothur \"# pre.cluster(fasta=$n,",
		"name="+INPUT+".shhh.trim.unique.good.filter.names,",
		"group="+INPUT+".shhh.good.groups, diffs=2,",
		"processors="+PROCESSORS+")\"; done"
	]))

	os.system(' '.join([
		"for n in",
		INPUT+".shhh.trim.unique.good.filter.unique.precluster.fasta; do mothur \"# chimera.uchime(fasta=$n,",
		"name="+INPUT+".shhh.trim.unique.good.filter.unique.precluster.names,",
		"group="+INPUT+".shhh.good.groups,",
		"processors="+PROCESSORS+")\"; done"
	]))

	os.system(' '.join([
		"for n in",
		INPUT+".shhh.trim.unique.good.filter.unique.precluster.uchime.accnos; do mothur \"# remove.seqs(accnos=$n,",
		"fasta="+INPUT+".shhh.trim.unique.good.filter.unique.precluster.fasta,",
		"name="+INPUT+".shhh.trim.unique.good.filter.unique.precluster.names,",
		"group="+INPUT+".shhh.good.groups)\"; done"
	]))

	os.system(' '.join([
		"for n in",
		INPUT+".shhh.trim.unique.good.filter.unique.precluster.pick.fasta; do mothur \"# summary.seqs(fasta=$n,",
		"name="+INPUT+".shhh.trim.unique.good.filter.unique.precluster.pick.names,",
		"processors="+PROCESSORS+")\"; done"
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

	os.system(' '.join([
		"for n in",
		INPUT+".shhh.trim.unique.good.filter.unique.precluster.pick.fasta; do mothur \"# remove.lineage(fasta=$n,",
		"name="+INPUT+".shhh.trim.unique.good.filter.unique.precluster.pick.names,",
		"group="+INPUT+".shhh.good.pick.groups,",
		"taxonomy="+INPUT+".shhh.trim.unique.good.filter.unique.precluster.pick.silva.wang.taxonomy,",
		"taxon=Mitochondria-Cyanobacteria_Chloroplast-Archaea-Eukarya-unknown)\"; done"
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

	##Rename Files to More Readable Names
	os.system(' '.join([
		"cp",
		INPUT+".shhh.trim.unique.good.filter.unique.precluster.pick.pick.names",
		NAME+"_final.names"
	]))
	
	os.system(' '.join([
		"cp",
		INPUT+".shhh.trim.unique.good.filter.unique.precluster.pick.pick.fasta",
		NAME+"_final.fasta"
	]))

	os.system(' '.join([
		"cp",
		INPUT+".shhh.good.pick.pick.groups",
		NAME+"_final.groups"
	]))

	os.system(' '.join([
		"cp",
		INPUT+".shhh.trim.unique.good.filter.unique.precluster.pick.silva.wang.pick.taxonomy",
		NAME+"_final_silva.taxonomy"
	]))

	os.system(' '.join([
		"cp",
		INPUT+".shhh.trim.unique.good.filter.unique.precluster.pick.pick.gg.wang.taxonomy",
		NAME+"_final_gg.taxonomy"
	]))

	#Make List File
	os.system(' '.join([
		"for n in",
		NAME+"_final.fasta; do mothur \"# dist.seqs(fasta=$n,",
		"cutoff=0.3,",
		"processors="+PROCESSORS+")\"; done"
	]))

	os.system(' '.join([
		"for n in",
		NAME+"_final.fasta; do mothur \"# cluster.split(fasta=$n,",
		"name="+NAME+"_final.names,",
		"taxonomy="+NAME+"_final_silva.taxonomy,",
		"splitmethod=classify, taxlevel=3,",
		"processors="+PROCESSORS+")\"; done"
	]))

	#Move all Final Files to Output directory
	os.system(' '.join([
		"mv",
		NAME+"_final.*",
		"./"+OUTPUT+"/"
	]))


	## Provide a version of the GG taxonomy file acceptable for importing into R
	## I do not use Silva for this b/c it does not have consistent delimiting according to taxonomic levels 
        os.system(' '.join([
        	"cp",
                "./"+OUTPUT+"/"+NAME+"_final.gg.taxonomy",
                "./"+OUTPUT+"/"+NAME+"_final.gg.R.taxonomy"
	]))

        os.system(' '.join([
        	"sed -i 's/     /;/g'",
                "./"+OUTPUT+"/"+NAME+"_final.gg.R.taxonomy"
	]))

        os.system(' '.join([
        	"sed -i 's/;//g'",
                "./"+OUTPUT+"/"+NAME+"_final.gg.R.taxonomy"
	]))


else:
	os.system(' '.join([
		"for n in",
		INPUT+".fasta; do mothur \"# summary.seqs(fasta=$n,",
		"name="+INPUT+".names,",
		"processors="+PROCESSORS+")\"; done"
	]))

	os.system(' '.join([
		"for n in",
		INPUT+".fasta; do mothur \"# trim.seqs(fasta=$n,",
		"oligos="+OLIGOS+".oligos,",
		"qfile="+INPUT+".qual,",
		"maxambig=0, maxhomop=8, bdiffs=1, pdiffs=2, qwindowaverage=35, qwindowsize=50, minlength=200,",
		"processors="+PROCESSORS+")\"; done"
	]))

	os.system(' '.join([
		"for n in",
		INPUT+".trim.fasta; do mothur \"# unique.seqs(fasta=$n,",
		"name="+INPUT+".names)\"; done"
	]))

	os.system(' '.join([
		"for n in",
		INPUT+".trim.unique.fasta; do mothur \"# align.seqs(fasta=$n,",
		"reference=~/Phylogenetic_Gene_Databases/silva.bacteria/silva.bacteria.fasta, flip=T,",
		"processors="+PROCESSORS+")\"; done"
	]))

	os.system(' '.join([
		"for n in",
		INPUT+".trim.unique.align; do mothur \"# summary.seqs(fasta=$n,",
		"name="+INPUT+".trim.names,",
		"processors="+PROCESSORS+")\"; done"
	]))

	#Get user feedback
	print "Based on the summary above, please choose a starting position and minimum length.\n"
	print "ENTER the starting position."
	start = raw_input()
	print "ENTER the minimum length."
	minlength = raw_input()

	os.system(' '.join([
		"for n in",
		INPUT+".trim.unique.align; do mothur \"# screen.seqs(fasta=$n,",
		"name="+INPUT+".trim.names,",
		"group="+INPUT+".groups,",
		"start="+start+",",
		"minlength="+minlength+",",
		"processors="+PROCESSORS+")\"; done"
	]))

	os.system(' '.join([
		"for n in",
		INPUT+".trim.unique.good.align; do mothur \"# filter.seqs(fasta=$n,",
		"vertical=T,",
		"processors="+PROCESSORS+")\"; done"
	]))

	os.system(' '.join([
		"for n in",
		INPUT+".trim.unique.good.filter.fasta; do mothur \"# unique.seqs(fasta=$n,",
		"name="+INPUT+".trim.good.names)\"; done"
	]))

	os.system(' '.join([
		"for n in",
		INPUT+".trim.unique.good.filter.unique.fasta; do mothur \"# pre.cluster(fasta=$n,",
		"name="+INPUT+".trim.unique.good.filter.names,",
		"group="+INPUT+".good.groups, diffs=2,",
		"processors="+PROCESSORS+")\"; done"
	]))

	os.system(' '.join([
		"for n in",
		INPUT+".trim.unique.good.filter.unique.precluster.fasta; do mothur \"# chimera.uchime(fasta=$n,",
		"name="+INPUT+".trim.unique.good.filter.unique.precluster.names,",
		"group="+INPUT+".good.groups,",
		"processors="+PROCESSORS+")\"; done"
	]))

	os.system(' '.join([
		"for n in",
		INPUT+".trim.unique.good.filter.unique.precluster.uchime.accnos; do mothur \"# remove.seqs(accnos=$n,",
		"fasta="+INPUT+".trim.unique.good.filter.unique.precluster.fasta,",
		"name="+INPUT+".trim.unique.good.filter.unique.precluster.names,",
		"group="+INPUT+".good.groups)\"; done"
	]))

	os.system(' '.join([
		"for n in",
		INPUT+".trim.unique.good.filter.unique.precluster.pick.fasta; do mothur \"# summary.seqs(fasta=$n,",
		"name="+INPUT+".trim.unique.good.filter.unique.precluster.pick.names,",
		"processors="+PROCESSORS+")\"; done"
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

	os.system(' '.join([
		"for n in",
		INPUT+".trim.unique.good.filter.unique.precluster.pick.fasta; do mothur \"# remove.lineage(fasta=$n,",
		"name="+INPUT+".trim.unique.good.filter.unique.precluster.pick.names,",
		"group="+INPUT+".good.pick.groups,",
		"taxonomy="+INPUT+".trim.unique.good.filter.unique.precluster.pick.silva.wang.taxonomy,",
		"taxon=Mitochondria-Cyanobacteria_Chloroplast-Archaea-Eukarya-unknown)\"; done"
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

	##Rename Files to More Readable Names
	os.system(' '.join([
		"cp",
		INPUT+".trim.unique.good.filter.unique.precluster.pick.pick.names",
		NAME+"_final.names"
	]))
	
	os.system(' '.join([
		"cp",
		INPUT+".trim.unique.good.filter.unique.precluster.pick.pick.fasta",
		NAME+"_final.fasta"
	]))

	os.system(' '.join([
		"cp",
		INPUT+".good.pick.pick.groups",
		NAME+"_final.groups"
	]))

	os.system(' '.join([
		"cp",
		INPUT+".trim.unique.good.filter.unique.precluster.pick.silva.wang.pick.taxonomy",
		NAME+"_final_silva.taxonomy"
	]))

	os.system(' '.join([
		"cp",
		INPUT+".trim.unique.good.filter.unique.precluster.pick.pick.gg.wang.taxonomy",
		NAME+"_final_gg.taxonomy"
	]))

	#Make List File
	os.system(' '.join([
		"for n in",
		NAME+"_final.fasta; do mothur \"# dist.seqs(fasta=$n,",
		"cutoff=0.3,",
		"processors="+PROCESSORS+")\"; done"
	]))

	os.system(' '.join([
		"for n in",
		NAME+"_final.fasta; do mothur \"# cluster.split(fasta=$n,",
		"name="+NAME+"_final.names,",
		"taxonomy="+NAME+"_final_silva.taxonomy,",
		"splitmethod=classify, taxlevel=3,",
		"processors="+PROCESSORS+")\"; done"
	]))

	#Move all Final Files to Output directory
	os.system(' '.join([
		"mv",
		NAME+"_final.*",
		"./"+OUTPUT+"/"
	]))

	## Provide a version of the GG taxonomy file acceptable for importing into R
	## I do not use Silva for this b/c it does not have consistent delimiting according to taxonomic levels 
        os.system(' '.join([
        	"cp",
                "./"+OUTPUT+"/"+NAME+"_final.gg.taxonomy",
                "./"+OUTPUT+"/"+NAME+"_final.gg.R.taxonomy"
	]))

        os.system(' '.join([
        	"sed -i 's/     /;/g'",
                "./"+OUTPUT+"/"+NAME+"_final.gg.R.taxonomy"
	]))

        os.system(' '.join([
        	"sed -i 's/;//g'",
                "./"+OUTPUT+"/"+NAME+"_final.gg.R.taxonomy"
	]))

end = timeit.default_timer()

print now - end

