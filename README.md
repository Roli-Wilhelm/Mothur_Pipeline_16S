Mothur_Pipeline_16S
===================

Standard Pipeline for Cleaning, Classifying and constructing OTUs for 16S rRNA gene phylogenetic libraries

Provides ALL necessary files to use the R script "Mothur_to_PhyloSeq_to_Differential_Abundance_Tools"
(in one case, a specially formatted taxonomy file)

This script is not meant to be a completely hands-off method for processing 16S rRNA data. It IS fully automated, but only in so far as it executes all the pre-set instructions. If your library requires different processing parameters, you will have to manually alter the script. It will do a decent job on your data as is, based on the fact multiple colleagues in my lab have used an nearly identical set of commands on large datasets (500+ libraries). 

Ths script is meant to offer a simple archetype and is very similar to the Schloss 454 SOP <http://www.mothur.org/wiki/454_SOP>. Further, there is a GUI version of Mothur available which will do a similarly automated process. I still find this faster, but take a look at it HERE: <http://www.mothur.org/wiki/Download_mothur#Graphical_Interfaces_for_Mothur>

Dependencies: 

- mothur (built with v.1.32.1)

- Green Genes "Mothur Formatted" Database (script expects it to be located in ~/Phylogenetic_Databases/GreenGenes_MothurFormatted/*)

- Silva Bacteria Database (script expects it to be located in ~/Phylogenetic_Databases/silva.bacteria/*)

This script will require mothur to be in your PATH. 
