# Run the docker 
sudo docker run -v $PWD:$PWD -w $PWD --rm -it dgarrimar/epigenomics_course

# clone the epigenomics repository 
git clone https://github.com/bborsari/epigenomics_uvic
# move to different directories inside the docker 
cd epigenomics_uvic
cd ATAC-seq 
ls


#Go to ENCODE (ENTEX) to look for our donant information and retrive the metadata
# (all the experiments regarding DNA accessibility for sigmoid colon and stomach) 
# Once downloaded we need to move the file to our docker 
cp /home/me/Downloads/files.txt /home/me/epigenomics_uvic/ATAC-seq/files.txt




#Move to folder ATAC-seq, and create folders to store bigBed data files and peaks #analyses files. Make sure the files are organized in a consistent way as done for ChIP-#seq.

# Using the first line (link directs us to our experiment) we download the files from the experiment 
./bin/download.metadata.sh "https://www.encodeproject.org/metadata/?type=Experiment&replicates.library.biosample.donor.uuid=d370683e-81e7-473f-8475-7716d027849b&status=released&assay_slims=DNA+accessibility&biosample_ontology.term_name=sigmoid+colon&biosample_ontology.term_name=stomach&assembly=GRCh38&assay_title=ATAC-seq" 


# Which fileds are we interested in from out metadata? 
head -1 metadata.tsv | awk 'BEGIN{FS=OFS="\t"}{for (i=1;i<=NF;i++){print $i, i}}'


# Organize a few the things
mkdir analyses
cd analyses
mkdir data
mkdir data/bigBed.files data/bigWig.files
ls # check that everything has been created

# Now that things are organized  we can download the BigBed and bigWig files
# we want the peaks from ATAC-seq experiments, so we need to extract info necessary:

awk '$8 == "ATAC-seq"' metadata.tsv |grep -F "bigBed_narrowPeak" |grep -F "pseudoreplicated_peaks" |grep -F "GRCh38" |awk 'BEGIN{FS=OFS="\t"}{print $1, $10, $22}' |sort -k2,2 -k1,1r | sort -k2,2 -u > bigBedATAC_peaks.txt

# move them 
mv bigBedATAC_peaks.txt /home/me/epigenomics_uvic/ATAC-seq/analyses/data/bigBed.files


# download them 
cut -f1 analyses/data/bigBed.files/bigBedATAC_peaks.txt | while read filename; do   wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"; done



# check 
for file_type in bigBed; do   ../bin/selectRows.sh <(cut -f1 analyses/"$file_type"ATAC*.txt) metadata.tsv | cut -f1,45 > data/"$file_type".files/md5sum.txt;    cat data/"$file_type".files/md5sum.txt |  while read filename original_md5sum; do      md5sum data/"$file_type".files/"$filename"."$file_type" |    awk -v filename="$filename" -v original_md5sum="$original_md5sum" 'BEGIN{FS=" "; OFS="\t"}{print filename, original_md5sum, $1}' ;   done > tmp ;   mv tmp data/"$file_type".files/md5sum.txt;    awk '$2!=$3' data/"$file_type".files/md5sum.txt;  done

# no output means that they coincide! So they are good!



# Convert bigBed files of ATAC-seq peaks to BED files with the bigBedToBed command:
cut -f1 analyses/bigBedATAC_peaks.txt |while read filename; do   bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$filename".bed; done


# Create an annotation folder for further steps:
mkdir annotation



# Download from a link the list of promoters ([-2 kb, +2 Kb] from TSS) of protein-coding genes. Store this file inside the annotation folder. (called gencode.v24.protein.coding.non.redundant.TSS.bed)
mv /home/me/Downloads/gencode.v24.protein.coding.non.redundant.TSS.bed /home/me/epigenomics_uvic/ATAC-seq/

mv gencode.v24.protein.coding.non.redundant.TSS.bed /home/me/epigenomics_uvic/ATAC-seq/annotation




# Retrieve genes of ATAC-seq peaks at the promoter region in each tissue.

cut -f-2 analyses/bigBedATAC_peaks.txt |\
while read filename tissue; do 
  bedtools intersect -a annotation/gencode.v24.protein.coding.non.redundant.TSS.bed -b data/bed.files/"$filename".bed -wb|\
  sort -u > analyses/genes.with.peaks."$tissue".ATAC.txt;
done



cd analyses 
ls 
# checking 
# files are there!


# How many genes for SIGMOID COLON:
wc -l genes.with.peaks.sigmoid_colon.ATAC.txt
		230540 

# How many genes for STOMACH:
wc -l genes.with.peaks.stomach.ATAC.txt
		206927 



# Downloaded the GENCODE reference files:
# but We need the same version used by ENCODE (v24) 
# named gencode.v24.primary_assembly.annotation.gtf.gz
mv /home/me/Downloads/gencode.v24.primary_assembly.annotation.gtf.gz /home/me/epigenomics_uvic/ATAC-seq/annotation



# Prepare a BED file with gene body coordinates of protein-coding genes
# Uncompress the gtf.gz file
gunzip annotation/gencode.v24.primary_assembly.annotation.gtf.gz



# Retrieve the protein coding genes for getting the gene coordinates 
awk '$3=="gene"' annotation/gencode.v24.primary_assembly.annotation.gtf |\
grep -F "protein_coding" |\
cut -d ";" -f1 |\
awk 'BEGIN{OFS="\t"}{print $1, $4, $5, $10, 0, $7, $10}' |\
sed 's/\"//g' |\
awk 'BEGIN{FS=OFS="\t"}$1!="chrM"{$2=($2-1); print $0}' > annotation/gencode.v24.protein.coding.gene.body.bed



# Check how it looks like
head annotation/gencode.v24.protein.coding.gene.body.bed



# Retrieve peaks that fall outside gene coordinates 
cut -f-2 analyses/bigBedATAC_peaks.txt |\
while read filename tissue; do 
  bedtools intersect -a annotation/gencode.v24.protein.coding.gene.body.bed -b data/bed.files/"$filename".bed -wb|\
  sort -u > analyses/genes.with.peaks."$tissue".ATAC_outside.bed
done


# How many for SIGMOID COLON:
wc -l genes.with.peaks.sigmoid_colon.ATAC_outside.bed
			77609 
			


# How many for STOMACH: 
wc -l genes.with.peaks.stomach.ATAC_outside.bed
			70606 
			




########################## SECOND PART (POINT 5) #############################
#Task 1: Create a folder regulatory_elements inside epigenomics_uvic. This will be the folder where you store all your subsequent results.
mkdir regulatory_elements



# Task 2: 

# Download the metadata for the epigenomics marks 
../bin/download.metadata.sh "https://www.encodeproject.org/metadata/?type=Experiment&replicates.library.biosample.donor.uuid=d370683e-81e7-473f-8475-7716d027849b&status=released&assembly=GRCh38&biosample_ontology.term_name=sigmoid+colon&biosample_ontology.term_name=stomach&assay_slims=DNA+binding"



# look for the marks of interest (H3K4me1 and H3K27ac)
grep 'H3K4me1\|H3K27ac' metadata.tsv |\
grep -F "bigBed_narrowPeak" |\
grep -F "pseudoreplicated_peaks" |\
grep -F "GRCh38" |\
awk 'BEGIN{FS=OFS="\t"}{print $1, $10, $22}' |\
sort -k2,2 -k1,1r |\
sort -k2,2 -u > analyses/bigBed.peaks.ids.txt




# download the files 
cut -f1 analyses/bigBed.peaks.ids.txt|\
while read filename; do
  wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
done


# convert them into a readable file (bed files)
mkdir data/bed.files
cut -f1 analyses/bigBed.peaks.ids.txt|\
while read filename; do
  bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$filename".bed
done



# We want to find those peaks (from H3K4me1 or H3K27ac) that coincide with the ones that layed outside the promoter region (obtained in the exercise before)
# move the files extracted from the previous exercise to our local environment 

mv /home/me/epigenomics_uvic/ATAC-seq/analyses/genes.with.peaks.sigmoid_colon.ATAC_outside.bed
mv /home/me/epigenomics_uvic/ATAC-seq/analyses/genes.with.peaks.stomach.ATAC_outside.bed


# Make both files mentiones before intersect:
cut -f-2 analyses/bigBed.peaks.ids.txt |\
while read filename tissue; do 
  bedtools intersect -a data/bed.files/"$filename".bed -b analyses/genes.with.peaks.stomach.ATAC_outside.bed -wa|\
  sort -u > analyses/peaks.analysis/genes.OUTSIDE.peaks."$tissue".me1ac.txt
done


# How many are they? 
wc -l genes.OUTSIDE.peaks.sigmoid_colon.me1ac.txt 
			17493


wc -l genes.OUTSIDE.peaks.stomach.me1ac.txt 
			21945 


### Task 3. Focus on chromosome 1 and keep the name of the original ATAC-seq peak and the start (5') coordinate of the region.

# One needs to extract the first column (chromosome name) the second column (start coordinate region) and the fourth column (peak name) from the files generated before.

# There fore I first, substract these columns of both files 
awk '{ print $1, $2, $4 }' genes.OUTSIDE.peaks.sigmoid_colon.me1ac.txt > sigmoidcolon_selected_columns_me1ac_outside.txt
awk '{ print $1, $2, $4 }' genes.OUTSIDE.peaks.stomach.me1ac.txt > stomach_selected_columns_me1ac_outside.txt

# Now I will substract the elements of each file that are in chromosome 1: 


awk '$1=="chr1"' sigmoidcolon_selected_columns_me1ac_outside.txt > sigmoidcolon_regulatory.elements.starts.tsv
awk '$1=="chr1"' stomach_selected_columns_me1ac_outside.txt > stomach_regulatory.elements.starts.tsv




# Task 4: Focus on protein-coding genes located on chromosome 1. From the BED file of gene body coordinates that you generated here, prepare a tab-separated file called gene.starts.tsv which will store the name of the gene in the first column, and the start coordinate of the gene on the second column (REMEMBER: for genes located on the minus strand, the start coordinate will be at the 3'). Use the command below as a starting point:

# First one is intersecting the annotation for protein coding genes with the bed files for each tissue. 
cut -f-2 analyses/bigBed.peaks.ids.txt |\
while read filename tissue; do 
  bedtools intersect -a annotation/gencode.v24.protein.coding.gene.body.bed -b data/bed.files/"$filename".bed -wb|\
  sort -u > analyses/input.chr1.with.peaks."$tissue".ATAC_outside.bed
done

# Then one is subsetting only the ones that are in chromosome one
awk '$1=="chr1"' input.chr1.with.peaks.sigmoid_colon.ATAC_outside.bed > sigmoidcolon_peaks.chr1.bed
awk '$1=="chr1"' input.chr1.with.peaks.stomach.ATAC_outside.bed > stomach_peaks.chr1.bed

# Storing the name of the gene (in the first column) and the start coordinate (in the second), but considering that for genes located in the minus (-) strand, the coordinate will not be in the third column but in the fourth. , because the start will be at 3' of the strand. 
awk 'BEGIN{FS=OFS="\t"}{if ($6=="+"){start=$2} else {start=$3}; print $4, start}' sigmoidcolon_peaks.chr1.bed > gene.starts.sigmoid_colon.tsv
awk 'BEGIN{FS=OFS="\t"}{if ($6=="+"){start=$2} else {start=$3}; print $4, start}' stomach_peaks.chr1.bed > gene.starts.stomach.tsv







# Task 5: Download or copy this python script inside the epigenomics_uvic/bin folder. Have a look at the help page of this script to understand how it works:

# One opened the link and save it (save page as) to the Desktop of my VM 
# Then move it to the bin folder inside epigenomics_uvic directory
mv get.distance.py /home/me/epigenomics_uvic/bin/
cd /home/me/epigenomics_uvic/bin/

# Used the command help page in order to understand how this script works: 
python ../bin/get.distance.py -h

# How the python script looks like after the changes applied on it:


			#!/usr/bin/env python


			#************
			# LIBRARIES *
			#************

			import sys
			from optparse import OptionParser


			#*****************
			# OPTION PARSING *
			#*****************

			parser = OptionParser()
			parser.add_option("-i", "--input", dest="input")
			parser.add_option("-s", "--start", dest="start")
			options, args = parser.parse_args()

			open_input = open(options.input)
			enhancer_start = int(options.start)


			#********
			# BEGIN *
			#********

			x=1000000 # set maximum distance to 1 Mb
			selectedGene="" # initialize the gene as empty
			selectedGeneStart=0 # initialize the start coordinate of the gene as empty

			for line in open_input.readlines(): # for each line in the input file
				gene, y = line.strip().split('\t') # split the line into two columns based on a tab 
				
				
				# define a variable called position that correspond to the integer of the start of the gene
				position= int(selectedGeneStart)
				
				# compute the absolute value of the difference between position and enhancer_start
				absolut_pos = abs(position - enhancer_start)
				
				# if this absolute value is lower than x
				if absolut_pos < x:
					# this value will now be your current x
					absolut_pos = x
					# save gene as selectedGene
					gene = selectedGene
					# save position as selectedGeneStart
					position = selectedGeneStart

			print "\t".join([selectedGene, str(selectedGeneStart), str(x)])
 



# Make sure that the Script is working

python bin/get.distance.py --input regulatory_elements/analyses/gene.starts.sigmoid_colon.tsv --start 980000
	ENSG00000156875.13	0	980000


python bin/get.distance.py --input regulatory_elements/analyses/gene.starts.stomach.tsv --start 980000
	ENSG00000156875.13	0	980000


# Task 6. For each regulatory element contained in the file regulatory.elements.starts.tsv, retrieve the closest gene and the distance to the closest gene using the python script you created above. Use the command below as a starting point:

# Python script modified: 

		#!/usr/bin/env python

		#************

		# LIBRARIES *

		#************

		import sys
		from optparse import OptionParser


		#*****************

		# OPTION PARSING *

		#*****************

		parser = OptionParser()
		parser.add_option("-g", "--gene", dest="gene")
		parser.add_option("-r", "--reg", dest="reg")
		options, args = parser.parse_args()

		open_gene = open(options.gene)
		open_reg = open(options.reg)

		#********

		# BEGIN *

		#********

		genes = {} # initialize an empty dictionary 

		# for GENES file
		for line in open_gene.readlines(): # for each line in the input file
			genes_l = line.split() # split the line into two columns based on a tab 
			genes[genes_l[0]] = int(genes_l[1]) # save results in a dictionary 

		regs = {} # empty dictionary 


		for line1 in open_reg.readlines(): # for each line in the input file
			peaks_l = line1.split() # split the line into two columns based on a tab 
			regs[peaks_l[2]] = int(peaks_l[1]) # save results in a dictionary 



		# iterate both dicts to get min distance
		results = {}
		for peaks, pos in regs.items():
			dist_min = 100000
			gen_value = ""

			for gen, pos1 in genes.items():
				dist = abs(pos - pos1)

				if dist <= dist_min:
					dist_min = dist
					gen_value = gen
					results[gen_value] = dist_min

		# print
		for g, v in results.items():
			print (g + '\t' + str(v))



# For calling the script and creating the resultant file in order to get the closest gene and its distance for each regulatory element, for each type, one will run the following command: 


# For sigmoid colon
cat regulatory_elements/analyses/peaks.analysis/sigmoidcolon_regulatory.elements.starts.tsv | while read element start; do     python ./bin/get.distance.py --gene regulatory_elements/analyses/gene.starts.sigmoid_colon.tsv --reg regulatory_elements/analyses/peaks.analysis/sigmoidcolon_regulatory.elements.starts.tsv; done > regulatoryElements_sigmoidcolon.genes.distances.tsv




# For stomach colon
cat regulatory_elements/analyses/peaks.analysis/stomach_regulatory.elements.starts.tsv | while read element start; do     python ./bin/get.distance.py --gene regulatory_elements/analyses/gene.starts.stomach.tsv --reg regulatory_elements/analyses/peaks.analysis/stomach_regulatory.elements.starts.tsv; done > regulatoryElements_stomach.genes.distances.tsv




# Task 7. Use R to compute the mean and the median of the distances stored in regulatoryElements.genes.distances.tsv.

# First we open R inside shell's environment
$R



# FOR SIGMOID COLON'S DATA 
# One is loading the data into this 'R' environment called before, and saving it in a object called 'data'.
data <- read.table("regulatoryElements_sigmoidcolon.genes.distances.tsv", sep="")

# As the data is composed of a column that corresponds to the gene and a column that corresponds to the distance, one is computing: (1) the mean and (2) the median running these codes: 
mean <- mean(data[,2])
mean 
	19368.98
	
med <- median(data[,2])
med

	2939
