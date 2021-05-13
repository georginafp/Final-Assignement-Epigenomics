# Final-Assignement-Epigenomics
Project based on bash, python and R that contains all the tasks entrusted on the Epigenomics subject. All the instructions can be found on the "Wiki" access, as well as in the README.txt file. 


##################################################################################
################## 4. EN TEx ATAC seq data: Downstream analysis ################## 
##################################################################################



## 1. Move to folder ATAC-seq, and create folders to store bigBed data files and peaks analyses files. Make sure the files are organized in a consistent way as done for ChIP-seq.
# Run the docker
sudo docker run -v $PWD:$PWD -w $PWD --rm -it dgarrimar/epigenomics_course

# Clone the epigenomics repository and move to different directories inside the docker
git clone https://github.com/bborsari/epigenomics_uvic
cd epigenomics_uvic
cd ATAC-seq ls

# Go to ENCODE (ENTEX) to look for our donant (ENCDO451RUA) information and retrive the metadata. 
# DNA accessibility for sigmoid colon and stomach 
# GRCh38 genome version
cp /home/me/Downloads/files.txt /home/me/epigenomics_uvic/ATAC-seq/files.txt

# Move to folder ATAC-seq 
# Create folders to store bigBed data files and peaks analyses files as done for ChIP-seq. 
# Download the files from the experiment. Finally, we subset the fields in which we are interested in (from our metadata file)
./bin/download.metadata.sh "https://www.encodeproject.org/metadata/?type=Experiment&replicates.library.biosample.donor.uuid=d370683e-81e7-473f-8475-7716d027849b&status=released&assay_slims=DNA+accessibility&biosample_ontology.term_name=sigmoid+colon&biosample_ontology.term_name=stomach&assembly=GRCh38&assay_title=ATAC-seq" 
head -1 metadata.tsv | awk 'BEGIN{FS=OFS="\t"}{for (i=1;i<=NF;i++){print $i, i}}'

# Organize a few the things
mkdir analyses
cd analyses
mkdir data
mkdir data/bigBed.files data/bigWig.files
ls



## 2. Retrieve from a newly generated metadata file ATAC-seq peaks (bigBed narrow, pseudoreplicated peaks, assembly GRCh38) for stomach and sigmoid_colon for the same donor used in the previous sections. 
# Hint: have a look at what we did here. Make sure your md5sum values coincide with the ones provided by ENCODE.
# Extract the necessary information
grep -F "bigBed_narrow" metadata.tsv | grep -F "pseudoreplicated_peaks" |\ 
grep -F "GRCh38" |awk 'BEGIN{FS=OFS="\t"}{print $1, $10, $22}' |\ 
sort -k2,2 -k1,1r | sort -k2,2 -u > bigBedATAC_peaks.txt

# Move and download them
mv bigBedATAC_peaks.txt /home/me/epigenomics_uvic/ATAC-seq/analyses/data/bigBed.files
cut -f1 analyses/data/bigBed.files/bigBedATAC_peaks.txt |\
 while read filename; do   wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"; |\
done

# Check the integrity of the downloaded files in order to verify their MD5 hash, whose is a kind of digital fingerprint of the file. 
# When no output is generated after the next command, means that they coincide so everything is fine.
../bin/selectRows.sh <(cut -f1 analyses/"$file_type"ATAC*.txt) metadata.tsv |\
cut -f1,45 > data/"$file_type".files/md5sum.txt cat data/"$file_type".files/md5sum.txt |\
while read filename original_md5sum; do md5sum data/"$file_type".files/"$filename"."$file_type" |\
awk -v filename="$filename" -v original_md5sum="$original_md5sum" 'BEGIN{FS=" "; OFS="\t"}{print filename, original_md5sum, $1}' |\
done > tmp mv tmp data/"$file_type".files/md5sum.txt awk '$2!=$3' data/"$file_type".files/md5sum.txt done

# Convert the files to BED files format with the bigBedToBed command
cut -f1 analyses/bigBedATAC_peaks.txt |\
while read filename; do   bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$filename".bed |\
done

# Create an annotation folder for further steps and download from a link the list of promoters ([-2 kb, +2 Kb] from TSS) of protein-coding genes. 
# Store this file inside the annotation folder (called gencode.v24.protein.coding.non.redundant.TSS.bed)
mkdir annotation
mv /home/me/Downloads/gencode.v24.protein.coding.non.redundant.TSS.bed /home/me/epigenomics_uvic/ATAC-seq/
mv gencode.v24.protein.coding.non.redundant.TSS.bed /home/me/epigenomics_uvic/ATAC-seq/annotation



## 3. For each tissue, run an intersection analysis using BEDTools report. Notice that the intersect function allows us to intersect two databases for retrieving those peaks from the first one that overlap the second one.
# Retrieve the number of peaks that intersect promoter regions
cat analyses/bigBedATAC_peaks.txt |\
while read filename tissue; do bedtools intersect -a data/bed.files/"$filename".bed -b annotation/gencode.v24.protein.coding.non.redundant.TSS.bed -u |\
cut -f7 | sort -u > analyses/genes.with.peaks."$tissue".ATAC.txt done

# Check that the files are there
cd analyses 
ls 

# Number of genes for SIGMOID COLON
wc -l genes.with.peaks.sigmoid_colon.ATAC.txt
		# 38071

# Number of genes for STOMACH
wc -l genes.with.peaks.stomach.ATAC.txt
		# 33169

# The number of peaks that fall outside gene coordinates (whole gene body, not just the promoter regions). 
# Hint: have a look at what we did here and here.
#  Downloaded the GENCODE reference files, but notice that one should point out that there is the need to be in the same version used by ENCODE (v24 in our case)
# named: named gencode.v24.primary_assembly.annotation.gtf.gz
mv /home/me/Downloads/gencode.v24.primary_assembly.annotation.gtf.gz /home/me/epigenomics_uvic/ATAC-seq/annotation

# Prepare a BED file with gene body coordinates of protein-coding genes, and uncompress the gtf.gz file
gunzip annotation/gencode.v24.primary_assembly.annotation.gtf.gz

# Retrieve the protein coding genes for getting the gene coordinates
awk '$3=="gene"' annotation/gencode.v24.primary_assembly.annotation.gtf |\
grep -F "protein_coding" |\
cut -d ";" -f1 |\
awk 'BEGIN{OFS="\t"}{print $1, $4, $5, $10, 0, $7, $10}' |\
sed 's/\"//g' |\
awk 'BEGIN{FS=OFS="\t"}$1!="chrM"{$2=($2-1); print $0}' > annotation/gencode.v24.protein.coding.gene.body.bed

head annotation/gencode.v24.protein.coding.gene.body.bed

# Retrieve peaks that fall outside gene coordinates
cat analyses/bigBedATAC_peaks.txt |\
while read filename tissue; do bedtools intersect -a data/bed.files/"$filename".bed -b annotation/gencode.v24.protein.coding.gene.body.bed -v |\
sort > analyses/genes.with.peaks."$tissue".ATAC_outside.bed

# Number of peaks for SIGMOID COLON
wc -l genes.with.peaks.sigmoid_colon.ATAC_outside.txt
		# 34537

# Number of peaks for for STOMACH
wc -l genes.with.peaks.stomach.ATAC_outside.txt
		# 37035



###################################################################
################## 5. Distal Regulatory Activity ################## 
###################################################################

## TASK 1. Create a folder regulatory_elements inside epigenomics_uvic. This will be the folder where you store all your subsequent results. 
## Plus, other directories that will be used for the analysis are also created below.
mkdir regulatory_elements
mkdir annotation
mkdir data
mkdir data/bed.files
mkdir data/bigBed.files
mkdir analyses
mkdir analyses/peaks.analyses



## TASK 2. Distal regulatory regions are usually found to be flanked by both H3K27ac and H3K4me1. 
## From your starting catalogue of open regions in each tissue, select those that overlap peaks of H3K27ac AND H3K4me1 in the corresponding tissue. 
## You will get a list of candidate distal regulatory elements for each tissue. How many are they?

# Download peaks and obtaining its bigBed files for H3K4me1 histone mark
grep -F H3K4me1 ../ChIP-seq/metadata.tsv | grep -F "bigBed_narrowPeak" |\
grep -F "pseudoreplicated_peaks" | grep -F "GRCh38" |\
awk 'BEGIN{FS=OFS="\t"}{print $1,$10,$22}' | sort -k2,2 -k1,1r |\
sort -k2,2 -u > analyses/bigBed.peaks.H3K4me1.txt

cut -f1 analyses/bigBed.peaks.H3K4me1.txt |\
while read filename; do wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed" done


# Download peaks and obtaining its bigBed files for H3K27ac histone mark
grep -F H3K27ac ../ChIP-seq/metadata.tsv | grep -F "bigBed_narrowPeak" |\
grep -F "pseudoreplicated_peaks" | grep -F "GRCh38" |\
awk 'BEGIN{FS=OFS="\t"}{print $1,$10,$22}' | sort -k2,2 -k1,1r |\
sort -k2,2 -u > analyses/bigBed.peaks.H3K27ac.txt

cut -f1 analyses/bigBed.peaks.H3K27ac.txt |\
while read filename; do wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed" done

# Check the integrity of the downloaded files in order to verify their MD5 hash, whose is a kind of digital fingerprint of the file. 
# When no output is generated after the next command, means that they coincide so everything is fine.
			# H3K4me1
for file_type in bigBed; do ../bin/selectRows.sh <(cut -f1 analyses/"$file_type".*.peaks.H3K4me1.txt) metadata.tsv |\
cut -f1,45 > data/"$file_type".files/md5sum.H3K4me1.txt |\ 
cat data/"$file_type".files/md5sum.H3K4me1.txt |\
while read filename original_md5sum; do md5sum data/"$file_type".files/"$filename"."$file_type" |\
awk -v filename="$filename" -v original_md5sum="$original_md5sum" 'BEGIN{FS=" "; OFS="\t"}{print filename, original_md5sum, $1}' done > tmp mv tmp data/"$file_type".files/md5sum.H3K4me1.txt awk '$2!=$3' data/"$file_type".files/md5sum.H3K4me1.txt done

head data/bigBed.files/md5sum.H3K4me1.txt

			# H3K27ac
for file_type in bigBed; do ../bin/selectRows.sh <(cut -f1 analyses/"$file_type".*.peaks.H3K27ac.txt) metadata.tsv |\
cut -f1,45 > data/"$file_type".files/md5sum.H3K27ac.txt |\ 
cat data/"$file_type".files/md5sum.H3K27ac.txt |\
while read filename original_md5sum; do md5sum data/"$file_type".files/"$filename"."$file_type" |\
awk -v filename="$filename" -v original_md5sum="$original_md5sum" 'BEGIN{FS=" "; OFS="\t"}{print filename, original_md5sum, $1}' done > tmp mv tmp data/"$file_type".files/md5sum.H3K27ac.txt awk '$2!=$3' data/"$file_type".files/md5sum.H3K27ac.txt done

head data/bigBed.files/md5sum.H3K27ac.txt

#  Convert bigBed files to bed files to be readable and treat and run an intersection

			# H3K4me1
cut -f1 analyses/bigBed.peaks.H3K4me1.txt |\
while read filename; do bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$filename".H3K4me1.bed done

			# H3K27ac
cut -f1 analyses/bigBed.peaks.H3K27ac.txt |\
while read filename; do bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$filename".H3K27ac.bed done


# Retrieve those peaks (from H3K4me1 or H3K27ac) that coincide with the ones that layed outside the promoter region.

			# sigmoid colon: 23130
bedtools intersect -a ../ATAC-seq/analyses/genes.with.peaks.sigmoid_colon.ATAC_outside.txt -b data/bed.files/ENCFF724ZOF.bed data/bed.files/ENCFF872UHN.bed -u |\
sort -u > analyses/peaks.analyses/genes.peaks.overlap.sigmoid_colon.txt

			# stomach: 18028
bedtools intersect -a ../ATAC-seq/analyses/genes.with.peaks.stomach.ATAC_outside.txt -b data/bed.files/ENCFF844XRN.bed data/bed.files/ENCFF977LBD.bed -u |\
sort -u > analyses/peaks.analyses/genes.peaks.overlap.stomach.txt



## TASK 3. Focus on regulatory elements that are located on chromosome 1, and generate a file regulatory.elements.starts.tsv that contains the name of the regulatory region 
## (i.e. the name of the original ATAC-seq peak) and the start (5') coordinate of the region.

# extract the first, second and the fourth column
			# sigmoid colon: 10734
grep -F “chr1” analyses/peaks.analyses/genes.peaks.overlap.sigmoid_colon.txt |\
awk 'BEGIN{FS=OFS="\t"} {print $4,$2}’ > analyses/peaks.analyses/regulatory.elements.starts.sigmoid_colon.tsv

wc -l regulatory.elements.starts.sigmoid_colon.tsv

			# stomach: 8602
grep -F “chr1” analyses/peaks.analyses/genes.peaks.overlap.stomach.txt |\
awk 'BEGIN{FS=OFS="\t"} {print $4,$2}’ > analyses/ peaks.analyses/regulatory.elements.starts.stomach.tsv

wc -l regulatory.elements.starts.stomach.tsv




TASK 4. Focus on protein-coding genes located on chromosome 1. From the BED file of gene body coordinates that you generated here 
## prepare a tab-separated file called gene.starts.tsv which will store the name of the gene in the first column, and the start coordinate of the gene on the second column 
## (REMEMBER: for genes located on the minus strand, the start coordinate will be at the 3'). Use the command below as a starting point

awk 'BEGIN{FS=OFS="\t"} $1==”chr1” {if ($6=="+"){start=$2} else {start=$3}; print $4, start}' ../ChIP-seq/annotation/gencode.v24.protein.coding.gene.body.bed > analyses/gene.starts.tsv
wc -l analyses/gene.starts.tsv
			# 2047



## TASK 5. Download or copy this python script inside the epigenomics_uvic/bin folder. Have a look at the help page of this script to understand how it works
# One opened the link and save it and then move it to the bin folder inside epigenomics_uvic directory
mv get.distance.py /home/me/epigenomics_uvic/bin/
cd /home/me/epigenomics_uvic/bin/

# help command
python ../bin/get.distance.py -h

# modified script
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
			position= int(y)
 				
			# compute the absolute value of the difference between position and enhancer_start
			absolut_pos = abs(position - enhancer_start)
 				
			# if this absolute value is lower than x
			if absolut_pos < x:
				# this value will now be your current x
				x = absolut_pos 
				# save gene as selectedGene
				selectedGene = gene
				# save position as selectedGeneStart
				selectedGeneStart = position

		print "\t".join([selectedGene, str(selectedGeneStart), str(x)])


# Make sure it works
python ../bin/get.distance.py --input analyses/gene.starts.tsv --start 980000

		# ENSG00000187642.9 982093 2093



## TASK 6. For each regulatory element contained in the file regulatory.elements.starts.tsv, retrieve the closest gene and the distance to the closest gene using the python script you created above. 
## Use the command below as a starting point
for tissue in stomach sigmoid_colon; do cat analyses/peaks.analyses/regulatory.elements.starts.$tissue.tsv |\
while read element start; do python ../bin/get.distance.py --input analyses/gene.starts.tsv --start $start done > analyses/peaks.analyses/regulatoryElements.genes.distances.$tissue.tsv done



## TASK 7. Use R to compute the mean and the median of the distances stored in regulatoryElements.genes.distances.tsv
# Open R
R

			# sigmoid colon
data_sigmoidcolon <- read.table("analyses/peaks.analyses/regulatoryElements.genes.distances.sigmoid_colon.tsv", sep="\t")

# MEAN
mean <- mean(data_sigmoidcolon[,2])
mean 
			# 155141.2
# MEDIAN
median <- median(data_sigmoidcolon[,2])
median
			# 51573


			# stomach
data_stomach <- read.table(""analyses/peaks.analyses/regulatoryElements.genes.distances.stomach.tsv", sep="\t")

# MEAN
mean <- mean(data_stomach[,2])
mean 
			# 142698.2
# MEDIAN
median <- median(data_stomach[,2])
median
			# 47592



