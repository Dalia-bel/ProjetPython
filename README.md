README

#Project Description

This script analyses SAM files, a common bioinformatics format used to store read alignments from sequencing data. It includes functions to count and evaluate mapped reads, assess mapping quality and   visualise data. The tool can handle large datasets efficiently and provides insight into the quality and distribution of sequencing reads.





#Prerequisites

To run this script, make sure you have

Python 3 installed on your system.

The matplotlib library is installed (pip install matplotlib).

Tkinter library for Python (sudo apt-get install python3-tk)





#How to use

Run the script in a terminal with the following command

python3 project.py <input_file.sam> [min_quality].

<input_file.sam>: The SAM file to analyse.

[min_quality] (optional): The minimum mapping quality (MAPQ) to consider during the analysis. The default is 0.

For example, on your terminal, type the following command:

python3 projet.py mapping.sam 20

(where 20 is the minimum mapping quality)





#Functionalities

1. Count Mapped Reads
2. 
Function: countMappedReads(sam_file)

Counts the total number of mapped reads in the given SAM file.

Mapped reads are identified by their SAM flags.

Output: The total number of mapped reads: <count>


2. Counting reads by flags
3. 
Function: countFlags(sam_file, minQ=0)

Analyses paired reads to count instances where:

Both reads in a pair are mapped.

Only one read in a pair is mapped.

Excludes reads with mapping quality below minQ.

Output:

Comment les reads et paires de reads sont mappés:

Pairs with both reads mapped: <count>

Pairs with a single read mapped : <count>


3. Analysing read distribution over reference segments
   
Function: countReadsBySegment(sam_file, segment_size=10000)

Groups reads by their start positions along the reference chromosomes, divided into segments (default size: 10,000 bp).
Summarises read counts per segment.


Output:

Distribution of reads per segment:

Chromosome Segment_Start Reads_Count

<chromosome> <segment_start> <count> 


4. Visualisation of mapping quality (MAPQ) values
   
Function: countReadsByMapQ(sam_file, minQ=0)

Groups reads by MAPQ scores into bins of 10 (e.g. 0-9, 10-19).

Generates a histogram to visualise the distribution of MAPQ scores.

Output (Console):

Distribution des scores de qualité de mapping (MAPQ) :

Tranche MAPQ <range>: <count> read

Output (Graphic):

A histogram with:

X-axis: MAPQ score ranges.

Y axis: Number of reads in each range.




#Example output

For an input SAM file with min_quality=20, the script will output

#Q1 Total number of mapped reads:350015

#Q2 Number of paired reads with both reads mapped: 154941

Number pf paired reads with one read mapped:0

#Q3 Distribution of reads across chromosome segments.

#Q4 A histogram visualising the distribution of MAPQ scores.

For both questions 3 and 4, look for answers in your terminal




#Notes

Graphical backend: The script uses the TkAgg backend for Matplotlib. Make sure your Python environment supports this




