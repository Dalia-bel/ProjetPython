import sys    # sys to retrieve our command at the end of the code
import matplotlib  # we import the matplotlib library to create graphs
matplotlib.use('TkAgg')   # use a stable graphical backend
import matplotlib.pyplot as plt  # we import the simplified pyplot interface to easily plot graphs

# Q1 count the number of mapped reads:

def countMappedReads(sam_file):  
    mapped_count = 0  # we initialize our counter to 0

    with open(sam_file, 'r') as file:  # we open the sam_file and place it in the parameter file; 'with as' is used to open the file in the variable file, and 'r' is for reading the file without modification (if we had used 'w', we could write to the file)
        for line in file:
            if line.startswith('@'):  # we ignore any line starting with '@' (the header)
                continue 
            
            columns = line.strip().split('\t') 
            # we remove leading and trailing spaces from the string using strip
            # the column becomes a list of strings split at each tab, as SAM is a tab-delimited file
            flag = int(columns[1])  # we convert the string in the flag column to an integer to process the bits
            
            if (flag & 4) == 0:  # the read is mapped if the 3rd bit from the right is 0; flag=4 indicates unmapped since, in binary, the 3rd bit is 1
                mapped_count += 1  # we add 1 to the counter

    print(f"Total number of mapped reads: {mapped_count}")  # f-string to interpret the value inside {}

    

#Q2 count the number of raeds per flag: translate the flag here and summarise only those that make sense
def countFlags(sam_file, minQ=0):
   # we create a dictionary to count the cases where the two pairs of reads are mapped or only one of the two is mapped
    paired_status = {
        "both_mapped": 0,   # we create a key for the two reads in the pair that are mapped and initialise the key value to 0
        "one_mapped": 0,    # we create a key for a single read in the pair and initialise the key value to 0
    }

    read_pairs = {}  # create an empty dictionary to store read flags by their QNAME because a pair of reads have the same Qname
    with open(sam_file, 'r') as file: #we open the sam_file file and place it in the file/ parameter with as c to open the file in the file variable and the ‘r’ is used to read the file if we had put w we could write it as follows
        for line in file:
            if line.startswith('@'): 
                continue  #here it reads line by line from the c=sam file and each time it finds an @ it ignores the line because it is the leading lines that start with an @.
            
            columns = line.strip().split('\t')  #we remove the spaces at the beginning and end with strip of the string and the column becomes a list of strings cut at each tab since SAM is a tabulated file.
            qname = columns[0]  
            flag = int(columns[1])  # transform the FLAG into an integer

            mapq = int(columns[4])  # we indicate that the mapping quality is in the 5th column (index 4)
            if mapq < minQ: # skip la ligne si la qualité est trop faible
                continue

            
            if qname not in read_pairs: # here we'll group together the flags of reads with the same QNAME (foward and reverse) and create a condition to check whether the qname already exists as a key in the read_pairs dictionary/ If it doesn't exist, this will be the first time we've encountered this read.
                read_pairs[qname] = [] # in the read_pairs dictionary, create an empty list associated with the key that is the qname encountered
            read_pairs[qname].append(flag) # append can be used to add the flag at the end of the list, which has the qname encountered as its key.

    # Analyse pairs of reads
    for qname, flags in read_pairs.items(): # for each qname, check whether the reads are even (if the reverse and foward exist) and ignore cases where there are no two reads for a pair
        if len(flags) != 2: #if the number of flags for each qname is not equal to 2, we ignore this qname 
            continue
        
         # We check whether the two reads in the pair with the same Qname are mapped by creating two variables and performing a Boolean operation between the read flag and bit 4.       
        first_mapped = (flags[0] & 0x4)==0  # the first read is mapped if the first read flag & 0x4 == 0 
        second_mapped = (flags[1] & 0x4)==0  # The second read is mapped if the second read flag & 0x4 == 0

        if first_mapped and second_mapped: #here we'll create a condition to find out if the two reads are mapped and add a +1 cursor when this is the case.
            paired_status["both_mapped"] += 1 #if the booleene condition returns true for both mapped reads, we add plus 1 to the cursor of the paired status dictionary key.
        elif first_mapped or second_mapped: # if not, if the Boolean condition returns true, only one of the reads is mapped, we add a cursor to the one mapped key of the paired status dictionary. 
            paired_status["one_mapped"] += 1

    # Show results
    print("Comment les reads et paires de reads sont mappés :")
    print(f"Paires avec les deux reads mappés : {paired_status['both_mapped']}") # f so that it takes the value between the {} and indicates that it should display the contents of the dictionary keys in string format.
    print(f"Paires avec un seul read mappé : {paired_status['one_mapped']}")

#Q3: Where are the reads mapped? Is the alignment homogeneous along the reference sequence? 
def countReadsBySegment(sam_file, segment_size=10000): #We create a function to count reads whose start position is on a chromosome segment of predefined size (here by default 10000bp). This function takes as input the SAM file and a default segment size of 10000 bases.
    chromosome_segments = {}  # here we create an empty dicitonaire to store the 10000 bp segments of the chromosome and the number of reads associated with each segment 


    with open(sam_file, 'r') as file: # open the sam file in read mode using the ‘r’ key
        for line in file: #for each line of the file it is going to analyse one by one 
            if line.startswith('@'):  #a condition is created if the line begins with an @.

                continue #it ignores it, so we'll ignore the leading lines that start with @.

            columns = line.strip().split('\t') # then remove any spaces at the beginning and end of the lines using strip and owith split divide the line into columns using the tab character as a separator
            chromosome = columns[2]  # we indicate that the name of the chromosome is in the 3rd column (index 2)
            start_position = int(columns[3])  # here we indicate that the start position of the read is in column 2 with index 3 and we transform the number into an integer with ‘int’ because in the file it is a character string.
            segment_index = start_position // segment_size # Calculates the index of the segment to which the read belongs by dividing the start position by the segment size

        
            if chromosome not in chromosome_segments: # a condition is used to check that the chromosome does not already exist in the dictionary 
                chromosome_segments[chromosome] = {}  # if not, add an entry for this chromosome with an empty dictionary for its segments           
            if segment_index not in chromosome_segments[chromosome]: #a second condition is created to check whether the segment index exists; if it does not exist, the counter is initialized in the next line. 
                chromosome_segments[chromosome][segment_index] = 0  # set the reads counter for this segment to 0
            chromosome_segments[chromosome][segment_index] += 1  # counter adds +1

    print("Distribution des reads par segment :") # the results of the reads are displayed in 10000 bp segments 
    print("Chromosome\tSegment_Start\tReads_Count") # a header is displayed with the following columns: Chromosome, Start of segment, Number of reads
    for chromosome, segments in chromosome_segments.items(): #for each chromosome in the chromosome_segments dictionary
        for segment_index, count in sorted(segments.items()): #sort segments by index using the ‘sorted’ function to display them in order
            segment_start = segment_index * segment_size     #here it will count the start position of the segment by multiplying the segment index x10000pb
            print(f"{chromosome}\t{segment_start}\t{count}") # display the results for each chromosome, the start position and the number of reads associated with this segment


# Q4: Count the number of reads per mapQ score range: here you need to make a diagram to see what this represents and do it in 10s.
def countReadsByMapQ(sam_file, minQ=0):  # Create a function that takes the SAM file as input and counts the number of reads for each MAPQ score.
    mapq_counts = {} #create an empty dictionary

    with open(sam_file, 'r') as file: #open the sam file and use the ‘r’ to read it as a file 
        for line in file: # create a ‘for’ loop where, for each line in the file 
            if line.startswith('@'): #create a condition if the line starts with an @ ignore it (leading lines)
                continue

            columns = line.strip().split('\t') # then remove any spaces at the beginning and end of the lines using strip and owith split divide the line into columns using the tab character as a separator
            mapq = int(columns[4])  # we indicate that the mapping quality is in the 5th column (index 4)

            if mapq < minQ: # skip la ligne si la qualité est trop faible
                continue

            # Calculation of the slice of 10 for each MAPQ score
            bin_range = (mapq // 10) * 10  # Round MAPQ score down to the nearest 10

            if bin_range not in mapq_counts: #create a condition if the mapq encountered is not in the map_counts dictionary
                mapq_counts[bin_range] = 0 # We initialise the counter for this MAPQ score to 0
            mapq_counts[bin_range] += 1 # add +1 to the counter

    print("Distribution des scores de qualité de mapping (MAPQ) :") #we display the results of the mapping quality score distribution 
    for bin_range, count in sorted(mapq_counts.items()): # using ‘sorted’, sort MAPQ scores in ascending order
        print(f"Tranche MAPQ {bin_range}-{bin_range + 9}: {count} reads") # display the score and the number of reads for each score and ‘f’ takes the value counted between the {}.

# Creating the histogram
    bins = sorted(mapq_counts.keys())  # we recover the sorted slices
    counts = [mapq_counts[bin] for bin in bins]  # read counts are retrieved for each slice

    # Creating the graphic
    plt.bar(bins, counts, width=8)  # Creation of a histogram with bars 8 wide
    plt.xlabel('Scores MAPQ (tranches de 10)')  # X-axis legend
    plt.ylabel('Nombre de reads')  # Y axis legend
    plt.title('Histogramme des scores MAPQ par tranches de 10')  # Graphic title
    plt.show()  # Graph display


# python3 projet.py mapping.sam   ->   sys.argv[0] = projet.py   sys.argv[1] = mapping.sam

input_sam = sys.argv[1]  # retrieve the name of the SAM file as an argument to the command
min_quality = 0 if len(sys.argv )<= 2 else int(sys.argv[2])
countMappedReads(input_sam) # Here we call the countMappedReads function, passing it the sam file (input_sam) as an argument.
countFlags(input_sam, min_quality) # Here we call the countFlags function, passing it the sam file (input_sam) as an argument.
countReadsBySegment(input_sam) # Here we call the coutReadsBySegments function, passing it the sam file (input_sam) as an argument.
countReadsByMapQ(input_sam, min_quality) # Here we call the countReadsByMapq function, passing it the sam file (input_sam) as an argument.