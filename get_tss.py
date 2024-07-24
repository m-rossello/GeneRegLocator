#!/usr/bin/env python3

"""
get_tss.py

Author: Maria Rossello
Date created: March 2024

Description:
    This script parses gene annotations from a GFF3 file, identifies transcription start sites (TSS), 
    converts them to BED format, and saves or prints the output.

Usage:
    python get_tss.py -i input.gff3 [-o output.bed]

Arguments:
    -i, --input_file: Path to the input GFF3 file containing gene annotations.
    -o, --output_file: (Optional) Path to the output BED file. If not provided, the output will be printed to stdout.

Requirements:
    - Python 3.9 or later
"""


#########################################################################################################
# PROGRAM ARGUMENTS
#########################################################################################################

# Importing necessary modules
import argparse
import sys

# Parsing command-line arguments
parser = argparse.ArgumentParser(prog='get_tss.py', 
                                 description='Find transcription start sites from gene models in a GFF3 file')

# Adding command-line arguments
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

# Defining required and optional arguments
required.add_argument('-i', '--input_file',
                      help= 'Path to the input GFF3 file. The file must be in GFF3 format containing gene models to extract transcription start sites.', 
                      type=argparse.FileType('r'))

optional.add_argument('-o','--output_file', 
                      help='Path to the output GFF file. If not specified, output will be written to stdout.')

# Parsing arguments
args = parser.parse_args()

# Checking if input file is provided
if args.input_file is None:
    parser.print_help()
    sys.exit(1)

#########################################################################################################
# FUNCTIONS
#########################################################################################################

def convert_to_bed(gff_list):
    """
    Convert GFF entries to BED format.

    Args:
        gff_list (list): A list of GFF entries.

    Returns:
        list: A list of entries converted to BED format.
    """

    bed_list = []

    for entry in gff_list:
        # Checking if the GFF entry has sufficient fields
        if len(entry) >= 9:
            chrom = entry[0]
            # Converting start position to BED format (0-based)
            start = str(int(entry[3]) - 1)  
            # Converting end position to BED format (0-based)
            end = str(int(entry[4]) - 1) 
            # Extracting gene name from GFF attributes
            name = entry[8].split(';')[0].split('=')[1].strip()  
            score = entry[5]
            strand = entry[6]

            # Creating BED entry
            out_print = [str(x) for x in [chrom, start, end, name, score, strand]]
            bed_list.append('\t'.join(out_print))

    return bed_list


#########################################################################################################
# PROCESS DATA
## This section processes the input data to identify the transcription start sites (TSS) of genes.
#########################################################################################################

tss_list = []  # Initialize a list to store TSS information

# Iterate over each line in the input GFF file
for line in args.input_file:
    if not line.startswith("#"):
        columns = line.split('\t')  # Split the line into columns using tab delimiter
    
        if columns[2] == "gene":  # Check if the entry is a gene
        
            # Modify the entry to represent a TSS
            columns[2] = "TSS"
            start = columns[3]
            end = columns[4]
            strand = columns[6]
    
            # Adjust start and end positions based on strand
            if strand == "+":  # For positive strand genes
                columns[4] = int(start) + 1  # Set TSS position
                tss_list.append(columns)  # Add TSS entry to the list
    
            elif strand == "-":  # For negative strand genes
                columns[3] = int(end) - 1  # Set TSS position
                tss_list.append(columns)  # Add TSS entry to the list

# Convert TSS list to BED format
bed_file = convert_to_bed(tss_list)

#--------------------------------------------------------------------------------------------------------
#save/print the output
#--------------------------------------------------------------------------------------------------------

try:
    if args.output_file:  # If output file path is provided

        with open(args.output_file, 'w') as f:

            # Write each BED line to the output file
            for bed_line in bed_file:
                f.write(bed_line + '\n')

    else:  # If no output file path is provided, print to stdout
        sys.stdout.write('\n'.join(bed_file) + "\n")

except BrokenPipeError:
    # Handle BrokenPipeError
    pass



