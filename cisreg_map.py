#!/usr/bin/env python3

"""
cisreg_map.py

Author: Maria Rossello
Date created: March 2024

Description:
This script processes genomic data to define regulatory zones (promoter, proximal, and gene body) from an annotation file 
and associates peaks with genes and regulatory zones. It generates two output tables: one associating peaks with genes 
and another associating peaks with regulatory zones. Additionally, it provides a summary of how the regulatory zones 
were defined.

Usage:
python cisreg_map.py -t <tss_bed_file> -g <genes_bed_file> -p <peaks_bed_file> -l <len_file> [-o <output_file>]

Arguments:
    -t, --tss_bed <tss_bed_file>: Path to the transcription start sites (TSS) BED file.
    -g, --genes_bed <genes_bed_file>: Path to the genes BED file.
    -p, --peaks_bed <peaks_bed_file>: Path to the peaks BED file.
    -l, --len_file <len_file>: Path to the tab-delimited file specifying the length of each chromosome or scaffold.
    -o, --output <output_file>: Path and suffix of the output files. Default is generated based on input files.

Requirements:
    - Python 3.9 or later
"""


#########################################################################################################
# PROGRAM ARGUMENTS
#########################################################################################################

# Importing necessary modules
import argparse
import os 
from pybedtools import BedTool
import sys
import pandas as pd
import warnings
from datetime import datetime
warnings.filterwarnings("ignore", category=RuntimeWarning)

parser = argparse.ArgumentParser(prog='cisreg_map.py',
                                 description='''This script processes genomic data to define regulatory zones 
                                 (promoter, proximal, and gene body) from an annotation file and associates peaks 
                                 with genes and regulatory zones. It generates two output tables: one associating 
                                 peaks with genes and another associating peaks with regulatory zones. Additionally, 
                                 it provides a summary of how the regulatory zones were defined.''',
                                 formatter_class=argparse.RawTextHelpFormatter)

required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

required.add_argument('-t', '--tss_bed',
                      help='Path to the transcription start sites (TSS) BED file.',
                      type=argparse.FileType('r'))

required.add_argument('-g', '--genes_bed',
                      help='Path to the genes BED file.',
                      type=argparse.FileType('r'))

required.add_argument('-p', '--peaks_bed',
                      help='Path to the peaks BED file.',
                      type=argparse.FileType('r'))

required.add_argument('-l', '--len_file',
                      help='Path to the tab-delimited file specifying the length of each chromosome or scaffold.',
                      type=argparse.FileType('r'))

optional.add_argument('-o', '--output',
                      help='Path and suffix of the output files.',
                       type=str, default=None)

args = parser.parse_args()

# Set default output filename if not provided
if args.output is None:
    output_basename = "ATAC"
    args.output = output_basename + "_reg_zones"

# Check if any of the required arguments are not provided
if args.tss_bed is None or args.genes_bed is None or args.peaks_bed is None or args.len_file is None:
    print("ERROR: Some required arguments are not supplied.")
    parser.print_help()
    sys.exit(1)


#########################################################################################################
# FUNCTIONS
#########################################################################################################

#--------------------------------------------------------------------------------------------------------
# File Handling Functions
#--------------------------------------------------------------------------------------------------------

def get_file_extension(file_path):
    """
    Get the extension of a file.

    Args:
        file_path: Path to the file.

    Returns:
        Lowercase file extension (e.g., ".bed").
    """

    _, extension = os.path.splitext(file_path)
    return extension.lower()

def get_col_num(file_path, delimiter='\t'):
    """
    Get the number of columns in a file.

    Args:
        file_path: Path to the file.
        delimiter: Delimiter used in the file. Default is tab.

    Returns:
        Number of columns in the file.
    """

    with open(file_path, 'r') as file:
        first_line = file.readline().rstrip('\n')
        num_columns = len(first_line.split(delimiter))
    return int(num_columns)

def open_bed(input_file):
    """
    Open and process a BED file.

    Args:
        input_file: File object of the BED file.

    Returns:
        Sorted BED file as a pybedtools object.
    """

    extension = get_file_extension(input_file.name)
    col_num=get_col_num(input_file.name, delimiter='\t')

    if extension.startswith(".bed"):
        #we have a bed file

        if col_num-1 < 3:
            print("ERROR: Bed files does not have at least 3 tab delimited fields")
            sys.exit(1)
        
        elif col_num-1 <= 6:
            bed_file = BedTool(input_file)
            bed_sorted=sort_bed_file(bed_file)

        else:
            bed_file = BedTool(input_file).cut(range(6))
            bed_sorted=sort_bed_file(bed_file)
        


    else:
        #file not recognized. Kill the script
        print("ERROR: File not recognized. Use a .bed file with the correct extension")
        sys.exit(1)

    return bed_sorted

#--------------------------------------------------------------------------------------------------------
# Utility Functions
#--------------------------------------------------------------------------------------------------------

def from_pybedtools_object_to_list(pybedtools_object):
    """
    Convert a pybedtools object to a list of BED format strings.

    Args:
        pybedtools_object: A pybedtools object.

    Returns:
        List of BED format strings.
    """

    str_bed = str(pybedtools_object)
    lines = str_bed.strip().split("\n")
    return lines

def sort_bed_file(bed_file):
    """
    Sort a BED file.

    Args:
        bed_file: A list of BED format strings or a pybedtools object.

    Returns:
        Sorted BED file as a list of strings.
    """
    
    if isinstance(bed_file, list):
        bed_sorted = sorted(bed_file, key=lambda x: (x.split('\t')[0], int(x.split('\t')[1]), int(x.split('\t')[2])))
    
    elif isinstance(bed_file, BedTool):
        str_bed=from_pybedtools_object_to_list(bed_file)
        bed_sorted_tmp = sorted(str_bed, key=lambda x: (x.split('\t')[0], int(x.split('\t')[1]), int(x.split('\t')[2])))
        bed_sorted = BedTool('\n'.join(bed_sorted_tmp), from_string=True)
    return bed_sorted

#--------------------------------------------------------------------------------------------------------
# Annotation Functions
#--------------------------------------------------------------------------------------------------------

def add_type_to_name(bed_list, type_of_zone):
    """
    Add a type of regulatory zone to the name field of each BED entry.

    Args:
        bed_list: List of BED format strings.
        type_of_zone: Type of regulatory zone (e.g., promoter, proximal, gene body).

    Returns:
        Modified list of BED format strings with type added to the name field.
    """

    modified_lines = []

    for line in bed_list:

        fields = line.split()
        fields[3] += str("_"+type_of_zone)
        modified_line = "\t".join(fields)
        modified_lines.append(modified_line)

    return modified_lines

def remove_type_to_name(bed_list):
    """
    Remove the type of regulatory zone from the name field of each BED entry.

    Args:
        bed_list: List of BED format strings.

    Returns:
        Modified list of BED format strings with type removed from the name field.
    """

    modified_lines = []

    for line in bed_list:

        fields = line.split()
        fields[3] = '_'.join(fields[3].split('_')[:-1])
        modified_line = "\t".join(fields)
        modified_lines.append(modified_line)

    return modified_lines

#--------------------------------------------------------------------------------------------------------
# Dictionary Functions
#--------------------------------------------------------------------------------------------------------

def remove_duplicate_values_inplace(dictionary):
    """
    Remove duplicate values from lists in a dictionary.

    Args:
        dictionary: Dictionary with lists as values.

    Returns:
        None. The function modifies the dictionary in place.
    """

    for key, value_list in dictionary.items():
        dictionary[key] = list(set(value_list))

#--------------------------------------------------------------------------------------------------------
# Cis-Regulatory Region Functions
#--------------------------------------------------------------------------------------------------------

def make_cis_reg_region(initial_bed, start, end, type_of_zone):
    """
    Generate cis-regulatory regions.

    Args:
        initial_bed: Initial BED file as a pybedtools object.
        start: Length of the region to include upstream of the feature.
        end: Length of the region to include downstream of the feature.
        type_of_zone: Type of regulatory zone (promoter, proximal, gene body).

    Returns:
        BED file of the generated cis-regulatory regions as a pybedtools object.
    """

    zone_bed=initial_bed.slop(g=args.len_file.name,
                        l=start, 
                        r=end,
                        s=True)

    zone_bed_sorted=sort_bed_file(zone_bed)
    zone_output_list=from_pybedtools_object_to_list(zone_bed_sorted)

    zone_list=add_type_to_name(zone_output_list, type_of_zone)

    zone_bed=BedTool(zone_list)
    zone_bed.saveas(str(args.output + "_" + type_of_zone + ".bed"))

    return zone_bed

#--------------------------------------------------------------------------------------------------------
# Cis-Regulatory Region Correction Functions
#--------------------------------------------------------------------------------------------------------

def get_overlapping_with_tss_list(zone_bed, tss_bed):
    """
    Get a list of cis-regulatory regions overlapping with transcription start sites (TSS).

    Args:
        zone_bed: BED file of cis-regulatory regions as a pybedtools object.
        tss_bed: BED file of TSS as a pybedtools object.

    Returns:
        List of overlapping cis-regulatory regions.
    """
 
    #find overlapping regions in the file
    closer_output=zone_bed.closest(tss_bed, 
                                        D='ref', 
                                        N=True,
                                        s=True,
                                        t='first',
                                        iu=True)

    closer_output_sorted=sort_bed_file(closer_output)
    closer_output_list=from_pybedtools_object_to_list(closer_output_sorted)
    
    #the zones that are overlapping but not with itself
    overlapping_zones_list=[]
    for line in closer_output_list:
        columns = line.split('\t')
        columns[3] = '_'.join(columns[3].split('_')[:-1])
        if columns[12] == '0' and columns[3] != columns[9]:
            overlapping_zones_list.append(line)

    return overlapping_zones_list

def make_dict_correct_tss_overlap(overlapping_zones_list):
    """
    Create a dictionary of corrected coordinates for overlapping cis-regulatory regions.

    Args:
        overlapping_zones_list: List of overlapping cis-regulatory regions.

    Returns:
        Dictionary containing corrected start and end coordinates for overlapping zones.
    """

    tss_overlap_dict = {}
    
    #correct the overlapping zones and save it in a dict
    for line in overlapping_zones_list:

        columns = line.split('\t')
 
        if columns[5]=="+" and int(columns[1]) < int(columns[7]) and int(columns[2]) != int(columns[8]):
            #(+) modify previous gene


            #modify end to avoid overlap
            name=columns[3]
            start = int(columns[8]) + 1
            end=int(columns[2])

            tss_overlap_dict[name] = [start,end]


        elif columns[5]=="-" and int(columns[1]) > int(columns[7]):
            #(-) modify next gene

            #modify start to avoid overlap
            name=columns[3]
            start=int(columns[1])
            end = int(columns[7]) - 1

            tss_overlap_dict[name] = [start,end]

    return tss_overlap_dict

def correct_tss_overlap(zone_bed, tss_bed):
    """
    Correct overlapping cis-regulatory regions with transcription start sites (TSS).

    Args:
        zone_bed: BED file of cis-regulatory regions as a pybedtools object.
        tss_bed: BED file of TSS as a pybedtools object.

    Returns:
        Corrected BED file of cis-regulatory regions as a pybedtools object.
    """

    overlapping_zones_list=get_overlapping_with_tss_list(zone_bed , tss_bed)

    dict_of_coord=make_dict_correct_tss_overlap(overlapping_zones_list)

    corrected_bed_list=[]
    
    bedfile_list=from_pybedtools_object_to_list(zone_bed)

    #correct the overlapping zones and save it in a dict
    for line in bedfile_list:

        columns = line.split('\t')
        
        name=columns[3]

        if name in dict_of_coord.keys():
            columns[1]=str(dict_of_coord[name][0])
            columns[2]=str(dict_of_coord[name][1])
        
        corrected_bed_list.append('\t'.join(columns))
    
    corrected_bedfile=BedTool(corrected_bed_list)
    
    return corrected_bedfile

#--------------------------------------------------------------------------------------------------------
# Peak-to-Regulatory-Region Association Functions
#--------------------------------------------------------------------------------------------------------

def make_peak2gene_to_df(intersect_output_list):
    """
    Convert the output of intersecting peaks with regulatory regions into a DataFrame containing peak-to-gene associations.

    Args:
        intersect_output_list: List of intersecting peaks with regulatory regions.

    Returns:
        DataFrame containing peak-to-gene associations.
    """

    genes_dict={}

    for line in intersect_output_list:
        columns = line.split('\t')

        peak=columns[3]
        gene='_'.join(columns[7].split('_')[:-1])       

        if peak not in genes_dict:
            genes_dict[peak] = [gene]

        else:
             genes_dict[peak].append(gene)
    
    remove_duplicate_values_inplace(genes_dict)

    zones_data = {key: ', '.join(value) for key, value in genes_dict.items()}
    df = pd.DataFrame(zones_data.items(), columns=['peak', 'genes'])

    return df

def update_peak2zone_dictionary(intersect_output_list, dictionary):
    """
    Update the peak-to-zone dictionary with new associations from intersecting peaks with regulatory regions.

    Args:
        intersect_output_list: List of intersecting peaks with regulatory regions.
        dictionary: Dictionary containing peak-to-zone associations.

    Returns:
        Updated dictionary containing peak-to-zone associations.
    """

    for line in intersect_output_list:
        columns = line.split('\t')
        
        peak=columns[3]
        zone=columns[7].split('_')[-1].strip()
        
        if peak not in dictionary.keys():
            dictionary[peak] = [zone]
    
    remove_duplicate_values_inplace(dictionary)

    return dictionary

def associate_peak_to_zone(all_peaks_bed, zone_bed):
    """
    Associate peaks with regulatory zones.

    Args:
        all_peaks_bed: BED file of all peaks as a pybedtools object.
        zone_bed: BED file of regulatory zones as a pybedtools object.

    Returns:
        List of peaks associated with regulatory zones.
    """

    zone_peaks=all_peaks_bed.intersect(zone_bed,
                                       wa=True,
                                       wb=True,
                                       f=0.6)

    zone_peaks_sorted=sort_bed_file(zone_peaks)
    zone_peaks_list=from_pybedtools_object_to_list(zone_peaks_sorted)
    return zone_peaks_list

#########################################################################################################
# PROCESS DATA
#########################################################################################################

#-------------------------------------------------------------------------------
# Open files and define important variables
#-------------------------------------------------------------------------------

# Open files
tss_bed = open_bed(args.tss_bed)
all_peaks_bed = open_bed(args.peaks_bed)
genes_bed = open_bed(args.genes_bed)


#-------------------------------------------------------------------------------
# Create DataFrames and dictionaries to store peak associations
#-------------------------------------------------------------------------------

# Get al peaks
peak_values = []        
for line in from_pybedtools_object_to_list(all_peaks_bed):
    columns = line.split('\t')
    peak_value = columns[3]
    peak_values.append(peak_value)

# Create DataFrame to store peak-to-gene associations
peak2gene_df = pd.DataFrame()
peak2gene_df['peak'] = peak_values

# Create DataFrame to store peak-to-zone associations
peak2zone_df = pd.DataFrame()
peak2zone_df['peak'] = peak_values

peak2zone_dict={}


#-------------------------------------------------------------------------------
# Work with the promoter zone
#-------------------------------------------------------------------------------

# Promoter is defined as 1000 bases upstream and 500 bases downstream of the TSS
promoter_start = 1000
promoter_end = 500

# Get the promoter cis reg (output bed is saved)
promoter_bed = make_cis_reg_region(tss_bed, 
                                start=promoter_start,
                                end=promoter_end,
                                type_of_zone="promoter")

# Get a list of peaks associated with the promoter
peaks_in_promoter_list = associate_peak_to_zone(all_peaks_bed, promoter_bed)

# Add the genes associated with each peak to peak2gene_df
promoter_df = make_peak2gene_to_df(peaks_in_promoter_list)
promoter_df.rename(columns={'genes': 'promoter'}, inplace=True)
peak2gene_df = pd.merge(peak2gene_df, promoter_df, on='peak', how='left')

# Add the zone associated with each peak to peak2zone_dict
update_peak2zone_dictionary(peaks_in_promoter_list, peak2zone_dict)


#-------------------------------------------------------------------------------
# Work with the proximal region
#-------------------------------------------------------------------------------

# Proximal region is defined as 4000 bases away from the TSS but not overlapping with the promoter
proximal_start = 5000
proximal_end = (promoter_start)*(-1)


# Get the proximal cis reg (output bed is saved)
proximal_bed = make_cis_reg_region(tss_bed, 
                                start=proximal_start,
                                end=proximal_end,
                                type_of_zone="proximal")


# Correct overlapping regions with neighboring TSS
proximal_bed_fixed = correct_tss_overlap(proximal_bed , tss_bed)

# Get a list of peaks associated with the proximal region
peaks_in_proximal_list=associate_peak_to_zone(all_peaks_bed, proximal_bed_fixed)

# Add the genes associated with each peak to peak2gene_df
proximal_df = make_peak2gene_to_df(peaks_in_proximal_list)
proximal_df.rename(columns={'genes': 'proximal'}, inplace=True)
peak2gene_df = pd.merge(peak2gene_df, proximal_df, on='peak', how='left')

# Add the zone associated with each peak to peak2zone_dict
update_peak2zone_dictionary(peaks_in_proximal_list, peak2zone_dict)


#-------------------------------------------------------------------------------
# Work with the gene body region
#-------------------------------------------------------------------------------

# Gene body is defined as the gene zone that is not promoter
gene_body_start = (promoter_end) * (-1)
gene_body_end = 0

# Ensure that the gene body region is not exceeding the boundaries
genes_bed_filt = genes_bed.filter(lambda x: len(x) > int(gene_body_start**(-1)))

# Get the gene body cis reg (output bed is saved)
gene_body_bed = make_cis_reg_region(genes_bed_filt,
                                  start=gene_body_start,
                                  end=gene_body_end,
                                  type_of_zone="genebody")

# Get a list of peaks associated with the gene body region
peaks_in_gene_body_list = associate_peak_to_zone(all_peaks_bed, gene_body_bed)

# Add the genes associated with each peak to peak2gene_df
gene_body_df = make_peak2gene_to_df(peaks_in_gene_body_list)
gene_body_df.rename(columns={'genes': 'gene_body'}, inplace=True)
peak2gene_df = pd.merge(peak2gene_df, gene_body_df, on='peak', how='left')

# Add the zone associated with each peak to peak2zone_dict
update_peak2zone_dictionary(peaks_in_gene_body_list, peak2zone_dict)


#-------------------------------------------------------------------------------
# Save the output tables
#-------------------------------------------------------------------------------

# Save the peak2gene dataframe
peak2gene_df.to_csv(str(args.output + "_" + 'peaks_and_genes.tsv'), sep='\t', index=False)

# Save the peak2gene dictionary
peak2zone_data = {key: ', '.join(value) for key, value in peak2zone_dict.items()}
peak2zone_data_df = pd.DataFrame(peak2zone_data.items(), columns=['peak', 'zone'])

peak2zone_df = pd.merge(peak2zone_df, peak2zone_data_df, on='peak', how='left')
peak2zone_df_final = peak2zone_df.fillna("distal")

peak2zone_df_final.to_csv(str(args.output + "_"  +'peaks_and_zones.tsv'), sep='\t', index=False)


#-------------------------------------------------------------------------------
# Finishing script
#-------------------------------------------------------------------------------

# A table summarizing the results
print("This is how the regulatory zones were defined")
table_header = "{:<12} {:<10} {:<10} {:<10}".format("Type", "From", "Start", "End")
promoter_row = "{:<12} {:<10} {:<10} {:<10}".format("Promoter", "tss", promoter_start, promoter_end)
proximal_row = "{:<12} {:<10} {:<10} {:<10}".format("Proximal", "tss", proximal_start, proximal_end)
gene_row = "{:<12} {:<10} {:<10} {:<10}".format("Gene body", "gene", gene_body_start, gene_body_end)

print(table_header)
print("-" * len(table_header))  # Separator line
print(promoter_row)
print(proximal_row)
print(gene_row)

print("\nYour files are saved as: " + str(args.output))

current_time = datetime.now()
time_str = current_time.strftime("%Y-%m-%d %H:%M:%S")

print("Finished at: " + time_str +
      "\n\nHappy analyzing.\nBYE! (>ᴗ•)❀")
