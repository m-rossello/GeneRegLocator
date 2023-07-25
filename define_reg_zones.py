#!/usr/bin/env python3

# File name: define_reg_zones.py
# Author: Maria Rossello
# Date created: 17/07/2023
# Python Version: 3.9.16

import sys
import argparse
from pybedtools import BedTool
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)

# program arguments
parser = argparse.ArgumentParser(prog='define_reg_zones.py',
                                 description='''Define Regulatory Zones is a program that extends the coordinates of a bed file containing TSS coordinates to define regulatory zones. 
                                 Users can specify the extension length and optionally filter the extended regions based on additional coordinates. 
                                 The program outputs a bed file with the extended coordinates and associated information for each region.''')

required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

required.add_argument('-a', '--file_to_transform',
                      help='Bed file containing the regulatory zones, typically TSS coordinates.',
                      type=argparse.FileType('r'))

required.add_argument('-s', '--start_add',
                      help='The number of base pairs to operate from the start coordinate.'
                           'To subtract, provide a negative number; to add, give a positive number. Integer.',
                      type=int)

required.add_argument('-e', '--end_add',
                      help='The number of base pairs to operate to the end coordinate.'
                           'To subtract, provide a negative number; to add, give a positive number. Integer.',
                      type=int)

required.add_argument('-l', '--len_file',
                      help='Tab-delimited file specifying the length of each chromosome or scaffold.',
                      type=argparse.FileType('r'))

optional.add_argument('-f', '--filtering_zone',
                    help='Bed file with coordinates for additional filtering. '
                         'The program will stop the zone when it encounters the start of this feature.'
                         'Must have the same features that ther re in file_to_transform',
                    type=argparse.FileType('r'))

optional.add_argument('-o', '--output',
                    help='Path and name of the output bed file',
                    type=str)

args = parser.parse_args()

if args.file_to_transform is None or args.start_add is None or args.end_add is None:
    parser.print_help()
    sys.exit(1)

#some basic functions that we will be using
def from_pybedtools_object_to_list(pybedtools_object):
    str_bed = str(pybedtools_object)
    lines = str_bed.strip().split("\n")
    return lines

def sort_bed_file(bed_file, from_pybedtools=False):
    if not from_pybedtools:
        bed_sorted = sorted(bed_file, key=lambda x: (x.split('\t')[0], int(x.split('\t')[1]), int(x.split('\t')[2])))
    else:
        str_bed=from_pybedtools_object_to_list(bed_file)
        bed_sorted_tmp = sorted(str_bed, key=lambda x: (x.split('\t')[0], int(x.split('\t')[1]), int(x.split('\t')[2])))
        bed_sorted = BedTool('\n'.join(bed_sorted_tmp), from_string=True)
    return bed_sorted

#functions to import files
def open_bed(bed_file):
    bed_file = BedTool(bed_file).cut(range(6))
    bed_sorted=sort_bed_file(bed_file, from_pybedtools=True)
    return bed_sorted

def scaffold_dict(scaffold_length_file):
    scaffold_lengths_dict = {}
    for line in scaffold_length_file:
        scaffold, length = line.strip().split('\t')
        scaffold_lengths_dict[scaffold] = int(length)
    return scaffold_lengths_dict

#functions to perform trimming
def trim_upstream_overlapping (bed1, bed2, scaffold_lengths):

    """Modifies 5' genomic intervals in bed1 based on overlaps with intervals in bed2, using scaffold lengths for trimming."""

    trim_up = []

    for i in range(len(bed1)):
        line1 = bed1[i]
        fields1 = line1.strip().split('\t')
        scaffold = fields1[0]

        if scaffold in scaffold_lengths:
            #get scaffold length
            scaffold_length = scaffold_lengths[scaffold]

        if i!= len(bed1) and i+1 < len(bed1):
            #not last line
            line_up=bed1[i+1]
            fields_up = line_up.strip().split('\t')

            if fields1[0] == fields_up[0] and int(fields1[2]) > int(fields_up[1]):
                #just overlaps
                line2 = bed2[i+1]
                fields2 = line2.strip().split('\t')

                if fields1[0] == fields2[0]:
                    #just in the same same scaffold
                    start_trim = max(min(int(fields2[1]), scaffold_length),int(fields1[1])+1)
                    fields1[2] = str(start_trim)
                    trim_up.append('\t'.join(fields1))

                else:
                    #in different scaffolds do not trim
                    trim_up.append('\t'.join(fields1))

            else:
                #no overlap
                trim_up.append('\t'.join(fields1))
        else:
            #last line
            trim_up.append('\t'.join(fields1))

    #sort the trimmed file
    sorted_trim_up = sort_bed_file(trim_up)

    return sorted_trim_up

def trim_downstream_overlapping (bed1, bed2, scaffold_lengths):

    """Modifies 3' genomic intervals in bed1 based on overlaps with intervals in bed2, using scaffold lengths for trimming."""

    trim_down = []
    for i in range(len(bed1)):
        line1 = bed1[i]
        fields1 = line1.strip().split('\t')
        scaffold = fields1[0]

        if scaffold in scaffold_lengths:
            #get scaffold length
            scaffold_length = scaffold_lengths[scaffold]
        
        if i != 0 and i+1 < len(bed1):
            #not first line
            line_down=bed1[i-1]
            fields_down = line_down.strip().split('\t')

            if fields1[0] == fields_down[0] and int(fields1[1]) < int(fields_down[2]):
                #just overlaps
                line2 = bed2[i-1]
                fields2 = line2.strip().split('\t')
                if fields1[0] == fields2[0]:
                    #just in the same same scaffold
                    start_trim = min(int(fields2[2]), scaffold_length, int(fields1[2])-1)
                    fields1[1] = str(start_trim)
                    trim_down.append('\t'.join(fields1))

                else:
                    #in different scaffolds do not trim
                    trim_down.append('\t'.join(fields1))

            else:
                #no overlap
                trim_down.append('\t'.join(fields1))
        else:
            #first line
            trim_down.append('\t'.join(fields1))

    #sort the trimmed file
    sorted_trim_down = sort_bed_file(trim_down)

    return sorted_trim_down

# functions to parse outputs
def check_for_closer_overlaps(closer_output):
    for item in closer_output:
        columns = item.split('\t')
        if len(columns) >= 13 and columns[12] == '0':
            return True
    return False

def filter_closer_output(closer_output):
    closer_output_sorted=sort_bed_file(closer_output, from_pybedtools=True)
    closer_output_list=from_pybedtools_object_to_list(closer_output_sorted)
    
    if check_for_closer_overlaps(closer_output_list):
        closer_output_filt=[]
        seen_combinations = {}

        #filter
        for line in closer_output_list:
            columns = line.strip().split('\t')
            if columns[12] == '0' and int(columns[1]) <= int(columns[7]) and int(columns[2]) <= int(columns[8]):
                #get overlapping genes downstream
                key = tuple(columns[0:2])
                if key not in seen_combinations:
                    #get one entry per gene
                    seen_combinations[key] = 1
                    closer_output_filt.append(line)
        closer_output_mod= []

        #modify end
        for line in closer_output_filt:
            columns = line.strip().split('\t')
            trim = int(columns[7]) - 1
            columns[2] = str(trim)
            modified_line = '\t'.join(columns)
            closer_output_mod.append(modified_line)

        #get bed6
        closer_output_bed = BedTool(closer_output_mod).cut(range(6))

        #sort
        closer_filtered_sorted=sort_bed_file(closer_output_bed, from_pybedtools=True)

    else:
        #if there are no overlaps just get bed6
        closer_filtered_sorted = BedTool(closer_output_list).cut(range(6))
    
    return closer_filtered_sorted

def replace_matching_lines(bed1_lines, bed2_lines):
    """Replace matching lines in BED1 with lines from BED2."""
    gene_ids_bed2 = {line.split()[3]: line for line in bed2_lines}
    new_bed1_lines = []
    for line in bed1_lines:
        gene_id = line.split()[3]
        if gene_id in gene_ids_bed2:
            new_bed1_lines.append(gene_ids_bed2[gene_id])
        else:
            new_bed1_lines.append(line)
    return new_bed1_lines

#open files
initial_bed = open_bed(args.file_to_transform)

if args.filtering_zone:
    filtering_zones = open_bed(args.filtering_zone)
else:
    filtering_zones=initial_bed

scaffold_lengths=scaffold_dict(args.len_file)

#extend the file to obtain the regulatory zone
start_add_corrected=args.start_add*(-1)
broad_zones=initial_bed.slop(g=args.len_file.name, 
                         l=start_add_corrected, 
                         r=int(args.end_add))

# eliminate multiple overlapping
trim_1=trim_upstream_overlapping (from_pybedtools_object_to_list(broad_zones),
                                  from_pybedtools_object_to_list(filtering_zones), 
                                  scaffold_lengths)

trim_2=trim_downstream_overlapping (trim_1,
                                    from_pybedtools_object_to_list(filtering_zones), 
                                    scaffold_lengths)
trimmed_bed = '\n'.join(trim_2)
trimmed_zones=BedTool(trimmed_bed, from_string=True)

#identify single overlaps and reduce their coordinates
near_zones=trimmed_zones.closest(trimmed_zones, D='ref', N=True)
corrected_overlapping_zones=filter_closer_output(near_zones)

final_zones=replace_matching_lines(trim_2, from_pybedtools_object_to_list(corrected_overlapping_zones))
final_zones_sorted = sort_bed_file(final_zones)

#save/print the output
try:
    if args.output:
        BedTool(final_zones_sorted).saveas(args.output)
    else:
        bed_zones = '\n'.join(final_zones_sorted)
        sys.stdout.write(bed_zones)
except BrokenPipeError:
    pass
