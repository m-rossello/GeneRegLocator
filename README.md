# GeneRegLocator

Streamline the assignment of ATAC peaks to cis-regulatory regions, facilitating links between these regions and genes.

## Installation

This repository contains Python functions designed to work with Python 3 (tested with Python 3.9, 3.10, and 3.11). To utilize these functions, ensure you have the `pybedtools` library installed. For installation options, please refer to the official documentation at [pybedtools](http://daler.github.io/pybedtools/main.html).

Stay tuned for the forthcoming conda environment, which will further simplify the setup process.

## Files Description

### `get_tss.py`
The `get_tss.py` script is essential for the initial processing of cis-reg annotations. It extracts transcription start sites (TSS) from a provided GFF3 file, and converts this information into BED format for further analysis. This script supports output customization through command-line arguments, allowing users to either save the results to a file or print them directly to the console. This utility is critical for preparing the data needed for downstream regulatory mapping.

Usage:
```bash
python get_tss.py -i your_input_file.gff3 -o your_output_file.bed
```

### `cisreg_map.py`

The `cisreg_map.py` script is a comprehensive tool for mapping ATAC peaks to specific putative cis-regulatory regions within a genome.  It categorizes genomic regions into promoters, proximal areas, and gene bodies based on definitions specified:

- Promoter: Defined as 1000 bases upstream and 500 bases downstream of the TSS.
- Proximal Region: Defined as the area extending 4000 bases from the TSS but not overlapping with the promoter.
- Gene Body: Defined as the part of the gene that is not included in the promoter.

This script utilizes input files detailing TSS, genes, and peaks in BED format, alongside a length file. It not only outputs tables linking peaks with genes and putative cis-regulatory regions but also includes a detailed summary of how each putative cis-regulatory region was defined, offering insights into the genomic landscape of regulation.

Usage:
```bash
python cisreg_map.py -t your_tss_file.bed -g your_genes_file.bed -p your_peaks_file.bed -l your_length_file.txt -o your_output_file_prefix
```
## Contact
For more information or support, please contact Maria Rossello at [mariarossello@ub.edu](mailto:mariarossello@ub.edu).
