# PHCP
PHCP is a modified version of ChromoPainter that takes pseudo-haplotypes as input and hence allows analysis of ancient DNA.  

PHCP "paints" a pseudo-haplotype as an imperfect mosaic of phased "donor"  haplotypes and calculates the expected length that the pseudo-haplotype copies from  each pair of donor haplotypes.   

For detailed information, please see …  

## Installation    
PHCP was implemented in R. Hence, to use it you need installed R and the R package "dplyr".

## Usage
Rscript phcp.R --target <pseudo_haplotype_ID> --modern <prefix_of_phased_data>  
--ancient <prefix_of_pseudo_haplotype_data> --chr <chromosome_number> --out <output_suffix>

To reduce running time, PHCP runs on each chromosome separately, so analyses of different chromosomes can be run in parallel.

## Inputs
PHCP requires the target pseudo-haplotypes, as well as phased haplotypes that will serve as donors.
Pseudo-haplotypes should be provided in two files:

1. Genomic data. The file requires the following columns (in this order): chromosome number, base pair (BP), allele 1, allele 2, SNP ID, genetic distance in Morgans, and then one column for each pseudo-haplotype. Pseudo-haplotypes columns should include the values 0,1 and NA for missing reads. The file can include pseudo-haplotypes that will not be analyzed.
For example, this is the representation of two SNPs and three individuals:  
1 752566 G A rs3094315 0.0201296 1 0 1  
1 842013 T G rs7419119 0.0225176 0 NA 1
The file name should be of the form <prefix_ancient>.ph  

2. Sample IDs and populations. The file should contain two columns, ID in the first column and population in the second column. The file should include all samples from the genomic data file, in the same order.
For example, here is a list of three samples from two different populations:
ID_1 popX
ID_2 popX
ID_3 popY
The file name should be of the form <prefix>.ids

Phased data should include two files:

1.	Genomic data in haps format (which is the output of several phasing packages such as SHAPEIT or EAGLE). The file does not have to include only the donors for the analysis and can include additional haplotypes. The file name should be of the form <prefix_modern>.haps

2.	Sample IDs and populations. The file should contain two columns, ID in the first column and population in the second column (as in the pseudo-haplotypes files). The file should include all samples from the haps file, in the same order. The file name should be of the form <prefix_modern>.ids  

Only SNPs that are identical in both modern and ancient data will be used in the analysis, i.e., SNPs need to have the same physical position (chromosome and BP) and the same alleles and allele order (such that 0 and 1 refer to the same alleles in the two datasets). SNPs that are not identical will be removed from the data.

## Required Parameters

--target &nbsp; &nbsp; &nbsp; ID of a single pseudo-haplotype to analyze  

--chr &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Chromosome number to analyze  

--ancient &nbsp; &nbsp; The prefixes of the pseudo-haplotypes files.  
&nbsp; &nbsp;  &nbsp; &nbsp; &nbsp; &nbsp;  &nbsp; &nbsp; &nbsp; &nbsp; If they have the same prefix, only one prefix can be provided.  
&nbsp; &nbsp;  &nbsp; &nbsp; &nbsp; &nbsp;  &nbsp; &nbsp; &nbsp; &nbsp; If they have different prefixes, two prefixes are required: the first for the genomic data and the second for the  
&nbsp; &nbsp;  &nbsp; &nbsp; &nbsp; &nbsp;  &nbsp; &nbsp; &nbsp; &nbsp; samples list. The genomic data file name should be in the form of <prefix_ancient>.ph and the name of the file  
&nbsp; &nbsp;  &nbsp; &nbsp; &nbsp; &nbsp;  &nbsp; &nbsp; &nbsp; &nbsp; with the IDs should be of the form <prefix_ancient>.ids.  

--modern &nbsp; &nbsp;The prefixes of the phased data files. If they have the same prefix, only one prefix can be provided.  
&nbsp; &nbsp;  &nbsp; &nbsp; &nbsp; &nbsp;  &nbsp; &nbsp; &nbsp; &nbsp; If they have different prefixes, two prefixes are required: the first for the genomic data and the second for the  
&nbsp; &nbsp;  &nbsp; &nbsp; &nbsp; &nbsp;  &nbsp; &nbsp; &nbsp; &nbsp; samples list. The genomic data file name should be in the form of <prefix_modern>.haps and the name of the  
&nbsp; &nbsp;  &nbsp; &nbsp; &nbsp; &nbsp;  &nbsp; &nbsp; &nbsp; &nbsp; file with ids list should be of the form <prefix_modern>.ids  

--out &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Output suffix. The output files names would be of the form: chrom<chr>_target_length_<out>  
  
## Optional parameters
--donors &nbsp; &nbsp; A file with list of donor populations, each population in a separate line.  
&nbsp; &nbsp;  &nbsp; &nbsp; &nbsp; &nbsp;  &nbsp; &nbsp; &nbsp; &nbsp;Default: all samples in <prefix_modern>.haps.  

--theta &nbsp; &nbsp; &nbsp; &nbsp;The mutation rate for ChromoPainter. Default: 0.00138488.  

--N &nbsp; &nbsp;  &nbsp; &nbsp;  &nbsp; &nbsp; &nbsp;The Ne parameter of ChromoPainter (related to the effective population size). Default: 64.5698.

## Output
The file chrom<chr>_target_length_<out> will contain the expected length of genetic material (in Morgan) copied from each ordered pair of donor haplotypes. It is in a matrix format, where the value in the row i and column j (x<sub>i,j</sub>) is the length of the genetic material copied from haplotypes i and j. The matrix is symmetric so x<sub>i,j</sub>=x<sub>i,j</sub>.  

## Example
For pseudo-haplotype sample called ID_1, pseudo-haplotypes data file called ancient_data.ph with a list of pseudo-haplotypes in ancient_data.ids, and phased data file called phased_data.haps with samples list in phased_data.ids, the command line for analyzing chromosome 1 would be as follows:  

Rscript phcp.R --target ID_1 --modern phased_data --ancient ancient_data --chr 1 --out example

  
