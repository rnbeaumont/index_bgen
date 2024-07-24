# index_bgen
Program to index and extract SNPs from bgen files

Usage:

To efficiently extract SNPs, index files for all bgen files will first need to be created.

To index files
./index_bgen -bgen <bgen_file> -index <index_file_name>

To extract SNPs
./index_bgen -bgen <bgen_file> -index <index_file_name> -sample <sample_file> -snps <SNPs_file> -chr <chr> -out <output_file>
  
  where <SNPS_file> is a tab delimited file containing SNPs to extract, one per line. Fields should be chr, pos, and 2 alleles,
  < chr > is chromosome number, <output_file> is the desired filename of the output file. In dosage format (default) the output is the   dosage of the  second allele listed in the output file (not necessarily the same order as in the <SNPs_file>)
  
optional arguments are:
-gen - output will be in gen format rather than SNP dosages
-rsid - SNPs will be named with the contents of the rsid field in the bgen file
