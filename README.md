The data and genotypes made available were obtained from the project:
Gough mouse project
Gray MM, Parmenter MD, Hogan CA, Ford I, Cuthbert RJ, Ryan PG, Broman
KW, Payseur BA (2015) Genetics of rapid and extreme size evolution in
island mice. Genetics 201:213-228

General description:
This Python script offers functionalities for reading, processing and analyzing genetic data.
It includes several functions for data quality control, generating different types of matrices and formatting the data for use with the BLUPF90 software focused on implementing non-additive genetic effects.
For this program, data and matrices must be provided with a renumbering from 1 to n (number of animals) in such a way that they are co-formable and compatible.
This script, in addition to renumbering and ordering operations, creates the various relationship matrices based on markers and vectorizes them to use them in BLUPF90 programs.

Key Features:
Quality control of genetic data, including call rate test, minimum allele frequency (MAF), detection of monomorphic markers and Hardy-Weinberg equilibrium test.
Generation of relationship matrices based on SNP markers, including additive, dominance, additive-additive, dominance-dominance and additive-dominance matrices.
Data formatting for use with the BLUPF90 software, forming two columns.

Instructions for use:
Users can run the script by providing the paths to the geno and data files as command line arguments. They can also choose between different types of genetic matrices to be generated and specify the name of the output file.
It is not necessary to change the code, just follow the instructions.

Dependencies:
This script requires the pandas and numpy libraries to run correctly. Make sure you have them installed before use.

Usage Examples:
To run the script, you can follow the following steps:

Open a terminal or command prompt.
Navigate to the directory where the script is located.
Run the python command script_name.py
load the geno.csv data.csv files
If you want the ordered file to be saved, then give it a name for output_file, otherwise leave it blank.

Additional Notes:
This script was developed for processing specific genetic data and requires that the genotypes are in tabular format and both files contain an ID column in this way
If you have any questions or encounter any problems using the script, please feel free to contact me.

  mateus.g.santos@ufv.br
  I'll be happy if you use and reference it.

