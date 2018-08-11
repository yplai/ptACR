# ptACR
ptACR: A statistical method to identify recombination in bacterial genomes based on SNP incompatibility

Source code is written in Python and the instructions are described as follows.

1. Version History

Vserion 1.0.0, initial release.


2. Requirements:

Python 2.7+ https://www.python.org

Numpy 1.2.1+ https://www.scipy.org/Download


3. Instructions:

By default, the number of permutation times is set as 10,000 and the cut-off value for empirical p-value is 0.05, adjusted by the Bonferroni correction. The parameters can be modified in the function of "find_bk()". The input file is a multiple sequence alignment of polymorphic sites (without gaps) in PHYLIP format, where the number of taxa and the length of the alignment is in the first line and the ID and the sequence of each taxon is shown in a line (see example file). One can analyze the entire alignment or a local region of interest by specifying the starting and ending positions, then, assign the sliding window size for identifying the recombination breakpoints. The ptACR allows up to five sliding window sizes at a time. 


python ptACR.py input_file starting_site ending_site sliding_window_size1 sliding_window_size2 ...


4. Data:

An example file is provided to test the execution of the script. The file (example.sites) consists of 19 taxa with 65244 polymorphic sites.

python ptACR.py example.sites 0 65243 125 250 500 ...


5. Copyright:

Authors: Yi-Pin Lai and Thomas R. Ioerger.

ptACR is a free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation (see http://www.gnu.org/licenses/).
