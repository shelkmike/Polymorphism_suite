This is a set of two scripts for the calculation of <b>polymorphism</b> from <b>multiple sequence alignments</b>. You give a script a multiple alignment of DNA, RNA or proteins as input, and it outputs various information about polymorphism. For details on the different metrics of polymorphism that the scripts calculate, see Table 1 in the article "The plastid genome of the non-photosynthetic plant <i>Rhopalocnemis phalloides</i> is one of the most polymorphic genomes known" (https://www.biorxiv.org/content/10.1101/2025.10.06.680833v1.full.pdf).<br>

The two scripts are as follows:<br>
1) <i>calculate_polymorphism_using_entire_alignment.py</i> — calculates polymorphism for the <b>entire</b> alignment.<br>
2) <i>calculate_polymorphism_in_windows.py</i> — splits the input alignment into <b>windows</b> (for example, 100 base pairs with a step of 1 base pair), then calculates the polymorphism for each window, and outputs all the results into a table.<br><br>

## Installation
1) Install the Python libraries <i>natsort</i> and <i>dendropy</i>. This can be done, for example, with the following command:<br>
`pip3 install --upgrade --user natsort dendropy`<br>
2) Download the archive with the latest release from https://github.com/shelkmike/Polymorphism_suite/releases and unpack it.<br><br>
That's it, nothing else needs to be done.<br>

## How to run the scripts
Both scripts use positional arguments.<br>
### 1) calculate_polymorphism_using_entire_alignment.py
This script has <b>3 positional arguments</b>:<br>
1) Path to the input multiple alignment in FASTA format.<br>
2) Number of threads.<br>
3) Path to the output folder. If it does not exist, the script will create it.<br><br>
Example:<br>
`python3 calculate_polymorphism_using_entire_alignment.py some_alignment.fasta 10 Output_folder`<br>

### 2) calculate_polymorphism_in_windows.py
This script has <b>7 positional arguments</b>:<br>
1) Path to the input multiple alignment in FASTA format.<br>
2) Window size.<br>
3) Window step.<br>
4) Circular or non-circular sequence. In the case of circular (e.g., alignment of an entire mitochondrial genome), the script will also consider windows that cross the beginning and end of the alignment. In the case of non-circular, the edge windows will be used truncated. The values should be "circular" or "non-circular" respectively.<br>
5) Should the center of the first window be at position 1 or at position (window_length / 2). In the second case, the first window will not cross the edge of the genome (if the previous option is "circular") and will not be shortened (if the previous option is "non-circular"). The value must be "start_from_1" or "start_from_window_center".<br>
6) Number of threads.<br>
7) Path to the output folder. If it does not exist, the script will create it. <br><br>
Example:<br>
`python3 calculate_polymorphism_in_windows.py some_alignment.fasta 100 1 circular start_from_1 10 Output_folder` <br><br>

## Output format
### 1) calculate_polymorphism_using_entire_alignment.py
The output folder will contain <b>CSV files</b> with the following tables:<br>
1) percent_identity_table.csv — A table in which percent identities for all pairs of sequences are recorded.<br>
2) MLD_table.csv — A table in which maximum likelihood phylogenetic distances for all pairs of sequences are recorded.<br>
3) p-distance_table.csv — A table in which p-distances for all pairs of sequences are recorded.<br><br>
Also, the output folder will have a file metrics_of_polymorphism.txt, which contains values of MPPI, MPMLD, and pi, that are essentially the arithmetic mean values for the three tables mentioned above, respectively.<br>

### 2) calculate_polymorphism_in_windows.py
The output folder will contain a <b>single CSV-file</b> in which the MPPI, MPMLD, and pi values are recorded for each position of the genome (the position is counted from 1, not 0) for the window whose center is at that position. If the window length is even, the rightmost of the two central nucleotides is considered the center of the window.<br><br>
## Questions and answers:
1) Can these scripts perform polymorphism analysis based on pairwise alignment rather than multiple alignment?<br>
Yes.<br>

2) If there is a capital letter in one sequence and the same lowercase letter in another sequence in the same column, will this be considered a mismatch?<br>
No.<br>

3) The version of IQ-TREE that comes with the package does not work. What should I do?<br>
These scripts come with IQ-TREE, which is used to calculate phylogenetic distances and, accordingly, MPMLD. If the file attached to these scripts is not working, you can compile IQ-TREE1 (not IQ-TREE2 or IQ-TREE3) yourself and replace the executable file in the Additional folder. IQ-TREE1 can be downloaded from there: https://iqtree.github.io/release/v1.6.12<br>

4) If there is a gap in two sequences in a given column in a multiple alignment, will this be considered a match or a mismatch in a pairwise comparison?<br>
In pairwise comparison, columns in which both sequences have gaps are not taken into account at all. That is, for example, if the alignment length is 1000 nucleotides, and both sequences have gaps in 10 positions, then the pairwise comparison will be made only on 990 positions.<br>

5) Where can I read about Polymorphism suite in more detail?<br>
See the section "Methods of polymorphism calculation" in the article "The plastid genome of the non-photosynthetic plant <i>Rhopalocnemis phalloides</i> is one of the most polymorphic genomes known" (https://www.biorxiv.org/content/10.1101/2025.10.06.680833v1.full.pdf).
