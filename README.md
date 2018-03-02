# fishermp
FisherMP: Motif Prediction using Fisher's exact test with openMP Parallel design

Option1: You can downlad the zip file and unzip it in your Linux computer. 
Option2: You can directly use to command "g++ -fopenmp fishermp.cpp -o fishermp" to compile the cpp file "fishermp.cpp".

For example, if given a fasta file, you can run the program with the following commands:
./fishermp input_file.fasta   > output.motifs

******
USAGE:
******
fishermp <dataset> [OPTIONS]  > OutputFile

<dataset>	file containing DNA sequences in FASTA format

OPTIONS:
-b		a background data file in FASTA format(default=directly produced by the program itself)
-m		minimum size of binding sites to find(default=5)
-M		Maximum size of binding sites to find(default=10)
-n		Number of motifs to find(default=10)
-t		number of threads to call(default=6 i.e. =M-m+1)
