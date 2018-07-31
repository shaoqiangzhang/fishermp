# fishermp
FisherMP: Motif Prediction using Fisher's exact test with openMP Parallel design

Option1: You can downlad the zip file and unzip it in your Linux computer. 

Option2: You can directly use to command "g++ -fopenmp fishermp.cpp -o fishermp" to compile the cpp file "fishermp.cpp".


If your computer is Linux based on X86_64, you can directly run the "./fishermp" in the terminal.

*******************************************
                  USAGE:
*******************************************
./fishermp <dataset> [OPTIONS]  > OutputFile

<dataset>	file containing DNA sequences in FASTA format

OPTIONS:
-b		a background data file in FASTA format(default=directly produced by the program itself)
-m		minimum size of binding sites to find(default=5)
-M		Maximum size of binding sites to find(default=10)
-n		Number of motifs to find(default=10)
-t		number of threads to call(default=6 i.e. =M-m+1)

*******************************************

You can run the example to test the prgram: 

./fishermp Klf1.fna -b Klf1.negative.seq -m 5 -M 9 -n 5 -t 10  > Klf1.motifs

or use the default settings:

./fishermp Klf1.fna > Klf1.motifs
