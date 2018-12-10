~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
https://github.com/shaoqiangzhang/fishermp

vesion: 1.1 (Dec 10, 2018) 
 
If FisherMP is having problems, call Shaoqiang Zhang( zhangshaoqiang@tjnu.edu.cn)
 
===============================================
Installation:

Before compiling the program, please first ensure that 
your computer has the openMP library (version>=4.0) 
(http://openmp.org) and GNU compiler (GCC >5.0). 

Please run "./GNUcompile" to compile the program. 
or type “g++ -fopenmp fishermp_current.cpp -o fishermp” in the terminal to compile the proram.

If your computer is Linux based on X86_64, 
you can directly run the executable file "./fishermp" in the terminal.

You can can also copy the execuable file "fishermp" to the "$USER/bin" directory. 

====================================================

USAGE:

./fishermp <dataset> [optional arguments]  > Output_File

<dataset>				file containing DNA sequences in FASTA format
[-b <file>]				a background data file in FASTA format (default=directly produced by the program itself)
[-m <min_size>]			minimum motif length (default=5)
[-M <max_size>]			Maximum motif length (default=10)
[-c <cooperative_size>]	define the number of cooperate motifs(default=2, you can select 1,2,3 for single motif, pair motifs, or triple motifs)
[-n <num>]				Number of motifs to find(default=all possible motifs)
[-t <threads>]			number of threads to call(default=6 i.e. =M-m+1)

=======================================================
Example:

If you have two fasta files (one is the foreground file, the other is the background file),
and you want to call 10 threads to compute, the motif lengths are set from 5 to 9, the command is as follows:

./fishermp foregroundfile.fasta -b backgroundfile.fasta -m 5 -M 9 -t 10  > output.file

