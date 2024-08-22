**Introduction:**

 Following the paper "Succint de Bruijn Graphs" by Bowe et al. (2012), we implement a way to represent de Bruijn graphs of k-mers in a DNA sequence with $m$ edges using $4m+o(m)$ bits, and show how this space can be compressed even further. Using this data structure we can find the outdegree and indegree of a node and the outgoing edge with a given label in constant time. The incoming edge with a given label is found in $\mathcal{O}(k\log\sigma)$ time, while the label of a node and membership queries can instead be answered in $\mathcal{O}(k)$ time.


**Usage:**

The file `boss.cpp` contains the implementation of the class that represents the BOSS data structure, while the file `main.cpp` is a simple example of how to use it. It takes as arguments the path of the input file and the value of k, then simply builds the data structure and prints all the information. It should be compiled with

```g++ -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib 'main.cpp' -o main -lsdsl -ldivsufsort -ldivsufsort64```

following the instructions from the [SDSL library](https://github.com/simongog/sdsl-lite) page. Then, we can execute it using the example file `kmers.txt`:

```./main kmers.txt 3```

Any other input file can be used, but it must follow the specifications described in the [report](Implementation_report.pdf). In the case of FASTQ files, 'fastq_preprocess.cpp' already takes care of the preprocessing, selecting only the lines with the DNA sequences. It takes as arguments the path to the FASTQ file, the number of lines to include and the value of k.

The file ```memory usage.cpp``` is the one we used to measure the performances (the results are summarized in the [report](Implementation_report.pdf)). Finally, the file ```timer.cpp``` contains a scope based timer by used to show visually in Chrome tracing the timing of the program.


 
