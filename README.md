**Introduction:**

 Following the paper "Succint de Bruijn Graphs" by Bowe et al. (2012), we implement a way to represent de Bruijn graphs of k-mers in a DNA sequence with $m$ edges using $4m+o(m)$ bits, and show how this space can be compressed even further. Using this data structure we can find the outdegree and indegree of a node and the outgoing edge with a given label in constant time. The incoming edge with a given label is found in $\mathcal{O}(k\log\sigma)$ time, while the label of a node and membership queries can instead be answered in $\mathcal{O}(k)$ time.


**Usage:**

boss.cpp contains the class
main.cpp is a simple example of how to use it
compile with g++ -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib 'main.cpp' -o main -lsdsl -ldivsufsort -ldivsufsort64 following sdsl
then ./main kmers.txt 3
memory usage is the file used to measure performance, the results are reported also in the report
timer is a scope based timer by used to benchmark the code


**Report:**

Append pdf [Some title here](FILE_NAME.pdf)
 
