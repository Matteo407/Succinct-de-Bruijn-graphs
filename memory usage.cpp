#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cmath>

#include "boss.cpp"

using namespace std::chrono;

void process_fastq(int k, int lines = 2e6) {
    std::cout << "starting preprocessing with lines = " << lines << '\n';
    std::ifstream input_file1("input.fastq");
    std::ofstream output_file("temp.txt");

    std::string line;
    int line_number = 0;

    while (std::getline(input_file1, line)) {
        if (line_number % 4 == 1 && line_number < lines) {
            reverse(line.begin(), line.end());
            output_file << line + '$';
        }

        line_number++;
    }

    output_file << std::string(k-1, '$');

    input_file1.close();
    output_file.close();

    std::cout << "finished preprocessing\n";
}

int main(int argc, char *argv[]) {
    int k = 31;

    std::ofstream data_file("measures.txt");
    data_file << "m W L F tot s" << std::endl;

    for (int i = 100; i < 2500000; i+= std::pow(10, std::floor(log10(i)))) { 

        process_fastq(k, i);

        auto start = high_resolution_clock::now();
        BOSS* boss = new BOSS("temp.txt", k);
        auto stop = high_resolution_clock::now();

        std::cout << "done construction\n";

        auto duration = duration_cast<microseconds>(stop - start);

        std::cout << i << " " << boss->m << " " << duration.count()/1e6 << std::endl << std::endl;
        data_file << boss->m << " " << sdsl::size_in_bytes(boss->W)*8 << " " << sdsl::size_in_bytes(boss->L)*8 << " " << sdsl::size_in_bytes(boss->F)*8 << " " << boss->memory_usage_in_bits() << " " << duration.count() << std::endl;

        // clear memory
        delete boss;
    }

    data_file.close();
}
