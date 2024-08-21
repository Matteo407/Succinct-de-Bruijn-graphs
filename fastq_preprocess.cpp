#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cmath>
#include <chrono>

using namespace std::chrono;

int main(int argc, char *argv[]) {
    if (argc != 3) {
        std::cout << "Usage: " << argv[0] << "<input_file_path> <num lines> <k>" << std::endl;
        return 1;
    }

    std::string path_to_input_file = argv[1];
    int lines = 4*atoi(argv[2]);
    int k = atoi(argv[3]);

    std::cout << "starting preprocessing with lines = " << lines << '\n';
    std::ifstream input_file(path_to_input_file);
    std::ofstream output_file("temp.txt");

    std::string line;
    int line_number = 0;

    while (std::getline(input_file, line)) {
        if (line_number % 4 == 1 && line_number < lines) {
            // std::cout << line_number << " " << line << '\n';

            reverse(line.begin(), line.end());
            output_file << line + '$';
        }

        line_number++;
    }

    output_file << std::string(k-1, '$');

    input_file.close();
    output_file.close();

    std::cout << "finished preprocessing\n";

    return 0;
}
