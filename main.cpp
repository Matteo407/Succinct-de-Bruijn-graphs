#include <iostream>
#include "boss.cpp"


int main(int argc, char *argv[]) {
    if (argc != 3) {
        std::cout << "Usage: " << argv[0] << " <input file> <k>" << std::endl;
        return 1;
    }

    // Read the input 
    std::string path_to_input_file = argv[1];
    int k = atoi(argv[2]);                     // Lenght of k-mers

    // BOSS structure
    BOSS* boss = new BOSS(path_to_input_file, k);

    // print the BOSS index
    boss->print_graph(true);
    
    // clear memory
    delete boss;
}
