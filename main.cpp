#include "bedpe_builder.h"
#include "apa.h"
#include <iostream>
#include <string>
#include <fstream>
#include <stdexcept>
#include <cstdint>

// APA4 Aggregate Peak Analysis
// first generate a bedpe file of all potential loop locations from bed files
// then iterates through the counts of the hic slice file
// add up the total counts for all reads that overlap a loop
// output the total counts
void printUsage() {
    std::cout << "Usage: apa4 [-only-inter] <forward.bed> <reverse.bed> "
              << "<min_genome_dist> <max_genome_dist> <hic_slice_file> <output.txt> [window_size]\n"
              << "\tCreate potential loop locations using the anchors\n"
              << "\t\tDefault is intra-chromosomal features\n"
              << "\t\tUse -only-inter for inter-chromosomal features\n"
              << "\t\t<hic_slice_file> is the path to the HiC slice file\n"
              << "\t\t<output.txt> is the path to the output file\n"
              << "\t\t[window_size] optional window size around loop (default: 10)\n";
}

bool fileExists(const std::string& filename) {
    std::ifstream f(filename);
    return f.good();
}

int main(int argc, char* argv[]) {
    try {
        int argOffset = 1;
        bool isInter = false;
        int window_size = 10;  // default window size
        
        if (argc < 6) {
            printUsage();
            return 1;
        }
        
        if (std::string(argv[1]) == "-only-inter") {
            isInter = true;
            argOffset++;
        }
        
        if (argc < (6 + argOffset)) {  // Check we have enough args after option
            printUsage();
            return 1;
        }
        
        std::string forward_bed = argv[argOffset];
        std::string reverse_bed = argv[argOffset + 1];
        
        // Check input files exist
        if (!fileExists(forward_bed)) {
            throw std::runtime_error("Forward BED file not found: " + forward_bed);
        }
        if (!fileExists(reverse_bed)) {
            throw std::runtime_error("Reverse BED file not found: " + reverse_bed);
        }
        
        long min_dist = std::stol(argv[argOffset + 2]);
        long max_dist = std::stol(argv[argOffset + 3]);
        
        if (min_dist < 0 || max_dist < min_dist) {
            throw std::runtime_error("Invalid distance parameters: min_dist must be >= 0 and max_dist must be >= min_dist");
        }
        
        std::string slice_file = argv[argOffset + 4];
        if (!fileExists(slice_file)) {
            throw std::runtime_error("Slice file not found: " + slice_file);
        }
        
        std::string output_file = argv[argOffset + 5];
        
        // Optional window size parameter
        if (argc > (6 + argOffset)) {
            window_size = std::stoi(argv[argOffset + 6]);
            if (window_size <= 0) {
                throw std::runtime_error("Window size must be positive");
            }
        }
        
        // Generate BEDPE entries
        BedpeBuilder builder(forward_bed, reverse_bed, min_dist, max_dist, isInter);
        auto bedpe_entries = builder.buildBedpe();

        APAMatrix apaMatrix = processSliceFile(slice_file, bedpe_entries, output_file, 
                                             window_size, isInter, min_dist, max_dist);
        std::cout << "APA matrix saved to: " << output_file << std::endl;
        std::cout << "Center pixel value: " << apaMatrix.matrix[apaMatrix.width/2][apaMatrix.width/2] << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
} 