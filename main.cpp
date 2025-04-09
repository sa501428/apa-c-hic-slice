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
    std::cout << "Usage: apa4 <inter|intra> <min_genome_dist> <max_genome_dist> <window_size> "
              << "<hic_slice_file> <forward.bed> <reverse.bed> <output.txt>\n"
              << "\tCreate potential loop locations using the anchors\n"
              << "\t\t'inter' for inter-chromosomal features\n"
              << "\t\t'intra' for intra-chromosomal features\n"
              << "\t\t<min_genome_dist> minimum genomic distance for loops\n"
              << "\t\t<max_genome_dist> maximum genomic distance for loops\n"
              << "\t\t<window_size> window size around loop\n"
              << "\t\t<hic_slice_file> path to the HiC slice file\n"
              << "\t\t<forward.bed> path to forward BED file\n"
              << "\t\t<reverse.bed> path to reverse BED file\n"
              << "\t\t<output.txt> path to the output file\n";
}

bool fileExists(const std::string& filename) {
    std::ifstream f(filename);
    return f.good();
}

int main(int argc, char* argv[]) {
    try {
        if (argc != 9) {
            printUsage();
            return 1;
        }

        // Parse arguments in new order
        std::string mode = argv[1];
        bool isInter;
        if (mode == "inter") {
            isInter = true;
        } else if (mode == "intra") {
            isInter = false;
        } else {
            throw std::runtime_error("First argument must be either 'inter' or 'intra'");
        }

        long min_dist = std::stol(argv[2]);
        long max_dist = std::stol(argv[3]);
        int window_size = std::stoi(argv[4]);
        std::string slice_file = argv[5];
        std::string forward_bed = argv[6];
        std::string reverse_bed = argv[7];
        std::string output_file = argv[8];

        // Validate parameters
        if (min_dist < 0 || max_dist < min_dist) {
            throw std::runtime_error("Invalid distance parameters: min_dist must be >= 0 and max_dist must be >= min_dist");
        }
        if (window_size <= 0) {
            throw std::runtime_error("Window size must be positive");
        }

        // Check input files exist
        if (!fileExists(forward_bed)) {
            throw std::runtime_error("Forward BED file not found: " + forward_bed);
        }
        if (!fileExists(reverse_bed)) {
            throw std::runtime_error("Reverse BED file not found: " + reverse_bed);
        }
        if (!fileExists(slice_file)) {
            throw std::runtime_error("Slice file not found: " + slice_file);
        }

        std::cout << "Loading BED files and generating BEDPE entries..." << std::endl;
        BedpeBuilder builder(forward_bed, reverse_bed, min_dist, max_dist, isInter);
        auto bedpe_entries = builder.buildBedpe();
        std::cout << "Generated " << bedpe_entries.size() << " BEDPE entries" << std::endl;

        std::cout << "Processing slice file: " << slice_file << std::endl;
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