#include "bedpe_builder.h"
#include <iostream>

void printUsage() {
    std::cout << "Usage: bedpe_builder [-both-intra-inter|-only-inter] <forward.bed> <reverse.bed> "
              << "<min_genome_dist> <max_genome_dist> <output.bedpe>\n"
              << "\tCreate potential loop locations using the anchors\n"
              << "\t\tDefault is intra-chromosomal only features\n"
              << "\t\tUse -only-inter for inter-chromosomal features only\n"
              << "\t\tUse -both-intra-inter for both types of features\n";
}

int main(int argc, char* argv[]) {
    try {
        int argOffset = 1;
        bool makeIntra = true;
        bool makeInter = false;
        
        if (argc < 6) {
            printUsage();
            return 1;
        }
        
        if (std::string(argv[1]) == "-both-intra-inter") {
            makeIntra = true;
            makeInter = true;
            argOffset++;
        } else if (std::string(argv[1]) == "-only-inter") {
            makeIntra = false;
            makeInter = true;
            argOffset++;
        }
        
        std::string forward_bed = argv[argOffset];
        std::string reverse_bed = argv[argOffset + 1];
        long min_dist = std::stol(argv[argOffset + 2]);
        long max_dist = std::stol(argv[argOffset + 3]);
        std::string output_file = argv[argOffset + 4];
        
        BedpeBuilder builder(forward_bed, reverse_bed, min_dist, max_dist, makeIntra, makeInter);
        builder.buildBedpe(output_file);
        
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
} 