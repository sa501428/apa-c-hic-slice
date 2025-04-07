#include "bedpe_builder.h"
#include "apa.h"
#include <iostream>
#include <string>

// APA4 Aggregate Peak Analysis
// first generate a bedpe file of all potential loop locations from bed files
// then iterates through the counts of the hic slice file
// add up the total counts for all reads that overlap a loop
// output the total counts
void printUsage() {
    std::cout << "Usage: apa4 [-both-intra-inter|-only-inter] <forward.bed> <reverse.bed> "
              << "<min_genome_dist> <max_genome_dist> <hic_slice_file> <output.txt>\n"
              << "\tCreate potential loop locations using the anchors\n"
              << "\t\tDefault is intra-chromosomal only features\n"
              << "\t\tUse -only-inter for inter-chromosomal features only\n"
              << "\t\tUse -both-intra-inter for both types of features\n"
              << "\t\t<hic_slice_file> is the path to the HiC slice file\n"
              << "\t\t<output.txt> is the path to the output file\n";
}

int main(int argc, char* argv[]) {
    try {
        int argOffset = 1;
        bool makeIntra = true;
        bool makeInter = false;
        
        if (argc < 7) {  // Need one more arg for output file
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
        
        if (min_dist < 0 || max_dist < min_dist) {
            throw std::runtime_error("Invalid distance parameters: min_dist must be >= 0 and max_dist must be >= min_dist");
        }
        
        std::string slice_file = argv[argOffset + 4];
        std::string output_file = argv[argOffset + 5];
        
        // Generate BEDPE entries
        BedpeBuilder builder(forward_bed, reverse_bed, min_dist, max_dist, makeIntra, makeInter);
        auto bedpe_entries = builder.buildBedpe();

        float totalCount = processSliceFile(slice_file, bedpe_entries);
        std::cout << "Total contact count in BEDPE regions: " << totalCount << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
} 