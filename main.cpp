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
              << "<hic_slice_file> [<forward.bed> <reverse.bed> <output.txt>]...\n"
              << "\tCreate potential loop locations using the anchors\n"
              << "\t\t'inter' for inter-chromosomal features\n"
              << "\t\t'intra' for intra-chromosomal features\n"
              << "\t\t<min_genome_dist> minimum genomic distance for loops\n"
              << "\t\t<max_genome_dist> maximum genomic distance for loops\n"
              << "\t\t<window_size> window size around loop\n"
              << "\t\t<hic_slice_file> path to the HiC slice file\n"
              << "\t\t<forward.bed> <reverse.bed> <output.txt> triplets (can have multiple)\n";
}

bool fileExists(const std::string& filename) {
    std::ifstream f(filename);
    return f.good();
}

struct BedpeSet {
    std::string forward_bed;
    std::string reverse_bed;
    std::string output_file;
};

int main(int argc, char* argv[]) {
    try {
        // Need at least 6 args (base args) + 3 (one bedpe set)
        if (argc < 9 || (argc - 6) % 3 != 0) {
            printUsage();
            return 1;
        }

        // Parse base arguments
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

        // Validate parameters
        if (min_dist < 0 || max_dist < min_dist) {
            throw std::runtime_error("Invalid distance parameters");
        }
        if (window_size <= 0) {
            throw std::runtime_error("Window size must be positive");
        }
        const int MAX_WINDOW_SIZE = 1000;
        if (window_size > MAX_WINDOW_SIZE) {
            throw std::runtime_error("Window size too large (max: " + 
                                   std::to_string(MAX_WINDOW_SIZE) + ")");
        }
        if (!fileExists(slice_file)) {
            throw std::runtime_error("Slice file not found: " + slice_file);
        }

        // Parse BEDPE sets
        std::vector<BedpeSet> bedpe_sets;
        for (int i = 6; i < argc; i += 3) {
            BedpeSet set = {
                argv[i],      // forward bed
                argv[i + 1],  // reverse bed
                argv[i + 2]   // output file
            };
            
            if (!fileExists(set.forward_bed)) {
                throw std::runtime_error("Forward BED file not found: " + set.forward_bed);
            }
            if (!fileExists(set.reverse_bed)) {
                throw std::runtime_error("Reverse BED file not found: " + set.reverse_bed);
            }
            
            bedpe_sets.push_back(set);
        }

        // Process each BEDPE set to generate entries
        std::cout << "Processing " << bedpe_sets.size() << " BEDPE sets..." << std::endl;
        std::vector<std::vector<BedpeEntry>> all_bedpe_entries(bedpe_sets.size());
        
        for (size_t i = 0; i < bedpe_sets.size(); i++) {
            const auto& set = bedpe_sets[i];
            std::cout << "Loading BED files: " << set.forward_bed << " and " << set.reverse_bed << std::endl;
            BedpeBuilder builder(set.forward_bed, set.reverse_bed, min_dist, max_dist, isInter);
            all_bedpe_entries[i] = builder.buildBedpe();
        }
        
        std::cout << "Processing slice file: " << slice_file << std::endl;
        auto matrices = processSliceFile(slice_file, all_bedpe_entries, 
                                       window_size, isInter, min_dist, max_dist);

        // Save all matrices
        for (size_t i = 0; i < matrices.size(); i++) {
            std::cout << "Saving matrix to: " << bedpe_sets[i].output_file << std::endl;
            matrices[i].save(bedpe_sets[i].output_file);
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
} 