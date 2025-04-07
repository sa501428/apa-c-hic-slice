#include "bedpe_builder.h"
#include "hic_slice_reader.h"
#include <iostream>
#include <zlib.h>
#include <stdexcept>
#include <cstring>
#include <cmath>
#include <map>


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

    std::vector<BedpeEntry> bedpe_entries;
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
        bedpe_entries = builder.buildBedpe(output_file);

        struct {
        int16_t chr1Key;
        int32_t binX;
        int16_t chr2Key;
        int32_t binY;
        float value;
        } compressedRecord;


        gzFile file = gzopen(filePath.c_str(), "rb");
        if (!file) {
            throw std::runtime_error("Could not open file: " + filePath);
        }

        // Read and verify magic string
        char magic[8];
        if (gzread(file, magic, 8) != 8 || strncmp(magic, "HICSLICE", 8) != 0) {
            throw std::runtime_error("Invalid file format: missing magic string");
        }

        int32_t resolution;
        int32_t numChromosomes; 

        // Read resolution
        if (gzread(file, &resolution, sizeof(int32_t)) != sizeof(int32_t)) {
            throw std::runtime_error("Failed to read resolution");
        }

        // Read number of chromosomes
        if (gzread(file, &numChromosomes, sizeof(int32_t)) != sizeof(int32_t)) {
            throw std::runtime_error("Failed to read chromosome count");
        }

        std::map<int16_t, std::string> chromosomeKeyToName;

        // Read chromosome mapping
        for (int i = 0; i < numChromosomes; i++) {
            int32_t nameLength;
            if (gzread(file, &nameLength, sizeof(int32_t)) != sizeof(int32_t)) {
                throw std::runtime_error("Failed to read chromosome name length");
            }

            std::vector<char> nameBuffer(nameLength + 1, 0);
            if (gzread(file, nameBuffer.data(), nameLength) != nameLength) {
                throw std::runtime_error("Failed to read chromosome name");
            }
            std::string chromosomeName(nameBuffer.data());

            int16_t key;
            if (gzread(file, &key, sizeof(int16_t)) != sizeof(int16_t)) {
                throw std::runtime_error("Failed to read chromosome key");
            }

            chromosomeKeyToName[key] = chromosomeName;
        }

        headerRead = true;

        while (gzread(file, &compressedRecord, sizeof(compressedRecord)) == sizeof(compressedRecord)) {
            SliceContactRecord record;
            record.chr1 = getChromosomeFromKey(compressedRecord.chr1Key);
            record.binX = compressedRecord.binX;
            record.chr2 = getChromosomeFromKey(compressedRecord.chr2Key);
            record.binY = compressedRecord.binY;
            record.value = compressedRecord.value;
        }

        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    // read the hic slice file
    HicSliceReader reader(hic_slice_file);
    int32_t resolution = reader.getResolution();

    // iterate through the bedpe entries
    for (const auto& bedpe_entry : bedpe_entries) {
    
} 