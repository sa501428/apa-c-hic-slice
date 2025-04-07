#include "bedpe_builder.h"
#include <iostream>
#include <zlib.h>
#include <stdexcept>
#include <cstring>
#include <cmath>
#include <map>
#include <string>
#include <vector>

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
        std::string slice_file = argv[argOffset + 4];
        std::string output_file = argv[argOffset + 5];
        
        // Generate BEDPE entries
        BedpeBuilder builder(forward_bed, reverse_bed, min_dist, max_dist, makeIntra, makeInter);
        auto bedpe_entries = builder.buildBedpe("");  // Empty string since we don't write to file

        // Open slice file
        gzFile file = gzopen(slice_file.c_str(), "rb");
        if (!file) {
            throw std::runtime_error("Could not open file: " + slice_file);
        }

        // Read and verify magic string
        char magic[8];
        if (gzread(file, magic, 8) != 8 || strncmp(magic, "HICSLICE", 8) != 0) {
            throw std::runtime_error("Invalid file format: missing magic string");
        }

        // Read resolution
        int32_t resolution;
        if (gzread(file, &resolution, sizeof(int32_t)) != sizeof(int32_t)) {
            throw std::runtime_error("Failed to read resolution");
        }

        // Read chromosome mapping
        int32_t numChromosomes;
        if (gzread(file, &numChromosomes, sizeof(int32_t)) != sizeof(int32_t)) {
            throw std::runtime_error("Failed to read chromosome count");
        }

        std::map<int16_t, std::string> chromosomeKeyToName;
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

        // Create map for faster BEDPE lookup
        std::map<std::pair<std::string, std::string>, std::vector<std::pair<std::pair<int32_t, int32_t>, std::pair<int32_t, int32_t>>>> regionMap;
        for (const auto& entry : bedpe_entries) {
            int32_t binStart1 = entry.start1 / resolution;
            int32_t binEnd1 = (entry.end1 + resolution - 1) / resolution;
            int32_t binStart2 = entry.start2 / resolution;
            int32_t binEnd2 = (entry.end2 + resolution - 1) / resolution;
            regionMap[{entry.chrom1, entry.chrom2}].push_back({{binStart1, binEnd1}, {binStart2, binEnd2}});
        }

        // Process contact records
        float totalCount = 0;
        struct {
            int16_t chr1Key;
            int32_t binX;
            int16_t chr2Key;
            int32_t binY;
            float value;
        } record;

        while (gzread(file, &record, sizeof(record)) == sizeof(record)) {
            if (std::isnan(record.value) || std::isinf(record.value) || record.value <= 0) {
                continue;
            }

            std::string chr1 = chromosomeKeyToName[record.chr1Key];
            std::string chr2 = chromosomeKeyToName[record.chr2Key];
            
            auto it = regionMap.find({chr1, chr2});
            if (it != regionMap.end()) {
                for (const auto& region : it->second) {
                    if (record.binX >= region.first.first && record.binX < region.first.second &&
                        record.binY >= region.second.first && record.binY < region.second.second) {
                        totalCount += record.value;
                        break;  // Count each contact only once
                    }
                }
            }
        }

        gzclose(file);
        std::cout << "Total contact count in BEDPE regions: " << totalCount << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
} 