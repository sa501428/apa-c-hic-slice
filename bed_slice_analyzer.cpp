#include "bedpe_builder.h"
#include "hic_slice_reader.h"
#include <iostream>
#include <map>

// Structure to hold a chromosome pair key (reusing from previous code)
struct ChromPair {
    std::string chr1, chr2;
    
    bool operator<(const ChromPair& other) const {
        if (chr1 != other.chr1) return chr1 < other.chr1;
        return chr2 < other.chr2;
    }
};

// Function to check if a bin overlaps with a region
bool isOverlapping(int32_t bin, int32_t binStart, int32_t binEnd) {
    return !(bin >= binEnd || bin < binStart);
}

void printUsage(const char* progname) {
    std::cout << "Usage: " << progname << " [-both-intra-inter|-only-inter] "
              << "<forward.bed> <reverse.bed> <min_genome_dist> <max_genome_dist> "
              << "<hicslice_file>\n";
}

int main(int argc, char* argv[]) {
    try {
        if (argc < 6) {
            printUsage(argv[0]);
            return 1;
        }

        int argOffset = 1;
        bool makeIntra = true;
        bool makeInter = false;
        
        if (std::string(argv[1]) == "-both-intra-inter") {
            makeIntra = true;
            makeInter = true;
            argOffset++;
        } else if (std::string(argv[1]) == "-only-inter") {
            makeIntra = false;
            makeInter = true;
            argOffset++;
        }

        // Parse command line arguments
        std::string forward_bed = argv[argOffset];
        std::string reverse_bed = argv[argOffset + 1];
        long min_dist = std::stol(argv[argOffset + 2]);
        long max_dist = std::stol(argv[argOffset + 3]);
        std::string slice_file = argv[argOffset + 4];

        // Initialize HiC slice reader to get resolution
        HicSliceReader reader(slice_file);
        int32_t resolution = reader.getResolution();

        // Create BEDPE builder and generate regions
        BedpeBuilder builder(forward_bed, reverse_bed, min_dist, max_dist, makeIntra, makeInter);
        
        // Create map to store binned regions
        std::map<ChromPair, std::vector<std::pair<int32_t, int32_t>>> bedpeRegions;

        // Process each BEDPE entry and convert to bin coordinates
        auto processBedpeEntry = [&](const BedpeEntry& entry) {
            ChromPair pair{entry.chrom1, entry.chrom2};
            
            // Convert to bin coordinates
            int32_t binStart1 = entry.start1 / resolution;
            int32_t binEnd1 = (entry.end1 / resolution) + 1;
            int32_t binStart2 = entry.start2 / resolution;
            int32_t binEnd2 = (entry.end2 / resolution) + 1;
            
            bedpeRegions[pair].push_back({binStart1, binEnd1});
            bedpeRegions[pair].push_back({binStart2, binEnd2});
        };

        // Generate and process BEDPE entries
        builder.generateBedpe(processBedpeEntry);

        float totalCount = 0;

        // Structure to read from slice file
        struct {
            int16_t chr1Key;
            int32_t binX;
            int16_t chr2Key;
            int32_t binY;
            float value;
        } record;

        // Process each record in the slice file
        while (gzread(reader.getFile(), &record, sizeof(record)) == sizeof(record)) {
            if (std::isnan(record.value) || std::isinf(record.value) || record.value <= 0) {
                continue;
            }

            std::string chr1 = reader.getChromosomeFromKey(record.chr1Key);
            std::string chr2 = reader.getChromosomeFromKey(record.chr2Key);
            
            ChromPair pair{chr1, chr2};
            auto it = bedpeRegions.find(pair);
            
            if (it != bedpeRegions.end()) {
                const auto& regions = it->second;
                for (size_t i = 0; i < regions.size(); i += 2) {
                    if (isOverlapping(record.binX, regions[i].first, regions[i].second) &&
                        isOverlapping(record.binY, regions[i+1].first, regions[i+1].second)) {
                        totalCount += record.value;
                        break;
                    }
                }
            }
        }

        std::cout << "Total contact count in regions: " << totalCount << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
} 