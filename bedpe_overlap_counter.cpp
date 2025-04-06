#include "hic_slice_reader.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <string>

// Structure to hold BEDPE regions (in bin coordinates)
struct BedpeRegion {
    int32_t binStart1, binEnd1;
    int32_t binStart2, binEnd2;
};

// Structure to hold a chromosome pair key
struct ChromPair {
    std::string chr1, chr2;
    
    bool operator<(const ChromPair& other) const {
        if (chr1 != other.chr1) return chr1 < other.chr1;
        return chr2 < other.chr2;
    }
};

// Function to check if a bin overlaps with a region (now in bin coordinates)
bool isOverlapping(int32_t bin, int32_t binStart, int32_t binEnd) {
    return !(bin >= binEnd || bin < binStart);
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <hicslice_file> <bedpe_file>" << std::endl;
        return 1;
    }

    try {
        // First get resolution from HiC file to convert BEDPE coordinates
        HicSliceReader reader(argv[1]);
        int32_t resolution = reader.getResolution();

        // Read BEDPE file into a map
        std::map<ChromPair, std::vector<BedpeRegion>> bedpeRegions;
        std::ifstream bedpeFile(argv[2]);
        std::string line;
        
        while (std::getline(bedpeFile, line)) {
            std::istringstream iss(line);
            std::string chr1, chr2;
            int32_t start1, end1, start2, end2;
            
            if (!(iss >> chr1 >> start1 >> end1 >> chr2 >> start2 >> end2)) {
                continue; // Skip malformed lines
            }
            
            // Convert genomic coordinates to bin coordinates
            int32_t binStart1 = start1 / resolution;
            int32_t binEnd1 = (end1 / resolution ) + 1; // Round up
            int32_t binStart2 = start2 / resolution;
            int32_t binEnd2 = (end2 / resolution) + 1; // Round up
            
            ChromPair pair{chr1, chr2};
            bedpeRegions[pair].push_back({binStart1, binEnd1, binStart2, binEnd2});
        }

        float totalCount = 0;

        // Structure to read from file
        struct {
            int16_t chr1Key;
            int32_t binX;
            int16_t chr2Key;
            int32_t binY;
            float value;
        } record;

        // Process each record
        while (gzread(reader.getFile(), &record, sizeof(record)) == sizeof(record)) {
            if (std::isnan(record.value) || std::isinf(record.value) || record.value <= 0) {
                continue; // Skip invalid values
            }

            std::string chr1 = reader.getChromosomeFromKey(record.chr1Key);
            std::string chr2 = reader.getChromosomeFromKey(record.chr2Key);
            
            ChromPair pair{chr1, chr2};
            auto it = bedpeRegions.find(pair);
            
            if (it != bedpeRegions.end()) {
                // Check each BEDPE region for this chromosome pair
                for (const auto& region : it->second) {
                    if (isOverlapping(record.binX, region.binStart1, region.binEnd1) &&
                        isOverlapping(record.binY, region.binStart2, region.binEnd2)) {
                        totalCount += record.value;
                        break; // Count each contact only once even if it overlaps multiple regions
                    }
                }
            }
        }

        std::cout << "Total contact count in BEDPE regions: " << totalCount << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
} 