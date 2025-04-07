#ifndef APA_H
#define APA_H

#include "bedpe_builder.h"
#include <string>
#include <vector>
#include <set>
#include <unordered_map>

// Structure to hold binned regions for faster lookup
struct BinRegion {
    int32_t start;
    int32_t end;
};

// Structure to hold possible bin positions for quick filtering
struct RegionsOfInterest {
    std::unordered_map<std::string, std::set<int32_t>> rowIndices;  // chrom -> set of binX
    std::unordered_map<std::string, std::set<int32_t>> colIndices;  // chrom -> set of binY
    int32_t resolution;
    int32_t window;

    RegionsOfInterest(const std::vector<BedpeEntry>& bedpe_entries, 
                     int32_t res, int32_t win) 
        : resolution(res), window(win) {
        for (const auto& entry : bedpe_entries) {
            // Convert BEDPE coordinates to bin positions
            int32_t bin1Start = entry.start1 / resolution;
            int32_t bin1End = (entry.end1 / resolution) + 1;  // +1 to include the full range
            int32_t bin2Start = entry.start2 / resolution;
            int32_t bin2End = (entry.end2 / resolution) + 1;  // +1 to include the full range
            
            // Calculate center positions
            int32_t centerX = (bin1Start + bin1End) / 2;
            int32_t centerY = (bin2Start + bin2End) / 2;
            
            // Add all possible bins within window of the loop center
            for (int32_t bin = centerX - win; bin <= centerX + win; bin++) {
                if (bin >= 0) rowIndices[entry.chrom1].insert(bin);
            }
            for (int32_t bin = centerY - win; bin <= centerY + win; bin++) {
                if (bin >= 0) colIndices[entry.chrom2].insert(bin);
            }
        }
    }

    bool probablyContainsRecord(const std::string& chr1, const std::string& chr2,
                              int32_t binX, int32_t binY) const {
        auto rowIt = rowIndices.find(chr1);
        auto colIt = colIndices.find(chr2);
        return rowIt != rowIndices.end() && colIt != colIndices.end() &&
               rowIt->second.count(binX) && colIt->second.count(binY);
    }
};

// Structure to hold APA matrix
struct APAMatrix {
    std::vector<std::vector<float>> matrix;
    int width;
    
    APAMatrix(int size) : width(size) {
        matrix.resize(size, std::vector<float>(size, 0.0f));
    }
    
    void add(int relX, int relY, float value) {
        if (relX >= 0 && relX < width && relY >= 0 && relY < width) {
            matrix[relX][relY] += value;
        }
    }

    // Save matrix to file (format determined by extension)
    void save(const std::string& filename) const;
};

APAMatrix processSliceFile(const std::string& slice_file, 
                         const std::vector<BedpeEntry>& bedpe_entries,
                         const std::string& output_file,
                         int window_size = 10);  // window size in bins

#endif 