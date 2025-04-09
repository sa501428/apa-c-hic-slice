#ifndef APA_H
#define APA_H

#include "bedpe_builder.h"
#include <string>
#include <vector>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <stdexcept>
#include <cstdint>  // For int32_t
#include <iostream>
#include <algorithm>
#include <map>

// Add after existing includes
#include <unordered_map>

// Add before the RegionsOfInterest struct
namespace {
    // Default chromosome sizes (in bp)
    const std::unordered_map<std::string, int64_t> DEFAULT_CHROM_SIZES = {
        {"chr1", 248956422}, {"chr2", 242193529}, {"chr3", 198295559},
        {"chr4", 190214555}, {"chr5", 181538259}, {"chr6", 170805979},
        {"chr7", 159345973}, {"chrX", 156040895}, {"chr8", 145138636},
        {"chr9", 138394717}, {"chr11", 135086622}, {"chr10", 133797422},
        {"chr12", 133275309}, {"chr13", 114364328}, {"chr14", 107043718},
        {"chr15", 101991189}, {"chr16", 90338345}, {"chr17", 83257441},
        {"chr18", 80373285}, {"chr20", 64444167}, {"chr19", 58617616},
        {"chrY", 57227415}, {"chr22", 50818468}, {"chr21", 46709983}
    };

    // Calculate number of bins for a chromosome at given resolution
    inline int32_t getChromBins(const std::string& chrom, int32_t resolution) {
        auto it = DEFAULT_CHROM_SIZES.find(chrom);
        if (it != DEFAULT_CHROM_SIZES.end()) {
            return (it->second / resolution) + 1;
        }
        return 20000000 / resolution; // Default fallback size
    }
}

// Structure to hold binned regions for faster lookup
struct BinRegion {
    int32_t start;
    int32_t end;
};

// Structure to hold possible bin positions for quick filtering
struct RegionsOfInterest {
    std::unordered_map<std::string, std::unordered_set<int32_t>> rowIndices;  // chrom -> set of binX
    std::unordered_map<std::string, std::unordered_set<int32_t>> colIndices;  // chrom -> set of binY
    int32_t resolution;
    int32_t window;
    bool isInter;  // true for inter-chromosomal, false for intra-chromosomal

    RegionsOfInterest(const std::vector<BedpeEntry>& bedpe_entries, 
                     int32_t res, int32_t win,
                     bool inter) 
        : resolution(res), window(win), isInter(inter) {
        // Pre-reserve space for better performance
        for (const auto& entry : bedpe_entries) {
            rowIndices[entry.chrom1].reserve(getChromBins(entry.chrom1, resolution));
            colIndices[entry.chrom2].reserve(getChromBins(entry.chrom2, resolution));
        }

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
        // Quick filter for inter/intra
        if (isInter && chr1 == chr2) return false;
        if (!isInter && chr1 != chr2) return false;
        
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
        if (size <= 0) {
            throw std::runtime_error("APAMatrix size must be positive");
        }
        matrix.resize(size, std::vector<float>(size, 0.0f));
    }
    
    void add(int relX, int relY, float value) {
        if (relX >= 0 && relX < width && relY >= 0 && relY < width) {
            matrix[relX][relY] += value;
        }
    }

    // Get average of non-zero values
    static float getAverage(const std::vector<float>& vec) {
        float sum = 0.0f;
        int count = 0;
        for (float val : vec) {
            if (val > 0) {
                sum += val;
                count++;
            }
        }
        return count > 0 ? sum / count : 1.0f;  // Return 1.0 instead of 0 to prevent division by zero
    }

    // Scale vector by its average
    static void scaleByAverage(std::vector<float>& vec) {
        float avg = getAverage(vec);
        if (avg > 0) {
            for (float& val : vec) {
                val /= avg;
            }
        }
    }

    // Normalize matrix using row and column sums
    void normalize(const std::vector<float>& rowSums, const std::vector<float>& colSums) {
        std::vector<std::vector<float>> normalized(width, std::vector<float>(width, 0.0f));
        
        for (int r = 0; r < width; ++r) {
            for (int c = 0; c < width; ++c) {
                float normVal = rowSums[r] * colSums[c];
                if (normVal > 0) {
                    normalized[r][c] = matrix[r][c] / normVal;
                }
                // If normVal is 0, leave normalized[r][c] as 0
            }
        }
        
        matrix = std::move(normalized);
    }

    void save(const std::string& filename) const;
};

// Structure to hold coverage vectors
struct CoverageVectors {
    std::unordered_map<std::string, std::vector<float>> vectors;  // chrom -> coverage vector
    int32_t resolution;
    
    CoverageVectors(int32_t res) : resolution(res) {}
    
    void add(const std::string& chrom, int32_t bin, float value) {
        auto& vec = vectors[chrom];
        if (vec.empty()) {
            vec.resize(getChromBins(chrom, resolution), 0.0f);
        }
        if (static_cast<size_t>(bin) >= vec.size()) {
            size_t new_size = getChromBins(chrom, resolution);
            vec.resize(new_size, 0.0f);
        }
        vec[bin] += value;
    }

    void addLocalSums(std::vector<float>& sums, const std::string& chrom, 
                     int32_t binStart) const {
        auto it = vectors.find(chrom);
        if (it != vectors.end()) {
            const auto& vec = it->second;
            for (size_t i = 0; i < sums.size(); i++) {
                int32_t bin = binStart + static_cast<int32_t>(i);
                if (bin >= 0 && static_cast<size_t>(bin) < vec.size()) {
                    sums[i] += vec[bin];
                }
            }
        }
    }
};

// Add this structure to help with loop lookup
struct ChromPair {
    std::string chrom1;
    std::string chrom2;
    
    bool operator<(const ChromPair& other) const {
        if (chrom1 != other.chrom1) return chrom1 < other.chrom1;
        return chrom2 < other.chrom2;
    }
};

// Structure to hold preprocessed loops for fast lookup
struct LoopIndex {
    static const int32_t BIN_GROUP_SIZE = 1000;  // Group size for binning
    std::map<ChromPair, std::map<int32_t, std::vector<const BedpeEntry*>>> loops;
    int32_t resolution;
    
    LoopIndex(const std::vector<BedpeEntry>& bedpe_entries, int32_t res) : resolution(res) {
        // Preprocess loops into the index
        for (const auto& loop : bedpe_entries) {
            ChromPair chrom_pair{loop.chrom1, loop.chrom2};
            
            // Calculate midpoint bin and group
            int32_t mid_bin = ((loop.start1 + loop.end1) / 2) / resolution;
            int32_t bin_group = mid_bin / BIN_GROUP_SIZE;
            
            // Store pointer to loop in appropriate group
            loops[chrom_pair][bin_group].push_back(&loop);
        }
    }
    
    // Get loops that might be relevant for a given contact
    std::vector<const BedpeEntry*> getNearbyLoops(const std::string& chr1, const std::string& chr2, 
                                                 int32_t binX) const {
        ChromPair chrom_pair{chr1, chr2};
        auto chrom_it = loops.find(chrom_pair);
        if (chrom_it == loops.end()) return {};
        
        // Calculate which bin group this contact belongs to
        int32_t bin_group = binX / BIN_GROUP_SIZE;
        
        std::vector<const BedpeEntry*> nearby_loops;
        nearby_loops.reserve(10);  // Pre-allocate space for 1M loops
        
        // Check the bin group and adjacent groups
        for (int32_t i = -1; i <= 1; i++) {
            auto group_it = chrom_it->second.find(bin_group + i);
            if (group_it != chrom_it->second.end()) {
                nearby_loops.insert(nearby_loops.end(), 
                                  group_it->second.begin(), 
                                  group_it->second.end());
            }
        }
        
        return nearby_loops;
    }
};

// Update function declaration
std::vector<APAMatrix> processSliceFile(
    const std::string& slice_file, 
    const std::vector<std::vector<BedpeEntry>>& all_bedpe_entries,
    int window_size = 10,
    bool isInter = false,
    long min_genome_dist = 0,
    long max_genome_dist = 0);

#endif 