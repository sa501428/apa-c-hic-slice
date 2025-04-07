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
    bool isInter;  // true for inter-chromosomal, false for intra-chromosomal

    RegionsOfInterest(const std::vector<BedpeEntry>& bedpe_entries, 
                     int32_t res, int32_t win,
                     bool inter) 
        : resolution(res), window(win), isInter(inter) {
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
        return count > 0 ? sum / count : 0.0f;
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
                if (normVal > 0.0f) {
                    normalized[r][c] = matrix[r][c] / normVal;
                }
            }
        }
        
        matrix = std::move(normalized);
    }

    void save(const std::string& filename) const;
};

// Structure to hold coverage vectors
struct CoverageVectors {
    std::unordered_map<std::string, std::vector<float>> vectors;  // chrom -> coverage vector
    static const int32_t INITIAL_SIZE = 20000000;  // 20Mb
    int32_t resolution;
    
    CoverageVectors(int32_t res) : resolution(res) {}
    
    void add(const std::string& chrom, int32_t bin, float value) {
        auto& vec = vectors[chrom];
        if (vec.empty()) {
            // Pre-allocate with reasonable size on first use
            vec.resize(INITIAL_SIZE / resolution, 0.0f);
        }
        if (bin >= vec.size()) {
            // If we need more space, double the size
            vec.resize(std::max(vec.size() * 2, static_cast<size_t>(bin * 2)), 0.0f);
        }
        vec[bin] += value;
    }

    // Get sums for a window around a bin position
    void addLocalSums(std::vector<float>& sums, const std::string& chrom, 
                     int32_t binStart, int32_t windowSize) const {
        auto it = vectors.find(chrom);
        if (it != vectors.end()) {
            const auto& vec = it->second;
            for (int i = 0; i < windowSize * 2 + 1; i++) {
                int32_t bin = binStart + i;
                if (bin >= 0 && bin < vec.size()) {
                    sums[i] += vec[bin];
                }
            }
        }
    }
};

APAMatrix processSliceFile(const std::string& slice_file, 
                         const std::vector<BedpeEntry>& bedpe_entries,
                         const std::string& output_file,
                         int window_size = 10,
                         bool isInter = false,
                         long min_genome_dist = 0,    // Add distance parameters
                         long max_genome_dist = 0);   // for intra filtering

#endif 