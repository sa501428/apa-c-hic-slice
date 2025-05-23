#ifndef APA_H
#define APA_H

#include "bedpe_builder.h"  // Must come first since it defines BedpeEntry
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
#include <sys/sysinfo.h>
#include <iomanip>

// Forward declarations
struct RegionsOfInterest;
struct LoopIndex;
struct APAMatrix;
struct CoverageVectors;

namespace detail {
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

    // Make all function definitions inline
    inline size_t estimateChromCoverageMemory(const std::string& chrom, int32_t resolution) {
        auto it = DEFAULT_CHROM_SIZES.find(chrom);
        if (it != DEFAULT_CHROM_SIZES.end()) {
            return (it->second / resolution + 1) * 4;  // 4 bytes per float
        }
        return 20000000 / resolution * 4; // Default size, 4 bytes per float
    }

    inline size_t estimateMemoryUsage(const std::vector<std::vector<BedpeEntry>>& bedpe_entries, 
                              int window_size) {
        size_t total_bedpes = 0;
        std::set<std::string> unique_chroms;
        
        // Count total BEDPEs and unique chromosomes
        for (const auto& entries : bedpe_entries) {
            total_bedpes += entries.size();
            for (const auto& entry : entries) {
                unique_chroms.insert(entry.chrom1);
                unique_chroms.insert(entry.chrom2);
            }
        }
        
        size_t peak_memory = 0;
        size_t current_memory = 0;
        
        // Phase 1: Initial loading and structure creation
        // BedpeEntry size: 2 strings (32 bytes each) + 4 longs (8 bytes each)
        current_memory += total_bedpes * (64 + 32);
        
        // RegionsOfInterest (2 maps of sets of bin indices)
        size_t roi_size = bedpe_entries.size() * // Number of ROIs
                         total_bedpes * // Entries per ROI
                         (window_size * 2 + 1) * 2 * // Indices per entry (both dimensions)
                         4; // Size of int32_t
        current_memory += roi_size;
        
        // LoopIndex structures with LoopInfo
        size_t loop_index_size = 0;
        // LoopInfo: 2 strings (32 bytes each) + 6 int32_t (4 bytes each)
        loop_index_size += total_bedpes * (64 + 24) * bedpe_entries.size();
        
        // ChromPair map structures
        size_t chrom_pairs = (unique_chroms.size() * unique_chroms.size()) / 4;
        loop_index_size += chrom_pairs * (
            64 + // ChromPair (2 strings, 32 bytes each)
            48 + // std::map overhead
            // Estimate bin groups per chrom pair: chr1 size = 248956422
            (248956422 / 1000000) * 24 // vector overhead, using 1M as bin group size (1000 * 1000)
        );
        current_memory += loop_index_size;
        
        peak_memory = std::max(peak_memory, current_memory);
        
        // Phase 2: After freeing BedpeEntries
        current_memory -= total_bedpes * (64 + 32);
        
        // APAMatrix (width x width float matrix for each BEDPE set)
        size_t matrix_width = window_size * 2 + 1;
        size_t matrix_size = matrix_width * matrix_width * 4; // 4 bytes per float
        current_memory += bedpe_entries.size() * matrix_size;
        
        // Coverage vectors (one float vector per chromosome)
        size_t coverage_size = 0;
        for (const auto& chrom : unique_chroms) {
            coverage_size += estimateChromCoverageMemory(chrom, 1000);
        }
        current_memory += coverage_size;
        
        peak_memory = std::max(peak_memory, current_memory);
        
        // Phase 3: After freeing ROI
        current_memory -= roi_size;
        
        // Row and column sums (float vector for each BEDPE set)
        size_t sums_size = bedpe_entries.size() * 2 * // Two vectors per BEDPE set
                          matrix_width * 4; // 4 bytes per float
        current_memory += sums_size;
        
        peak_memory = std::max(peak_memory, current_memory);

        // Add 30% overhead for STL containers and other overhead
        peak_memory = static_cast<size_t>(peak_memory * 1.3);
        
        return peak_memory;
    }

    inline void checkMemoryRequirements(const std::vector<std::vector<BedpeEntry>>& bedpe_entries,
                               int window_size) {
        size_t estimated_bytes = estimateMemoryUsage(bedpe_entries, window_size);
        
        // Get available system memory
        struct sysinfo si;
        if (sysinfo(&si) == 0) {
            uint64_t total_ram = si.totalram * si.mem_unit;
            uint64_t available_ram = si.freeram * si.mem_unit;
            
            double est_gb = estimated_bytes / (1024.0 * 1024.0 * 1024.0);
            double total_gb = total_ram / (1024.0 * 1024.0 * 1024.0);
            double avail_gb = available_ram / (1024.0 * 1024.0 * 1024.0);
            
            std::cout << "\nMemory Requirements:\n"
                      << "Estimated memory needed: " << std::fixed << std::setprecision(2) 
                      << est_gb << " GB\n"
                      << "Total system RAM: " << std::setprecision(2) << total_gb << " GB\n"
                      << "Available RAM: " << std::setprecision(2) << avail_gb << " GB\n\n";
            
            if (estimated_bytes > available_ram * 0.9) { // Leave 10% buffer
                throw std::runtime_error(
                    "Insufficient memory available. Need " + 
                    std::to_string(est_gb) + " GB but only " + 
                    std::to_string(avail_gb) + " GB available"
                );
            }
        } else {
            std::cerr << "Warning: Could not check system memory. Continuing without verification.\n";
        }
    }
}

// Define structures first
struct LoopInfo {
    std::string chrom1;
    std::string chrom2;
    int32_t start1;
    int32_t end1;
    int32_t start2;
    int32_t end2;
    
    LoopInfo(const BedpeEntry& entry) 
        : chrom1(entry.chrom1), chrom2(entry.chrom2),
          start1(entry.start1), end1(entry.end1),
          start2(entry.start2), end2(entry.end2) {}
};

struct ChromPair {
    std::string chrom1;
    std::string chrom2;
    
    bool operator<(const ChromPair& other) const {
        if (chrom1 != other.chrom1) return chrom1 < other.chrom1;
        return chrom2 < other.chrom2;
    }
};

// Structure to hold binned regions for faster lookup
struct BinRegion {
    int32_t start;
    int32_t end;
};

struct RegionsOfInterest {
    std::unordered_map<std::string, std::unordered_set<int32_t>> rowIndices;
    std::unordered_map<std::string, std::unordered_set<int32_t>> colIndices;
    int32_t resolution;
    int32_t window;
    bool isInter;

    RegionsOfInterest(const std::vector<BedpeEntry>& bedpe_entries, 
                     int32_t res, int32_t win,
                     bool inter) 
        : resolution(res), window(win), isInter(inter) {
        // Pre-reserve space for better performance
        for (const auto& entry : bedpe_entries) {
            rowIndices[entry.chrom1].reserve(detail::getChromBins(entry.chrom1, resolution));
            colIndices[entry.chrom2].reserve(detail::getChromBins(entry.chrom2, resolution));
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

struct LoopIndex {
    static const int32_t BIN_GROUP_SIZE = 1000;
    std::map<ChromPair, std::map<int32_t, std::vector<LoopInfo>>> loops;
    int32_t resolution;
    
    LoopIndex(const std::vector<BedpeEntry>& bedpe_entries, int32_t res) : resolution(res) {
        for (const auto& loop : bedpe_entries) {
            ChromPair chrom_pair{loop.chrom1, loop.chrom2};
            int32_t mid_bin = ((loop.start1 + loop.end1) / 2) / resolution;
            int32_t bin_group = mid_bin / BIN_GROUP_SIZE;
            loops[chrom_pair][bin_group].emplace_back(loop);
        }
    }
    
    std::vector<const LoopInfo*> getNearbyLoops(const std::string& chr1, const std::string& chr2, 
                                               int32_t binX) const {
        ChromPair chrom_pair{chr1, chr2};
        auto chrom_it = loops.find(chrom_pair);
        if (chrom_it == loops.end()) return {};
        
        int32_t bin_group = binX / BIN_GROUP_SIZE;
        std::vector<const LoopInfo*> nearby_loops;
        nearby_loops.reserve(10);
        
        for (int32_t i = -1; i <= 1; i++) {
            auto group_it = chrom_it->second.find(bin_group + i);
            if (group_it != chrom_it->second.end()) {
                for (const auto& loop : group_it->second) {
                    nearby_loops.push_back(&loop);
                }
            }
        }
        return nearby_loops;
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
            vec.resize(detail::getChromBins(chrom, resolution), 0.0f);
        }
        if (static_cast<size_t>(bin) >= vec.size()) {
            size_t new_size = detail::getChromBins(chrom, resolution);
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

// Update function declaration
std::vector<APAMatrix> processSliceFile(
    const std::string& slice_file, 
    const std::vector<std::vector<BedpeEntry>>& all_bedpe_entries,
    int window_size = 10,
    bool isInter = false,
    long min_genome_dist = 0,
    long max_genome_dist = 0);

#endif 