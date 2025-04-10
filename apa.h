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
#include "apa_matrix.h"

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
        
        // Phase 1: Initial data structures
        
        // LoopInfo size (2 int16_t + 2 int32_t = 12 bytes)
        size_t loop_info_size = total_bedpes * 12;
        current_memory += loop_info_size;
        
        // RegionsOfInterest - bin indices for each chromosome
        size_t matrix_width = window_size * 2 + 1;
        // Each bin is an int32_t (4 bytes) in an unordered_set
        size_t bins_per_chrom = matrix_width * 4;
        size_t roi_size = unique_chroms.size() * bins_per_chrom * 2; // For both row and col indices
        current_memory += roi_size;
        
        peak_memory = std::max(peak_memory, current_memory);
        
        // Phase 2: Analysis structures
        
        // APAMatrix (one per BEDPE set)
        size_t matrix_size = matrix_width * matrix_width * 4; // 4 bytes per float
        current_memory += bedpe_entries.size() * matrix_size;
        
        // Coverage vectors (one per chromosome)
        size_t avg_chrom_size = 100000000; // 100Mb average chromosome size
        size_t bins_per_coverage = avg_chrom_size / 1000; // 1kb resolution
        size_t coverage_size = unique_chroms.size() * bins_per_coverage * 4; // 4 bytes per float
        current_memory += coverage_size;
        
        peak_memory = std::max(peak_memory, current_memory);
        
        // Add 5% overhead for STL containers and memory fragmentation
        peak_memory = static_cast<size_t>(peak_memory * 1.05);
        
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
            
            if (estimated_bytes > available_ram) {
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
    int16_t chrom1Key;
    int16_t chrom2Key;
    int32_t gmid1;
    int32_t gmid2;
    
    LoopInfo(const BedpeEntry& entry, 
             const std::map<std::string, int16_t>& chromNameToKey) 
        : chrom1Key(chromNameToKey.at(entry.chrom1))
        , chrom2Key(chromNameToKey.at(entry.chrom2))
        , gmid1(entry.gmid1)
        , gmid2(entry.gmid2)
        {}
};

struct ChromPair {
    int16_t chrom1Key;
    int16_t chrom2Key;
    
    bool operator<(const ChromPair& other) const {
        if (chrom1Key != other.chrom1Key) return chrom1Key < other.chrom1Key;
        return chrom2Key < other.chrom2Key;
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
            int32_t centerX = entry.gmid1 / resolution;
            int32_t centerY = entry.gmid2 / resolution;
            
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
    
    LoopIndex(const std::vector<BedpeEntry>& bedpe_entries, 
              int32_t res,
              const std::map<std::string, int16_t>& chromNameToKey) 
        : resolution(res) {
        for (const auto& loop : bedpe_entries) {
            ChromPair chrom_pair{
                chromNameToKey.at(loop.chrom1),
                chromNameToKey.at(loop.chrom2)
            };
            int32_t mid_bin = loop.gmid1 / resolution;
            int32_t bin_group = mid_bin / BIN_GROUP_SIZE;
            loops[chrom_pair][bin_group].emplace_back(loop, chromNameToKey);
        }
    }
    
    std::vector<const LoopInfo*> getNearbyLoops(int16_t chr1Key, int16_t chr2Key, 
                                               int32_t binX) const {
        ChromPair chrom_pair{chr1Key, chr2Key};
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

// Structure to hold coverage vectors
struct CoverageVectors {
    std::unordered_map<int16_t, std::vector<float>> vectors;  // chromKey -> coverage vector
    int32_t resolution;
    
    CoverageVectors(int32_t res) : resolution(res) {}
    
    void add(int16_t chromKey, const std::string& chromName, int32_t bin, float value) {
        static const int32_t MAX_VECTOR_SIZE = 30000000;  // 30 million elements
        if (bin >= MAX_VECTOR_SIZE) {
            throw std::runtime_error("Bin index exceeds maximum allowed size");
        }
        auto& vec = vectors[chromKey];
        if (vec.empty()) {
            vec.resize(detail::getChromBins(chromName, resolution), 0.0f);
        }
        if (static_cast<size_t>(bin) >= vec.size()) {
            size_t new_size = detail::getChromBins(chromName, resolution);
            vec.resize(new_size, 0.0f);
        }
        vec[bin] += value;
    }

    void addLocalSums(std::vector<float>& sums, int16_t chromKey, int32_t binStart) const {
        auto it = vectors.find(chromKey);
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