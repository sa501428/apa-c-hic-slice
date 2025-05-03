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
#include <iomanip>
#include "apa_matrix.h"
#include <chrono>

#ifdef __linux__
#include <sys/sysinfo.h>
#endif

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

    inline double estimateMemoryBytes(uint64_t N_loops, uint32_t window_size, uint32_t N_sets) {
        // Estimate simplified peak memory usage
        double matrixW = 2.0 * window_size + 1.0;
        double perSet  = 4.0 * matrixW * matrixW    // APA matrix
                       + 8.0 * matrixW;             // row+col sums
        double total   = 12.0 * N_loops             // LoopInfo structs
                       + perSet * N_sets; 
        return total * 1.1;                        // +10% overhead
    }

    inline void checkMemoryRequirements(
        const std::vector<std::vector<BedpeEntry>>& bedpe_entries,
        int window_size,
        bool verbose = false)
    {
        // Count total loops and sets
        uint64_t total_loops = 0;
        for (const auto& entries : bedpe_entries) {
            total_loops += entries.size();
        }
        uint32_t num_sets = bedpe_entries.size();

        // Calculate memory estimate
        double need_bytes = estimateMemoryBytes(total_loops, window_size, num_sets);
        double need_gb = need_bytes / (1024.0 * 1024 * 1024);

        // Get memory information
        double total_gb = 0.0;
        double avail_gb = 0.0;

        // First try SLURM environment variable
        char* slurm_mem = std::getenv("SLURM_MEM_PER_NODE");
        if (slurm_mem != nullptr) {
            // SLURM_MEM_PER_NODE is in MB
            total_gb = std::stod(slurm_mem) / 1024.0;
            avail_gb = total_gb;  // On SLURM, we get the full node memory
        } else {
            // Try sysinfo() for non-SLURM environments
            #ifdef __linux__
            struct sysinfo si;
            if (sysinfo(&si) == 0) {
                total_gb = (si.totalram * si.mem_unit) / (1024.0 * 1024 * 1024);
                avail_gb = ((si.freeram + si.bufferram) * si.mem_unit) / (1024.0 * 1024 * 1024);
            } else {
                // Fallback to conservative estimate if sysinfo fails
                total_gb = 16.0;
                avail_gb = 14.0;
            }
            #else
            // Non-Linux systems get conservative estimate
            total_gb = 16.0;
            avail_gb = 14.0;
            #endif
        }

        std::cout << "\nMemory Requirements:\n"
                  << "  Needed:       " << std::fixed << std::setprecision(2) << need_gb  << " GB\n"
                  << "  Total system: " << std::setprecision(2) << total_gb << " GB\n"
                  << "  Available:    " << std::setprecision(2) << avail_gb << " GB\n\n";

        if (need_gb > avail_gb) {
            throw std::runtime_error(
                "Insufficient memory: need " +
                std::to_string(need_gb) + " GB but only " +
                std::to_string(avail_gb) + " GB available"
            );
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

    RegionsOfInterest(int32_t res, int32_t win, bool inter) 
        : resolution(res), window(win), isInter(inter) {}

    void clear() {
        rowIndices.clear();
        colIndices.clear();
    }

    void addEntries(const std::vector<BedpeEntry>& bedpe_entries) {
        // Pre-reserve space for better performance if this is the first addition
        if (rowIndices.empty()) {
            for (const auto& entry : bedpe_entries) {
                rowIndices[entry.chrom1].reserve(detail::getChromBins(entry.chrom1, resolution));
                colIndices[entry.chrom2].reserve(detail::getChromBins(entry.chrom2, resolution));
            }
        }

        for (const auto& entry : bedpe_entries) {
            // Convert BEDPE coordinates to bin positions
            int32_t centerX = entry.gmid1 / resolution;
            int32_t centerY = entry.gmid2 / resolution;
            
            // Add all possible bins within window of the loop center
            for (int32_t bin = centerX - window; bin <= centerX + window; bin++) {
                if (bin >= 0) rowIndices[entry.chrom1].insert(bin);
            }
            for (int32_t bin = centerY - window; bin <= centerY + window; bin++) {
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

    bool probablyContainsPartialRecord(const std::string& chr1, const std::string& chr2,
                               int32_t binX, int32_t binY) const {
        // Quick filter for inter/intra
        if (isInter && chr1 == chr2) return false;
        if (!isInter && chr1 != chr2) return false;
        
        // Check if either bin is in any region of interest
        auto rowIt = rowIndices.find(chr1);
        if (rowIt != rowIndices.end() && rowIt->second.count(binX)) {
            return true;
        }
        
        auto colIt = colIndices.find(chr2);
        if (colIt != colIndices.end() && colIt->second.count(binY)) {
            return true;
        }
        
        return false;
    }
};

struct LoopIndex {
    const int32_t BIN_GROUP_SIZE;
    // Map of (chrom1,chrom2) -> (binX,binY) -> loops
    std::map<ChromPair, std::map<std::pair<int32_t, int32_t>, std::vector<LoopInfo>>> loops;
    int32_t resolution;
    
    LoopIndex(const std::vector<BedpeEntry>& bedpe_entries, 
              int32_t res,
              const std::map<std::string, int16_t>& chromNameToKey,
              int32_t window_size) 
        : BIN_GROUP_SIZE(3 * window_size), resolution(res) {
        for (const auto& loop : bedpe_entries) {
            ChromPair chrom_pair{
                chromNameToKey.at(loop.chrom1),
                chromNameToKey.at(loop.chrom2)
            };
            int32_t binX = loop.gmid1 / resolution;
            int32_t binY = loop.gmid2 / resolution;
            int32_t bin_group_x = binX / BIN_GROUP_SIZE;
            int32_t bin_group_y = binY / BIN_GROUP_SIZE;
            loops[chrom_pair][{bin_group_x, bin_group_y}].emplace_back(loop, chromNameToKey);
        }
    }
    
    std::vector<const LoopInfo*> getNearbyLoops(int16_t chr1Key, int16_t chr2Key, 
                                               int32_t binX, int32_t binY) const {
        ChromPair chrom_pair{chr1Key, chr2Key};
        auto chrom_it = loops.find(chrom_pair);
        if (chrom_it == loops.end()) return {};
        
        int32_t bin_group_x = binX / BIN_GROUP_SIZE;
        int32_t bin_group_y = binY / BIN_GROUP_SIZE;
        std::vector<const LoopInfo*> nearby_loops;
        nearby_loops.reserve(10);
        
        // Check 3x3 grid of bin groups around the target
        for (int32_t i = -1; i <= 1; i++) {
            for (int32_t j = -1; j <= 1; j++) {
                auto group_it = chrom_it->second.find({bin_group_x + i, bin_group_y + j});
                if (group_it != chrom_it->second.end()) {
                    for (const auto& loop : group_it->second) {
                        nearby_loops.push_back(&loop);
                    }
                }
            }
        }
        return nearby_loops;
    }
};

// Structure to hold coverage vectors
struct CoverageVectors {
    // Change to sparse representation using unordered_map
    std::unordered_map<int16_t, std::unordered_map<int32_t, float>> vectors;  // chromKey -> (bin -> coverage)
    int32_t resolution;
    
    CoverageVectors(int32_t res) : resolution(res) {}
    
    void add(int16_t chromKey, int32_t bin, float value) {
        static const int32_t MAX_VECTOR_SIZE = 30000000;  // 30 million elements
        if (bin >= MAX_VECTOR_SIZE) {
            throw std::runtime_error("Bin index exceeds maximum allowed size");
        }
        // Only store non-zero values
        if (value > 0) {
            vectors[chromKey][bin] += value;
        }
    }

    void addLocalSums(std::vector<float>& sums, int16_t chromKey, int32_t binStart) const {
        auto it = vectors.find(chromKey);
        if (it != vectors.end()) {
            const auto& sparse_vec = it->second;
            for (size_t i = 0; i < sums.size(); i++) {
                int32_t bin = binStart + static_cast<int32_t>(i);
                auto bin_it = sparse_vec.find(bin);
                if (bin_it != sparse_vec.end()) {
                    sums[i] += bin_it->second;
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
    long max_genome_dist = 0,
    long job_id = 0,
    bool verbose = false);

// Helper function to print timestamp with ID
inline void printTimestamp(const std::string& message, long id, bool verbose = false) {
    if (!verbose) return;
    auto now = std::chrono::system_clock::now();
    auto time = std::chrono::system_clock::to_time_t(now);
    std::cout << message << " [ID: " << id << "] (" << std::ctime(&time) << ")" << std::endl;
}

#endif 