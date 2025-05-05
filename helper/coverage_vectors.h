#ifndef COVERAGE_VECTORS_H
#define COVERAGE_VECTORS_H

#include <cstdint>
#include <unordered_map>
#include <vector>
#include <stdexcept>

// CoverageVectors accumulates contact coverage per chromosomal bin.
// Uses sparse storage: vectors[chrKey][bin] holds the summed values.
class CoverageVectors {
public:
    // resolution unused here but kept for potential future use
    explicit CoverageVectors(int32_t resolution)
        : resolution_(resolution) {
        if (resolution_ <= 0) {
            throw std::invalid_argument("Resolution must be positive");
        }
    }

    // Add a contact value to a given chromosome key and bin index
    void add(int16_t chromKey, int32_t bin, float value) {
        // Only store non-zero values
        if (value > 0) {
            auto& chrom_map = vectors_[chromKey];
            auto it = chrom_map.find(bin);
            if (it != chrom_map.end()) {
                it->second += value;
            } else {
                chrom_map.emplace(bin, value);
            }
        }
    }

    // Access the internal coverage map
    const std::unordered_map<int16_t, std::unordered_map<int32_t, float>>& getVectors() const {
        return vectors_;
    }

    void addLocalSums(std::vector<float>& sums, int16_t chromKey, int32_t binStart) const {
        auto it = vectors_.find(chromKey);
        if (it != vectors_.end()) {
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

private:
    int32_t resolution_;
    std::unordered_map<int16_t, std::unordered_map<int32_t, float>> vectors_;
};

#endif // COVERAGE_VECTORS_H
