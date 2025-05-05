#ifndef COVERAGE_VECTORS_H
#define COVERAGE_VECTORS_H

#include <cstdint>
#include <map>
#include <vector>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <string>

// CoverageVectors accumulates contact coverage per chromosomal bin.
// Bins are created on demand: vectors[chrKey][bin] holds the summed values.
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
    void add(int16_t chrKey, int32_t bin, float value) {
        if (bin < 0) return;  // ignore negative bins
        auto& vec = vectors_[chrKey];
        if (bin >= static_cast<int32_t>(vec.size())) {
            vec.resize(bin + 1, 0.0f);
        }
        vec[bin] += value;
    }

    // Access the internal coverage map
    const std::map<int16_t, std::vector<float>>& getVectors() const {
        return vectors_;
    }

    // Read coverage vectors from a TSV file
    void readFromTSV(const std::string& filename, const std::map<std::string, int16_t>& chromNameToKey);

private:
    int32_t resolution_;
    std::map<int16_t, std::vector<float>> vectors_;
};

#endif // COVERAGE_VECTORS_H
