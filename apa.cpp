#include "apa.h"
#include <zlib.h>
#include <stdexcept>
#include <cstring>
#include <cmath>
#include <map>


float processSliceFile(const std::string& slice_file, 
                      const std::vector<BedpeEntry>& bedpe_entries) {
    // Open slice file
    gzFile file = gzopen(slice_file.c_str(), "rb");
    if (!file) {
        throw std::runtime_error("Could not open file: " + slice_file);
    }

    try {
        // Read and verify magic string
        char magic[8];
        if (gzread(file, magic, 8) != 8 || strncmp(magic, "HICSLICE", 8) != 0) {
            throw std::runtime_error("Invalid file format: missing magic string");
        }

        // Read resolution
        int32_t resolution;
        if (gzread(file, &resolution, sizeof(int32_t)) != sizeof(int32_t)) {
            throw std::runtime_error("Failed to read resolution");
        }

        // Read chromosome mapping
        int32_t numChromosomes;
        if (gzread(file, &numChromosomes, sizeof(int32_t)) != sizeof(int32_t)) {
            throw std::runtime_error("Failed to read chromosome count");
        }

        std::map<int16_t, std::string> chromosomeKeyToName;
        for (int i = 0; i < numChromosomes; i++) {
            int32_t nameLength;
            if (gzread(file, &nameLength, sizeof(int32_t)) != sizeof(int32_t)) {
                throw std::runtime_error("Failed to read chromosome name length");
            }

            std::vector<char> nameBuffer(nameLength + 1, 0);
            if (gzread(file, nameBuffer.data(), nameLength) != nameLength) {
                throw std::runtime_error("Failed to read chromosome name");
            }
            std::string chromosomeName(nameBuffer.data());

            int16_t key;
            if (gzread(file, &key, sizeof(int16_t)) != sizeof(int16_t)) {
                throw std::runtime_error("Failed to read chromosome key");
            }

            chromosomeKeyToName[key] = chromosomeName;
        }

        // Create map for faster BEDPE lookup
        using RegionPair = std::pair<BinRegion, BinRegion>;
        std::map<std::pair<std::string, std::string>, std::vector<RegionPair>> regionMap;
        
        // Convert BEDPE entries to bin coordinates and store in both orientations
        for (const auto& entry : bedpe_entries) {
            BinRegion region1{entry.start1 / resolution, (entry.end1 + resolution - 1) / resolution};
            BinRegion region2{entry.start2 / resolution, (entry.end2 + resolution - 1) / resolution};
            
            // Store in both chromosome orientations since Hi-C contacts can be in either order
            regionMap[{entry.chrom1, entry.chrom2}].push_back({region1, region2});
            if (entry.chrom1 != entry.chrom2) {
                regionMap[{entry.chrom2, entry.chrom1}].push_back({region2, region1});
            }
        }

        // Process contact records
        float totalCount = 0;
        struct {
            int16_t chr1Key;
            int32_t binX;
            int16_t chr2Key;
            int32_t binY;
            float value;
        } record;

        while (gzread(file, &record, sizeof(record)) == sizeof(record)) {
            if (std::isnan(record.value) || std::isinf(record.value) || record.value <= 0) {
                continue;
            }

            std::string chr1 = chromosomeKeyToName[record.chr1Key];
            std::string chr2 = chromosomeKeyToName[record.chr2Key];
            
            auto it = regionMap.find({chr1, chr2});
            if (it != regionMap.end()) {
                for (const auto& region : it->second) {
                    if (record.binX >= region.first.start && record.binX < region.first.end &&
                        record.binY >= region.second.start && record.binY < region.second.end) {
                        totalCount += record.value;
                        break;  // Count each contact only once
                    }
                }
            }
        }

        gzclose(file);
        return totalCount;

    } catch (...) {
        gzclose(file);
        throw;
    }
}
