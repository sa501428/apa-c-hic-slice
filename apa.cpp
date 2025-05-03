#include "apa.h"
#include <zlib.h>
#include <stdexcept>
#include <cstring>
#include <cmath>
#include <map>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <iostream>
#include "apa_matrix.h"
#include "vector_tools.h"


std::vector<APAMatrix> processSliceFile(
    const std::string& slice_file,
    const std::vector<std::vector<BedpeEntry>>& all_bedpe_entries_const,
    int window_size,
    bool isInter,
    long min_genome_dist,
    long max_genome_dist) {
    
    // Make a copy we can modify
    auto all_bedpe_entries = all_bedpe_entries_const;
    
    std::cout << "Opening slice file..." << std::endl;
    if (window_size <= 0) {
        throw std::runtime_error("Window size must be positive");
    }

    // Declare chromosome mapping here
    std::map<int16_t, std::string> chromosomeKeyToName;
    std::map<std::string, int16_t> chromNameToKey;

    // Try opening as uncompressed first
    FILE* raw_file = fopen(slice_file.c_str(), "rb");
    gzFile gz_file = nullptr;
    bool is_compressed = false;

    if (!raw_file) {
        // Try opening as compressed
        gz_file = gzopen(slice_file.c_str(), "rb");
        if (!gz_file) {
            throw std::runtime_error("Could not open file: " + slice_file);
        }
        is_compressed = true;
    }

    std::cout << "File opened..." << std::endl;

    try {
        // Read and verify magic string
        char magic[8];
        if (is_compressed) {
            if (gzread(gz_file, magic, 8) != 8 || strncmp(magic, "HICSLICE", 8) != 0) {
                throw std::runtime_error("Invalid file format: missing magic string");
            }
        } else {
            if (fread(magic, 1, 8, raw_file) != 8 || strncmp(magic, "HICSLICE", 8) != 0) {
                throw std::runtime_error("Invalid file format: missing magic string");
            }
        }

        // Read resolution
        int32_t resolution;
        if (is_compressed) {
            if (gzread(gz_file, &resolution, sizeof(int32_t)) != sizeof(int32_t)) {
                throw std::runtime_error("Failed to read resolution");
            }
        } else {
            if (fread(&resolution, sizeof(int32_t), 1, raw_file) != 1) {
                throw std::runtime_error("Failed to read resolution");
            }
        }

        if (resolution <= 0) {
            throw std::runtime_error("Invalid resolution in slice file");
        }

        std::cout << "Resolution is " << resolution << std::endl;

        // Perform memory check now that we have the resolution
        detail::checkMemoryRequirements(all_bedpe_entries, window_size, resolution);

        // Move chromosome mapping read to here, before creating data structures
        int32_t numChromosomes;
        if (is_compressed) {
            if (gzread(gz_file, &numChromosomes, sizeof(int32_t)) != sizeof(int32_t)) {
                throw std::runtime_error("Failed to read chromosome count");
            }
        } else {
            if (fread(&numChromosomes, sizeof(int32_t), 1, raw_file) != 1) {
                throw std::runtime_error("Failed to read chromosome count");
            }
        }

        if (numChromosomes <= 0) {
            throw std::runtime_error("Invalid number of chromosomes in slice file");
        }

        std::cout << "Number of chromosomes: " << numChromosomes << std::endl;

        for (int i = 0; i < numChromosomes; i++) {
            int32_t nameLength;
            if (is_compressed) {
                if (gzread(gz_file, &nameLength, sizeof(int32_t)) != sizeof(int32_t)) {
                    throw std::runtime_error("Failed to read chromosome name length");
                }
            } else {
                if (fread(&nameLength, sizeof(int32_t), 1, raw_file) != 1) {
                    throw std::runtime_error("Failed to read chromosome name length");
                }
            }

            std::vector<char> nameBuffer(nameLength + 1, 0);
            if (is_compressed) {
                if (gzread(gz_file, nameBuffer.data(), nameLength) != nameLength) {
                    throw std::runtime_error("Failed to read chromosome name");
                }
            } else {
                if (fread(nameBuffer.data(), 1, nameLength, raw_file) != static_cast<size_t>(nameLength)) {
                    throw std::runtime_error("Failed to read chromosome name");
                }
            }
            std::string chromosomeName(nameBuffer.data());

            int16_t key;
            if (is_compressed) {
                if (gzread(gz_file, &key, sizeof(int16_t)) != sizeof(int16_t)) {
                    throw std::runtime_error("Failed to read chromosome key");
                }
            } else {
                if (fread(&key, sizeof(int16_t), 1, raw_file) != 1) {
                    throw std::runtime_error("Failed to read chromosome key");
                }
            }

            chromosomeKeyToName[key] = chromosomeName;
            chromNameToKey[chromosomeName] = key;
            std::cout << "Read chromosome: " << chromosomeName << " (key=" << key << ")" << std::endl;
        }

        // Create a single RegionsOfInterest for all BEDPE sets
        RegionsOfInterest roi(resolution, window_size, isInter);
        
        // Add regions from all BEDPE sets to the single ROI
        for (const auto& entries : all_bedpe_entries) {
            roi.addEntries(entries);
        }

        // Create loop indices for each BEDPE set
        std::vector<LoopIndex> all_indices;
        all_indices.reserve(all_bedpe_entries.size());

        for (const auto& entries : all_bedpe_entries) {
            all_indices.emplace_back(entries, resolution, chromNameToKey, window_size);
        }

        // Clear original BedpeEntries as they're no longer needed
        all_bedpe_entries.clear();
        all_bedpe_entries.shrink_to_fit();

        // Create vectors to hold per-bedpe data structures
        size_t num_bedpes = all_indices.size();
        std::vector<APAMatrix> all_matrices;
        std::vector<std::vector<float>> all_rowSums(num_bedpes);
        std::vector<std::vector<float>> all_colSums(num_bedpes);

        // Initialize data structures for each BEDPE set
        all_matrices.reserve(num_bedpes);
        for (size_t i = 0; i < num_bedpes; i++) {
            all_matrices.push_back(APAMatrix(window_size * 2 + 1));
            all_rowSums[i].resize(window_size * 2 + 1, 0.0f);
            all_colSums[i].resize(window_size * 2 + 1, 0.0f);
        }

        std::cout << "Data structures initialized..." << std::endl;

        // Create single coverage vectors instance (shared across all BEDPEs)
        CoverageVectors coverage(resolution);

        // Single pass: process contacts for both coverage and APA
        struct {
            int16_t chr1Key;
            int32_t binX;
            int16_t chr2Key;
            int32_t binY;
            float value;
        } record;

        std::cout << "Processing contacts..." << std::endl;
        int64_t contact_count = 0;

        int32_t max_genome_dist_as_bin = (max_genome_dist / resolution) + (3 * window_size);
        int32_t min_genome_dist_as_bin = (min_genome_dist / resolution) - (3 * window_size);

        size_t record_size = sizeof(record);
        while ((is_compressed ? 
                (gzread(gz_file, &record, record_size) == static_cast<int>(record_size)) :
                (fread(&record, record_size, 1, raw_file) == 1))) {
            contact_count++;
            
            // Print first two records for debugging
            if (contact_count <= 2) {
                std::string chr1 = chromosomeKeyToName[record.chr1Key];
                std::string chr2 = chromosomeKeyToName[record.chr2Key];
                std::cout << "Contact " << contact_count << ": " 
                         << chr1 << ":" << record.binX << " - "
                         << chr2 << ":" << record.binY 
                         << " value=" << record.value << std::endl;
            }

            if (std::isnan(record.value) || std::isinf(record.value) || record.value <= 0) {
                continue;
            }

            if (isInter && record.chr1Key == record.chr2Key) continue;
            if (!isInter && record.chr1Key != record.chr2Key) continue;

            std::string chr1 = chromosomeKeyToName[record.chr1Key];
            std::string chr2 = chromosomeKeyToName[record.chr2Key];
            

            // Add to coverage vectors (after inter/intra but before any other filtering)
            coverage.add(record.chr1Key, chr1, record.binX, record.value);
            if (record.chr1Key != record.chr2Key || record.binX != record.binY) {  // Don't double count diagonal
                coverage.add(record.chr2Key, chr2, record.binY, record.value);
            }

            // For intra-chromosomal analysis, filter by distance from diagonal
            if (!isInter) {
                int32_t bin_distance = std::abs(record.binX - record.binY);
                if (bin_distance < min_genome_dist_as_bin || bin_distance > max_genome_dist_as_bin) {
                    continue;
                }
            }

            // Check against the single ROI
            if (roi.probablyContainsRecord(chr1, chr2, record.binX, record.binY)) {
                // Process for each BEDPE set
                for (size_t bedpe_idx = 0; bedpe_idx < all_indices.size(); bedpe_idx++) {
                    auto nearby_loops = all_indices[bedpe_idx].getNearbyLoops(record.chr1Key, record.chr2Key, 
                                                                           record.binX, record.binY);
                    for (const auto* loop : nearby_loops) {
                        int32_t loopCenterBinX = loop->gmid1 / resolution;
                        int32_t loopCenterBinY = loop->gmid2 / resolution;
                        if (std::abs(record.binX - loopCenterBinX) <= window_size &&
                            std::abs(record.binY - loopCenterBinY) <= window_size) {
                            
                            // Calculate relative position and add to matrix
                            int relX = record.binX - (loopCenterBinX - window_size);
                            int relY = record.binY - (loopCenterBinY - window_size);
                            all_matrices[bedpe_idx].add(relX, relY, record.value);
                        }
                    }
                }
            }
        }
        std::cout << std::endl;

        std::cout << "Finished processing " << contact_count << " contacts" << std::endl;
        
        // Free RegionsOfInterest as it's no longer needed for contact processing
        roi.clear();

        std::cout << "Calculating coverage normalization..." << std::endl;
        // After processing all contacts, normalize each matrix
        for (size_t bedpe_idx = 0; bedpe_idx < all_matrices.size(); bedpe_idx++) {
            // Calculate row and column sums for this matrix
            for (const auto& chrom_pair : all_indices[bedpe_idx].loops) {
                for (const auto& bin_group : chrom_pair.second) {
                    for (const auto& loop : bin_group.second) {
                        int32_t bin1Start = (loop.gmid1/resolution) - window_size;
                        int32_t bin2Start = (loop.gmid2/resolution) - window_size;
                        
                        coverage.addLocalSums(all_rowSums[bedpe_idx], loop.chrom1Key, bin1Start);
                        coverage.addLocalSums(all_colSums[bedpe_idx], loop.chrom2Key, bin2Start);
                    }
                }
            }

            // Scale sums and normalize matrix
            vector_tools::scaleByAverage(all_rowSums[bedpe_idx]);
            vector_tools::scaleByAverage(all_colSums[bedpe_idx]);
            all_matrices[bedpe_idx].normalize(all_rowSums[bedpe_idx], all_colSums[bedpe_idx]);
        }

        if (is_compressed) {
            gzclose(gz_file);
        } else {
            fclose(raw_file);
        }
        
        return all_matrices;

    } catch (const std::exception& e) {
        if (is_compressed && gz_file) {
            gzclose(gz_file);
        } else if (!is_compressed && raw_file) {
            fclose(raw_file);
        }
        throw;
    }
}
