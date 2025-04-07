#include "apa.h"
#include <zlib.h>
#include <stdexcept>
#include <cstring>
#include <cmath>
#include <map>
#include <fstream>
#include <iomanip>
#include <algorithm>

void APAMatrix::save(const std::string& filename) const {
    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("Cannot open output file: " + filename);
    }

    // Write as space-separated values with high precision
    out << std::fixed << std::setprecision(6);
    for (int i = 0; i < width; i++) {
        for (int j = 0; j < width; j++) {
            if (j > 0) out << " ";
            out << matrix[i][j];
        }
        out << "\n";
    }
}

APAMatrix processSliceFile(const std::string& slice_file, 
                         const std::vector<BedpeEntry>& bedpe_entries,
                         const std::string& output_file,
                         int window_size,
                         bool isInter,
                         long min_genome_dist,
                         long max_genome_dist) {
    if (window_size <= 0) {
        throw std::runtime_error("Window size must be positive");
    }

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

        if (resolution <= 0) {
            throw std::runtime_error("Invalid resolution in slice file");
        }

        // Read chromosome mapping
        int32_t numChromosomes;
        if (gzread(file, &numChromosomes, sizeof(int32_t)) != sizeof(int32_t)) {
            throw std::runtime_error("Failed to read chromosome count");
        }

        if (numChromosomes <= 0) {
            throw std::runtime_error("Invalid number of chromosomes in slice file");
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

        // Create coverage vectors and APA matrix
        CoverageVectors coverage(resolution);  // Pass resolution to constructor
        RegionsOfInterest roi(bedpe_entries, resolution, window_size, isInter);
        APAMatrix apaMatrix(window_size * 2 + 1);

        // Single pass: process contacts for both coverage and APA
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
            
            // Quick filter for inter/intra chromosomal contacts
            if (isInter && chr1 == chr2) continue;
            if (!isInter && chr1 != chr2) continue;

            // Add to coverage vectors (after inter/intra but before any other filtering)
            coverage.add(chr1, record.binX, record.value);
            if (chr1 != chr2 || record.binX != record.binY) {  // Don't double count diagonal
                coverage.add(chr2, record.binY, record.value);
            }

            // For intra-chromosomal analysis, filter by genomic distance
            if (!isInter) {
                // Convert bin distance to genomic distance
                int32_t bin_distance = std::abs(record.binX - record.binY);
                int64_t genomic_distance = static_cast<int64_t>(bin_distance) * resolution;
                
                // Add buffer of 3*window_size to both min and max distances
                int64_t buffer = static_cast<int64_t>(3 * window_size) * resolution;
                
                // Skip if outside the distance range (with buffer) used to generate BEDPE
                if (genomic_distance < (min_genome_dist - buffer) || 
                    genomic_distance > (max_genome_dist + buffer)) {
                    continue;
                }
            }

            // Process for APA matrix
            if (roi.probablyContainsRecord(chr1, chr2, record.binX, record.binY)) {
                // Check each loop region
                for (const auto& loop : bedpe_entries) {
                    if (loop.chrom1 == chr1 && loop.chrom2 == chr2) {
                        // Convert BEDPE coordinates to bin positions
                        int32_t bin1Start = loop.start1 / resolution;
                        int32_t bin1End = (loop.end1 / resolution) + 1;  // +1 to include the full range
                        int32_t bin2Start = loop.start2 / resolution;
                        int32_t bin2End = (loop.end2 / resolution) + 1;  // +1 to include the full range
                        
                        // Calculate center bin positions
                        int32_t loopCenterX = (bin1Start + bin1End) / 2;
                        int32_t loopCenterY = (bin2Start + bin2End) / 2;
                        
                        // Calculate upper-left corner of window
                        int32_t windowStartX = loopCenterX - window_size;
                        int32_t windowStartY = loopCenterY - window_size;
                        
                        // Calculate relative position from upper-left corner
                        int relX = record.binX - windowStartX;
                        int relY = record.binY - windowStartY;
                        
                        // Add to matrix if within bounds
                        apaMatrix.add(relX, relY, record.value);
                    }
                }
            }
        }

        gzclose(file);
        
        // After processing all contacts, normalize the matrix
        std::vector<float> rowSums(apaMatrix.width, 0.0f);
        std::vector<float> colSums(apaMatrix.width, 0.0f);

        // Calculate row and column sums from coverage vectors
        for (const auto& loop : bedpe_entries) {
            // Convert to bin coordinates
            int32_t bin1Start = ((loop.start1 + loop.end1) / 2) / resolution - window_size;
            int32_t bin2Start = ((loop.start2 + loop.end2) / 2) / resolution - window_size;
            
            // Add local sums from coverage vectors
            coverage.addLocalSums(rowSums, loop.chrom1, bin1Start);
            coverage.addLocalSums(colSums, loop.chrom2, bin2Start);
        }

        // Scale sums by their averages
        APAMatrix::scaleByAverage(rowSums);
        APAMatrix::scaleByAverage(colSums);

        // Normalize the matrix
        apaMatrix.normalize(rowSums, colSums);
        
        // Save the matrix to file
        apaMatrix.save(output_file);
        
        return apaMatrix;

    } catch (...) {
        gzclose(file);
        throw;
    }
}
