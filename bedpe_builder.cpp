#include "bedpe_builder.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <cctype>

BedpeBuilder::BedpeBuilder(const std::string& forward_bed, 
                          const std::string& reverse_bed,
                          long min_dist,
                          long max_dist,
                          bool isInter)
    : forward_bed_file(forward_bed)
    , reverse_bed_file(reverse_bed)
    , min_genome_dist(min_dist)
    , max_genome_dist(max_dist)
    , isInter(isInter) {}

std::map<std::string, std::vector<BedEntry>> BedpeBuilder::loadBedFile(const std::string& filename) {
    std::map<std::string, std::vector<BedEntry>> bed_data;
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open input file: " + filename);
    }
    std::string line;
    
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string chrom;
        long start, end;
        iss >> chrom >> start >> end;
        
        BedEntry entry;
        entry.chrom = chrom;
        entry.gmid = (start + end) / 2;
        
        bed_data[entry.chrom].push_back(entry);
    }
    
    // Sort entries by position
    for (auto& pair : bed_data) {
        std::sort(pair.second.begin(), pair.second.end(),
                 [](const BedEntry& a, const BedEntry& b) {
                     return a.gmid < b.gmid;
                 });
    }
    
    return bed_data;
}

std::vector<BedpeEntry> BedpeBuilder::generateIntraChromosomal(
    const std::string& chrom,
    const std::vector<BedEntry>& forwards,
    const std::vector<BedEntry>& reverses) {
    
    std::vector<BedpeEntry> results;
    
    for (const auto& forward : forwards) {
        for (const auto& reverse : reverses) {
            long dist = abs(reverse.gmid - forward.gmid);
            
            if (dist > min_genome_dist && dist <= max_genome_dist) {
                BedpeEntry bedpe;
                bedpe.chrom1 = chrom;
                bedpe.gmid1 = forward.gmid;
                bedpe.chrom2 = chrom;
                bedpe.gmid2 = reverse.gmid;
                results.push_back(bedpe);
            }
        }
    }
    
    return results;
}

static bool isStandardChromosome(const std::string& chrom) {
    // Check if it starts with "chr"
    if (chrom.substr(0, 3) != "chr") {
        return false;
    }
    
    // Get the part after "chr"
    std::string num = chrom.substr(3);
    
    // Check if it's non-empty and contains only digits
    return !num.empty() && std::all_of(num.begin(), num.end(), ::isdigit);
}

std::vector<BedpeEntry> BedpeBuilder::generateInterChromosomal(
    const std::string& chrom1,
    const std::string& chrom2,
    const std::vector<BedEntry>& forwards,
    const std::vector<BedEntry>& reverses) {
    
    std::vector<BedpeEntry> results;
    
    // Skip if same chromosome or non-standard chromosomes
    if (chrom1 == chrom2 || 
        !isStandardChromosome(chrom1) || 
        !isStandardChromosome(chrom2)) {
        return results;
    }

    // Extract chromosome numbers for comparison
    int chr1_num = std::stoi(chrom1.substr(3));
    int chr2_num = std::stoi(chrom2.substr(3));

    if(chr1_num > chr2_num) {
        return results;
    }

    for (const auto& first : forwards) {
        for (const auto& second : reverses) {
            BedpeEntry bedpe;
            bedpe.chrom1 = chrom1;
            bedpe.gmid1 = first.gmid;
            bedpe.chrom2 = chrom2;
            bedpe.gmid2 = second.gmid;
            results.push_back(bedpe);
        }
    }
    
    return results;
}

std::vector<BedpeEntry> BedpeBuilder::buildBedpe() {
    std::cout << "Loading forward BED file: " << forward_bed_file << std::endl;
    auto forward_data = loadBedFile(forward_bed_file);
    std::cout << "Loading reverse BED file: " << reverse_bed_file << std::endl;
    auto reverse_data = loadBedFile(reverse_bed_file);
    
    std::vector<BedpeEntry> all_results;
    
    if (isInter) {
        std::cout << "Generating inter-chromosomal BEDPE entries..." << std::endl;
        for (const auto& forward_pair : forward_data) {
            for (const auto& reverse_pair : reverse_data) {
                if (forward_pair.first != reverse_pair.first) {
                    auto results = generateInterChromosomal(forward_pair.first,
                                                          reverse_pair.first,
                                                          forward_pair.second,
                                                          reverse_pair.second);
                    all_results.insert(all_results.end(), results.begin(), results.end());
                }
            }
        }
    } else {
        std::cout << "Generating intra-chromosomal BEDPE entries..." << std::endl;
        for (const auto& forward_pair : forward_data) {
            const std::string& chrom = forward_pair.first;
            if (reverse_data.count(chrom) > 0) {
                auto results = generateIntraChromosomal(chrom, 
                                                      forward_pair.second,
                                                      reverse_data[chrom]);
                all_results.insert(all_results.end(), results.begin(), results.end());
            }
        }
    }
    
    std::cout << "Sorting and removing duplicates..." << std::endl;
    std::sort(all_results.begin(), all_results.end());
    auto last = std::unique(all_results.begin(), all_results.end());
    all_results.erase(last, all_results.end());

    std::cout << "Generated " << all_results.size() << " unique BEDPE entries" << std::endl;
    return all_results;
}

bool BedpeEntry::operator<(const BedpeEntry& other) const {
    if (chrom1 != other.chrom1) return chrom1 < other.chrom1;
    if (chrom2 != other.chrom2) return chrom2 < other.chrom2;
    if (gmid1 != other.gmid1) return gmid1 < other.gmid1;
    return gmid2 < other.gmid2;
}

bool BedpeEntry::operator==(const BedpeEntry& other) const {
    return chrom1 == other.chrom1 &&
           gmid1 == other.gmid1 &&
           chrom2 == other.chrom2 &&
           gmid2 == other.gmid2;
} 