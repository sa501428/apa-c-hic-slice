#ifndef BEDPE_BUILDER_H
#define BEDPE_BUILDER_H

#include <string>
#include <vector>
#include <map>

struct BedEntry {
    std::string chrom;
    long start;
    long end;
    std::string name;
    
    long getMid() const {
        return (start + end) / 2;
    }
};

struct BedpeEntry {
    std::string chrom1;
    long start1;
    long end1;
    std::string chrom2;
    long start2;
    long end2;
    std::string name = ".";
    std::string score = ".";
    std::string strand1 = ".";
    std::string strand2 = ".";
    std::string color = "0,0,0";

    // Helper function to normalize chromosome order (only for inter-chromosomal)
    void normalizeOrder() {
        // Only normalize if chromosomes are different
        if (chrom1 != chrom2) {
            // Extract chromosome numbers
            int chr1_num = std::stoi(chrom1.substr(3));  // skip "chr" prefix
            int chr2_num = std::stoi(chrom2.substr(3));

            if (chr1_num > chr2_num) {
                std::swap(chrom1, chrom2);
                std::swap(start1, start2);
                std::swap(end1, end2);
                std::swap(strand1, strand2);
            }
        }
    }

    // Modified comparison operator for ascending order with chr1 <= chr2 convention
    bool operator<(const BedpeEntry& other) const {
        if (chrom1 != other.chrom1) return chrom1 < other.chrom1;
        if (chrom2 != other.chrom2) return chrom2 < other.chrom2;
        if (start1 != other.start1) return start1 < other.start1;
        return start2 < other.start2;
    }

    // Add equality operator for std::unique
    bool operator==(const BedpeEntry& other) const {
        return chrom1 == other.chrom1 &&
               start1 == other.start1 &&
               end1 == other.end1 &&
               chrom2 == other.chrom2 &&
               start2 == other.start2 &&
               end2 == other.end2;
    }
};

class BedpeBuilder {
public:
    BedpeBuilder(const std::string& forward_bed, 
                 const std::string& reverse_bed,
                 long min_dist,
                 long max_dist,
                 bool make_intra = true,
                 bool make_inter = false);

    void buildBedpe(const std::string& output_file);

    template<typename Callback>
    void generateBedpe(Callback processEntry) {
        auto forward_data = loadBedFile(forward_bed_file);
        auto reverse_data = loadBedFile(reverse_bed_file);
        
        // Generate intra-chromosomal pairs
        if (make_intra) {
            for (const auto& forward_pair : forward_data) {
                const std::string& chrom = forward_pair.first;
                if (reverse_data.count(chrom) > 0) {
                    auto results = generateIntraChromosomal(chrom, 
                                                          forward_pair.second,
                                                          reverse_data[chrom]);
                    for (const auto& entry : results) {
                        processEntry(entry);
                    }
                }
            }
        }
        
        // Generate inter-chromosomal pairs
        if (make_inter) {
            for (const auto& forward_pair : forward_data) {
                for (const auto& reverse_pair : reverse_data) {
                    if (forward_pair.first != reverse_pair.first) {
                        auto results = generateInterChromosomal(forward_pair.first,
                                                              reverse_pair.first,
                                                              forward_pair.second,
                                                              reverse_pair.second);
                        for (const auto& entry : results) {
                            processEntry(entry);
                        }
                    }
                }
            }
        }
    }

private:
    std::string forward_bed_file;
    std::string reverse_bed_file;
    long min_genome_dist;
    long max_genome_dist;
    bool make_intra;
    bool make_inter;

    std::map<std::string, std::vector<BedEntry>> loadBedFile(const std::string& filename);
    std::vector<BedpeEntry> generateIntraChromosomal(const std::string& chrom,
                                                    const std::vector<BedEntry>& forwards,
                                                    const std::vector<BedEntry>& reverses);
    std::vector<BedpeEntry> generateInterChromosomal(const std::string& chrom1,
                                                    const std::string& chrom2,
                                                    const std::vector<BedEntry>& forwards,
                                                    const std::vector<BedEntry>& reverses);
};

#endif 