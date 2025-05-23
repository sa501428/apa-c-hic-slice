#ifndef BEDPE_BUILDER_H
#define BEDPE_BUILDER_H

#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <cstdint>  // For int32_t
#include <cctype>   // For isdigit

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

    bool operator<(const BedpeEntry& other) const;
    bool operator==(const BedpeEntry& other) const;
};

class BedpeBuilder {
public:
    BedpeBuilder(const std::string& forward_bed, 
                 const std::string& reverse_bed,
                 long min_dist,
                 long max_dist,
                 bool isInter);  // true for inter-chromosomal, false for intra-chromosomal

    std::vector<BedpeEntry> buildBedpe();

private:
    std::string forward_bed_file;
    std::string reverse_bed_file;
    long min_genome_dist;
    long max_genome_dist;
    bool isInter;

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