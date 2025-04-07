#ifndef APA_H
#define APA_H

#include "bedpe_builder.h"
#include <string>
#include <vector>

// Structure to hold binned regions for faster lookup
struct BinRegion {
    int32_t start;
    int32_t end;
};

float processSliceFile(const std::string& slice_file, 
                      const std::vector<BedpeEntry>& bedpe_entries);

#endif 