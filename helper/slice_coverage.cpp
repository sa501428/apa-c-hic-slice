#include <zlib.h>
#include <stdexcept>
#include <cstring>
#include <cmath>
#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "coverage_vectors.h"  

void calculateCoverageAndDump(
    const std::string& slice_file,
    const std::string& out_csv,
    bool verbose = false)
{
    // 1 Open slice file (gz or raw) ---
    FILE* raw = fopen(slice_file.c_str(), "rb");
    gzFile gz = nullptr;
    bool is_gz = false;
    if (!raw) {
        gz = gzopen(slice_file.c_str(), "rb");
        if (!gz) throw std::runtime_error("Could not open file: " + slice_file);
        is_gz = true;
    }

    // helper lambdas for read
    auto read_exact = [&](void* buf, size_t sz){
        if (is_gz) {
            if (gzread(gz, buf, sz) != (int)sz)
                throw std::runtime_error("Unexpected EOF");
        } else {
            if (fread(buf, 1, sz, raw) != sz)
                throw std::runtime_error("Unexpected EOF");
        }
    };

    // 2 Read & check magic ---
    char magic[8];
    read_exact(magic, 8);
    if (std::strncmp(magic, "HICSLICE", 8) != 0)
        throw std::runtime_error("Not a HICSLICE file");

    // 3 Read resolution ---
    int32_t resolution = 0;
    read_exact(&resolution, sizeof(resolution));
    if (resolution <= 0)
        throw std::runtime_error("Invalid resolution");

    if (verbose) std::cout << "Resolution = " << resolution << "\n";

    // 4 Read chromosome map ---
    int32_t nChr = 0;
    read_exact(&nChr, sizeof(nChr));
    if (nChr <= 0) throw std::runtime_error("No chromosomes listed");
    std::map<int16_t,std::string> key2name;
    for (int i = 0; i < nChr; ++i) {
        int32_t len = 0;
        read_exact(&len, sizeof(len));
        std::string name(len, '\0');
        read_exact(&name[0], len);
        int16_t key = 0;
        read_exact(&key, sizeof(key));
        key2name[key] = name;
        if (verbose) std::cout << "  chr " << key << " = " << name << "\n";
    }

    // 5 Init coverage accumulator ---
    CoverageVectors coverage(resolution);

    struct Record {
        int16_t chr1; int32_t bin1;
        int16_t chr2; int32_t bin2;
        float   val;
    } rec;

    // 6 Read all contact records ---
    size_t recSize = sizeof(rec);
    int64_t count = 0;
    while (true) {
        try {
            read_exact(&rec, recSize);
        } catch (...) {
            break;  // EOF
        }
        ++count;
        if (std::isnan(rec.val) || std::isinf(rec.val) || rec.val <= 0) 
            continue;

        // add coverage: but don't double‐count exact diagonal contacts
        if (rec.chr1 == rec.chr2 && rec.bin1 == rec.bin2) {
            coverage.add(rec.chr1, rec.bin1, rec.val);
        } else {
            coverage.add(rec.chr1, rec.bin1, rec.val);
            coverage.add(rec.chr2, rec.bin2, rec.val);
        }

        if (count % 10000000 == 0 && verbose) {
            std::cout << "Processed " << count << " records...\n";
        }
    }
    if (verbose) std::cout << "Total records: " << count << "\n";

    // 7 Write out CSV
    std::ofstream ofs(out_csv);
    if (!ofs) throw std::runtime_error("Could not open output: " + out_csv);
    ofs << "Chromosome\tBin\tCoverage\n";
    for (const auto& chr_pair : coverage.getVectors()) {
        auto it = key2name.find(chr_pair.first);
        std::string chrName = (it!=key2name.end() ? it->second : std::to_string(chr_pair.first));
        for (const auto& bin_pair : chr_pair.second) {
            ofs << chrName << '\t' << bin_pair.first << '\t' 
                << std::fixed << std::setprecision(3) 
                << bin_pair.second << "\n";
        }
    }
    ofs.close();

    if (is_gz) gzclose(gz);
    else       fclose(raw);

    if (verbose) std::cout << "Coverage written to " << out_csv << "\n";
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: calc_coverage <slice_file> <out_csv>\n";
        return 1;
    }
    try {
        calculateCoverageAndDump(argv[1], argv[2], true);
    } catch (std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
