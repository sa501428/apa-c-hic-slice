// Filters a HiC Slice file for contacts whose first and second bin fall within regions specified in a BED file.
// Usage: bed_slice_filter_SQ <input.slice> <regions.bed> <output.slice> [--gz]
//
// Compile with: g++ -std=c++11 -lz -o bed_slice_filter_SQ bed_slice_filter_squared.cpp

#include <zlib.h>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <iostream>
#include <sstream>

static const int MAGIC_LEN = 8;
static const char MAGIC[MAGIC_LEN+1] = "HICSLICE";

struct ContactRecord {
    int16_t chr1Key;
    int32_t binX;
    int16_t chr2Key;
    int32_t binY;
    float   value;
};

void usage(const char* prog) {
    std::cerr << "Usage: " << prog << " <input.slice> <regions.bed> <output.slice> [--gz]\n";
    std::exit(EXIT_FAILURE);
}

void gzReadOrDie(gzFile f, void* buf, size_t size) {
    int bytes = gzread(f, buf, size);
    if (bytes != (int)size) {
        std::cerr << "Error: unexpected EOF or read error in compressed input\n";
        std::exit(EXIT_FAILURE);
    }
}
void freadOrDie(FILE* f, void* buf, size_t size) {
    size_t got = fread(buf, 1, size, f);
    if (got != size) {
        std::cerr << "Error: unexpected EOF or read error in input\n";
        std::exit(EXIT_FAILURE);
    }
}
void gzWriteOrDie(gzFile f, const void* buf, size_t size) {
    int written = gzwrite(f, buf, size);
    if (written != (int)size) {
        std::cerr << "Error: write error in compressed output\n";
        std::exit(EXIT_FAILURE);
    }
}
void fwriteOrDie(FILE* f, const void* buf, size_t size) {
    size_t written = fwrite(buf, 1, size, f);
    if (written != size) {
        std::cerr << "Error: write error in output\n";
        std::exit(EXIT_FAILURE);
    }
}

int main(int argc, char** argv) {
    if (argc < 4 || argc > 5) usage(argv[0]);

    std::string inPath  = argv[1];
    std::string bedPath = argv[2];
    std::string outPath = argv[3];
    bool outCompressed  = (argc == 5 && std::string(argv[4]) == "--gz");

    FILE*  rawIn = fopen(inPath.c_str(), "rb");
    gzFile gzIn  = nullptr;
    bool   inCompressed = false;
    if (!rawIn) {
        gzIn = gzopen(inPath.c_str(), "rb");
        if (!gzIn) {
            std::cerr << "Error: could not open input " << inPath << "\n";
            return EXIT_FAILURE;
        }
        inCompressed = true;
    }

    char magicBuf[MAGIC_LEN];
    if (inCompressed) gzReadOrDie(gzIn, magicBuf, MAGIC_LEN);
    else              freadOrDie(rawIn, magicBuf, MAGIC_LEN);
    if (std::memcmp(magicBuf, MAGIC, MAGIC_LEN) != 0) {
        std::cerr << "Error: invalid slice file (magic mismatch)\n";
        return EXIT_FAILURE;
    }

    int32_t resolution = 0;
    if (inCompressed) gzReadOrDie(gzIn, &resolution, sizeof(resolution));
    else              freadOrDie(rawIn, &resolution, sizeof(resolution));
    if (resolution <= 0) {
        std::cerr << "Error: invalid resolution " << resolution << "\n";
        return EXIT_FAILURE;
    }

    int32_t numChroms = 0;
    if (inCompressed) gzReadOrDie(gzIn, &numChroms, sizeof(numChroms));
    else              freadOrDie(rawIn, &numChroms, sizeof(numChroms));
    if (numChroms <= 0) {
        std::cerr << "Error: invalid chromosome count " << numChroms << "\n";
        return EXIT_FAILURE;
    }

    std::unordered_map<int16_t,std::string> keyToName;
    std::unordered_map<std::string,int16_t> nameToKey;
    for (int i = 0; i < numChroms; i++) {
        int32_t nameLen = 0;
        if (inCompressed) gzReadOrDie(gzIn, &nameLen, sizeof(nameLen));
        else              freadOrDie(rawIn, &nameLen, sizeof(nameLen));

        std::vector<char> buf(nameLen+1, 0);
        if (inCompressed) gzReadOrDie(gzIn, buf.data(), nameLen);
        else              freadOrDie(rawIn, buf.data(), nameLen);

        int16_t key;
        if (inCompressed) gzReadOrDie(gzIn, &key, sizeof(key));
        else              freadOrDie(rawIn, &key, sizeof(key));

        std::string chrName(buf.data());
        keyToName[key]    = chrName;
        nameToKey[chrName] = key;
    }

    std::unordered_map<int16_t, std::unordered_set<int32_t>> bedBins;
    {
        std::ifstream bedIn(bedPath);
        if (!bedIn) {
            std::cerr << "Error: cannot open BED file " << bedPath << "\n";
            return EXIT_FAILURE;
        }
        std::string line;
        while (std::getline(bedIn, line)) {
            if (line.empty() || line[0]=='#') continue;
            std::istringstream iss(line);
            std::string chr;
            int32_t start, end;
            if (!(iss >> chr >> start >> end)) continue;
            auto it = nameToKey.find(chr);
            if (it == nameToKey.end()) {
                std::cerr << "Warning: unknown chr " << chr << " in BED, skipping\n";
                continue;
            }
            int16_t key = it->second;
            int32_t bstart = start / resolution;
            int32_t bend   = (end-1) / resolution;
            for (int32_t b = bstart; b <= bend; ++b) {
                bedBins[key].insert(b);
            }
        }
    }

    FILE*  rawOut = nullptr;
    gzFile gzOut  = nullptr;
    if (outCompressed) {
        gzOut = gzopen(outPath.c_str(), "wb");
        if (!gzOut) {
            std::cerr << "Error: could not open compressed output " << outPath << "\n";
            return EXIT_FAILURE;
        }
        gzWriteOrDie(gzOut, MAGIC, MAGIC_LEN);
        gzWriteOrDie(gzOut, &resolution, sizeof(resolution));
        gzWriteOrDie(gzOut, &numChroms, sizeof(numChroms));
        for (auto& kv : keyToName) {
            int32_t L = kv.second.size();
            gzWriteOrDie(gzOut, &L, sizeof(L));
            gzWriteOrDie(gzOut, kv.second.c_str(), L);
            gzWriteOrDie(gzOut, &kv.first, sizeof(kv.first));
        }
    } else {
        rawOut = fopen(outPath.c_str(), "wb");
        if (!rawOut) {
            std::cerr << "Error: could not open output " << outPath << "\n";
            return EXIT_FAILURE;
        }
        fwriteOrDie(rawOut, MAGIC, MAGIC_LEN);
        fwriteOrDie(rawOut, &resolution, sizeof(resolution));
        fwriteOrDie(rawOut, &numChroms, sizeof(numChroms));
        for (auto& kv : keyToName) {
            int32_t L = kv.second.size();
            fwriteOrDie(rawOut, &L, sizeof(L));
            fwriteOrDie(rawOut, kv.second.c_str(), L);
            fwriteOrDie(rawOut, &kv.first, sizeof(kv.first));
        }
    }

    ContactRecord rec{};
    while (true) {
        // read next record
        if (inCompressed) {
            int r = gzread(gzIn, &rec, sizeof(rec));
            if (r != (int)sizeof(rec)) break;
        } else {
            if (fread(&rec, sizeof(rec), 1, rawIn) != 1) break;
        }

        // lookup both ends
        auto it1 = bedBins.find(rec.chr1Key);
        auto it2 = bedBins.find(rec.chr2Key);
        bool in1 = (it1 != bedBins.end() && it1->second.count(rec.binX));
        bool in2 = (it2 != bedBins.end() && it2->second.count(rec.binY));

        // keep only if BOTH ends overlap
        if (!(in1 && in2)) 
            continue;

        // write surviving record
        if (outCompressed) gzWriteOrDie(gzOut, &rec, sizeof(rec));
        else              fwriteOrDie(rawOut, &rec, sizeof(rec));
    }

    if (inCompressed) gzclose(gzIn);
    else              fclose(rawIn);
    if (outCompressed) gzclose(gzOut);
    else              fclose(rawOut);

    return EXIT_SUCCESS;
}
