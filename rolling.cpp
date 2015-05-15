#include <string>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <stdint.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

namespace opt {
unsigned kmerLen = 28;
unsigned ibits = 8;
unsigned nhash = 1;
size_t fsize;
}

// cpOff is the offset for the complement base in the random seeds table
const int cpOff = -20;

// Table of random seeds for rolling hash function created by seedgen.cpp
static const uint64_t seedTab[256] = {
    0, 0, 0, 0, 0, 0, 0, 0, // 0..7
    0, 0, 0, 0, 0, 0, 0, 0, // 8..15
    0, 0, 0, 0, 0, 0, 0, 0, // 16..23
    0, 0, 0, 0, 0, 0, 0, 0, // 24..31
    0, 0, 0, 0, 0, 0, 0, 0, // 32..39
    0, 0, 0, 0, 0, 0x9be9b3d202074c81, 0, 0x90b45d39fb6da1fa, // 40..47
    0, 0, 0, 0x664ecc8f7cd05f14, 0, 0, 0, 0, // 48..55
    0, 0, 0, 0, 0, 0, 0, 0, // 56..63
    0x6d13226485bab26f, 0x6d13226485bab26f, 0, 0x664ecc8f7cd05f14, 0, 0, 0, 0x90b45d39fb6da1fa, // 64..71
    0, 0, 0, 0, 0, 0, 0, 0, // 72..79
    0, 0, 0, 0, 0x9be9b3d202074c81, 0, 0, 0, // 80..87
    0, 0, 0, 0, 0, 0, 0, 0, // 88..95
    0, 0, 0, 0, 0, 0, 0, 0, // 96..103
    0, 0, 0, 0, 0, 0, 0, 0, // 104..111
    0, 0, 0, 0, 0, 0, 0, 0, // 112..119
    0, 0, 0, 0, 0, 0, 0, 0, // 120..127
    0, 0, 0, 0, 0, 0, 0, 0, // 128..135
    0, 0, 0, 0, 0, 0, 0, 0, // 136..143
    0, 0, 0, 0, 0, 0, 0, 0, // 144..151
    0, 0, 0, 0, 0, 0, 0, 0, // 152..159
    0, 0, 0, 0, 0, 0, 0, 0, // 160..167
    0, 0, 0, 0, 0, 0, 0, 0, // 168..175
    0, 0, 0, 0, 0, 0, 0, 0, // 176..183
    0, 0, 0, 0, 0, 0, 0, 0, // 184..191
    0, 0, 0, 0, 0, 0, 0, 0, // 192..199
    0, 0, 0, 0, 0, 0, 0, 0, // 200..207
    0, 0, 0, 0, 0, 0, 0, 0, // 208..215
    0, 0, 0, 0, 0, 0, 0, 0, // 216..223
    0, 0, 0, 0, 0, 0, 0, 0, // 224..231
    0, 0, 0, 0, 0, 0, 0, 0, // 232..239
    0, 0, 0, 0, 0, 0, 0, 0, // 240..247
    0, 0, 0, 0, 0, 0, 0, 0  // 248..255
};

// Get total memory on the system
size_t getTotalSystemMemory() {
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    return pages * page_size;
}

// Get total number of k-mers in the input file
size_t getInfo(const char *aName, unsigned k) {
    string line;
    ifstream faFile(aName);

    getline(faFile, line);
    if (line[0]!='>') {
        cerr << "Target file is not in correct format!\n";
        exit(EXIT_FAILURE);
    }
    size_t totItm=0, uLen=0;
    while (getline(faFile, line)) {
        if (line[0] != '>')
            uLen += line.length();
        else {
            if (uLen>=k)
                totItm+=uLen-k+1;
            uLen = 0;
        }
    }
    if (uLen>=k)
        totItm+=uLen-k+1;

    cerr << "|totLen|=" << totItm << "\n";
    faFile.close();
    return totItm;
}

// Rotate "value" to the left by "shift" positions
static inline uint64_t rol(uint64_t value, int shift) {
    return (value << shift) | (value >> (64 - shift));
}

// Rotate "value" to the right by "shift" positions
static inline uint64_t ror(uint64_t value, int shift) {
    return (value >> shift) | (value << (64 - shift));
}

// Get forward-strand hash value of the base kmer, i.e. fhval(kmer_0)
static inline uint64_t getFhval(const string &kmerSeq) {
    uint64_t hVal=0;
    for(unsigned i=0; i<opt::kmerLen; i++)
        hVal ^= rol(seedTab[(unsigned)kmerSeq[i]], opt::kmerLen-1-i);
    return hVal;
}

// Get reverse-strand hash value of the base kmer, i.e. rhval(kmer_0)
static inline uint64_t getRhval(const string &kmerSeq) {
    uint64_t hVal=0;
    for(unsigned i=0; i<opt::kmerLen; i++)
        hVal ^= rol(seedTab[(unsigned)kmerSeq[i]+cpOff], i);
    return hVal;
}

// Parallel loading sequences into the lock-free Bloom filter by omp atomic (or //Built-in atomic )
void loadSeq(unsigned char *bloomFilter, string& seq) {
    if (seq.size() < opt::kmerLen)
        return;
    transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    string kmer = seq.substr(0,opt::kmerLen);

    uint64_t fhVal = getFhval(kmer);
    uint64_t rhVal = getRhval(kmer);
    uint64_t mhVal = (rhVal<fhVal)? rhVal: fhVal;
    size_t hashLoc = mhVal%opt::fsize;
    #pragma omp atomic
    bloomFilter[hashLoc/8] |= (1 << (7 - hashLoc % 8));
    //__sync_or_and_fetch(&bloomFilter[hashLoc/8], (1 << (7 - hashLoc % 8)));
    for (size_t i = 1; i < seq.size() - opt::kmerLen + 1; i++) {
        fhVal = rol(fhVal, 1) ^ rol(seedTab[(unsigned)seq[i-1]], opt::kmerLen) ^ seedTab[(unsigned)seq[i-1+opt::kmerLen]];
        rhVal = ror(rhVal, 1) ^ ror(seedTab[(unsigned)seq[i-1]+cpOff], 1) ^ rol(seedTab[(unsigned)seq[i-1+opt::kmerLen]+cpOff], opt::kmerLen-1);
        mhVal = (rhVal<fhVal)? rhVal: fhVal;
        hashLoc = mhVal%opt::fsize;
        #pragma omp atomic
        bloomFilter[hashLoc/8] |= (1 << (7 - hashLoc % 8));
        //__sync_or_and_fetch(&bloomFilter[hashLoc/8], (1 << (7 - hashLoc % 8)));
    }
}

// Get the hamming weight of x
int popCnt(unsigned char x) {
    return ((0x876543210 >>
             (((0x4332322132212110 >> ((x & 0xF) << 2)) & 0xF) << 2)) >>
            ((0x4332322132212110 >> (((x & 0xF0) >> 2)) & 0xF) << 2))
           & 0xf;
}

int main(int argc, const char* argv[]) {
    opt::fsize = opt::ibits*getInfo(argv[1], opt::kmerLen);

    unsigned char *bloomFilter = new unsigned char [(opt::fsize + 7)/8];
    for(size_t i=0; i<(opt::fsize + 7)/8; i++) bloomFilter[i]=0;
    ifstream uFile(argv[1]);

#ifdef _OPENMP
    double sTime = omp_get_wtime();
#else
    clock_t sTime = clock();
#endif

    bool good = true;
    #pragma omp parallel shared(bloomFilter,good)
    for(string line; good;) {
        #pragma omp critical(uFile)
        {
            good = getline(uFile, line);
            good = getline(uFile, line);
        }
        if(good)
            loadSeq(bloomFilter, line);
    }

#ifdef _OPENMP
    cerr << "Running time in sec: " << omp_get_wtime() - sTime << "\n";
#else
    cerr << "Running time in sec: " << (double)(clock() - sTime)/CLOCKS_PER_SEC << "\n";
#endif

    delete [] bloomFilter;
    uFile.close();

    return 0;
}
