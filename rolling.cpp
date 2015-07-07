#include <string>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <stdint.h>
#include <unistd.h>
#include "rolling.h"

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
