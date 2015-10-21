#include <string>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <stdint.h>
#include <unistd.h>
#include "rolling.h"
#include "city.h"
#include "xxhash.h"
//#include "gzstream.h"
#include "Uncompress.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

namespace opt {
unsigned kmerLen = 30;
unsigned ibits = 64;
unsigned nhash = 5;
size_t fsize;
size_t mfsize;
size_t fbits;
size_t totitm=0;
bool fastq=false;
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

    if (line[0]!='>' && line[0]!='@') {
        cerr << "Target file is not in correct format!\n";
        exit(EXIT_FAILURE);
    }
    size_t totItm=0, uLen=0;

    if (line[0]=='>') {
        while (getline(faFile, line)) {
            if (line[0] != '>')
                uLen += line.length();
            else {
                if (uLen>=k)
                    totItm+=uLen-k+1;
                uLen = 0;
            }
        }
        if (uLen>=k) totItm+=uLen-k+1;
    }
    else if(line[0]=='@') {
        opt::fastq=true;
        do {
            getline(faFile, line);
            uLen = line.length();
            if(uLen>=k) totItm+=uLen-k+1;
            getline(faFile, line);
            getline(faFile, line);
        } while (getline(faFile, line));
    }

    cerr << "|totLen|=" << totItm << "\n";
    faFile.close();

    size_t totSize;
    for(totSize=1; totSize<=totItm; totSize*=2);
        
    //cerr << "|totSiz|=" << totSize << "\n";
    
    return totItm;
    //return totSize;
}

void loadSeqr(unsigned char *bloomFilter, const string& seq) {
    if (seq.size() < opt::kmerLen) return;
    string kmer = seq.substr(0,opt::kmerLen);
    uint64_t fhVal, rhVal, mhVal;
    mhVal= initHashes(kmer.c_str(), opt::kmerLen, fhVal, rhVal);
    for (size_t i = 0; i < seq.size() - opt::kmerLen + 1; i++) {
        size_t hashLoc = mhVal%opt::fsize;
        //size_t hashLoc = mhVal&opt::fsize-1;
        //#pragma omp atomic
        //bloomFilter[hashLoc/8] |= (1 << (7 - hashLoc % 8));
        __sync_or_and_fetch(&bloomFilter[hashLoc/8], (1 << (7 - hashLoc % 8)));
	mhVal = rollHashesRight(fhVal, rhVal, seq[i], seq[i+opt::kmerLen], opt::kmerLen);
    }
}

void loadSeqrd(unsigned char *bloomFilter, string& seq) {
    if (seq.size() < opt::kmerLen) return;
    string kmer = seq.substr(0,opt::kmerLen);
    uint64_t fhVal, rhVal, mhVal;
    mhVal= initHashes(kmer.c_str(), opt::kmerLen, fhVal, rhVal);
    size_t hashLoc = mhVal%opt::fsize;
    #pragma omp atomic
    bloomFilter[hashLoc/8] |= (1 << (7 - hashLoc % 8));
    
    for (size_t i = 1; i < seq.size() - opt::kmerLen + 1; i++) {
        mhVal = rollHashesRight(fhVal, rhVal, seq[i-1], seq[i+opt::kmerLen-1], opt::kmerLen);
        hashLoc = mhVal%opt::fsize;
        #pragma omp atomic
        bloomFilter[hashLoc/8] |= (1 << (7 - hashLoc % 8));
    }
}

void loadSeqrM(unsigned char *bloomFilter, const string& seq) {
    if (seq.size() < opt::kmerLen) return;
    uint64_t h[opt::nhash];
    for(unsigned j=0; j<opt::nhash;j++)
        h[j] = j ^ (opt::kmerLen * varSeed);
    string kmer = seq.substr(0,opt::kmerLen);
    uint64_t fhVal, rhVal, mhVal;
    mhVal= initHashes(kmer.c_str(), opt::kmerLen, fhVal, rhVal);
    for (size_t i = 0; i < seq.size() - opt::kmerLen + 1; i++) {
        for (unsigned j=0; j<opt::nhash; j++) {
            size_t hashLoc = (mhVal^h[j])%opt::fsize;
            //#pragma omp atomic
            //bloomFilter[hashLoc/8] |= (1 << (7 - hashLoc % 8));
            __sync_or_and_fetch(&bloomFilter[hashLoc/8], (1 << (7 - hashLoc % 8)));
        }
        mhVal = rollHashesRight(fhVal, rhVal, seq[i], seq[i+opt::kmerLen], opt::kmerLen);
    }
}

static const unsigned char b2r[256] = {
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //0
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //1
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //2
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //3
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', //4   'A' 'C' 'G'
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', //5   'T'
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', //6   'a' 'c' 'g'
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', //7   't'
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //8
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //9
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //10
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //11
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //12
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //13
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //14
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //15
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'
};

void getCanon(std::string &bMer) {
    int p=0, hLen=(opt::kmerLen-1)/2;
    while (bMer[p] == b2r[(unsigned char)bMer[opt::kmerLen-1-p]]) {
        ++p;
        if(p>=hLen) break;
    }
    if (bMer[p] > b2r[(unsigned char)bMer[opt::kmerLen-1-p]]) {
        for (int lIndex = p, rIndex = opt::kmerLen-1-p; lIndex<=rIndex; ++lIndex,--rIndex) {
            char tmp = b2r[(unsigned char)bMer[rIndex]];
            bMer[rIndex] = b2r[(unsigned char)bMer[lIndex]];
            bMer[lIndex] = tmp;
        }
    }
}

void loadSeqc(unsigned char *bloomFilter, string& seq) {
    if (seq.size() < opt::kmerLen) return;
    for (size_t i = 0; i < seq.size() - opt::kmerLen + 1; i++) {
        string kmer = seq.substr(i, opt::kmerLen);
        getCanon(kmer);
        uint64_t mhVal = CityHash64(kmer.c_str(), kmer.length());
        size_t hashLoc = mhVal%opt::fsize;
        //size_t hashLoc = mhVal&opt::fsize-1;
#pragma omp atomic
        bloomFilter[hashLoc/8] |= (1 << (7 - hashLoc % 8));
        
    }
}

uint64_t MurmurHash64A ( const void * key, int len, unsigned int seed ) {
    const uint64_t m = 0xc6a4a7935bd1e995;
    const int r = 47;
    
    uint64_t h = seed ^ (len * m);
    
    const uint64_t * data = (const uint64_t *)key;
    const uint64_t * end = data + (len/8);
    
    while (data != end)
    {
        uint64_t k = *data++;
        
        k *= m;
        k ^= k >> r;
        k *= m;
        
        h ^= k;
        h *= m;
    }
    
    const unsigned char * data2 = (const unsigned char*)data;
    
    switch(len & 7)
    {
        case 7:
            h ^= uint64_t(data2[6]) << 48;
        case 6:
            h ^= uint64_t(data2[5]) << 40;
        case 5:
            h ^= uint64_t(data2[4]) << 32;
        case 4:
            h ^= uint64_t(data2[3]) << 24;
        case 3:
            h ^= uint64_t(data2[2]) << 16;
        case 2:
            h ^= uint64_t(data2[1]) << 8;
        case 1:
            h ^= uint64_t(data2[0]);
            h *= m;
    };
    
    h ^= h >> r;
    h *= m;
    h ^= h >> r;
    
    return h;
}

void loadSeqm(unsigned char *bloomFilter, string& seq) {
    if (seq.size() < opt::kmerLen)
        return;
    
    for (size_t i = 0; i < seq.size() - opt::kmerLen + 1; i++) {
        string kmer = seq.substr(i, opt::kmerLen);
        getCanon(kmer);
        uint64_t mhVal = MurmurHash64A(kmer.c_str(), kmer.length(), 0);
        size_t hashLoc = mhVal%opt::fsize;
        //size_t hashLoc = mhVal&opt::fsize-1;
#pragma omp atomic
        bloomFilter[hashLoc/8] |= (1 << (7 - hashLoc % 8));
    }
}

void loadSeqmM(unsigned char *bloomFilter, string& seq) {
    if (seq.size() < opt::kmerLen) return;
    for (size_t i = 0; i < seq.size() - opt::kmerLen + 1; i++) {
        string kmer = seq.substr(i, opt::kmerLen);
        getCanon(kmer);
        for (unsigned j=0; j<opt::nhash; j++) {
            uint64_t mhVal = MurmurHash64A(kmer.c_str(), kmer.length(), j);
            size_t hashLoc = mhVal%opt::fsize;
            #pragma omp atomic
            bloomFilter[hashLoc/8] |= (1 << (7 - hashLoc % 8));
        }
    }
}

void loadSeqcM(unsigned char *bloomFilter, string& seq) {
    if (seq.size() < opt::kmerLen) return;
    for (size_t i = 0; i < seq.size() - opt::kmerLen + 1; i++) {
        string kmer = seq.substr(i, opt::kmerLen);
        getCanon(kmer);
        for (unsigned j=0; j<opt::nhash; j++) {
            uint64_t mhVal = CityHash64WithSeed(kmer.c_str(), kmer.length(), j);
            size_t hashLoc = mhVal%opt::fsize;
            #pragma omp atomic
            bloomFilter[hashLoc/8] |= (1 << (7 - hashLoc % 8));
        }
    }
}

void loadSeqx(unsigned char *bloomFilter, const string& seq) {
    if (seq.size() < opt::kmerLen) return;
    for (size_t i = 0; i < seq.size() - opt::kmerLen + 1; i++) {
        string kmer = seq.substr(i, opt::kmerLen);
        getCanon(kmer);
        uint64_t mhVal = XXH64(kmer.c_str(), kmer.length(), 0);
        size_t hashLoc = mhVal%opt::fsize;
        //size_t hashLoc = mhVal&opt::fsize-1;
#pragma omp atomic
        bloomFilter[hashLoc/8] |= (1 << (7 - hashLoc % 8));
    }
}

void loadSeqxM(unsigned char *bloomFilter, string& seq) {
    if (seq.size() < opt::kmerLen) return;
    for (size_t i = 0; i < seq.size() - opt::kmerLen + 1; i++) {
        string kmer = seq.substr(i, opt::kmerLen);
        getCanon(kmer);
        for (unsigned j=0; j<opt::nhash; j++) {
            uint64_t mhVal = XXH64(kmer.c_str(), kmer.length(), j);
            size_t hashLoc = mhVal%opt::fsize;
            #pragma omp atomic
            bloomFilter[hashLoc/8] |= (1 << (7 - hashLoc % 8));
        }
    }
}

// Get the hamming weight of x
inline unsigned popCnt(unsigned char x) {
    return ((0x876543210 >>
             (((0x4332322132212110 >> ((x & 0xF) << 2)) & 0xF) << 2)) >>
            ((0x4332322132212110 >> (((x & 0xF0) >> 2)) & 0xF) << 2))
           & 0xf;
}

void printBin(uint64_t n) {
    int count=0,count1=0;
    while (++count<=64) {
        if (n & ((uint64_t)1<<63)) {
            printf("1 ");
            ++count1;
        }
        else
            printf("0 ");
        
        n <<= 1;
    }
    printf(" count1=%d\n",count1);
}

void loadBf(unsigned char *bloomFilter, const char* faqFile) {
    for(int fIndex = 3; fIndex <= 3; fIndex++) {
        ifstream uFile(faqFile);
        bool good = true;
        #pragma omp parallel shared(bloomFilter,good)
        for(string line, hline; good;) {
            #pragma omp critical(uFile)
            {
                good = getline(uFile, hline);
                good = getline(uFile, line);
                good = getline(uFile, hline);
                good = getline(uFile, hline);
            }
            if(good) loadSeqrM(bloomFilter, line);
        }
        uFile.close();
    }
}

bool bfContain(const unsigned char *bloomFilter, const string& seq){
    if (seq.size() < opt::kmerLen) return false;
    uint64_t h[opt::nhash];
    for(unsigned j=0; j<opt::nhash;j++) h[j] = j ^ (opt::kmerLen * varSeed);
    string kmer = seq.substr(0,opt::kmerLen);
    uint64_t fhVal, rhVal, mhVal;
    mhVal= initHashes(kmer.c_str(), opt::kmerLen, fhVal, rhVal);
    for (size_t i = 0; i < seq.size() - opt::kmerLen + 1; i++) {
        bool inBf=true;
        for (unsigned j=0; j<opt::nhash; j++) {
            size_t hashLoc = (mhVal^h[j])%opt::mfsize;
            if((bloomFilter[hashLoc/8] & (1 << (7 - hashLoc % 8)))== 0) {
                inBf=false;
                break;
            }
        }
        if(inBf) return true;
        mhVal = rollHashesRight(fhVal, rhVal, seq[i], seq[i+opt::kmerLen], opt::kmerLen);
    }
    return false;
}

void loadpBf(unsigned char *petBf, const unsigned char *mpetBf, const char* faqFile1, const char* faqFile2) {
    ifstream uFile1(faqFile1), uFile2(faqFile2);
    bool good = true;
    size_t totKmer=0;
    #pragma omp parallel shared(good,totKmer)
    for(string line1, hline1, line2, hline2; good && totKmer<638400000;) {
        #pragma omp critical(uFile)
        {
            good = getline(uFile1, hline1); good = getline(uFile2, hline2);
            good = getline(uFile1, line1); good = getline(uFile2, line2);
            good = getline(uFile1, hline1); good = getline(uFile2, hline2);
            good = getline(uFile1, hline1); good = getline(uFile2, hline2);
        }
        if(good) {
            if(bfContain(mpetBf,line1) || bfContain(mpetBf,line2)) {
                loadSeqrM(petBf, line1);
                loadSeqrM(petBf, line2);
                size_t tmpLen=0;
                if(line1.length()>=opt::kmerLen) tmpLen=line1.length()-opt::kmerLen+1;
                if(line2.length()>=opt::kmerLen) tmpLen+=line2.length()-opt::kmerLen+1;
                #pragma omp atomic
                totKmer += tmpLen;
            }
        }
    }
    cerr << totKmer << "\n";
    uFile1.close();
    uFile2.close();
}


size_t popBf(const unsigned char *bloomFilter) {
    size_t i, popBF=0;
    #pragma omp parallel for reduction(+:popBF)
    for(i=0; i<(opt::fsize + 7)/8; i++)
        popBF = popBF + popCnt(bloomFilter[i]);
    return popBF;
}

int main(int argc, const char* argv[]) {
    
    opt::nhash= atoi(argv[1]);
    string hashName(argv[2]);
    //opt::fsize = opt::ibits*getInfo(argv[1], opt::kmerLen);
    //opt::fsize = opt::ibits*5998515280;

    double sTime = omp_get_wtime();
    //opt::fsize = opt::mfsize =  opt::ibits*638400000;
    opt::fsize = 80857600000;
    unsigned char *mpetBf = new unsigned char [(opt::fsize + 7)/8];
    for(size_t i=0; i<(opt::fsize + 7)/8; i++) mpetBf[i]=0;
    loadBf(mpetBf, argv[3]);
    cerr << "|popBF|=" << popBf(mpetBf) << " ";
    cerr << setprecision(4) << fixed << omp_get_wtime() - sTime << "\n";
    
    /*sTime = omp_get_wtime();
    //opt::fsize = opt::ibits*638400000;
    
    unsigned char *petBf = new unsigned char [(opt::fsize + 7)/8];
    for(size_t i=0; i<(opt::fsize + 7)/8; i++) petBf[i]=0;
    //loadpBf(petBf, mpetBf, argv[4], argv[5]);
    //loadpBf(petBf, mpetBf, argv[6], argv[7]);
    loadBf(petBf, argv[3]);
    cerr << "|popBF|=" << popBf(petBf) << " ";
    cerr << setprecision(4) << fixed << omp_get_wtime() - sTime << "\n";*/
    
    delete [] mpetBf;
    //delete [] petBf;
    return 0;
}
