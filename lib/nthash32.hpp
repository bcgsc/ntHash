/*
 *
 * nthash.hpp
 * Author: Hamid Mohamadi
 * Genome Sciences Centre,
 * British Columbia Cancer Agency
 */

#ifndef NT_HASH32_H
#define NT_HASH32_H

#include <stdint.h>

// offset for the complement base in the random seeds table
const int cpOff = -20;

// shift for gerenerating multiple hash values
const int multiShift = 13;

// seed for gerenerating multiple hash values
static const uint32_t multiSeed = 0xa4547c9b;

// 64-bit random seeds corresponding to bases and their complements
static const uint32_t seedA = 0x72b2c14e;
static const uint32_t seedC = 0xb43a3ac1;
static const uint32_t seedG = 0x0d4ddc33;
static const uint32_t seedT = 0xcbc527bc;
static const uint32_t seedN = 0x00000000;

static const uint32_t vecA[32] = {
    0x72b2c14e,0xe565829c,0xcacb0539,0x95960a73,0x2b2c14e7,0x565829ce,0xacb0539c,0x5960a739,
    0xb2c14e72,0x65829ce5,0xcb0539ca,0x960a7395,0x2c14e72b,0x5829ce56,0xb0539cac,0x60a73959,
    0xc14e72b2,0x829ce565,0x0539cacb,0x0a739596,0x14e72b2c,0x29ce5658,0x539cacb0,0xa7395960,
    0x4e72b2c1,0x9ce56582,0x39cacb05,0x7395960a,0xe72b2c14,0xce565829,0x9cacb053,0x0395960a7
};

static const uint32_t vecC[32] = {
    0xb43a3ac1,0x68747583,0xd0e8eb06,0xa1d1d60d,0x43a3ac1b,0x87475836,0x0e8eb06d,0x1d1d60da,
    0x3a3ac1b4,0x74758368,0xe8eb06d0,0xd1d60da1,0xa3ac1b43,0x47583687,0x8eb06d0e,0x1d60da1d,
    0x3ac1b43a,0x75836874,0xeb06d0e8,0xd60da1d1,0xac1b43a3,0x58368747,0xb06d0e8e,0x60da1d1d,
    0xc1b43a3a,0x83687475,0x06d0e8eb,0x0da1d1d6,0x1b43a3ac,0x36874758,0x6d0e8eb0,0xda1d1d60
};

static const uint32_t vecG[32] = {
    0x0d4ddc33,0x1a9bb866,0x353770cc,0x6a6ee198,0xd4ddc330,0xa9bb8661,0x53770cc3,0xa6ee1986,
    0x4ddc330d,0x9bb8661a,0x3770cc35,0x6ee1986a,0xddc330d4,0xbb8661a9,0x770cc353,0xee1986a6,
    0xdc330d4d,0xb8661a9b,0x70cc3537,0xe1986a6e,0xc330d4dd,0x8661a9bb,0x0cc35377,0x1986a6ee,
    0x330d4ddc,0x661a9bb8,0xcc353770,0x986a6ee1,0x30d4ddc3,0x61a9bb86,0xc353770c,0x86a6ee19
};

static const uint32_t vecT[32] = {
    0xcbc527bc,0x978a4f79,0x2f149ef3,0x5e293de6,0xbc527bcc,0x78a4f799,0xf149ef32,0xe293de65,
    0xc527bccb,0x8a4f7997,0x149ef32f,0x293de65e,0x527bccbc,0xa4f79978,0x49ef32f1,0x93de65e2,
    0x27bccbc5,0x4f79978a,0x9ef32f14,0x3de65e29,0x7bccbc52,0xf79978a4,0xef32f149,0xde65e293,
    0xbccbc527,0x79978a4f,0xf32f149e,0xe65e293d,0xccbc527b,0x9978a4f7,0x32f149ef,0x65e293de
};

static const uint32_t vecN[32] = {
    seedN,seedN,seedN,seedN,seedN,seedN,seedN,seedN,
    seedN,seedN,seedN,seedN,seedN,seedN,seedN,seedN,
    seedN,seedN,seedN,seedN,seedN,seedN,seedN,seedN,
    seedN,seedN,seedN,seedN,seedN,seedN,seedN,seedN
};

static const uint32_t *msTab[256] = {
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 0..7
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 8..15
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 16..23
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 24..31
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 32..39
    vecN, vecN, vecN, vecN, vecN, vecT, vecN, vecG, // 40..47
    vecN, vecN, vecN, vecC, vecN, vecN, vecN, vecN, // 48..55
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 56..63
    vecA, vecA, vecN, vecC, vecN, vecN, vecN, vecG, // 64..71
    vecN, vecN, vecN, vecN, vecN, vecT, vecN, vecG, // 72..79
    vecN, vecN, vecN, vecC, vecT, vecN, vecN, vecN, // 80..87
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 88..95
    vecA, vecA, vecN, vecC, vecN, vecN, vecN, vecG, // 96..103
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 104..111
    vecN, vecN, vecN, vecN, vecT, vecN, vecN, vecN, // 112..119
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 120..127
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 128..135
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 136..143
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 144..151
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 152..159
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 160..167
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 168..175
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 176..183
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 184..191
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 192..199
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 200..207
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 208..215
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 216..223
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 224..231
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 232..239
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 240..247
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN  // 248..255
};

static const uint32_t seedTab[256] = {
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 0..7
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 8..15
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 16..23
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 24..31
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 32..39
    seedN, seedN, seedN, seedN, seedN, seedT, seedN, seedG, // 40..47
    seedN, seedN, seedN, seedC, seedN, seedN, seedN, seedN, // 48..55
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 56..63
    seedA, seedA, seedN, seedC, seedN, seedN, seedN, seedG, // 64..71
    seedN, seedN, seedN, seedN, seedN, seedT, seedN, seedG, // 72..79
    seedN, seedN, seedN, seedC, seedT, seedN, seedN, seedN, // 80..87
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 88..95
    seedA, seedA, seedN, seedC, seedN, seedN, seedN, seedG, // 96..103
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 104..111
    seedN, seedN, seedN, seedN, seedT, seedN, seedN, seedN, // 112..119
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 120..127
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 128..135
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 136..143
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 144..151
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 152..159
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 160..167
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 168..175
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 176..183
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 184..191
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 192..199
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 200..207
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 208..215
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 216..223
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 224..231
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 232..239
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 240..247
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN  // 248..255
};

// rotate "v" to the left by "s" positions
inline uint32_t rol(const uint32_t v, const int s) {
    return (v << s) | (v >> (32 - s));
}

// rotate "v" to the right by "s" positions
inline uint32_t ror(const uint32_t v, const int s) {
    return (v >> s) | (v << (32 - s));
}

// forward-strand hash value of the base kmer, i.e. fhval(kmer_0)
inline uint32_t getFhval(const char * kmerSeq, const unsigned k) {
    uint32_t hVal=0;
    for(unsigned i=0; i<k; i++)
        hVal ^= rol(seedTab[(unsigned char)kmerSeq[i]], k-1-i);
    return hVal;
}

// reverse-strand hash value of the base kmer, i.e. rhval(kmer_0)
inline uint32_t getRhval(const char * kmerSeq, const unsigned k) {
    uint32_t hVal=0;
    for(unsigned i=0; i<k; i++)
        hVal ^= rol(seedTab[(unsigned char)kmerSeq[i]+cpOff], i);
    return hVal;
}

// ntHash basic function, i.e. ntBase
inline uint32_t NT32(const char * kmerSeq, const unsigned k) {
    return getFhval(kmerSeq, k);
}

// ntHash for sliding k-mers
inline uint32_t NT32(const uint32_t fhVal, const unsigned char charOut, const unsigned char charIn, const unsigned k) {
    return(rol(fhVal, 1) ^ rol(seedTab[charOut], k) ^ seedTab[charIn]);
}

// ntHash with seeding option
inline uint32_t NT32(const char * kmerSeq, const unsigned k, const unsigned seed) {
    uint32_t hVal=0;
    for(unsigned i=0; i<k; i++)
        hVal ^= rol(seedTab[(unsigned char)kmerSeq[i]], k-1-i);
    hVal *= seed ^ k * multiSeed;
    hVal ^= hVal >> multiShift;
    return hVal;
}

// canonical ntHash, base function
inline uint32_t NTC32(const char * kmerSeq, const unsigned k) {
    return (getFhval(kmerSeq, k) ^ getRhval(kmerSeq, k));
}

// canonical ntHash
inline uint32_t NTC32(const char * kmerSeq, const unsigned k, uint32_t& fhVal, uint32_t& rhVal) {
    fhVal = getFhval(kmerSeq, k);
    rhVal = getRhval(kmerSeq, k);
    return (rhVal^fhVal);
}

// canonical ntHash for sliding k-mers
inline uint32_t NTC32(uint32_t& fhVal, uint32_t& rhVal, const unsigned char charOut, const unsigned char charIn, const unsigned k) {
    fhVal = rol(fhVal, 1) ^ rol(seedTab[charOut], k) ^ seedTab[charIn];
    rhVal = ror(rhVal, 1) ^ ror(seedTab[charOut+cpOff], 1) ^ rol(seedTab[charIn+cpOff], k-1);
    return rhVal ^ fhVal;
}


// canonical ntHash with seeding option
inline uint32_t NTC32(const char * kmerSeq, const unsigned k, const unsigned seed) {
    uint32_t hVal=0;
    for(unsigned i=0; i<k; i++)
        hVal ^= rol(seedTab[(unsigned char)kmerSeq[i]], k-1-i);
    hVal *= seed ^ k * multiSeed;
    hVal ^= hVal >> multiShift;
    return hVal;
}

// ntHash using pre-computed seed matrix (msTab)
inline uint32_t NTP32(const char * kmerSeq, const unsigned k) {
    uint32_t hVal=0;
    for(unsigned i=0; i<k; i++)
        hVal ^= msTab[(unsigned char)kmerSeq[i]][(k-1-i)%32];
    return hVal;
}

// ntHash for sliding k-mers using pre-computed seed matrix (msTab)
inline uint32_t NTP32(const uint32_t fhVal, const unsigned char charOut, const unsigned char charIn, const unsigned k) {
    return(rol(fhVal, 1) ^ msTab[charOut][k%32] ^ msTab[charIn][0]);
}

// ntHash with seeding option using pre-computed seed matrix (msTab)
inline uint32_t NTP32(const char * kmerSeq, const unsigned k, const unsigned seed) {
    uint32_t hVal=0;
    for(unsigned i=0; i<k; i++)
        hVal ^= msTab[(unsigned char)kmerSeq[i]][(k-1-i)%32];
    hVal *= seed ^ k * multiSeed;
    hVal ^= hVal >> multiShift;
    return hVal;
}

// canonical ntHash using pre-computed seed matrix (msTab)
inline uint32_t NTPC32(const char * kmerSeq, const unsigned k) {
    uint32_t fhVal=0, rhVal=0;
    for(unsigned i=0; i<k; i++) {
        fhVal ^= msTab[(unsigned char)kmerSeq[i]][(k-1-i)%32];
        rhVal ^= msTab[(unsigned char)kmerSeq[i]+cpOff][i%32];
    }
    return (rhVal^fhVal);
}

// canonical ntHash using pre-computed seed matrix (msTab)
inline uint32_t NTPC32(const char * kmerSeq, const unsigned k, uint32_t& fhVal, uint32_t& rhVal) {
    fhVal=0, rhVal=0;
    for(unsigned i=0; i<k; i++) {
        fhVal ^= msTab[(unsigned char)kmerSeq[i]][(k-1-i)%32];
        rhVal ^= msTab[(unsigned char)kmerSeq[i]+cpOff][i%32];
    }
    return (rhVal^fhVal);
}

// canonical ntHash for sliding k-mers using pre-computed seed matrix (msTab)
inline uint32_t NTPC32(uint32_t& fhVal, uint32_t& rhVal, const unsigned char charOut, const unsigned char charIn, const unsigned k) {
    fhVal = rol(fhVal, 1) ^ rol(seedTab[charOut], k) ^ seedTab[charIn];
    rhVal = ror(rhVal, 1) ^ ror(seedTab[charOut+cpOff], 1) ^ rol(seedTab[charIn+cpOff], k-1);
    return (rhVal^fhVal);
}

// canonical ntHash with seeding option using pre-computed seed matrix (msTab)
inline uint32_t NTPC32(const char * kmerSeq, const unsigned k, const unsigned seed) {
    uint32_t fhVal=0, rhVal=0, hVal=0;
    for(unsigned i=0; i<k; i++) {
        fhVal ^= msTab[(unsigned char)kmerSeq[i]][(k-1-i)%32];
        rhVal ^= msTab[(unsigned char)kmerSeq[i]+cpOff][i%32];
    }
    hVal = rhVal^fhVal;
    hVal *= seed ^ k * multiSeed;
    hVal ^= hVal >> multiShift;
    return hVal;
}

// multi-hash version of ntHash
void NTM32(const char * kmerSeq, const unsigned k, const unsigned m, uint32_t *hVal) {
    uint32_t bVal=0, tVal=0;
    hVal[0] = bVal = NTP32(kmerSeq, k);
    for(unsigned i=1; i<m; i++) {
        tVal = bVal * (i ^ k * multiSeed);
        tVal ^= tVal >> multiShift;
        hVal[i] =  tVal;
    }
}

// multi-hash version of ntHash for sliding k-mers
void NTM32(const unsigned char charOut, const unsigned char charIn, const unsigned k, const unsigned m, uint32_t *hVal) {
    uint32_t bVal=0, tVal=0;
    hVal[0] = bVal = rol(hVal[0], 1) ^ msTab[charOut][k%32] ^ msTab[charIn][0];
    for(unsigned i=1; i<m; i++) {
        tVal = bVal * (i ^ k * multiSeed);
        tVal ^= tVal >> multiShift;
        hVal[i] =  tVal;
    }
}

// canonical multi-hash version of ntHash
void NTMC32(const char * kmerSeq, const unsigned k, const unsigned m, uint32_t *hVal, uint32_t& fhVal, uint32_t& rhVal) {
    uint32_t bVal=0, tVal=0;
    hVal[0] = bVal = NTPC32(kmerSeq, k, fhVal, rhVal);
    for(unsigned i=1; i<m; i++) {
        tVal = bVal * (i ^ k * multiSeed);
        tVal ^= tVal >> multiShift;
        hVal[i] =  tVal;
    }
}

// canonical multi-hash version of ntHash for sliding k-mers
void NTMC32(uint32_t& fhVal, uint32_t& rhVal, const unsigned char charOut, const unsigned char charIn, const unsigned k, const unsigned m, uint32_t *hVal) {
    uint32_t bVal=0, tVal=0;
    hVal[0] = bVal = NTPC32(fhVal, rhVal, charOut, charIn, k);
    for(unsigned i=1; i<m; i++) {
        tVal = bVal * (i ^ k * multiSeed);
        tVal ^= tVal >> multiShift;
        hVal[i] =  tVal;
    }
}

#endif
