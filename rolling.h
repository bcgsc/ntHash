#ifndef ROLLING_HASH_H
#define ROLLING_HASH_H

#include <string>
#include <stdint.h>

// cpOff is the offset for the complement base in the random seeds table
static const int cpOff = -20;

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

// Rotate "value" to the left by "shift" positions
static inline uint64_t rol(uint64_t value, int shift) {
    return (value << shift) | (value >> (64 - shift));
}

// Rotate "value" to the right by "shift" positions
static inline uint64_t ror(uint64_t value, int shift) {
    return (value >> shift) | (value << (64 - shift));
}

// Get forward-strand hash value of the base kmer, i.e. fhval(kmer_0)
static inline uint64_t getFhval(const std::string &kmerSeq) {
    uint64_t hVal=0;
    unsigned k = kmerSeq.length();
    for(unsigned i=0; i<k; i++)
        hVal ^= rol(seedTab[(unsigned)kmerSeq[i]], k-1-i);
    return hVal;
}

// Get reverse-strand hash value of the base kmer, i.e. rhval(kmer_0)
static inline uint64_t getRhval(const std::string &kmerSeq) {
    uint64_t hVal=0;
    unsigned k = kmerSeq.length();
    for(unsigned i=0; i<k; i++)
        hVal ^= rol(seedTab[(unsigned)kmerSeq[i]+cpOff], i);
    return hVal;
}

static inline uint64_t initHashes(const std::string &kmerSeq,
    uint64_t& fhVal, uint64_t& rhVal)
{
    fhVal = getFhval(kmerSeq);
    rhVal = getRhval(kmerSeq);
    return (rhVal<fhVal)? rhVal : fhVal;
}

static inline uint64_t rollHashesRight(uint64_t& fhVal, uint64_t& rhVal,
    unsigned char charOut, unsigned char charIn, unsigned k)
{
    fhVal = rol(fhVal, 1) ^ rol(seedTab[charOut], k) ^ seedTab[charIn];
    rhVal = ror(rhVal, 1) ^ ror(seedTab[charOut+cpOff], 1) ^ rol(seedTab[charIn+cpOff], k-1);
    return (rhVal<fhVal)? rhVal : fhVal;
}

#endif
