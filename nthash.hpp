/*
 *
 * nthash.hpp
 * Author: Hamid Mohamadi
 * Genome Sciences Centre,
 * British Columbia Cancer Agency
 */

#ifndef NT_HASH_H
#define NT_HASH_H

#include <stdint.h>

// offset for the complement base in the random seeds table
const int cpOff = -20;

// shift for gerenerating multiple hash values
const int multiShift = 27;

// seed for gerenerating multiple hash values
static const uint64_t multiSeed = 0x90b45d39fb6da1fa;

// 64-bit random seeds corresponding to bases and their complements
static const uint64_t seedA = 0x359d7bbf07501017;
static const uint64_t seedC = 0x4326d1e1c820effb;
static const uint64_t seedG = 0xec612e52b7efe400;
static const uint64_t seedT = 0x9ada840c789f1bec;
static const uint64_t seedN = 0x0000000000000000;

static const uint64_t vecA[64] = {
    0x359d7bbf07501017,0x6b3af77e0ea0202e,0xd675eefc1d40405c,0xacebddf83a8080b9,0x59d7bbf075010173,0xb3af77e0ea0202e6,0x675eefc1d40405cd,0xcebddf83a8080b9a,
    0x9d7bbf0750101735,0x3af77e0ea0202e6b,0x75eefc1d40405cd6,0xebddf83a8080b9ac,0xd7bbf07501017359,0xaf77e0ea0202e6b3,0x5eefc1d40405cd67,0xbddf83a8080b9ace,
    0x7bbf07501017359d,0xf77e0ea0202e6b3a,0xeefc1d40405cd675,0xddf83a8080b9aceb,0xbbf07501017359d7,0x77e0ea0202e6b3af,0xefc1d40405cd675e,0xdf83a8080b9acebd,
    0xbf07501017359d7b,0x7e0ea0202e6b3af7,0xfc1d40405cd675ee,0xf83a8080b9acebdd,0xf07501017359d7bb,0xe0ea0202e6b3af77,0xc1d40405cd675eef,0x83a8080b9acebddf,
    0x07501017359d7bbf,0x0ea0202e6b3af77e,0x1d40405cd675eefc,0x3a8080b9acebddf8,0x7501017359d7bbf0,0xea0202e6b3af77e0,0xd40405cd675eefc1,0xa8080b9acebddf83,
    0x501017359d7bbf07,0xa0202e6b3af77e0e,0x40405cd675eefc1d,0x8080b9acebddf83a,0x01017359d7bbf075,0x0202e6b3af77e0ea,0x0405cd675eefc1d4,0x080b9acebddf83a8,
    0x1017359d7bbf0750,0x202e6b3af77e0ea0,0x405cd675eefc1d40,0x80b9acebddf83a80,0x017359d7bbf07501,0x02e6b3af77e0ea02,0x05cd675eefc1d404,0x0b9acebddf83a808,
    0x17359d7bbf075010,0x2e6b3af77e0ea020,0x5cd675eefc1d4040,0xb9acebddf83a8080,0x7359d7bbf0750101,0xe6b3af77e0ea0202,0xcd675eefc1d40405,0x9acebddf83a8080b
};

static const uint64_t vecC[64] = {
    0x4326d1e1c820effb,0x864da3c39041dff6,0x0c9b47872083bfed,0x19368f0e41077fda,0x326d1e1c820effb4,0x64da3c39041dff68,0xc9b47872083bfed0,0x9368f0e41077fda1,
    0x26d1e1c820effb43,0x4da3c39041dff686,0x9b47872083bfed0c,0x368f0e41077fda19,0x6d1e1c820effb432,0xda3c39041dff6864,0xb47872083bfed0c9,0x68f0e41077fda193,
    0xd1e1c820effb4326,0xa3c39041dff6864d,0x47872083bfed0c9b,0x8f0e41077fda1936,0x1e1c820effb4326d,0x3c39041dff6864da,0x7872083bfed0c9b4,0xf0e41077fda19368,
    0xe1c820effb4326d1,0xc39041dff6864da3,0x872083bfed0c9b47,0x0e41077fda19368f,0x1c820effb4326d1e,0x39041dff6864da3c,0x72083bfed0c9b478,0xe41077fda19368f0,
    0xc820effb4326d1e1,0x9041dff6864da3c3,0x2083bfed0c9b4787,0x41077fda19368f0e,0x820effb4326d1e1c,0x041dff6864da3c39,0x083bfed0c9b47872,0x1077fda19368f0e4,
    0x20effb4326d1e1c8,0x41dff6864da3c390,0x83bfed0c9b478720,0x077fda19368f0e41,0x0effb4326d1e1c82,0x1dff6864da3c3904,0x3bfed0c9b4787208,0x77fda19368f0e410,
    0xeffb4326d1e1c820,0xdff6864da3c39041,0xbfed0c9b47872083,0x7fda19368f0e4107,0xffb4326d1e1c820e,0xff6864da3c39041d,0xfed0c9b47872083b,0xfda19368f0e41077,
    0xfb4326d1e1c820ef,0xf6864da3c39041df,0xed0c9b47872083bf,0xda19368f0e41077f,0xb4326d1e1c820eff,0x6864da3c39041dff,0xd0c9b47872083bfe,0xa19368f0e41077fd
};

static const uint64_t vecG[64] = {
    0xec612e52b7efe400,0xd8c25ca56fdfc801,0xb184b94adfbf9003,0x63097295bf7f2007,0xc612e52b7efe400e,0x8c25ca56fdfc801d,0x184b94adfbf9003b,0x3097295bf7f20076,
    0x612e52b7efe400ec,0xc25ca56fdfc801d8,0x84b94adfbf9003b1,0x097295bf7f200763,0x12e52b7efe400ec6,0x25ca56fdfc801d8c,0x4b94adfbf9003b18,0x97295bf7f2007630,
    0x2e52b7efe400ec61,0x5ca56fdfc801d8c2,0xb94adfbf9003b184,0x7295bf7f20076309,0xe52b7efe400ec612,0xca56fdfc801d8c25,0x94adfbf9003b184b,0x295bf7f200763097,
    0x52b7efe400ec612e,0xa56fdfc801d8c25c,0x4adfbf9003b184b9,0x95bf7f2007630972,0x2b7efe400ec612e5,0x56fdfc801d8c25ca,0xadfbf9003b184b94,0x5bf7f20076309729,
    0xb7efe400ec612e52,0x6fdfc801d8c25ca5,0xdfbf9003b184b94a,0xbf7f200763097295,0x7efe400ec612e52b,0xfdfc801d8c25ca56,0xfbf9003b184b94ad,0xf7f200763097295b,
    0xefe400ec612e52b7,0xdfc801d8c25ca56f,0xbf9003b184b94adf,0x7f200763097295bf,0xfe400ec612e52b7e,0xfc801d8c25ca56fd,0xf9003b184b94adfb,0xf200763097295bf7,
    0xe400ec612e52b7ef,0xc801d8c25ca56fdf,0x9003b184b94adfbf,0x200763097295bf7f,0x400ec612e52b7efe,0x801d8c25ca56fdfc,0x003b184b94adfbf9,0x00763097295bf7f2,
    0x00ec612e52b7efe4,0x01d8c25ca56fdfc8,0x03b184b94adfbf90,0x0763097295bf7f20,0x0ec612e52b7efe40,0x1d8c25ca56fdfc80,0x3b184b94adfbf900,0x763097295bf7f200
};

static const uint64_t vecT[64] = {
    0x9ada840c789f1bec,0x35b50818f13e37d9,0x6b6a1031e27c6fb2,0xd6d42063c4f8df64,0xada840c789f1bec9,0x5b50818f13e37d93,0xb6a1031e27c6fb26,0x6d42063c4f8df64d,
    0xda840c789f1bec9a,0xb50818f13e37d935,0x6a1031e27c6fb26b,0xd42063c4f8df64d6,0xa840c789f1bec9ad,0x50818f13e37d935b,0xa1031e27c6fb26b6,0x42063c4f8df64d6d,
    0x840c789f1bec9ada,0x0818f13e37d935b5,0x1031e27c6fb26b6a,0x2063c4f8df64d6d4,0x40c789f1bec9ada8,0x818f13e37d935b50,0x031e27c6fb26b6a1,0x063c4f8df64d6d42,
    0x0c789f1bec9ada84,0x18f13e37d935b508,0x31e27c6fb26b6a10,0x63c4f8df64d6d420,0xc789f1bec9ada840,0x8f13e37d935b5081,0x1e27c6fb26b6a103,0x3c4f8df64d6d4206,
    0x789f1bec9ada840c,0xf13e37d935b50818,0xe27c6fb26b6a1031,0xc4f8df64d6d42063,0x89f1bec9ada840c7,0x13e37d935b50818f,0x27c6fb26b6a1031e,0x4f8df64d6d42063c,
    0x9f1bec9ada840c78,0x3e37d935b50818f1,0x7c6fb26b6a1031e2,0xf8df64d6d42063c4,0xf1bec9ada840c789,0xe37d935b50818f13,0xc6fb26b6a1031e27,0x8df64d6d42063c4f,
    0x1bec9ada840c789f,0x37d935b50818f13e,0x6fb26b6a1031e27c,0xdf64d6d42063c4f8,0xbec9ada840c789f1,0x7d935b50818f13e3,0xfb26b6a1031e27c6,0xf64d6d42063c4f8d,
    0xec9ada840c789f1b,0xd935b50818f13e37,0xb26b6a1031e27c6f,0x64d6d42063c4f8df,0xc9ada840c789f1be,0x935b50818f13e37d,0x26b6a1031e27c6fb,0x4d6d42063c4f8df6
};

static const uint64_t vecN[64] = {
    seedN,seedN,seedN,seedN,seedN,seedN,seedN,seedN,
    seedN,seedN,seedN,seedN,seedN,seedN,seedN,seedN,
    seedN,seedN,seedN,seedN,seedN,seedN,seedN,seedN,
    seedN,seedN,seedN,seedN,seedN,seedN,seedN,seedN,
    seedN,seedN,seedN,seedN,seedN,seedN,seedN,seedN,
    seedN,seedN,seedN,seedN,seedN,seedN,seedN,seedN,
    seedN,seedN,seedN,seedN,seedN,seedN,seedN,seedN,
    seedN,seedN,seedN,seedN,seedN,seedN,seedN,seedN
};

static const uint64_t *msTab[256] = {
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

static const uint64_t seedTab[256] = {
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
inline uint64_t rol(const uint64_t v, const int s) {
    return (v << s) | (v >> (64 - s));
}

// rotate "v" to the right by "s" positions
inline uint64_t ror(const uint64_t v, const int s) {
    return (v >> s) | (v << (64 - s));
}

// forward-strand hash value of the base kmer, i.e. fhval(kmer_0)
inline uint64_t getFhval(const char * kmerSeq, const unsigned k) {
    uint64_t hVal=0;
    for(unsigned i=0; i<k; i++)
        hVal ^= rol(seedTab[(unsigned char)kmerSeq[i]], k-1-i);
    return hVal;
}

// reverse-strand hash value of the base kmer, i.e. rhval(kmer_0)
inline uint64_t getRhval(const char * kmerSeq, const unsigned k) {
    uint64_t hVal=0;
    for(unsigned i=0; i<k; i++)
        hVal ^= rol(seedTab[(unsigned char)kmerSeq[i]+cpOff], i);
    return hVal;
}

// ntHash basic function, i.e. ntBase
inline uint64_t NT64(const char * kmerSeq, const unsigned k) {
    return getFhval(kmerSeq, k);
}

// ntHash for sliding k-mers
inline uint64_t NT64(const uint64_t fhVal, const unsigned char charOut, const unsigned char charIn, const unsigned k) {
    return(rol(fhVal, 1) ^ rol(seedTab[charOut], k) ^ seedTab[charIn]);
}

// ntHash with seeding option
inline uint64_t NT64(const char * kmerSeq, const unsigned k, const unsigned seed) {
    uint64_t hVal=0;
    for(unsigned i=0; i<k; i++)
        hVal ^= rol(seedTab[(unsigned char)kmerSeq[i]], k-1-i);
    hVal *= seed ^ k * multiSeed;
    hVal ^= hVal >> multiShift;
    return hVal;
}

// canonical ntHash, base function
inline uint64_t NTC64(const char * kmerSeq, const unsigned k) {
    return (getFhval(kmerSeq, k) ^ getRhval(kmerSeq, k));
}

// canonical ntHash
inline uint64_t NTC64(const char * kmerSeq, const unsigned k, uint64_t& fhVal, uint64_t& rhVal) {
    fhVal = getFhval(kmerSeq, k);
    rhVal = getRhval(kmerSeq, k);
    return (rhVal^fhVal);
}

// canonical ntHash for sliding k-mers
inline uint64_t NTC64(uint64_t& fhVal, uint64_t& rhVal, const unsigned char charOut, const unsigned char charIn, const unsigned k) {
    fhVal = ror(fhVal, 1) ^ ror(seedTab[charOut], 1) ^ rol(seedTab[charIn], k-1);
    rhVal = rol(rhVal, 1) ^ rol(seedTab[charOut+cpOff], k) ^ seedTab[charIn+cpOff];
    return rhVal ^ fhVal;
}

// canonical ntHash with seeding option
inline uint64_t NTC64(const char * kmerSeq, const unsigned k, const unsigned seed) {
    uint64_t hVal=0;
    for(unsigned i=0; i<k; i++)
        hVal ^= rol(seedTab[(unsigned char)kmerSeq[i]], k-1-i);
    hVal *= seed ^ k * multiSeed;
    hVal ^= hVal >> multiShift;
    return hVal;
}

// ntHash using pre-computed seed matrix (msTab)
inline uint64_t NTP64(const char * kmerSeq, const unsigned k) {
    uint64_t hVal=0;
    for(unsigned i=0; i<k; i++)
        hVal ^= msTab[(unsigned char)kmerSeq[i]][(k-1-i)%64];
    return hVal;
}

// ntHash for sliding k-mers using pre-computed seed matrix (msTab)
inline uint64_t NTP64(const uint64_t fhVal, const unsigned char charOut, const unsigned char charIn, const unsigned k) {
    return(rol(fhVal, 1) ^ msTab[charOut][k%64] ^ msTab[charIn][0]);
}

// ntHash with seeding option using pre-computed seed matrix (msTab)
inline uint64_t NTP64(const char * kmerSeq, const unsigned k, const unsigned seed) {
    uint64_t hVal=0;
    for(unsigned i=0; i<k; i++)
        hVal ^= msTab[(unsigned char)kmerSeq[i]][(k-1-i)%64];
    hVal *= seed ^ k * multiSeed;
    hVal ^= hVal >> multiShift;
    return hVal;
}

// canonical ntHash using pre-computed seed matrix (msTab)
inline uint64_t NTPC64(const char * kmerSeq, const unsigned k) {
    uint64_t fhVal=0, rhVal=0;
    for(unsigned i=0; i<k; i++) {
        fhVal ^= msTab[(unsigned char)kmerSeq[i]][(k-1-i)%64];
        rhVal ^= msTab[(unsigned char)kmerSeq[i]+cpOff][i%64];
    }
    return (rhVal^fhVal);
}

// canonical ntHash using pre-computed seed matrix (msTab)
inline uint64_t NTPC64(const char * kmerSeq, const unsigned k, uint64_t& fhVal, uint64_t& rhVal) {
    fhVal=0, rhVal=0;
    for(unsigned i=0; i<k; i++) {
        fhVal ^= msTab[(unsigned char)kmerSeq[i]][(k-1-i)%64];
        rhVal ^= msTab[(unsigned char)kmerSeq[i]+cpOff][i%64];
    }
    return (rhVal^fhVal);
}

// canonical ntHash for sliding k-mers using pre-computed seed matrix (msTab)
inline uint64_t NTPC64(uint64_t& fhVal, uint64_t& rhVal, const unsigned char charOut, const unsigned char charIn, const unsigned k) {
    fhVal = rol(fhVal, 1) ^ rol(seedTab[charOut], k) ^ seedTab[charIn];
    rhVal = ror(rhVal, 1) ^ ror(seedTab[charOut+cpOff], 1) ^ rol(seedTab[charIn+cpOff], k-1);
    return (rhVal^fhVal);
}

// canonical ntHash with seeding option using pre-computed seed matrix (msTab)
inline uint64_t NTPC64(const char * kmerSeq, const unsigned k, const unsigned seed) {
    uint64_t fhVal=0, rhVal=0, hVal=0;
    for(unsigned i=0; i<k; i++) {
        fhVal ^= msTab[(unsigned char)kmerSeq[i]][(k-1-i)%64];
        rhVal ^= msTab[(unsigned char)kmerSeq[i]+cpOff][i%64];
    }
    hVal = rhVal^fhVal;
    hVal *= seed ^ k * multiSeed;
    hVal ^= hVal >> multiShift;
    return hVal;
}

// multi-hash version of ntHash
void NTM64(const char * kmerSeq, const unsigned k, const unsigned m, uint64_t *hVal) {
    uint64_t bVal=0, tVal=0;
    hVal[0] = bVal = NTP64(kmerSeq, k);
    for(unsigned i=1; i<m; i++) {
        tVal = bVal * (i ^ k * multiSeed);
        tVal ^= tVal >> multiShift;
        hVal[i] =  tVal;
    }
}

// multi-hash version of ntHash for sliding k-mers
void NTM64(const unsigned char charOut, const unsigned char charIn, const unsigned k, const unsigned m, uint64_t *hVal) {
    uint64_t bVal=0, tVal=0;
    hVal[0] = bVal = rol(hVal[0], 1) ^ msTab[charOut][k%64] ^ msTab[charIn][0];
    for(unsigned i=1; i<m; i++) {
        tVal = bVal * (i ^ k * multiSeed);
        tVal ^= tVal >> multiShift;
        hVal[i] =  tVal;
    }
}

// canonical multi-hash version of ntHash
void NTMC64(const char * kmerSeq, const unsigned k, const unsigned m, uint64_t *hVal, uint64_t& fhVal, uint64_t& rhVal) {
    uint64_t bVal=0, tVal=0;
    hVal[0] = bVal = NTPC64(kmerSeq, k, fhVal, rhVal);
    for(unsigned i=1; i<m; i++) {
        tVal = bVal * (i ^ k * multiSeed);
        tVal ^= tVal >> multiShift;
        hVal[i] =  tVal;
    }
}

// canonical multi-hash version of ntHash for sliding k-mers
void NTMC64(uint64_t& fhVal, uint64_t& rhVal, const unsigned char charOut, const unsigned char charIn, const unsigned k, const unsigned m, uint64_t *hVal) {
    uint64_t bVal=0, tVal=0;
    hVal[0] = bVal = NTPC64(fhVal, rhVal, charOut, charIn, k);
    for(unsigned i=1; i<m; i++) {
        tVal = bVal * (i ^ k * multiSeed);
        tVal ^= tVal >> multiShift;
        hVal[i] =  tVal;
    }
}

#endif
