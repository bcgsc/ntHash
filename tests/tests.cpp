#include "nthash/nthash.hpp"

#include <chrono>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <queue>
#include <random>
#include <stack>
#include <string>

#define PRINT_TEST_NAME(TEST_NAME)                                             \
  std::cerr << __FILE__ << ": Testing " << TEST_NAME << std::endl;

#define TEST_ASSERT(x)                                                         \
  if (!(x)) {                                                                  \
    std::cerr << __FILE__ ":" << __LINE__ << ":TEST_ASSERT: " #x " is false"   \
              << std::endl;                                                    \
    std::exit(EXIT_FAILURE);                                                   \
  }

#define TEST_ASSERT_RELATIONAL(x, y, op)                                       \
  if (!((x)op(y))) {                                                           \
    std::cerr << __FILE__ ":" << __LINE__                                      \
              << ":TEST_ASSERT_RELATIONAL: " #x " " #op " " #y << '\n'         \
              << #x " = " << x << '\n'                                         \
              << #y " = " << y << std::endl;                                   \
    std::exit(EXIT_FAILURE);                                                   \
  }

#define TEST_ASSERT_EQ(x, y) TEST_ASSERT_RELATIONAL(x, y, ==)
#define TEST_ASSERT_NE(x, y) TEST_ASSERT_RELATIONAL(x, y, !=)
#define TEST_ASSERT_GE(x, y) TEST_ASSERT_RELATIONAL(x, y, >=)
#define TEST_ASSERT_GT(x, y) TEST_ASSERT_RELATIONAL(x, y, >)
#define TEST_ASSERT_LE(x, y) TEST_ASSERT_RELATIONAL(x, y, <=)
#define TEST_ASSERT_LT(x, y) TEST_ASSERT_RELATIONAL(x, y, <)

#define TEST_ASSERT_ARRAY_EQ(x, y, size)                                       \
  for (unsigned i = 0; i < size; i++) {                                        \
    TEST_ASSERT_EQ(x[i], y[i])                                                 \
  }

int
main()
{

  {
    PRINT_TEST_NAME("k-mer hash values")

    std::string seq = "ACATGCATGCA";
    const unsigned k = 5;
    const unsigned h = 3;

    const std::vector<std::array<uint64_t, h>> hashes = {
      { 0xf59ecb45f0e22b9c, 0x4969c33ac240c129, 0x688d616f0d7e08c3 },
      { 0x38cc00f940aebdae, 0xab7e1b110e086fc6, 0x11a1818bcfdd553 },
      { 0x603a48c5a11c794a, 0xe66016e61816b9c4, 0xc5b13cb146996ffe }
    };

    nthash::NtHash nthash(seq, h, k);
    nthash::BlindNtHash ntblind(seq.substr(0, k), h, k);

    for (const auto& h_vals : hashes) {
      nthash.roll();
      TEST_ASSERT_ARRAY_EQ(h_vals, nthash.hashes(), h);
      ntblind.roll(seq[ntblind.get_pos() + 1]);
      TEST_ASSERT_ARRAY_EQ(h_vals, ntblind.hashes(), h);
    }
  }

  {
    PRINT_TEST_NAME("k-mer rolling")

    std::string seq = "AGTCAGTC";
    unsigned h = 3;
    unsigned k = 4;

    nthash::NtHash nthash(seq, h, k);
    std::vector<uint64_t*> hashes;

    while (nthash.roll()) {
      uint64_t* h_vals = new uint64_t[h];
      std::copy(nthash.hashes(), nthash.hashes() + h, h_vals);
      hashes.push_back(h_vals);
    }

    TEST_ASSERT_EQ(hashes.size(), seq.length() - k + 1);

    // check same hash value for identical k-mers (first and last)
    TEST_ASSERT_ARRAY_EQ(hashes[0], hashes[hashes.size() - 1], h);
  }

  {
    PRINT_TEST_NAME("k-mer rolling vs ntbase hash values")

    std::string seq = "ACGTACACTGGACTGAGTCT";

    nthash::NtHash nthash(seq, 3, seq.size() - 2);
    /* 18-mers of kmer*/
    std::string kmer1 = "ACGTACACTGGACTGAGT";
    std::string kmer2 = "CGTACACTGGACTGAGTC";
    std::string kmer3 = "GTACACTGGACTGAGTCT";

    nthash::NtHash nthash_vector[] = {
      nthash::NtHash(kmer1, nthash.get_hash_num(), kmer1.size()),
      nthash::NtHash(kmer2, nthash.get_hash_num(), kmer2.size()),
      nthash::NtHash(kmer3, nthash.get_hash_num(), kmer3.size())
    };

    size_t i;
    for (i = 0; nthash.roll() && nthash_vector[i].roll(); ++i) {
      for (size_t j = 0; j < nthash.get_hash_num(); ++j) {
        TEST_ASSERT_EQ(nthash.hashes()[j], nthash_vector[i].hashes()[j]);
      }
    }
    TEST_ASSERT_EQ(i, 3);
  }

  {
    PRINT_TEST_NAME("canonical hashing")

    std::string seq_f = "ACGTACACTGGACTGAGTCT";
    std::string seq_r = "AGACTCAGTCCAGTGTACGT";
    unsigned h = 3;

    nthash::NtHash nthash_f(seq_f, h, seq_f.size());
    nthash::NtHash nthash_r(seq_r, h, seq_r.size());

    nthash_f.roll();
    nthash_r.roll();
    TEST_ASSERT_EQ(nthash_f.get_hash_num(), nthash_r.get_hash_num())
    TEST_ASSERT_ARRAY_EQ(nthash_f.hashes(), nthash_r.hashes(), h)
  }

  {
    PRINT_TEST_NAME("k-mer back rolling")

    std::string seq = "ACTAGCTG";
    unsigned h = 3;
    unsigned k = 5;

    nthash::NtHash nthash(seq, h, k);
    std::stack<uint64_t*> hashes;

    while (nthash.roll()) {
      uint64_t* h_vals = new uint64_t[h];
      std::copy(nthash.hashes(), nthash.hashes() + h, h_vals);
      hashes.push(h_vals);
    }

    TEST_ASSERT_EQ(hashes.size(), seq.length() - k + 1)

    do {
      TEST_ASSERT_ARRAY_EQ(nthash.hashes(), hashes.top(), h)
      hashes.pop();
    } while (nthash.roll_back());
  }

  {
    PRINT_TEST_NAME("k-mer peeking")

    std::string seq = "ACTGATCAG";
    unsigned h = 3;
    unsigned k = 6;

    nthash::NtHash nthash(seq, h, k);
    nthash.roll();

    size_t steps = 3;
    while (steps--) {
      nthash.peek();
      uint64_t* h_peek = new uint64_t[h];
      std::copy(nthash.hashes(), nthash.hashes() + h, h_peek);
      nthash.peek(seq[nthash.get_pos() + k]);
      TEST_ASSERT_ARRAY_EQ(nthash.hashes(), h_peek, h);
      nthash.roll();
      TEST_ASSERT_ARRAY_EQ(nthash.hashes(), h_peek, h);
    }
  }

  {
    PRINT_TEST_NAME("skipping Ns")

    std::string seq = "ACGTACACTGGACTGAGTCT";
    std::string seq_with_ns = seq;

    TEST_ASSERT_GE(seq_with_ns.size(), 10)
    seq_with_ns[seq_with_ns.size() / 2] = 'N';
    seq_with_ns[seq_with_ns.size() / 2 + 1] = 'N';
    unsigned k = (seq.size() - 2) / 2 - 1;
    nthash::NtHash nthash(seq_with_ns, 3, k);

    std::vector<uint64_t> positions;
    for (size_t i = 0; i < seq_with_ns.size() / 2 - k + 1; i++) {
      positions.push_back(i);
    }
    for (size_t i = seq_with_ns.size() / 2 + 2; i < seq_with_ns.size() - k + 1;
         i++) {
      positions.push_back(i);
    }

    size_t i = 0;
    while (nthash.roll()) {
      TEST_ASSERT_EQ(nthash.get_pos(), positions[i])
      i++;
    }
    TEST_ASSERT_EQ(positions.size(), i)
  }

  {
    PRINT_TEST_NAME("base substitution")

    std::string seq = "ACGTACACTGGACTGAGTCT";
    std::string sub = "ACGCGCACTGGACTGAGTCT";

    nthash::NtHash nthash(seq, 3, seq.size());
    nthash::NtHash nthash_subbed(sub, 3, sub.size());

    nthash.roll();
    nthash.sub({ 3, 4 }, { 'C', 'G' });
    nthash_subbed.roll();
    TEST_ASSERT_EQ(nthash.get_hash_num(), nthash_subbed.get_hash_num());
    size_t i;
    for (i = 0; i < nthash.get_hash_num(); ++i) {
      TEST_ASSERT_EQ(nthash.hashes()[i], nthash_subbed.hashes()[i]);
    }
    TEST_ASSERT_EQ(i, 3);
  }

  {
    std::cerr << "Testing RNA" << std::endl;

    std::string seq = "ACGTACACTGGACTGAGTCT";
    nthash::NtHash dna_nthash(seq, 3, 20);

    std::string rna_seq = "ACGUACACUGGACUGAGUCU";
    nthash::NtHash rna_nthash(rna_seq, 3, 20);

    dna_nthash.roll();
    rna_nthash.roll();
    size_t i;
    for (i = 0; i < dna_nthash.get_hash_num(); ++i) {
      TEST_ASSERT_EQ(dna_nthash.hashes()[i], rna_nthash.hashes()[i]);
    }
    TEST_ASSERT_EQ(i, 3);
  }

  {
    PRINT_TEST_NAME("block parsing");

    std::vector<nthash::SpacedSeedBlocks> blocks_out;
    std::vector<nthash::SpacedSeedMonomers> monos_out;

    std::vector<std::string> seeds = { "1101001001011", "11011" };
    std::vector<nthash::SpacedSeedBlocks> blocks_true = {
      { { 0, 2 }, { 11, 13 } }, { { 0, 5 } }
    };
    std::vector<nthash::SpacedSeedMonomers> monos_true = { { 3, 6, 9 }, { 2 } };

    nthash::parse_seeds(seeds, blocks_out, monos_out);

    for (unsigned i_seed = 0; i_seed < seeds.size(); i_seed++) {
      TEST_ASSERT_EQ(blocks_out[i_seed].size(), blocks_true[i_seed].size());
      for (unsigned i_block = 0; i_block < blocks_out[i_seed].size();
           i_block++) {
        TEST_ASSERT_EQ(blocks_out[i_seed][i_block][0],
                       blocks_true[i_seed][i_block][0]);
        TEST_ASSERT_EQ(blocks_out[i_seed][i_block][1],
                       blocks_true[i_seed][i_block][1]);
      }
      TEST_ASSERT_EQ(monos_out[i_seed].size(), monos_true[i_seed].size());
      for (unsigned i = 0; i < monos_true[i_seed].size(); i++) {
        TEST_ASSERT_EQ(monos_out[i_seed][i], monos_true[i_seed][i]);
      }
    }
  }

  {
    PRINT_TEST_NAME("spaced seed hash values")

    std::string seq = "ACATGCATGCA";
    std::vector<std::string> seeds = { "11100111" };
    const unsigned k = seeds[0].length();
    const unsigned h = 3;

    const std::vector<std::array<uint64_t, h>> hashes = {
      { 0x10be4904ad8de5d, 0x3e29e4f4c991628c, 0x3f35c984b13feb20 },
      { 0x8200a7aa3eaf17c8, 0x344198402f4c2a9c, 0xb6423fe62e69c40c },
      { 0x3ce8adcbeaa56532, 0x162e91a4dbedbf11, 0x53173f786a031f45 }
    };

    nthash::SeedNtHash nthash(seq, seeds, h, k);

    for (const auto& h_vals : hashes) {
      nthash.roll();
      TEST_ASSERT_ARRAY_EQ(h_vals, nthash.hashes(), h);
    }
  }

  {
    PRINT_TEST_NAME("spaced seeds")

    std::string seq = "ACGTACACTGGACTGAGTCT";
    std::vector<std::string> seeds = { "111110000000011111",
                                       "111111100001111111" };

    /* Point mutations of k-mer */
    std::string seqM1 = "ACGTACACTTGACTGAGTCT";
    std::string seqM2 = "ACGTACACTGTACTGAGTCT";
    std::string seqM3 = "ACGTACACTGCACTGAGTCT";

    unsigned k = seq.size() - 2;
    TEST_ASSERT_EQ(k, seeds[0].size());
    TEST_ASSERT_EQ(k, seeds[1].size());

    nthash::SeedNtHash seed_nthash(seq, seeds, 2, k);
    nthash::SeedNtHash seed_nthashM1(seqM1, seeds, 2, k);
    nthash::SeedNtHash seed_nthashM2(seqM2, seeds, 2, k);
    nthash::SeedNtHash seed_nthashM3(seqM3, seeds, 2, k);

    std::vector<std::vector<uint64_t>> hashes;

    TEST_ASSERT_EQ(seed_nthash.get_hash_num(), seeds.size() * 2);

    size_t steps = 0;
    for (; seed_nthash.roll(); steps++) {
      TEST_ASSERT_EQ(seed_nthashM1.roll(), true);
      TEST_ASSERT_EQ(seed_nthashM2.roll(), true);
      TEST_ASSERT_EQ(seed_nthashM3.roll(), true);

      const std::string seq_sub = seq.substr(steps, k);
      const std::string seqM1_sub = seqM1.substr(steps, k);
      const std::string seqM2_sub = seqM2.substr(steps, k);
      const std::string seqM3_sub = seqM3.substr(steps, k);
      nthash::SeedNtHash seed_nthash_base(seq_sub, seeds, 2, k);
      nthash::SeedNtHash seed_nthashM1_base(seqM1_sub, seeds, 2, k);
      nthash::SeedNtHash seed_nthashM2_base(seqM2_sub, seeds, 2, k);
      nthash::SeedNtHash seed_nthashM3_base(seqM3_sub, seeds, 2, k);

      TEST_ASSERT_EQ(seed_nthash_base.roll(), true);
      TEST_ASSERT_EQ(seed_nthashM1_base.roll(), true);
      TEST_ASSERT_EQ(seed_nthashM2_base.roll(), true);
      TEST_ASSERT_EQ(seed_nthashM3_base.roll(), true);

      hashes.push_back({});
      for (size_t i = 0; i < seed_nthash.get_hash_num(); i++) {
        const auto hval = seed_nthash.hashes()[i];
        TEST_ASSERT_EQ(hval, seed_nthashM1.hashes()[i]);
        TEST_ASSERT_EQ(hval, seed_nthashM2.hashes()[i]);
        TEST_ASSERT_EQ(hval, seed_nthashM3.hashes()[i]);
        TEST_ASSERT_EQ(hval, seed_nthashM1_base.hashes()[i]);
        TEST_ASSERT_EQ(hval, seed_nthashM2_base.hashes()[i]);
        TEST_ASSERT_EQ(hval, seed_nthashM3_base.hashes()[i]);
        hashes.back().push_back(hval);
      }

      if (seed_nthash.get_pos() > 0) {
        seed_nthash.peek_back();
        for (size_t i = 0; i < seed_nthash.get_hash_num(); i++) {
          TEST_ASSERT_EQ(seed_nthash.hashes()[i], hashes[hashes.size() - 2][i]);
        }
        seed_nthash.peek_back(seq[seed_nthash.get_pos() - 1]);
        for (size_t i = 0; i < seed_nthash.get_hash_num(); i++) {
          TEST_ASSERT_EQ(seed_nthash.hashes()[i], hashes[hashes.size() - 2][i]);
        }
      }
    }
    TEST_ASSERT_EQ(seed_nthashM1.roll(), false);
    TEST_ASSERT_EQ(seed_nthashM2.roll(), false);
    TEST_ASSERT_EQ(seed_nthashM3.roll(), false);
    TEST_ASSERT_EQ(steps, seq.size() - k + 1);
  }

  {
    PRINT_TEST_NAME("spaced seed back roll")

    std::string seq = "ACTAGCTG";
    std::string seed = "110011";
    unsigned h = 3;
    unsigned k = seed.length();

    nthash::SeedNtHash nthash(seq, { seed }, h, k);
    std::stack<uint64_t*> hashes;

    while (nthash.roll()) {
      uint64_t* h_vals = new uint64_t[h];
      std::copy(nthash.hashes(), nthash.hashes() + h, h_vals);
      hashes.push(h_vals);
    }

    TEST_ASSERT_EQ(hashes.size(), seq.length() - k + 1)

    do {
      TEST_ASSERT_ARRAY_EQ(nthash.hashes(), hashes.top(), h)
      hashes.pop();
    } while (nthash.roll_back());
  }

  {
    std::cerr << "Testing rolling spaced seeds vs masked hashing" << std::endl;

    std::string seq = "AACGTGACTACTGACTAGCTAGCTAGCTGATCGT";
    std::vector<std::string> seeds = { "111111111101111111111",
                                       "110111010010010111011" };
    const unsigned k = seeds[0].length();
    const unsigned h = 2;

    std::queue<uint64_t> masked_hashes;
    for (unsigned i = 0; i <= seq.size() - k; i++) {
      for (std::string seed : seeds) {
        uint64_t h_vals[h];
        uint64_t fk = nthash::ntf64(seq.data() + i, k);
        uint64_t rk = nthash::ntr64(seq.data() + i, k);
        uint64_t hs = nthash::mask_hash(fk, rk, seed.data(), seq.data() + i, k);
        nthash::nte64(hs, k, h, h_vals);
        for (unsigned i = 0; i < h; i++) {
          masked_hashes.push(h_vals[i]);
        }
      }
    }

    std::queue<uint64_t> rolled_hashes;
    nthash::SeedNtHash roller(seq, seeds, h, k);
    unsigned num_hashes = roller.get_hash_num_per_seed() * seeds.size();
    while (roller.roll()) {
      for (unsigned i = 0; i < num_hashes; i++) {
        rolled_hashes.push(roller.hashes()[i]);
      }
    }

    TEST_ASSERT_EQ(masked_hashes.size(), rolled_hashes.size());
    while (masked_hashes.size() > 0) {
      uint64_t masked_hash = masked_hashes.front();
      uint64_t rolled_hash = rolled_hashes.front();
      TEST_ASSERT_EQ(masked_hash, rolled_hash);
      masked_hashes.pop();
      rolled_hashes.pop();
    }
  }

  {
    std::cerr << "Testing reset function" << std::endl;

    std::string seq1 = "ACATAAGT";
    std::string seed = "1001";
    unsigned k = seed.length();
    std::queue<uint64_t> hashes1, hashes2;

    nthash::NtHash kmer_hash(seq1, 1, k);
    while (kmer_hash.roll()) {
      hashes1.push(kmer_hash.hashes()[0]);
    }
    kmer_hash.change_seq(seq1, 0);
    while (kmer_hash.roll()) {
      hashes2.push(kmer_hash.hashes()[0]);
    }
    TEST_ASSERT_EQ(hashes1.size(), hashes2.size());
    while (!hashes1.empty()) {
      TEST_ASSERT_EQ(hashes1.front(), hashes2.front());
      hashes1.pop();
      hashes2.pop();
    }

    nthash::SeedNtHash seed_hash(seq1, { seed }, 1, k);
    while (seed_hash.roll()) {
      hashes1.push(seed_hash.hashes()[0]);
    }
    seed_hash.change_seq(seq1);
    while (seed_hash.roll()) {
      hashes2.push(seed_hash.hashes()[0]);
    }
    TEST_ASSERT_EQ(hashes1.size(), hashes2.size());
    while (!hashes1.empty()) {
      TEST_ASSERT_EQ(hashes1.front(), hashes2.front());
      hashes1.pop();
      hashes2.pop();
    }
  }

  {
    PRINT_TEST_NAME("canonical hashing in spaced seeds")

    std::string seq_fwd = "CACTCGGCCACACACACACACACACACCCTCACACACACAAAACGCACAC";
    std::string seq_rev = "GTGTGCGTTTTGTGTGTGTGAGGGTGTGTGTGTGTGTGTGTGGCCGAGTG";

    std::vector<std::string> seeds = {
      "11011000001100101101011000011010110100110000011011",
      "01010000101001110100111011011100101110010100001010",
      "11100000100111010111000100100011101011100100000111",
      "01111000011000111101000011000010111100011000011110",
      "00111000011000111101000011000010111100011000011100",
      "00000000000000000000000011000000000000000000000000",
      "11111111111111111111111100111111111111111111111111",
      "11111111111111111111111111111111111111111111111111",
    };

    const unsigned h = 4;

    nthash::SeedNtHash nthash1(seq_fwd, seeds, h, seeds[0].length());
    nthash::SeedNtHash nthash2(seq_rev, seeds, h, seeds[0].length());

    while (nthash1.roll() && nthash2.roll()) {
      TEST_ASSERT_ARRAY_EQ(nthash1.hashes(), nthash2.hashes(), h);
    }
  }

  return 0;
}