#include "internal.hpp"
#include "nthash/nthash.hpp"

namespace {

using nthash::canonical;
using nthash::CP_OFF;
using nthash::MULTISEED;
using nthash::MULTISHIFT;
using nthash::raise_error;
using nthash::raise_warning;
using nthash::SEED_N;
using nthash::srol;
using nthash::srol_table;
using nthash::sror;
using nthash::typedefs::SpacedSeedBlocks;
using nthash::typedefs::SpacedSeedMonomers;

void
get_blocks(const std::vector<std::string>& seed_strings,
           std::vector<SpacedSeedBlocks>& blocks,
           std::vector<SpacedSeedMonomers>& monomers)
{
  for (const auto& seed_string : seed_strings) {
    const char pad = seed_string[seed_string.length() - 1] == '1' ? '0' : '1';
    const std::string padded_string = seed_string + pad;
    SpacedSeedBlocks care_blocks, ignore_blocks;
    std::vector<unsigned> care_monos, ignore_monos;
    unsigned i_start = 0;
    bool is_care_block = padded_string[0] == '1';
    for (unsigned pos = 0; pos < padded_string.length(); pos++) {
      if (is_care_block && padded_string[pos] == '0') {
        if (pos - i_start == 1) {
          care_monos.push_back(i_start);
        } else {
          const std::array<unsigned, 2> block{ { i_start, pos } };
          care_blocks.push_back(block);
        }
        i_start = pos;
        is_care_block = false;
      } else if (!is_care_block && padded_string[pos] == '1') {
        if (pos - i_start == 1) {
          ignore_monos.push_back(i_start);
        } else {
          const std::array<unsigned, 2> block{ { i_start, pos } };
          ignore_blocks.push_back(block);
        }
        i_start = pos;
        is_care_block = true;
      }
    }
    const unsigned num_cares = care_blocks.size() * 2 + care_monos.size();
    const unsigned num_ignores =
      ignore_blocks.size() * 2 + ignore_monos.size() + 2;
    if (num_ignores < num_cares) {
      const unsigned string_end = seed_string.length();
      const std::array<unsigned, 2> block{ { 0, string_end } };
      ignore_blocks.push_back(block);
      blocks.push_back(ignore_blocks);
      monomers.push_back(ignore_monos);
    } else {
      blocks.push_back(care_blocks);
      monomers.push_back(care_monos);
    }
  }
}

void
parsed_seeds_to_blocks(const std::vector<std::vector<unsigned>>& seeds,
                       unsigned k,
                       std::vector<SpacedSeedBlocks>& blocks,
                       std::vector<SpacedSeedMonomers>& monomers)
{
  std::vector<std::string> seed_strings;
  for (const auto& seed : seeds) {
    std::string seed_string(k, '1');
    for (const auto& i : seed) {
      seed_string[i] = '0';
    }
    seed_strings.push_back(seed_string);
  }
  get_blocks(seed_strings, blocks, monomers);
}

void
check_seeds(const std::vector<std::string>& seeds, unsigned k)
{
  for (const auto& seed : seeds) {

    if (seed.length() != k) {
      raise_error("SeedNtHash",
                  "Spaced seed string length (" +
                    std::to_string(seed.length()) +
                    ") not equal to k=" + std::to_string(k) + " in " + seed);
    }
    const std::string reversed(seed.rbegin(), seed.rend());
    if (seed != reversed) {
      raise_warning(
        "SeedNtHash",
        "Seed " + seed +
          " is not symmetric, reverse-complement hashing will be inconsistent");
    }
  }
}

/**
 * Generate multiple hash values for the input spaced seeds and first k-mer.
 *
 * @param kmer_seq Array of characters representing the k-mer.
 * @param seed_seq Array of SpacedSeed objects representing the seeds' blocks.
 * @param monomers List of the positions that represent blocks of size one for
 * each seed.
 * @param k k-mer size.
 * @param m Number of spaced seeds.
 * @param m2 Number of hashes per seed.
 * @param fh_nomonos Container for the forward hash values before including the
 * size-one blocks.
 * @param rh_nomonos Container for the reverse hash values before including the
 * size-one blocks.
 * @param fh_val Container for the forward hash values after including the
 * size-one blocks.
 * @param rh_val Container for the reverse hash values after including the
 * size-one blocks.
 * @param loc_n Location of the first unknown character in the first sequence.
 * @param h_val Array of size m * m2 for storing the output hash values.
 *
 * @return true if all the care positions of the first k-mer are valid,
 * otherwise false.
 */
inline bool
ntmsm64(const char* kmer_seq,
        const std::vector<SpacedSeedBlocks>& seeds_blocks,
        const std::vector<SpacedSeedMonomers>& seeds_monomers,
        unsigned k,
        unsigned m,
        unsigned m2,
        uint64_t* fh_nomonos,
        uint64_t* rh_nomonos,
        uint64_t* fh_val,
        uint64_t* rh_val,
        unsigned& loc_n,
        uint64_t* h_val)
{
  unsigned i_base;
  uint64_t fh_seed, rh_seed;
  for (unsigned i_seed = 0; i_seed < m; i_seed++) {
    fh_seed = 0;
    rh_seed = 0;
    for (const auto& block : seeds_blocks[i_seed]) {
      for (unsigned pos = block[0]; pos < block[1]; pos++) {
        if (kmer_seq[pos] == SEED_N) {
          loc_n = pos;
          return false;
        }
        fh_seed ^= srol_table((unsigned char)kmer_seq[pos], k - 1 - pos);
        rh_seed ^= srol_table((unsigned char)kmer_seq[pos] & CP_OFF, pos);
      }
    }
    fh_nomonos[i_seed] = fh_seed;
    rh_nomonos[i_seed] = rh_seed;
    for (const auto& pos : seeds_monomers[i_seed]) {
      fh_seed ^= srol_table((unsigned char)kmer_seq[pos], k - 1 - pos);
      rh_seed ^= srol_table((unsigned char)kmer_seq[pos] & CP_OFF, pos);
    }
    fh_val[i_seed] = fh_seed;
    rh_val[i_seed] = rh_seed;
    i_base = i_seed * m2;
    h_val[i_base] = canonical(fh_seed, rh_seed);
    for (unsigned i_hash = 1; i_hash < m2; i_hash++) {
      h_val[i_base + i_hash] = h_val[i_base] * (i_hash ^ k * MULTISEED);
      h_val[i_base + i_hash] ^= h_val[i_base + i_hash] >> MULTISHIFT;
    }
  }
  return true;
}

#define NTMSM64(ROL_HANDLING, IN_HANDLING, OUT_HANDLING, ROR_HANDLING)         \
  unsigned char char_out, char_in;                                             \
  uint64_t fh_seed, rh_seed;                                                   \
  unsigned i_out, i_in, i_base;                                                \
  for (unsigned i_seed = 0; i_seed < m; i_seed++) {                            \
    ROL_HANDLING /* NOLINT(bugprone-macro-parentheses) */                      \
      for (const auto& block : seeds_blocks[i_seed])                           \
    {                                                                          \
      IN_HANDLING                                                              \
      OUT_HANDLING                                                             \
      fh_seed ^= srol_table(char_out, k - i_out);                              \
      fh_seed ^= srol_table(char_in, k - i_in);                                \
      rh_seed ^= srol_table(char_out & CP_OFF, i_out);                         \
      rh_seed ^= srol_table(char_in & CP_OFF, i_in);                           \
    }                                                                          \
    ROR_HANDLING /* NOLINT(bugprone-macro-parentheses) */                      \
      fh_nomonos[i_seed] = fh_seed;                                            \
    rh_nomonos[i_seed] = rh_seed;                                              \
    for (const auto& pos : seeds_monomers[i_seed]) {                           \
      fh_seed ^= srol_table((unsigned char)kmer_seq[pos + 1], k - 1 - pos);    \
      rh_seed ^= srol_table((unsigned char)kmer_seq[pos + 1] & CP_OFF, pos);   \
    }                                                                          \
    fh_val[i_seed] = fh_seed;                                                  \
    rh_val[i_seed] = rh_seed;                                                  \
    i_base = i_seed * m2;                                                      \
    h_val[i_base] = canonical(fh_seed, rh_seed);                               \
    for (unsigned i_hash = 1; i_hash < m2; i_hash++) {                         \
      h_val[i_base + i_hash] = h_val[i_base] * (i_hash ^ k * MULTISEED);       \
      h_val[i_base + i_hash] ^= h_val[i_base + i_hash] >> MULTISHIFT;          \
    }                                                                          \
  }

/**
 * Generate multiple hash values for the input spaced seeds and the next
 * k-mer by performing a forward roll operation.
 *
 * @param kmer_seq Array of characters representing the previous k-mer.
 * @param seed_seq Array of SpacedSeed objects representing the seeds' blocks.
 * @param monomers List of the positions that represent blocks of size one for
 * each seed.
 * @param k k-mer size.
 * @param m Number of spaced seeds.
 * @param m2 Number of hashes per seed.
 * @param fh_nomonos Previous forward hash values before including the size-one
 * blocks.
 * @param rh_nomonos Previous reverse hash values before including the size-one
 * blocks.
 * @param fh_val Previous forward hash values after including the size-one
 * blocks.
 * @param rh_val Previous reverse hash values after including the size-one
 * blocks.
 * @param h_val Array of size m * m2 for storing the output hash values.
 */
inline void
ntmsm64(const char* kmer_seq,
        const std::vector<SpacedSeedBlocks>& seeds_blocks,
        const std::vector<SpacedSeedMonomers>& seeds_monomers,
        unsigned k,
        unsigned m,
        unsigned m2,
        uint64_t* fh_nomonos,
        uint64_t* rh_nomonos,
        uint64_t* fh_val,
        uint64_t* rh_val,
        uint64_t* h_val)
{
  NTMSM64(fh_seed = srol(fh_nomonos[i_seed]); rh_seed = rh_nomonos[i_seed];
          , i_in = block[1];
          char_in = (unsigned char)kmer_seq[i_in];
          , i_out = block[0];
          char_out = (unsigned char)kmer_seq[i_out];
          , rh_seed = sror(rh_seed);)
}

inline void
ntmsm64(const std::deque<char>& kmer_seq,
        const std::vector<SpacedSeedBlocks>& seeds_blocks,
        const std::vector<SpacedSeedMonomers>& seeds_monomers,
        unsigned k,
        unsigned m,
        unsigned m2,
        uint64_t* fh_nomonos,
        uint64_t* rh_nomonos,
        uint64_t* fh_val,
        uint64_t* rh_val,
        uint64_t* h_val)
{
  NTMSM64(fh_seed = srol(fh_nomonos[i_seed]); rh_seed = rh_nomonos[i_seed];
          , i_in = block[1];
          char_in = (unsigned char)kmer_seq[i_in];
          , i_out = block[0];
          char_out = (unsigned char)kmer_seq[i_out];
          , rh_seed = sror(rh_seed);)
}

/**
 * Generate multiple hash values for the input spaced seeds and the next
 * k-mer by performing a backward roll operation.
 *
 * @param kmer_seq Array of characters representing the previous k-mer.
 * @param seed_seq Array of SpacedSeed objects representing the seeds' blocks.
 * @param monomers List of the positions that represent blocks of size one for
 * each seed.
 * @param k k-mer size.
 * @param m Number of spaced seeds.
 * @param m2 Number of hashes per seed.
 * @param fh_nomonos Previous forward hash values before including the size-one
 * blocks.
 * @param rh_nomonos Previous reverse hash values before including the size-one
 * blocks.
 * @param fh_val Previous forward hash values after including the size-one
 * blocks.
 * @param rh_val Previous reverse hash values after including the size-one
 * blocks.
 * @param h_val Array of size m * m2 for storing the output hash values.
 */
inline void
ntmsm64l(const char* kmer_seq,
         const std::vector<SpacedSeedBlocks>& seeds_blocks,
         const std::vector<SpacedSeedMonomers>& seeds_monomers,
         unsigned k,
         unsigned m,
         unsigned m2,
         uint64_t* fh_nomonos,
         uint64_t* rh_nomonos,
         uint64_t* fh_val,
         uint64_t* rh_val,
         uint64_t* h_val)
{
  NTMSM64(fh_seed = fh_nomonos[i_seed]; rh_seed = srol(rh_nomonos[i_seed]);
          , i_in = block[0];
          char_in = (unsigned char)kmer_seq[i_in];
          , i_out = block[1];
          char_out = (unsigned char)kmer_seq[i_out];
          , fh_seed = sror(fh_seed);)
}

inline void
ntmsm64l(const std::deque<char>& kmer_seq,
         const std::vector<SpacedSeedBlocks>& seeds_blocks,
         const std::vector<SpacedSeedMonomers>& seeds_monomers,
         unsigned k,
         unsigned m,
         unsigned m2,
         uint64_t* fh_nomonos,
         uint64_t* rh_nomonos,
         uint64_t* fh_val,
         uint64_t* rh_val,
         uint64_t* h_val)
{
  NTMSM64(fh_seed = fh_nomonos[i_seed]; rh_seed = srol(rh_nomonos[i_seed]);
          , i_in = block[0];
          char_in = (unsigned char)kmer_seq[i_in];
          , i_out = block[1];
          char_out = (unsigned char)kmer_seq[i_out];
          , fh_seed = sror(fh_seed);)
}

/**
 * Generate multiple hash values for the input spaced seeds and the next
 * k-mer by performing a forward peek operation.
 *
 * @param kmer_seq Array of characters representing the previous k-mer.
 * @param seed_seq Array of SpacedSeed objects representing the seeds' blocks.
 * @param monomers List of the positions that represent blocks of size one for
 * each seed.
 * @param k k-mer size.
 * @param m Number of spaced seeds.
 * @param m2 Number of hashes per seed.
 * @param fh_nomonos Previous forward hash values before including the size-one
 * blocks.
 * @param rh_nomonos Previous reverse hash values before including the size-one
 * blocks.
 * @param fh_val Previous forward hash values after including the size-one
 * blocks.
 * @param rh_val Previous reverse hash values after including the size-one
 * blocks.
 * @param h_val Array of size m * m2 for storing the output hash values.
 */
inline void
ntmsm64(const char* kmer_seq,
        char in,
        const std::vector<SpacedSeedBlocks>& seeds_blocks,
        const std::vector<SpacedSeedMonomers>& seeds_monomers,
        unsigned k,
        unsigned m,
        unsigned m2,
        uint64_t* fh_nomonos,
        uint64_t* rh_nomonos,
        uint64_t* fh_val,
        uint64_t* rh_val,
        uint64_t* h_val)
{
  NTMSM64(
    fh_seed = srol(fh_nomonos[i_seed]); rh_seed = rh_nomonos[i_seed];
    , i_in = block[1];
    if (i_in > k - 1) { char_in = in; } else {
      char_in = (unsigned char)kmer_seq[i_in];
    },
    i_out = block[0];
    char_out = (unsigned char)kmer_seq[i_out];
    , rh_seed = sror(rh_seed);)
}

/**
 * Generate multiple hash values for the input spaced seeds and the next
 * k-mer by performing a backwards peek operation.
 *
 * @param kmer_seq Array of characters representing the previous k-mer.
 * @param seed_seq Array of SpacedSeed objects representing the seeds' blocks.
 * @param monomers List of the positions that represent blocks of size one for
 * each seed.
 * @param k k-mer size.
 * @param m Number of spaced seeds.
 * @param m2 Number of hashes per seed.
 * @param fh_nomonos Previous forward hash values before including the size-one
 * blocks.
 * @param rh_nomonos Previous reverse hash values before including the size-one
 * blocks.
 * @param fh_val Previous forward hash values after including the size-one
 * blocks.
 * @param rh_val Previous reverse hash values after including the size-one
 * blocks.
 * @param h_val Array of size m * m2 for storing the output hash values.
 */
inline void
ntmsm64l(const char* kmer_seq,
         char in,
         const std::vector<SpacedSeedBlocks>& seeds_blocks,
         const std::vector<SpacedSeedMonomers>& seeds_monomers,
         unsigned k,
         unsigned m,
         unsigned m2,
         uint64_t* fh_nomonos,
         uint64_t* rh_nomonos,
         uint64_t* fh_val,
         uint64_t* rh_val,
         uint64_t* h_val)
{
  NTMSM64(
    fh_seed = fh_nomonos[i_seed]; rh_seed = srol(rh_nomonos[i_seed]);
    , i_in = block[0];
    if (i_in > k - 1) { char_in = in; } else {
      char_in = (unsigned char)kmer_seq[i_in];
    },
    i_out = block[1];
    char_out = (unsigned char)kmer_seq[i_out];
    , fh_seed = sror(fh_seed);)
}

} // namespace

namespace nthash {

std::vector<std::vector<unsigned>>
parse_seeds(const std::vector<std::string>& seed_strings)
{
  std::vector<std::vector<unsigned>> seed_set;
  for (const auto& seed_string : seed_strings) {
    std::vector<unsigned> seed;
    size_t pos = 0;
    for (const auto& c : seed_string) {
      if (c != '1') {
        seed.push_back(pos);
      }
      ++pos;
    }
    seed_set.push_back(seed);
  }
  return seed_set;
}

SeedNtHash::SeedNtHash(const char* seq,
                       size_t seq_len,
                       const std::vector<std::string>& seeds,
                       typedefs::NUM_HASHES_TYPE num_hashes_per_seed,
                       typedefs::K_TYPE k,
                       size_t pos)
  : seq(seq, seq_len)
  , num_hashes_per_seed(num_hashes_per_seed)
  , k(k)
  , pos(pos)
  , initialized(false)
  , fwd_hash_nomonos(new uint64_t[seeds.size()])
  , rev_hash_nomonos(new uint64_t[seeds.size()])
  , fwd_hash(new uint64_t[seeds.size()])
  , rev_hash(new uint64_t[seeds.size()])
  , hash_arr(new uint64_t[num_hashes_per_seed * seeds.size()])
{
  check_seeds(seeds, k);
  if (seeds[0].size() != k) {
    raise_error("SeedNtHash", "k should be equal to seed string lengths");
  }
  get_blocks(seeds, blocks, monomers);
}

SeedNtHash::SeedNtHash(const char* seq,
                       size_t seq_len,
                       const std::vector<std::vector<unsigned>>& seeds,
                       typedefs::NUM_HASHES_TYPE num_hashes_per_seed,
                       typedefs::K_TYPE k,
                       size_t pos)
  : seq(seq, seq_len)
  , num_hashes_per_seed(num_hashes_per_seed)
  , k(k)
  , pos(pos)
  , initialized(false)
  , fwd_hash_nomonos(new uint64_t[seeds.size()])
  , rev_hash_nomonos(new uint64_t[seeds.size()])
  , fwd_hash(new uint64_t[seeds.size()])
  , rev_hash(new uint64_t[seeds.size()])
  , hash_arr(new uint64_t[num_hashes_per_seed * seeds.size()])
{
  parsed_seeds_to_blocks(seeds, k, blocks, monomers);
}

bool
SeedNtHash::init()
{
  unsigned pos_n = 0;
  while (pos < seq.size() - k + 1 && !ntmsm64(seq.data() + pos,
                                              blocks,
                                              monomers,
                                              k,
                                              blocks.size(),
                                              num_hashes_per_seed,
                                              fwd_hash_nomonos.get(),
                                              rev_hash_nomonos.get(),
                                              fwd_hash.get(),
                                              rev_hash.get(),
                                              pos_n,
                                              hash_arr.get())) {
    pos += pos_n + 1;
  }
  if (pos > seq.size() - k) {
    return false;
  }
  initialized = true;
  return true;
}

bool
SeedNtHash::roll()
{
  if (!initialized) {
    return init();
  }
  if (pos >= seq.size() - k) {
    return false;
  }
  if (SEED_TAB[(unsigned char)seq[pos + k]] == SEED_N) {
    pos += k;
    return init();
  }
  ntmsm64(seq.data() + pos,
          blocks,
          monomers,
          k,
          blocks.size(),
          num_hashes_per_seed,
          fwd_hash_nomonos.get(),
          rev_hash_nomonos.get(),
          fwd_hash.get(),
          rev_hash.get(),
          hash_arr.get());
  ++pos;
  return true;
}

bool
SeedNtHash::roll_back()
{
  if (!initialized) {
    return init();
  }
  if (pos == 0) {
    return false;
  }
  if (SEED_TAB[(unsigned char)seq[pos - 1]] == SEED_N && pos >= k) {
    pos -= k;
    return init();
  }
  if (SEED_TAB[(unsigned char)seq[pos - 1]] == SEED_N) {
    return false;
  }
  ntmsm64l(seq.data() + pos - 1,
           blocks,
           monomers,
           k,
           blocks.size(),
           num_hashes_per_seed,
           fwd_hash_nomonos.get(),
           rev_hash_nomonos.get(),
           fwd_hash.get(),
           rev_hash.get(),
           hash_arr.get());
  --pos;
  return true;
}

bool
SeedNtHash::peek()
{
  if (pos >= seq.size() - k) {
    return false;
  }
  return peek(seq[pos + k]);
}

bool
SeedNtHash::peek(char char_in)
{
  if (!initialized) {
    return init();
  }
  const std::unique_ptr<uint64_t[]> fwd_hash_nomonos_cpy(
    new uint64_t[blocks.size()]);
  const std::unique_ptr<uint64_t[]> rev_hash_nomonos_cpy(
    new uint64_t[blocks.size()]);
  const std::unique_ptr<uint64_t[]> fwd_hash_cpy(new uint64_t[blocks.size()]);
  const std::unique_ptr<uint64_t[]> rev_hash_cpy(new uint64_t[blocks.size()]);
  std::memcpy(fwd_hash_nomonos_cpy.get(),
              fwd_hash_nomonos.get(),
              blocks.size() * sizeof(uint64_t));
  std::memcpy(rev_hash_nomonos_cpy.get(),
              rev_hash_nomonos.get(),
              blocks.size() * sizeof(uint64_t));
  std::memcpy(
    fwd_hash_cpy.get(), fwd_hash.get(), blocks.size() * sizeof(uint64_t));
  std::memcpy(
    rev_hash_cpy.get(), rev_hash.get(), blocks.size() * sizeof(uint64_t));
  ntmsm64(seq.data() + pos,
          char_in,
          blocks,
          monomers,
          k,
          blocks.size(),
          num_hashes_per_seed,
          fwd_hash_nomonos_cpy.get(),
          rev_hash_nomonos_cpy.get(),
          fwd_hash_cpy.get(),
          rev_hash_cpy.get(),
          hash_arr.get());
  return true;
}

bool
SeedNtHash::peek_back()
{
  if (pos == 0) {
    return false;
  }
  return peek_back(seq[pos - 1]);
}

bool
SeedNtHash::peek_back(char char_in)
{
  if (!initialized) {
    return init();
  }
  const std::unique_ptr<uint64_t[]> fwd_hash_nomonos_cpy(
    new uint64_t[blocks.size()]);
  const std::unique_ptr<uint64_t[]> rev_hash_nomonos_cpy(
    new uint64_t[blocks.size()]);
  const std::unique_ptr<uint64_t[]> fwd_hash_cpy(new uint64_t[blocks.size()]);
  const std::unique_ptr<uint64_t[]> rev_hash_cpy(new uint64_t[blocks.size()]);
  std::memcpy(fwd_hash_nomonos_cpy.get(),
              fwd_hash_nomonos.get(),
              blocks.size() * sizeof(uint64_t));
  std::memcpy(rev_hash_nomonos_cpy.get(),
              rev_hash_nomonos.get(),
              blocks.size() * sizeof(uint64_t));
  std::memcpy(
    fwd_hash_cpy.get(), fwd_hash.get(), blocks.size() * sizeof(uint64_t));
  std::memcpy(
    rev_hash_cpy.get(), rev_hash.get(), blocks.size() * sizeof(uint64_t));
  ntmsm64l(seq.data() + pos - 1,
           char_in,
           blocks,
           monomers,
           k,
           blocks.size(),
           num_hashes_per_seed,
           fwd_hash_nomonos_cpy.get(),
           rev_hash_nomonos_cpy.get(),
           fwd_hash_cpy.get(),
           rev_hash_cpy.get(),
           hash_arr.get());
  return true;
}

BlindSeedNtHash::BlindSeedNtHash(const char* seq,
                                 const std::vector<std::string>& seeds,
                                 typedefs::NUM_HASHES_TYPE num_hashes_per_seed,
                                 typedefs::K_TYPE k,
                                 ssize_t pos)
  : seq(seq + pos, seq + pos + k)
  , num_hashes_per_seed(num_hashes_per_seed)
  , k(k)
  , pos(pos)
  , fwd_hash_nomonos(new uint64_t[seeds.size()])
  , rev_hash_nomonos(new uint64_t[seeds.size()])
  , fwd_hash(new uint64_t[seeds.size()])
  , rev_hash(new uint64_t[seeds.size()])
  , hash_arr(new uint64_t[num_hashes_per_seed * seeds.size()])
{
  check_seeds(seeds, k);
  get_blocks(seeds, blocks, monomers);
  unsigned pos_n = 0;
  ntmsm64(seq + pos,
          blocks,
          monomers,
          k,
          blocks.size(),
          num_hashes_per_seed,
          fwd_hash_nomonos.get(),
          rev_hash_nomonos.get(),
          fwd_hash.get(),
          rev_hash.get(),
          pos_n,
          hash_arr.get());
}

void
BlindSeedNtHash::roll(char char_in)
{
  seq.push_back(char_in);
  ntmsm64(seq,
          blocks,
          monomers,
          k,
          blocks.size(),
          num_hashes_per_seed,
          fwd_hash_nomonos.get(),
          rev_hash_nomonos.get(),
          fwd_hash.get(),
          rev_hash.get(),
          hash_arr.get());
  seq.pop_front();
  ++pos;
}

void
BlindSeedNtHash::roll_back(char char_in)
{
  seq.push_front(char_in);
  ntmsm64l(seq,
           blocks,
           monomers,
           k,
           blocks.size(),
           num_hashes_per_seed,
           fwd_hash_nomonos.get(),
           rev_hash_nomonos.get(),
           fwd_hash.get(),
           rev_hash.get(),
           hash_arr.get());
  seq.pop_back();
  --pos;
}

} // namespace nthash