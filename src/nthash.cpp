#include "nthash/nthash.hpp"

namespace {

inline void
raise_warning(const std::string& class_name, const std::string& msg)
{
  std::cerr << "[ntHash::" << class_name << "] \33[33mWARNING: \33[0m" << msg
            << std::endl;
}

inline void
raise_error(const std::string& class_name, const std::string& msg)
{
  std::cerr << "[ntHash::" << class_name << "] \33[31mERROR: \33[0m" << msg
            << std::endl;
  std::exit(1); // NOLINT(concurrency-mt-unsafe)
}

}
namespace nthash {

NtHash::NtHash(const char* seq,
               size_t seq_len,
               unsigned hash_num,
               unsigned k,
               size_t pos)
  : seq(seq)
  , seq_len(seq_len)
  , hash_num(hash_num)
  , k(k)
  , pos(pos)
  , initialized(false)
  , hashes_array(new uint64_t[hash_num])
{
  if (k == 0) {
    raise_error("NtHash", "k must be greater than 0");
  }
  if (k > NTHASH_K_MAX) {
    raise_error("NtHash",
                "passed k value (" + std::to_string(k) +
                  ") is larger than allowed (" + std::to_string(NTHASH_K_MAX) +
                  ")");
  }
  if (hash_num > NTHASH_HASH_NUM_MAX) {
    raise_error("NtHash",
                "passed number of hashes (" + std::to_string(hash_num) +
                  ") is larger than allowed (" +
                  std::to_string(NTHASH_HASH_NUM_MAX) + ")");
  }
  if (seq_len < k) {
    raise_error("NtHash",
                "sequence length (" + std::to_string(seq_len) +
                  ") is smaller than k (" + std::to_string(k) + ")");
  }
  if (pos >= seq_len) {
    raise_error("NtHash",
                "passed position (" + std::to_string(pos) +
                  ") is larger than sequence length (" +
                  std::to_string(seq_len) + ")");
  }
}

NtHash::NtHash(const std::string& seq,
               unsigned hash_num,
               unsigned k,
               size_t pos)
  : NtHash(seq.c_str(), seq.size(), hash_num, k, pos)
{}

NtHash::NtHash(const NtHash& nthash)
  : seq(nthash.seq)
  , seq_len(nthash.seq_len)
  , hash_num(nthash.hash_num)
  , k(nthash.k)
  , pos(nthash.pos)
  , initialized(nthash.initialized)
  , hashes_array(new uint64_t[hash_num])
  , forward_hash(nthash.forward_hash)
  , reverse_hash(nthash.reverse_hash)
{
  std::memcpy(
    hashes_array.get(), nthash.hashes_array.get(), hash_num * sizeof(uint64_t));
}

BlindNtHash::BlindNtHash(const char* seq,
                         size_t seq_len,
                         unsigned hash_num,
                         unsigned k,
                         size_t pos)
  : seq(new char[seq_len])
  , seq_len(seq_len)
  , hash_num(hash_num)
  , k(k)
  , pos(pos)
  , initialized(false)
  , hashes_array(new uint64_t[hash_num])
{
  if (k > NTHASH_K_MAX) {
    raise_error("BlindNtHash",
                "passed k value (" + std::to_string(k) +
                  ") is larger than allowed (" + std::to_string(NTHASH_K_MAX) +
                  ")");
  }
  if (hash_num > NTHASH_HASH_NUM_MAX) {
    raise_error("BlindNtHash",
                "passed number of hashes (" + std::to_string(hash_num) +
                  ") is larger than allowed (" +
                  std::to_string(NTHASH_HASH_NUM_MAX) + ")");
  }
  std::memcpy(this->seq.get(), seq, seq_len);
}

BlindNtHash::BlindNtHash(const std::string& seq,
                         unsigned hash_num,
                         unsigned k,
                         size_t pos)
  : BlindNtHash(seq.c_str(), seq.size(), hash_num, k, pos)
{}

BlindNtHash::BlindNtHash(const BlindNtHash& nthash)
  : seq(new char[nthash.seq_len])
  , seq_len(nthash.seq_len)
  , hash_num(nthash.hash_num)
  , k(nthash.k)
  , pos(nthash.pos)
  , initialized(nthash.initialized)
  , hashes_array(new uint64_t[hash_num])
  , forward_hash(nthash.forward_hash)
  , reverse_hash(nthash.reverse_hash)
{
  std::memcpy(this->seq.get(), nthash.seq.get(), nthash.seq_len);
  std::memcpy(
    hashes_array.get(), nthash.hashes_array.get(), hash_num * sizeof(uint64_t));
}

SeedNtHash::SeedNtHash(const char* seq,
                       size_t seq_len,
                       const std::vector<SpacedSeed>& seeds,
                       unsigned hash_num_per_seed,
                       unsigned k,
                       size_t pos)
  : nthash(seq, seq_len, seeds.size() * hash_num_per_seed, k, pos)
  , hash_num_per_seed(hash_num_per_seed)
  , fh_no_monomers(new uint64_t[seeds.size()])
  , rh_no_monomers(new uint64_t[seeds.size()])
  , forward_hash(new uint64_t[seeds.size()])
  , reverse_hash(new uint64_t[seeds.size()])
{
  parsed_seeds_to_blocks(seeds, k, blocks, monomers);
}

SeedNtHash::SeedNtHash(const std::string& seq,
                       const std::vector<SpacedSeed>& seeds,
                       unsigned hash_num_per_seed,
                       unsigned k,
                       size_t pos)
  : nthash(seq, seeds.size() * hash_num_per_seed, k, pos)
  , hash_num_per_seed(hash_num_per_seed)
  , fh_no_monomers(new uint64_t[seeds.size()])
  , rh_no_monomers(new uint64_t[seeds.size()])
  , forward_hash(new uint64_t[seeds.size()])
  , reverse_hash(new uint64_t[seeds.size()])
{
  parsed_seeds_to_blocks(seeds, k, blocks, monomers);
}

SeedNtHash::SeedNtHash(const char* seq,
                       size_t seq_len,
                       const std::vector<std::string>& seeds,
                       unsigned hash_num_per_seed,
                       unsigned k,
                       size_t pos)
  : nthash(seq, seq_len, seeds.size() * hash_num_per_seed, k, pos)
  , hash_num_per_seed(hash_num_per_seed)
  , fh_no_monomers(new uint64_t[seeds.size()])
  , rh_no_monomers(new uint64_t[seeds.size()])
  , forward_hash(new uint64_t[seeds.size()])
  , reverse_hash(new uint64_t[seeds.size()])
{
  check_seeds(seeds, k);
  parse_seeds(seeds, blocks, monomers);
}

SeedNtHash::SeedNtHash(const std::string& seq,
                       const std::vector<std::string>& seeds,
                       unsigned hash_num_per_seed,
                       unsigned k,
                       size_t pos)
  : nthash(seq, seeds.size() * hash_num_per_seed, k, pos)
  , hash_num_per_seed(hash_num_per_seed)
  , fh_no_monomers(new uint64_t[seeds.size()])
  , rh_no_monomers(new uint64_t[seeds.size()])
  , forward_hash(new uint64_t[seeds.size()])
  , reverse_hash(new uint64_t[seeds.size()])
{
  check_seeds(seeds, k);
  parse_seeds(seeds, blocks, monomers);
}

SeedNtHash::SeedNtHash(const SeedNtHash& seed_nthash)
  : nthash(seed_nthash.nthash)
  , hash_num_per_seed(seed_nthash.hash_num_per_seed)
  , blocks(seed_nthash.blocks)
  , monomers(seed_nthash.monomers)
  , fh_no_monomers(new uint64_t[seed_nthash.blocks.size()])
  , rh_no_monomers(new uint64_t[seed_nthash.blocks.size()])
  , forward_hash(new uint64_t[seed_nthash.blocks.size()])
  , reverse_hash(new uint64_t[seed_nthash.blocks.size()])
{
  std::memcpy(fh_no_monomers.get(),
              seed_nthash.fh_no_monomers.get(),
              seed_nthash.blocks.size() * sizeof(uint64_t));
  std::memcpy(rh_no_monomers.get(),
              seed_nthash.rh_no_monomers.get(),
              seed_nthash.blocks.size() * sizeof(uint64_t));
  std::memcpy(forward_hash.get(),
              seed_nthash.forward_hash.get(),
              seed_nthash.blocks.size() * sizeof(uint64_t));
  std::memcpy(reverse_hash.get(),
              seed_nthash.reverse_hash.get(),
              seed_nthash.blocks.size() * sizeof(uint64_t));
}

BlindSeedNtHash::BlindSeedNtHash(const char* seq,
                                 const std::vector<std::string>& seeds,
                                 unsigned hash_num_per_seed,
                                 unsigned k,
                                 size_t pos)
  : seq(seq, seq + k)
  , hash_num_per_seed(hash_num_per_seed)
  , k(k)
  , pos(pos)
  , initialized(true)
  , fh_no_monomers(new uint64_t[seeds.size()])
  , rh_no_monomers(new uint64_t[seeds.size()])
  , forward_hash(new uint64_t[seeds.size()])
  , reverse_hash(new uint64_t[seeds.size()])
  , hashes_array(new uint64_t[hash_num_per_seed * seeds.size()])
{
  check_seeds(seeds, k);
  parse_seeds(seeds, blocks, monomers);
  init();
}

BlindSeedNtHash::BlindSeedNtHash(const std::string& seq,
                                 const std::vector<std::string>& seeds,
                                 unsigned hash_num_per_seed,
                                 unsigned k,
                                 size_t pos)
  : BlindSeedNtHash(seq.data(), seeds, hash_num_per_seed, k, pos)
{}

BlindSeedNtHash::BlindSeedNtHash(BlindSeedNtHash& seed_nthash)
  : seq(seed_nthash.seq)
  , hash_num_per_seed(seed_nthash.hash_num_per_seed)
  , k(seed_nthash.k)
  , pos(seed_nthash.pos)
  , initialized(seed_nthash.initialized)
  , blocks(seed_nthash.blocks)
  , monomers(seed_nthash.monomers)
  , fh_no_monomers(new uint64_t[seed_nthash.blocks.size()])
  , rh_no_monomers(new uint64_t[seed_nthash.blocks.size()])
  , forward_hash(new uint64_t[seed_nthash.blocks.size()])
  , reverse_hash(new uint64_t[seed_nthash.blocks.size()])
  , hashes_array(new uint64_t[hash_num_per_seed * seed_nthash.blocks.size()])
{
  std::memcpy(fh_no_monomers.get(),
              seed_nthash.fh_no_monomers.get(),
              seed_nthash.blocks.size() * sizeof(uint64_t));
  std::memcpy(rh_no_monomers.get(),
              seed_nthash.rh_no_monomers.get(),
              seed_nthash.blocks.size() * sizeof(uint64_t));
  std::memcpy(forward_hash.get(),
              seed_nthash.forward_hash.get(),
              seed_nthash.blocks.size() * sizeof(uint64_t));
  std::memcpy(reverse_hash.get(),
              seed_nthash.reverse_hash.get(),
              seed_nthash.blocks.size() * sizeof(uint64_t));
  std::memcpy(hashes_array.get(),
              seed_nthash.hashes_array.get(),
              hash_num_per_seed * seed_nthash.blocks.size() * sizeof(uint64_t));
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
    std::string reversed(seed.rbegin(), seed.rend());
    if (seed != reversed) {
      raise_warning(
        "SeedNtHash",
        "Seed " + seed +
          " is not symmetric, reverse-complement hashing will be inconsistent");
    }
  }
}

std::vector<SpacedSeed>
parse_seeds(const std::vector<std::string>& seed_strings)
{
  std::vector<SpacedSeed> seed_set;
  for (const auto& seed_string : seed_strings) {
    SpacedSeed seed;
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

void
parse_seeds(const std::vector<std::string>& seed_strings,
            std::vector<SpacedSeedBlocks>& out_blocks,
            std::vector<SpacedSeedMonomers>& out_monomers)
{
  for (const auto& seed_string : seed_strings) {
    char pad = seed_string[seed_string.length() - 1] == '1' ? '0' : '1';
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
          std::array<unsigned, 2> block{ { i_start, pos } };
          care_blocks.push_back(block);
        }
        i_start = pos;
        is_care_block = false;
      } else if (!is_care_block && padded_string[pos] == '1') {
        if (pos - i_start == 1) {
          ignore_monos.push_back(i_start);
        } else {
          std::array<unsigned, 2> block{ { i_start, pos } };
          ignore_blocks.push_back(block);
        }
        i_start = pos;
        is_care_block = true;
      }
    }
    unsigned num_cares = care_blocks.size() * 2 + care_monos.size();
    unsigned num_ignores = ignore_blocks.size() * 2 + ignore_monos.size() + 2;
    if (num_ignores < num_cares) {
      unsigned string_end = seed_string.length();
      std::array<unsigned, 2> block{ { 0, string_end } };
      ignore_blocks.push_back(block);
      out_blocks.push_back(ignore_blocks);
      out_monomers.push_back(ignore_monos);
    } else {
      out_blocks.push_back(care_blocks);
      out_monomers.push_back(care_monos);
    }
  }
}

void
parsed_seeds_to_blocks(const std::vector<SpacedSeed>& seeds,
                       unsigned k,
                       std::vector<SpacedSeedBlocks>& out_blocks,
                       std::vector<SpacedSeedMonomers>& out_monomers)
{
  std::vector<std::string> seed_strings;
  for (const SpacedSeed& seed : seeds) {
    std::string seed_string(k, '1');
    for (const auto& i : seed) {
      seed_string[i] = '0';
    }
    seed_strings.push_back(seed_string);
  }
  parse_seeds(seed_strings, out_blocks, out_monomers);
}

void
NtHash::sub(const std::vector<unsigned>& positions,
            const std::vector<unsigned char>& new_bases)
{
  sub_hash(forward_hash,
           reverse_hash,
           seq + pos,
           positions,
           new_bases,
           get_k(),
           get_hash_num(),
           hashes_array.get());
}

void
BlindNtHash::sub(const std::vector<unsigned>& positions,
                 const std::vector<unsigned char>& new_bases)
{
  sub_hash(forward_hash,
           reverse_hash,
           seq.get() + pos,
           positions,
           new_bases,
           get_k(),
           get_hash_num(),
           hashes_array.get());
}

} // namespace nthash