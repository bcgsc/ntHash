#include "nthash/nthash_lowlevel.hpp"

namespace nthash {

uint64_t
ntf64(const char* kmer_seq, const unsigned k)
{
  uint64_t h_val = 0;
  for (unsigned i = 0; i < k / 4; i++) {
    h_val = srol(h_val, 4);
    uint8_t curr_offset = 4 * i;
    uint8_t tetramer_loc =
      64 * CONVERT_TAB[(unsigned char)kmer_seq[curr_offset]] +     // NOLINT
      16 * CONVERT_TAB[(unsigned char)kmer_seq[curr_offset + 1]] + // NOLINT
      4 * CONVERT_TAB[(unsigned char)kmer_seq[curr_offset + 2]] +
      CONVERT_TAB[(unsigned char)kmer_seq[curr_offset + 3]];
    h_val ^= TETRAMER_TAB[tetramer_loc];
  }
  unsigned remainder = k % 4;
  h_val = srol(h_val, remainder);
  if (remainder == 3) {
    uint8_t trimer_loc =
      16 * CONVERT_TAB[(unsigned char)kmer_seq[k - 3]] + // NOLINT
      4 * CONVERT_TAB[(unsigned char)kmer_seq[k - 2]] +
      CONVERT_TAB[(unsigned char)kmer_seq[k - 1]];
    h_val ^= TRIMER_TAB[trimer_loc];
  } else if (remainder == 2) {
    uint8_t dimer_loc = 4 * CONVERT_TAB[(unsigned char)kmer_seq[k - 2]] +
                        CONVERT_TAB[(unsigned char)kmer_seq[k - 1]];
    h_val ^= DIMER_TAB[dimer_loc];
  } else if (remainder == 1) {
    h_val ^= SEED_TAB[(unsigned char)kmer_seq[k - 1]];
  }
  return h_val;
}

uint64_t
ntr64(const char* kmer_seq, const unsigned k)
{
  uint64_t h_val = 0;
  unsigned remainder = k % 4;
  if (remainder == 3) {
    uint8_t trimer_loc =
      16 * RC_CONVERT_TAB[(unsigned char)kmer_seq[k - 1]] + // NOLINT
      4 * RC_CONVERT_TAB[(unsigned char)kmer_seq[k - 2]] +
      RC_CONVERT_TAB[(unsigned char)kmer_seq[k - 3]];
    h_val ^= TRIMER_TAB[trimer_loc];
  } else if (remainder == 2) {
    uint8_t dimer_loc = 4 * RC_CONVERT_TAB[(unsigned char)kmer_seq[k - 1]] +
                        RC_CONVERT_TAB[(unsigned char)kmer_seq[k - 2]];
    h_val ^= DIMER_TAB[dimer_loc];
  } else if (remainder == 1) {
    h_val ^= SEED_TAB[(unsigned char)kmer_seq[k - 1] & CP_OFF];
  }
  for (unsigned i = 0; i < k / 4; i++) {
    h_val = srol(h_val, 4);
    uint8_t curr_offset = 4 * (k / 4 - i) - 1;
    uint8_t tetramer_loc =
      64 * RC_CONVERT_TAB[(unsigned char)kmer_seq[curr_offset]] +     // NOLINT
      16 * RC_CONVERT_TAB[(unsigned char)kmer_seq[curr_offset - 1]] + // NOLINT
      4 * RC_CONVERT_TAB[(unsigned char)kmer_seq[curr_offset - 2]] +
      RC_CONVERT_TAB[(unsigned char)kmer_seq[curr_offset - 3]];
    h_val ^= TETRAMER_TAB[tetramer_loc];
  }
  return h_val;
}

uint64_t
ntf64(const uint64_t fh_val,
      const unsigned k,
      const unsigned char char_out,
      const unsigned char char_in)
{
  uint64_t h_val = srol(fh_val);
  h_val ^= SEED_TAB[char_in];
  h_val ^= MS_TAB(char_out, k);
  return h_val;
}

uint64_t
ntr64(const uint64_t rh_val,
      const unsigned k,
      const unsigned char char_out,
      const unsigned char char_in)
{
  uint64_t h_val = rh_val ^ MS_TAB(char_in & CP_OFF, k);
  h_val ^= SEED_TAB[char_out & CP_OFF];
  h_val = sror(h_val);
  return h_val;
}

uint64_t
ntc64(const char* kmer_seq, const unsigned k)
{
  uint64_t fh_val = 0, rh_val = 0;
  fh_val = ntf64(kmer_seq, k);
  rh_val = ntr64(kmer_seq, k);
  return canonical(fh_val, rh_val);
}

uint64_t
ntc64(const char* kmer_seq,
      const unsigned k,
      uint64_t& fh_val,
      uint64_t& rh_val)
{
  fh_val = ntf64(kmer_seq, k);
  rh_val = ntr64(kmer_seq, k);
  return canonical(fh_val, rh_val);
}

uint64_t
ntc64(const unsigned char char_out,
      const unsigned char char_in,
      const unsigned k,
      uint64_t& fh_val,
      uint64_t& rh_val)
{
  fh_val = ntf64(fh_val, k, char_out, char_in);
  rh_val = ntr64(rh_val, k, char_out, char_in);
  return canonical(fh_val, rh_val);
}

uint64_t
ntf64l(const uint64_t rh_val,
       const unsigned k,
       const unsigned char char_out,
       const unsigned char char_in)
{
  uint64_t h_val = rh_val ^ MS_TAB(char_in, k);
  h_val ^= SEED_TAB[char_out];
  h_val = sror(h_val);
  return h_val;
}

uint64_t
ntr64l(const uint64_t fh_val,
       const unsigned k,
       const unsigned char char_out,
       const unsigned char char_in)
{
  uint64_t h_val = srol(fh_val);
  h_val ^= SEED_TAB[char_in & CP_OFF];
  h_val ^= MS_TAB(char_out & CP_OFF, k);
  return h_val;
}

uint64_t
ntc64l(const unsigned char char_out,
       const unsigned char char_in,
       const unsigned k,
       uint64_t& fh_val,
       uint64_t& rh_val)
{
  fh_val = ntf64l(fh_val, k, char_out, char_in);
  rh_val = ntr64l(rh_val, k, char_out, char_in);
  return canonical(fh_val, rh_val);
}

void
nte64(const uint64_t bh_val,
      const unsigned k,
      const unsigned h,
      uint64_t* h_val)
{
  uint64_t t_val;
  h_val[0] = bh_val;
  for (unsigned i = 1; i < h; i++) {
    t_val = bh_val * (i ^ k * MULTISEED);
    t_val ^= t_val >> MULTISHIFT;
    h_val[i] = t_val;
  }
}

void
ntmc64(const char* kmer_seq,
       const unsigned k,
       const unsigned m,
       uint64_t* h_val)
{
  uint64_t b_val = ntc64(kmer_seq, k);
  nte64(b_val, k, m, h_val);
}

void
ntmc64(const char* kmer_seq,
       const unsigned k,
       const unsigned m,
       uint64_t& fh_val,
       uint64_t& rh_val,
       uint64_t* h_val)
{
  uint64_t b_val = ntc64(kmer_seq, k, fh_val, rh_val);
  nte64(b_val, k, m, h_val);
}

void
ntmc64(const unsigned char char_out,
       const unsigned char char_in,
       const unsigned k,
       const unsigned m,
       uint64_t& fh_val,
       uint64_t& rh_val,
       uint64_t* h_val)
{
  uint64_t b_val = ntc64(char_out, char_in, k, fh_val, rh_val);
  nte64(b_val, k, m, h_val);
}

void
ntmc64l(const unsigned char char_out,
        const unsigned char char_in,
        const unsigned k,
        const unsigned m,
        uint64_t& fh_val,
        uint64_t& rh_val,
        uint64_t* h_val)
{
  uint64_t b_val = ntc64l(char_out, char_in, k, fh_val, rh_val);
  nte64(b_val, k, m, h_val);
}

bool
ntc64(const char* kmer_seq, const unsigned k, uint64_t& h_val, unsigned& loc_n)
{
  h_val = 0;
  loc_n = 0;
  uint64_t fh_val = 0, rh_val = 0;
  for (int i = int(k - 1); i >= 0; i--) {
    if (SEED_TAB[(unsigned char)kmer_seq[i]] == SEED_N) {
      loc_n = i;
      return false;
    }
    fh_val = srol(fh_val);
    fh_val ^= SEED_TAB[(unsigned char)kmer_seq[k - 1 - i]];

    rh_val = srol(rh_val);
    rh_val ^= SEED_TAB[(unsigned char)kmer_seq[i] & CP_OFF];
  }
  h_val = canonical(fh_val, rh_val);
  return true;
}

bool
ntmc64(const char* kmer_seq,
       const unsigned k,
       const unsigned m,
       unsigned& loc_n,
       uint64_t* h_val)
{
  uint64_t b_val = 0, fh_val = 0, rh_val = 0;
  loc_n = 0;
  for (int i = int(k - 1); i >= 0; i--) {
    if (SEED_TAB[(unsigned char)kmer_seq[i]] == SEED_N) {
      loc_n = i;
      return false;
    }
    fh_val = srol(fh_val);
    fh_val ^= SEED_TAB[(unsigned char)kmer_seq[k - 1 - i]];

    rh_val = srol(rh_val);
    rh_val ^= SEED_TAB[(unsigned char)kmer_seq[i] & CP_OFF];
  }
  b_val = canonical(fh_val, rh_val);
  nte64(b_val, k, m, h_val);
  return true;
}

bool
ntc64(const char* kmer_seq,
      const unsigned k,
      uint64_t& fh_val,
      uint64_t& rh_val,
      uint64_t& h_val,
      unsigned& loc_n)
{
  h_val = fh_val = rh_val = 0;
  loc_n = 0;
  for (int i = int(k - 1); i >= 0; i--) {
    if (SEED_TAB[(unsigned char)kmer_seq[i]] == SEED_N) {
      loc_n = i;
      return false;
    }
    fh_val = srol(fh_val);
    fh_val ^= SEED_TAB[(unsigned char)kmer_seq[k - 1 - i]];

    rh_val = srol(rh_val);
    rh_val ^= SEED_TAB[(unsigned char)kmer_seq[i] & CP_OFF];
  }
  h_val = canonical(fh_val, rh_val);
  return true;
}

bool
ntmc64(const char* kmer_seq,
       const unsigned k,
       const unsigned m,
       uint64_t& fh_val,
       uint64_t& rh_val,
       unsigned& loc_n,
       uint64_t* h_val)
{
  fh_val = rh_val = 0;
  uint64_t b_val = 0;
  loc_n = 0;
  for (int i = int(k - 1); i >= 0; i--) {
    if (SEED_TAB[(unsigned char)kmer_seq[i]] == SEED_N) {
      loc_n = i;
      return false;
    }
    fh_val = srol(fh_val);
    fh_val ^= SEED_TAB[(unsigned char)kmer_seq[k - 1 - i]];

    rh_val = srol(rh_val);
    rh_val ^= SEED_TAB[(unsigned char)kmer_seq[i] & CP_OFF];
  }
  b_val = canonical(fh_val, rh_val);
  nte64(b_val, k, m, h_val);
  return true;
}

bool
ntmc64(const char* kmer_seq,
       const unsigned k,
       const unsigned m,
       uint64_t& fh_val,
       uint64_t& rh_val,
       unsigned& loc_n,
       uint64_t* h_val,
       bool& h_stn)
{
  fh_val = rh_val = 0;
  uint64_t b_val = 0;
  loc_n = 0;
  for (int i = int(k - 1); i >= 0; i--) {
    if (SEED_TAB[(unsigned char)kmer_seq[i]] == SEED_N) {
      loc_n = i;
      return false;
    }
    fh_val = srol(fh_val);
    fh_val ^= SEED_TAB[(unsigned char)kmer_seq[k - 1 - i]];

    rh_val = srol(rh_val);
    rh_val ^= SEED_TAB[(unsigned char)kmer_seq[i] & CP_OFF];
  }
  h_stn = rh_val < fh_val;
  b_val = canonical(fh_val, rh_val);
  nte64(b_val, k, m, h_val);
  return true;
}

void
ntmc64(const unsigned char char_out,
       const unsigned char char_in,
       const unsigned k,
       const unsigned m,
       uint64_t& fh_val,
       uint64_t& rh_val,
       uint64_t* h_val,
       bool& h_stn)
{
  uint64_t b_val = ntc64(char_out, char_in, k, fh_val, rh_val);
  h_stn = rh_val < fh_val;
  nte64(b_val, k, m, h_val);
}

uint64_t
mask_hash(uint64_t& fk_val,
          uint64_t& rk_val,
          const char* seed_seq,
          const char* kmer_seq,
          const unsigned k)
{
  uint64_t fs_val = fk_val, rs_val = rk_val;
  for (unsigned i = 0; i < k; i++) {
    if (seed_seq[i] != '1') {
      fs_val ^= MS_TAB((unsigned char)kmer_seq[i], k - 1 - i);
      rs_val ^= MS_TAB((unsigned char)kmer_seq[i] & CP_OFF, i);
    }
  }
  return canonical(fs_val, rs_val);
}

void
sub_hash(uint64_t fh_val,
         uint64_t rh_val,
         const char* kmer_seq,
         const std::vector<unsigned>& positions,
         const std::vector<unsigned char>& new_bases,
         const unsigned k,
         const unsigned m,
         uint64_t* h_val)
{
  uint64_t b_val = 0;

  for (size_t i = 0; i < positions.size(); i++) {
    const auto pos = positions[i];
    const auto new_base = new_bases[i];

    fh_val ^= MS_TAB((unsigned char)kmer_seq[pos], k - 1 - pos);
    fh_val ^= MS_TAB(new_base, k - 1 - pos);

    rh_val ^= MS_TAB((unsigned char)kmer_seq[pos] & CP_OFF, pos);
    rh_val ^= MS_TAB(new_base & CP_OFF, pos);
  }

  b_val = canonical(fh_val, rh_val);
  nte64(b_val, k, m, h_val);
}

bool
ntmsm64(const char* kmer_seq,
        const std::vector<SpacedSeedBlocks>& seeds_blocks,
        const std::vector<SpacedSeedMonomers>& seeds_monomers,
        const unsigned k,
        const unsigned m,
        const unsigned m2,
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
        fh_seed ^= MS_TAB((unsigned char)kmer_seq[pos], k - 1 - pos);
        rh_seed ^= MS_TAB((unsigned char)kmer_seq[pos] & CP_OFF, pos);
      }
    }
    fh_nomonos[i_seed] = fh_seed;
    rh_nomonos[i_seed] = rh_seed;
    for (unsigned pos : seeds_monomers[i_seed]) {
      fh_seed ^= MS_TAB((unsigned char)kmer_seq[pos], k - 1 - pos);
      rh_seed ^= MS_TAB((unsigned char)kmer_seq[pos] & CP_OFF, pos);
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

void
ntmsm64(const char* kmer_seq,
        const std::vector<SpacedSeedBlocks>& seeds_blocks,
        const std::vector<SpacedSeedMonomers>& seeds_monomers,
        const unsigned k,
        const unsigned m,
        const unsigned m2,
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

void
ntmsm64l(const char* kmer_seq,
         const std::vector<SpacedSeedBlocks>& seeds_blocks,
         const std::vector<SpacedSeedMonomers>& seeds_monomers,
         const unsigned k,
         const unsigned m,
         const unsigned m2,
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

void
ntmsm64(const char* kmer_seq,
        const char in,
        const std::vector<SpacedSeedBlocks>& seeds_blocks,
        const std::vector<SpacedSeedMonomers>& seeds_monomers,
        const unsigned k,
        const unsigned m,
        const unsigned m2,
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

void
ntmsm64l(const char* kmer_seq,
         const char in,
         const std::vector<SpacedSeedBlocks>& seeds_blocks,
         const std::vector<SpacedSeedMonomers>& seeds_monomers,
         const unsigned k,
         const unsigned m,
         const unsigned m2,
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

} // namespace nthash