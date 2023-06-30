#include "internal.hpp"
#include "nthash/nthash.hpp"

namespace {

using namespace nthash;

/**
 * Check the current k-mer for non ACGTU's
 * @param seq C array containing the sequence's characters
 * @param k k-mer size
 * @return `true` if any of the first k characters is not an ACGTU, `false`
 * otherwise
 */
inline bool
is_invalid_kmer(const char* seq, unsigned k, size_t& pos_n)
{
  for (int i = k - 1; i >= 0; i--) {
    if (SEED_TAB[(unsigned char)seq[i]] == SEED_N) {
      pos_n = i;
      return true;
    }
  }
  return false;
}

/**
 * Generate the forward-strand hash value of the first k-mer in the sequence.
 * @param seq C array containing the sequence's characters
 * @param k k-mer size
 * @return Hash value of k-mer_0
 */
inline uint64_t
base_forward_hash(const char* seq, unsigned k)
{
  uint64_t h_val = 0;
  for (unsigned i = 0; i < k / 4; i++) {
    h_val = srol(h_val, 4);
    uint8_t curr_offset = 4 * i;
    uint8_t tetramer_loc =
      64 * CONVERT_TAB[(unsigned char)seq[curr_offset]] +     // NOLINT
      16 * CONVERT_TAB[(unsigned char)seq[curr_offset + 1]] + // NOLINT
      4 * CONVERT_TAB[(unsigned char)seq[curr_offset + 2]] +
      CONVERT_TAB[(unsigned char)seq[curr_offset + 3]];
    h_val ^= TETRAMER_TAB[tetramer_loc];
  }
  unsigned remainder = k % 4;
  h_val = srol(h_val, remainder);
  if (remainder == 3) {
    uint8_t trimer_loc = 16 * CONVERT_TAB[(unsigned char)seq[k - 3]] + // NOLINT
                         4 * CONVERT_TAB[(unsigned char)seq[k - 2]] +
                         CONVERT_TAB[(unsigned char)seq[k - 1]];
    h_val ^= TRIMER_TAB[trimer_loc];
  } else if (remainder == 2) {
    uint8_t dimer_loc = 4 * CONVERT_TAB[(unsigned char)seq[k - 2]] +
                        CONVERT_TAB[(unsigned char)seq[k - 1]];
    h_val ^= DIMER_TAB[dimer_loc];
  } else if (remainder == 1) {
    h_val ^= SEED_TAB[(unsigned char)seq[k - 1]];
  }
  return h_val;
}

/**
 * Perform a roll operation on the forward strand by removing char_out and
 * including char_in.
 * @param fh_val Previous hash value computed for the sequence
 * @param k k-mer size
 * @param char_out Character to be removed
 * @param char_in Character to be included
 * @return Rolled forward hash value
 */
inline uint64_t
next_forward_hash(uint64_t fh_val,
                  unsigned k,
                  unsigned char char_out,
                  unsigned char char_in)
{
  uint64_t h_val = srol(fh_val);
  h_val ^= SEED_TAB[char_in];
  h_val ^= MS_TAB(char_out, k);
  return h_val;
}

/**
 * Perform a roll back operation on the forward strand.
 * @param fh_val Previous hash value computed for the sequence
 * @param k k-mer size
 * @param char_out Character to be removed
 * @param char_in Character to be included
 * @return Forward hash value rolled back
 */
inline uint64_t
prev_forward_hash(uint64_t fh_val,
                  unsigned k,
                  unsigned char char_out,
                  unsigned char char_in)
{
  uint64_t h_val = fh_val ^ MS_TAB(char_in, k);
  h_val ^= SEED_TAB[char_out];
  h_val = sror(h_val);
  return h_val;
}

/**
 * Generate a hash value for the reverse-complement of the first k-mer in the
 * sequence.
 * @param seq C array containing the sequence's characters
 * @param k k-mer size
 * @return Hash value of the reverse-complement of k-mer_0
 */
inline uint64_t
base_reverse_hash(const char* seq, unsigned k)
{
  uint64_t h_val = 0;
  unsigned remainder = k % 4;
  if (remainder == 3) {
    uint8_t trimer_loc =
      16 * RC_CONVERT_TAB[(unsigned char)seq[k - 1]] + // NOLINT
      4 * RC_CONVERT_TAB[(unsigned char)seq[k - 2]] +
      RC_CONVERT_TAB[(unsigned char)seq[k - 3]];
    h_val ^= TRIMER_TAB[trimer_loc];
  } else if (remainder == 2) {
    uint8_t dimer_loc = 4 * RC_CONVERT_TAB[(unsigned char)seq[k - 1]] +
                        RC_CONVERT_TAB[(unsigned char)seq[k - 2]];
    h_val ^= DIMER_TAB[dimer_loc];
  } else if (remainder == 1) {
    h_val ^= SEED_TAB[(unsigned char)seq[k - 1] & CP_OFF];
  }
  for (unsigned i = 0; i < k / 4; i++) {
    h_val = srol(h_val, 4);
    uint8_t curr_offset = 4 * (k / 4 - i) - 1;
    uint8_t tetramer_loc =
      64 * RC_CONVERT_TAB[(unsigned char)seq[curr_offset]] +     // NOLINT
      16 * RC_CONVERT_TAB[(unsigned char)seq[curr_offset - 1]] + // NOLINT
      4 * RC_CONVERT_TAB[(unsigned char)seq[curr_offset - 2]] +
      RC_CONVERT_TAB[(unsigned char)seq[curr_offset - 3]];
    h_val ^= TETRAMER_TAB[tetramer_loc];
  }
  return h_val;
}

/**
 * Perform a roll operation on the reverse-complement by removing char_out and
 * including char_in.
 * @param rh_val Previous reverse-complement hash value computed for the
 * sequence
 * @param k k-mer size
 * @param char_out Character to be removed
 * @param char_in Character to be included
 * @return Rolled hash value for the reverse-complement
 */
inline uint64_t
next_reverse_hash(uint64_t rh_val,
                  unsigned k,
                  unsigned char char_out,
                  unsigned char char_in)
{
  uint64_t h_val = rh_val ^ MS_TAB(char_in & CP_OFF, k);
  h_val ^= SEED_TAB[char_out & CP_OFF];
  h_val = sror(h_val);
  return h_val;
}

/**
 * Perform a roll back operation on the reverse strand.
 * @param rh_val Previous hash value computed for the sequence
 * @param k k-mer size
 * @param char_out Character to be removed
 * @param char_in Character to be included
 * @return Reverse hash value rolled back
 */
inline uint64_t
prev_reverse_hash(uint64_t rh_val,
                  unsigned k,
                  unsigned char char_out,
                  unsigned char char_in)
{
  uint64_t h_val = srol(rh_val);
  h_val ^= SEED_TAB[char_in & CP_OFF];
  h_val ^= MS_TAB(char_out & CP_OFF, k);
  return h_val;
}

} // namespace

namespace nthash {

NtHash::NtHash(const char* seq,
               size_t seq_len,
               typedefs::NUM_HASHES_TYPE num_hashes,
               typedefs::K_TYPE k,
               size_t pos)
  : seq(seq, seq_len)
  , num_hashes(num_hashes)
  , k(k)
  , pos(pos)
  , initialized(false)
  , hash_arr(new uint64_t[num_hashes])
{
  if (k == 0) {
    raise_error("NtHash", "k must be greater than 0");
  }
  if (this->seq.size() < k) {
    raise_error("NtHash",
                "sequence length (" + std::to_string(this->seq.size()) +
                  ") is smaller than k (" + std::to_string(k) + ")");
  }
  if (pos > this->seq.size() - k) {
    raise_error("NtHash",
                "passed position (" + std::to_string(pos) +
                  ") is larger than sequence length (" +
                  std::to_string(this->seq.size()) + ")");
  }
}

bool
NtHash::init()
{
  size_t pos_n = 0;
  while (pos <= seq.size() - k + 1 &&
         is_invalid_kmer(seq.data() + pos, k, pos_n)) {
    pos += pos_n + 1;
  }
  if (pos > seq.size() - k) {
    return false;
  }
  fwd_hash = base_forward_hash(seq.data() + pos, k);
  rev_hash = base_reverse_hash(seq.data() + pos, k);
  extend_hashes(fwd_hash, rev_hash, k, num_hashes, hash_arr.get());
  initialized = true;
  return true;
}

bool
NtHash::roll()
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
  fwd_hash = next_forward_hash(fwd_hash, k, seq[pos], seq[pos + k]);
  rev_hash = next_reverse_hash(rev_hash, k, seq[pos], seq[pos + k]);
  extend_hashes(fwd_hash, rev_hash, k, num_hashes, hash_arr.get());
  ++pos;
  return true;
}

bool
NtHash::roll_back()
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
  } else if (SEED_TAB[(unsigned char)seq[pos - 1]] == SEED_N) {
    return false;
  }
  fwd_hash = prev_forward_hash(fwd_hash, k, seq[pos + k - 1], seq[pos - 1]);
  rev_hash = prev_reverse_hash(rev_hash, k, seq[pos + k - 1], seq[pos - 1]);
  extend_hashes(fwd_hash, rev_hash, k, num_hashes, hash_arr.get());
  --pos;
  return true;
}

bool
NtHash::peek()
{
  if (pos >= seq.size() - k) {
    return false;
  }
  return peek(seq[pos + k]);
}

bool
NtHash::peek(char char_in)
{
  if (!initialized) {
    return init();
  }
  if (SEED_TAB[(unsigned char)char_in] == SEED_N) {
    return false;
  }
  uint64_t fwd = next_forward_hash(fwd_hash, k, seq[pos], char_in);
  uint64_t rev = next_reverse_hash(rev_hash, k, seq[pos], char_in);
  extend_hashes(fwd, rev, k, num_hashes, hash_arr.get());
  return true;
}

bool
NtHash::peek_back()
{
  if (pos == 0) {
    return false;
  }
  return peek_back(seq[pos - 1]);
}

bool
NtHash::peek_back(char char_in)
{
  if (!initialized) {
    return init();
  }
  if (SEED_TAB[(unsigned char)char_in] == SEED_N) {
    return false;
  }
  uint64_t fwd = prev_forward_hash(fwd_hash, k, seq[pos + k - 1], char_in);
  uint64_t rev = prev_reverse_hash(rev_hash, k, seq[pos + k - 1], char_in);
  extend_hashes(fwd, rev, k, num_hashes, hash_arr.get());
  return true;
}

BlindNtHash::BlindNtHash(const char* seq,
                         typedefs::NUM_HASHES_TYPE num_hashes,
                         typedefs::K_TYPE k,
                         size_t pos)
  : seq(seq, seq + k)
  , num_hashes(num_hashes)
  , pos(pos)
  , hash_arr(new uint64_t[num_hashes])
{
  if (k == 0) {
    raise_error("BlindNtHash", "k must be greater than 0");
  }
  fwd_hash = base_forward_hash(seq, k);
  rev_hash = base_reverse_hash(seq, k);
  extend_hashes(fwd_hash, rev_hash, k, num_hashes, hash_arr.get());
}

void
BlindNtHash::roll(char char_in)
{
  fwd_hash = next_forward_hash(fwd_hash, seq.size(), seq.front(), char_in);
  rev_hash = next_reverse_hash(rev_hash, seq.size(), seq.front(), char_in);
  extend_hashes(fwd_hash, rev_hash, seq.size(), num_hashes, hash_arr.get());
  seq.pop_front();
  seq.push_back(char_in);
  ++pos;
}

void
BlindNtHash::roll_back(char char_in)
{
  fwd_hash = prev_forward_hash(fwd_hash, seq.size(), seq.back(), char_in);
  rev_hash = prev_reverse_hash(rev_hash, seq.size(), seq.back(), char_in);
  extend_hashes(fwd_hash, rev_hash, seq.size(), num_hashes, hash_arr.get());
  seq.pop_back();
  seq.push_front(char_in);
  --pos;
}

void
BlindNtHash::peek(char char_in)
{
  uint64_t fwd = next_forward_hash(fwd_hash, seq.size(), seq.front(), char_in);
  uint64_t rev = next_reverse_hash(rev_hash, seq.size(), seq.front(), char_in);
  extend_hashes(fwd, rev, seq.size(), num_hashes, hash_arr.get());
}

void
BlindNtHash::peek_back(char char_in)
{
  uint64_t fwd = prev_forward_hash(fwd_hash, seq.size(), seq.back(), char_in);
  uint64_t rev = prev_reverse_hash(rev_hash, seq.size(), seq.back(), char_in);
  extend_hashes(fwd, rev, seq.size(), num_hashes, hash_arr.get());
}

}