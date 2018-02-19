#ifndef ROLLING_HASH_ITERATOR_H
#define ROLLING_HASH_ITERATOR_H 1

#include <cstring>
#include <vector>
#include <cassert>
#include <limits>
#include <string>
#include "RollingHash.h"

/**
 * Permitted characters in k-mers. All k-mers containing
 * other characters will be skipped.
 */
#define ACGT_CHARS "acgtACGT"

//TODO: implement with string iterator rather than m_pos

/**
 * Iterate over hash values for k-mers in a
 * given DNA sequence.
 *
 * This implementation uses a rolling hash
 * function to efficiently calculate
 * hash values for successive k-mers.
 */
class RollingHashIterator
{
private:

	/**
	 * Advance iterator right to the next valid k-mer.
	 */
	void next()
	{
		if (m_seq.length() < m_k) {
			m_pos = std::numeric_limits<std::size_t>::max();
			return;
		}

		while(m_pos < m_seq.length() - m_k + 1) {
			/* skip over k-mers with non-ACGT chars */
			if (m_nextInvalidChar - m_pos < m_k) {
				m_pos = m_nextInvalidChar + 1;
				m_nextInvalidChar = m_pos + strspn(m_seq.c_str() + m_pos, ACGT_CHARS);
				m_rollNextHash = false;
			} else {
				/* we are positioned at the next valid k-mer */
				if (!m_rollNextHash) {
					/* we don't have hash values for the
					 * preceding k-mer, so we must compute
					 * the hash values from scratch */
					m_rollingHash.reset(m_seq.begin() + m_pos);
					m_rollNextHash = true;
				} else {
					/* compute new hash values based on
					 * hash values of preceding k-mer */
					assert(m_pos > 0);
					m_rollingHash.rollRight(m_seq.begin() + m_pos - 1);
				}
				return;
			}
		}
		/* there are no more valid k-mers */
		m_pos = std::numeric_limits<std::size_t>::max();
	}

public:

	/**
	 * Default constructor. Creates an iterator pointing to
	 * the end of the iterator range.
	 */
	RollingHashIterator() :
			m_k(0), m_rollingHash(m_k, vector<vector<unsigned> >()), m_rollNextHash(
					false), m_start(true), m_pos(
					std::numeric_limits<std::size_t>::max()), m_nextInvalidChar(
					0) {
	}

	/**
	 * Constructor.
	 * @param seq DNA sequence to be hashed
	 * @param k k-mer size
	 * @param numHashes number of hash values to compute
	 * for each k-mer
	 */
	RollingHashIterator(const std::string& seq, unsigned k,
			const vector<vector<unsigned> > &ssVals) :
			m_seq(seq), m_k(k), m_rollingHash(m_k, ssVals), m_rollNextHash(
					false), m_start(true), m_pos(0) {
		//TURN OFF INVALID CHAR CHECKING
//		m_nextInvalidChar = strspn(m_seq.c_str(), ACGT_CHARS);
		next();
	}

	/** get reference to hash values for current k-mer */
	const uint64_t* operator*() const
	{
		assert(m_pos + m_k <= m_seq.length());
		return m_rollingHash.getHash();
	}

	/** test equality with another iterator */
	bool operator==(const RollingHashIterator& it) const
	{
		return m_pos == it.m_pos;
	}

	/** test inequality with another iterator */
	bool operator!=(const RollingHashIterator& it) const
	{
		return !(*this == it);
	}

	/** pre-increment operator */
	RollingHashIterator& operator++()
	{
		++m_pos;
		next();
		return *this;
	}

	/** post-increment operator */
	RollingHashIterator operator++(int)
	{
		RollingHashIterator &it = *this;
		++*this;
		return it;
	}

	/** increments and get reference to hash values for current k-mer
	 *	Returns empty hash when after last element
	 **/
	const uint64_t* getNext() {
		if (m_pos >= m_seq.length() - m_k) {
			return NULL;
		}
		if (m_start) {
			m_start = false;
			return m_rollingHash.getHash();
		}
		++m_pos;
		next();
		return m_rollingHash.getHash();
	}

	/** iterator pointing to one past last element */
	static const RollingHashIterator end() {
		return RollingHashIterator();
	}

	/** return k-mer at current position */
	std::string kmer() const
	{
		return std::string(m_seq, m_pos, m_k);
	}

	/** return position of current k-mer */
    size_t pos() const{
    	return m_pos;
    }

private:

	/** DNA sequence being hashed */
	const std::string m_seq;
	/** k-mer size */
	unsigned m_k;
	/** number of hash values to compute for each k-mer */
//	unsigned m_numHashes;
	/** internal state for rolling hash */
	RollingHash m_rollingHash;
	/** true whenever we can "roll" the hash values for
	 * the current k-mer to compute the hash values for the
	 * next k-mer */
	bool m_rollNextHash;
	/** completely new start of string used for getNext() **/
	bool m_start;
//	/** position of current k-mer */
	size_t m_pos;
	/** position of next non-ACGT char */
	size_t m_nextInvalidChar;
};

#endif
