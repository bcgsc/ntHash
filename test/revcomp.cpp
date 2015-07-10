#include <string>
#include <iostream>
#include "rolling.h"

using namespace std;

static int g_test_num = 1;

static inline bool test(const string& desc, bool passed)
{
	cerr << "TEST " << g_test_num << ": "
		<< desc << "... " <<
		(passed ? "PASS" : "FAIL")
		<< endl;
	g_test_num++;

	return passed;
}

int main(int argc, char** argv)
{
	bool all_passed = true;

	/*
	 * test sequence:      GCAATGT
	 * reverse complement: ACATTGC
	 * k-mer size = 6
	 */

	const unsigned k = 6;
	string kmer1("GCAATG");
	string kmer1_rc("CATTGC");
	string kmer2("CAATGT");
	string kmer2_rc("ACATTG");
	
	uint64_t fwdHash, rcHash;
	uint64_t kmer1_hash = initHashes(kmer1, fwdHash, rcHash);
	uint64_t kmer1_rc_hash = initHashes(kmer1_rc, fwdHash, rcHash);

	if (!test("hash(kmer) == hash(rc(kmer))",
		kmer1_hash == kmer1_rc_hash))
		all_passed = false;

	uint64_t rolled_kmer1_hash = initHashes(kmer1, fwdHash, rcHash);
	rolled_kmer1_hash = rollHashesRight(fwdHash, rcHash, 'G', 'T', k);
	uint64_t kmer2_hash = initHashes(kmer2, fwdHash, rcHash);

	if (!test("rollRight(hash(kmer1)) == hash(kmer2)",
		rolled_kmer1_hash == kmer2_hash))
		all_passed = false;
	
	uint64_t hash_1 = initHashes(kmer1, fwdHash, rcHash);
	uint64_t hash_2 = rollHashesRight(fwdHash, rcHash, 'G', 'T', k);
	uint64_t rc_hash_1 = initHashes(kmer2_rc, fwdHash, rcHash);
	uint64_t rc_hash_2 = rollHashesRight(fwdHash, rcHash, 'A', 'C', k);
	
	if (!test("seq and rc(seq) hash values agree",
		hash_1 == rc_hash_2 && hash_2 == rc_hash_1))
		all_passed = false;

	return !all_passed;
}
