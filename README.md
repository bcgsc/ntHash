ntHash 
=
ntHash is a recursive hash function for hashing all possible k-mers in a DNA/RNA sequence. 

## Build the test suite
Run:
```
$ make
```
The nttest suite has the options for *runtime* and *uniformity* tests. 

## Runtime test
For the runtime test the program has the following options:
```
./nttest [OPTIONS] ... [FILE]
```
Parameters:
  * `-k`,  `--kmer=SIZE`: the length of k-mer used for runtime test hashing `[50]`
  * `-h`,  `--hash=SIZE`: the number of generated hashes for each k-mer `[1]`
  * `-f` , `--fastq`: to run the test on a fastq file
  * `FILE`: is the input fasta or fastq file

For example to evaluate the runtime of different hash methods on the test file `reads.fa` in DATA/ folder for k-mer length `50`, run:
```
$ nttest -k50 reads.fa 
```

# Uniformity test
For the uniformity test using the Bloom filter data structure the program has the following options:
```
./nttest --uniformity [OPTIONS] ... [REF_FILE] [QUERY_FILE]
```

Parameters:
  * `-q`, `--qnum=SIZE`: number of queries in query file
  * `-l`, `--qlen=SIZE`: length of reads in query file
  * `-t`, `--tnum=SIZE`: number of sequences in reference file
  * `-g`, `--tlen=SIZE`: length of reference sequence
  * `-i`, `--input`: generate random query and reference files
  * `-j`, `threads=SIZE`: number of threads to run uniformity test `[1]`
  * `REF_FILE`: the reference file name
  * `QUERY_FILE`: the query file name

For example, to evaluate the uniformity of different hash methods using the Bloom filter data structure on randomly generated data sets with following options:
  * `100` genes of length `5,000,000bp` as reference in file `genes.fa`
  * `4,000,000` reads of length `250bp` as query in file `reads.fa`
  * `12` threads

run:
```
$ ./nttest --uniformity --input -q4000000 -l250 -t100 -g5000000 -j12 genes.fa reads.fa 
```

## Code samples
To hash all k-mers of length `k` in a given sequence `seq`:
```bash
    string kmer = seq.substr(0, k);
    uint64_t hVal=0;
    hVal = NTP64(kmer.c_str(), k); // initial hash value
    ...
    for (size_t i = 0; i < seq.length() - k; i++) 
    {
        hVal = NTP64(hVal, seq[i], seq[i+k], k); // consecutive hash values
        ...
    }
```
To canonical hash all k-mers of length `k` in a given sequence `seq`:
```bash
    string kmer = seq.substr(0, k);
    uint64_t hVal, fhVal=0, rhVal=0; // canonical, forward, and reverse-strand hash values
    hVal = NTC64(kmer.c_str(), k, fhVal, rhVal); // initial hash value
    ...
    for (size_t i = 0; i < seq.length() - k; i++) 
    {
        hVal = NTC64(fhVal, rhVal, seq[i], seq[i+k], k); // consecutive hash values
        ...
    }
```
To multi-hash with `h` hash values all k-mers of length `k` in a given sequence `seq`: 
```bash
    string kmer = seq.substr(0, k);
    uint64_t hVec[h];
    NTM64(kmer.c_str(), k, h, hVec); // initial hash vector
    ...
    for (size_t i = 0; i < seq.length() - k; i++) 
    {
        NTM64(seq[i], seq[i+k], k, h, hVec); // consecutive hash vectors
        ...
    }
```
