To build the test suite (nttest) do:
```
$ make
```
The nttest suite has the options for runtime test and uniformity test. 

For the runtime test the program has the following options:
```
./nttest [OPTIONS] ... [FILE]
```
Description:
* -k,  --kmer=SIZE : the length of k-mer used for runtime test hashing [50].
* -h,  --hash=SIZE : the number of generated hashes for each k-mer [1].
* -f ,--fastq: to run the test on fastq file.
* FILE: is the fasta or fastq file of sequences.

For example to evaluate the runtime of different hash methods on a fasta file reads.fa for k-mer length 50, run:
```
$ ./nttest -k50 reads.fa 
```
For a fastq file reads.fq, just add --fastq option, run:
```
$ ./nttest -k50 --fastq reads.fq 
```

For the uniformity test using the Bloom filter data structure the program has the following options:
```
./nttest --uniformity [OPTIONS] ... [REF_FILE] [QUERY_FILE]
```

Description:
* -q, --qnum=SIZE: number of queries in query file.
* -l, --qlen=SIZE: length of reads in query file.
* -t, --tnum=SIZE: number of sequences in reference file.
* -g, --tlen=SIZE: length of reference sequence.
* -i, --input: generate random query and reference files.
* -j, threads=SIZE: number of threads to run uniformity test [1].
* REF_FILE: the reference file name.
* QUERY_FILE: the query file name.

For example, to evaluate the uniformity of different hash methods using the Bloom filter data structure on randomly generated data sets with following options:
* 100 genes of length 5,000,000bp as reference in file genes.fa
* 4,000,000 reads of length 250bp as query in file reads.fa
* 12 threads
run:
```
$ ./nttest --uniformity --input -q4000000 -l250 -t100 -g5000000 -j12 genes.fa reads.fa 
```
