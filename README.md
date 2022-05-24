[![Release](https://img.shields.io/github/release/bcgsc/ntHash.svg)](https://github.com/bcgsc/ntHash/releases)
[![Downloads](https://img.shields.io/github/downloads/bcgsc/ntHash/total?logo=github)](https://github.com/bcgsc/ntHash/archive/master.zip)
[![Issues](https://img.shields.io/github/issues/bcgsc/ntHash.svg)](https://github.com/bcgsc/ntHash/issues)

![Logo](nthash-logo.png)

ntHash 
=
ntHash is a recursive hash function for hashing all possible k-mers in a DNA/RNA sequence.

# Installation & Usage

## btllib (Recommended)

ntHash is included in [btllib](https://github.com/bcgsc/btllib), which can be installed using Conda:

```shell
conda install -c bioconda btllib
```

This installs other classes such as Bloom filters that can be used with ntHash.  Also, ntHash can be imported in Python projects via btllib's wrappers. ntHash's code is automatically synced to btllib on each release. Refer to [btllib's readme](https://github.com/bcgsc/btllib/blob/master/README.md) for usage instructions.

## Manual Compilation

To use ntHash without installing btllib from conda, download the repo and generate a cmake buildsystem in an arbitrary directory (e.g. `release`), by running the following command in the project's root:

```shell
cmake -S . -B release
```

Then, build the project and its dependencies using:

```shell
cmake --build release --target all
```

The project's executable (`./nthash`), static library (`libnthash.a`), and tests (`nthash_tests`) will be generated alongside `btllib` in the `release` folder.

Tests can be run using `ctest`:

```shell
cd release && ctest
```

### Executable

After building the project, use the output binary file to generate hash values from input data and store the results on disk:

```
Usage: nthash [options] files 

Positional arguments:
files           Input sequence files [required]

Optional arguments:
-v --version    prints version information and exits [default: false]
-k              k-mer size [required]
-o              Output file (for -f collect) or directory path [required]
-f              Output file organization (create files containing hashes for each 'file', 'record', or 'collect' all hashes into a single file [default: "file"]
-h              Number of hashes per k-mer [default: 1]
-s              Input spaced seed patterns separated by commas (e.g. 1110111,11011011). Performs k-mer hashing if no value provided.
--long          Optimize file reader for long sequences (>5kbp) [default: false]
--binary        Output hashes in binary files (otherwise plain text) [default: false]
--verbose       Print progress to stdout [default: false]
```

For example, given two input files `1.fa` and `2.fa`, two hash values for each 64-mer can be saved to a binary file called `out.bin` by running this command in `release`:

```shell
./nthash -k 64 -h 2 -o out.bin -f collect --binary --verbose 1.fa 2.fa
```

In the plain text format (tab separated values), the rows consist of the hashes of k-mers/seeds in the same order seen in the input sequences. For binary output, these values are dumped without any delimiters. E.g. in the case of `out.bin` generated above, the first and second hashes of the first k-mer are stored in the beginning of the file, followed by the first and second hashes of the second k-mer.

### Static Library Usage

To use ntHash in a C++ project:
+ Link the code with `libnthash.a` (i.e. pass `-L path/to/nthash/release -l nthash` to the compiler).
+ Add the `include` directory (pass `-I path/to/nthash/include` to the compiler).
+ Repeat for btllib (add flags: `-L path/to/nthash/release/btllib -l btllib -I path/to/nthash/vendor/btllib/include`)
+ Import ntHash in the code using `#include <nthash/nthash.hpp>`.
+ Access ntHash classes from the `nthash` namespace.

# Examples

Generally, the `nthash::NtHash` and `nthash::SeedNtHash` classes are used for hashing sequences:

```C++
nthash::NtHash nth("TGACTGATCGAGTCGTACTAG", 1, 5);  // 1 hash per 5-mer
while (nth.roll()) {
    // use nth.hashes() for canonical hashes
    //     nth.get_forward_hash() for forward strand hashes
    //     nth.get_reverse_hash() for reverse strand hashes
}
```

```C++
std::vector<std::string> seeds = {"10001", "11011"};
nthash::SeedNtHash nth("TGACTGATCGAGTCGTACTAG", seeds, 3, 5);
while (nth.roll()) {
    // nth.hashes()[0] = "T###T"'s first hash
    // nth.hashes()[1] = "T###T"'s second hash
    // nth.hashes()[2] = "T###T"'s third hash
    // nth.hashes()[3] = "TG#CT"'s first hash
}
```

Refer to docs for more information.


Publications
============

## [ntHash](http://bioinformatics.oxfordjournals.org/content/early/2016/08/01/bioinformatics.btw397)

Hamid Mohamadi, Justin Chu, Benjamin P Vandervalk, and Inanc Birol.
**ntHash: recursive nucleotide hashing**.
*Bioinformatics* (2016) 32 (22): 3492-3494.
[doi:10.1093/bioinformatics/btw397 ](http://dx.doi.org/10.1093/bioinformatics/btw397)


# Acknowledgements

+ ntHash: [Hamid Mohamadi](https://github.com/mohamadi)
+ ntHash2: [Parham Kazemi](https://github.com/parham-k)
+ btllib: [Vladimir Nikolic](https://github.com/vlad0x00)
+ [argparse for C++](https://github.com/p-ranav/argparse)
