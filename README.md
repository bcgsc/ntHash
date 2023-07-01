[![Release](https://img.shields.io/github/release/bcgsc/ntHash.svg)](https://github.com/bcgsc/ntHash/releases)
[![Downloads](https://img.shields.io/github/downloads/bcgsc/ntHash/total?logo=github)](https://github.com/bcgsc/ntHash/archive/master.zip)
[![Issues](https://img.shields.io/github/issues/bcgsc/ntHash.svg)](https://github.com/bcgsc/ntHash/issues)

![Logo](nthash-logo.png)

ntHash 
=
ntHash is a recursive hash function for hashing all possible k-mers in a DNA/RNA sequence.

# Installation & Usage

Make sure [Meson](https://mesonbuild.com/) is installed on the system.

Download the repo (either from the releases section or close using `git clone https://github.com/bcgsc/ntHash`). Setup meson in an arbitrary directory (e.g. `build`), by running the following command in the project's root:

```shell
meson setup build
```

Then, build the project and its dependencies using:

```shell
cd build && ninja
```

The project's static library (`libnthash.a`) and tests (`nthash`) will be generated alongside in the `build` folder.

Tests can be run using `ninja` in the `build` directory:

```shell
ninja test
```

### Static Library Usage

To use ntHash in a C++ project:
+ Link the code with `libnthash.a` (i.e. pass `-Lpath/to/nthash/build -lnthash` to the compiler).
+ Add the `include` directory (pass `-I path/to/nthash/include` to the compiler).
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

## [ntHash2](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btac564/6674501)
Parham Kazemi, Johnathan Wong, Vladimir Nikolić, Hamid Mohamadi, René L Warren, Inanç Birol, ntHash2: recursive spaced seed hashing for nucleotide sequences, Bioinformatics, 2022;, btac564, [https://doi.org/10.1093/bioinformatics/btac564](https://doi.org/10.1093/bioinformatics/btac564)

## [ntHash](http://bioinformatics.oxfordjournals.org/content/early/2016/08/01/bioinformatics.btw397)

Hamid Mohamadi, Justin Chu, Benjamin P Vandervalk, and Inanc Birol.
**ntHash: recursive nucleotide hashing**.
*Bioinformatics* (2016) 32 (22): 3492-3494.
[doi:10.1093/bioinformatics/btw397 ](http://dx.doi.org/10.1093/bioinformatics/btw397)


# Acknowledgements

+ ntHash: [Hamid Mohamadi](https://github.com/mohamadi)
+ ntHash2: [Parham Kazemi](https://github.com/parham-k)
