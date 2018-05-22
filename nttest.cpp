#include <iostream>
#include <string>
#include <vector>
#include "ntHashIterator.hpp"

using namespace std;

int main()
{
    /* test sequence */
    std::string seq = "GAGTGTCAAACATTCAGACAACAGCAGGGGTGCTCTGGAATCCTATGTGAGGAACAAACATTCAGGCCACAGTAG";
    
    /* h is the number of hashes for each k-mer */
    unsigned h = 4;

    /* k is the length of k-mer */
    unsigned k = 42;

    ntHashIterator itr(seq, h, k);
    
    while (itr != itr.end()) {
        std::cout << (*itr)[0] << "\t" << (*itr)[1] <<"\t" << (*itr)[2] << "\t" <<(*itr)[3] <<"\t" << itr.strand() << std::endl;
        ++itr;
    }
    return 0;
}