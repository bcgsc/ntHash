#include <iostream>
#include <string>
#include <vector>
#include "ntHashIterator.hpp"

using namespace std;

int main(int argc, const char* argv[])
{
    /* test sequence */
    std::string seq = "GAGTGTCAAACATTCAGACAACAGCAGGGGTGCTCTGGAATCCTATGTGAGGAACAAACATTCAGGCCACAGTAG";
    
    /* k is the k-mer length */
    unsigned k = 42;
    
    /* h is the number of hashes for each k-mer */
    unsigned h = 4;
    
    std::vector<std::string> seedSeq(h);
    seedSeq[0]="110001100111000001110100110110100010011101";
    seedSeq[1]="000001111100101110111000001011010100011110";
    seedSeq[2]="011110001010110100000111011101001111100000";
    seedSeq[3]="101110010001011011001011100000111001100011";
    
    //cerr << seedSeq[3] << endl;
//    ,"110011110011101","111110110011101","100111110011101"};
    
    ntHashIterator ssitr(seq, seedSeq, h, k);
    
    while (ssitr != ssitr.end()) {
        std::cout << (*ssitr)[0] << "\t" << (*ssitr)[1] <<"\t" << (*ssitr)[2] << "\t" <<(*ssitr)[3] <<"\t" << std::endl;
        --ssitr;
    }
    
    std::cout << endl << endl;
    
    uint64_t* hVal = new uint64_t[4];
    for(unsigned i=0; i<seq.length()-k+1;i++) {
        string kmer = seq.substr(i,k);
        NTMS64_test(kmer.data(), seedSeq, k, h, hVal);
        std::cout << hVal[0] << "\t" << hVal[1] <<"\t" << hVal[2] << "\t" << hVal[3] <<"\t" << std::endl;
    }
        
    
    return 0;
}