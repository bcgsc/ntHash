#include <iostream>
#include <string>
#include <vector>
#include "stHashIterator.hpp"

using namespace std;

int main()
{
    /* test sequence */
    std::string seq = "GAGTGTCAAACATTCAGACAACAGCAGGGGTGCTCTGGAATCCTATGTGAGGAACAAACATTCAGGCCACAGTAG";
    
   
    /* h is the number of hashes for each k-mer */
    unsigned h = 4;
    
    std::vector<std::string> seedSeq(h);
    seedSeq[0]="110001100111000001110100110110100010011101";
    seedSeq[1]="000001111100101110111000001011010100011110";
    seedSeq[2]="011110001010110100000111011101001111100000";
    seedSeq[3]="101110010001011011001011100000111001100011";
    
    //cerr << seedSeq[3] << endl;
//    ,"110011110011101","111110110011101","100111110011101"};
    
    stHashIterator ssitr(seq, seedSeq, h, seedSeq[0].size());
    
    while (ssitr != ssitr.end()) {
        std::cout << (*ssitr)[0] << "\t" << (*ssitr)[1] <<"\t" << (*ssitr)[2] << "\t" <<(*ssitr)[3] <<"\t" << std::endl;
        std::cout << (ssitr.strand())[0] << "\t" << (ssitr.strand())[1] <<"\t" << (ssitr.strand())[2] << "\t" <<(ssitr.strand())[3] <<"\t" << std::endl;

        ++ssitr;
    }
    
    std::cout << endl << endl;
    
    
    return 0;
}