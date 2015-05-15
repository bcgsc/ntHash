#include <stdio.h>
#include <iostream>
#include <stdlib.h>

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

void printBin(uint64_t n) {
    int count=0,count1=0;
    while (++count<=64) {
        if (n & ((uint64_t)1<<63)) {
            printf("1 ");
            ++count1;
        }
        else
            printf("0 ");

        n <<= 1;
    }
    printf(" count1=%d\n",count1);
}

int main(int argc, const char* argv[]) {
    srand(time(NULL));
    const int seedNum=4;
    int hashSeed[seedNum][64];
    for (int i=0; i<seedNum; i++)
        for (int j=0; j<64; j++)
            hashSeed[i][j]=0;

    for (int j=0; j<64; j++) {
        int ranVec[seedNum];
        for (int i=0; i<seedNum; i++)
            ranVec[i]=i;
        for (int i=0; i<seedNum/2; i++) {
            int ranInd = rand() % seedNum;
            int tmp = ranVec[seedNum-i-1];
            ranVec[seedNum-i-1] = ranVec[ranInd];
            ranVec[ranInd] = tmp;
        }
        for (int i=0; i<seedNum/2; i++)
            hashSeed[ranVec[i]][j] = 1;
    }

    for (int i=0; i<seedNum; i++) {
        for (int j=0; j<64; j++)
            std::cout << hashSeed[i][j] << " ";
        std::cout << "\n";
    }

    for (int i=0; i<seedNum; i++) {
        uint64_t hSeed=0;
        for (int j=0; j<64; j++) {
            if(hashSeed[i][j]==1)
                hSeed = (uint64_t) (hSeed << 1 | 1);
            else
                hSeed = (uint64_t) (hSeed << 1 | 0);
        }
        printf("%" PRIu64 "\n", hSeed);
        printBin(hSeed);

    }
    return 0;
}
