/*
 *
 * test.hpp
 * Author: Justin Chu
 * Genome Sciences Centre,
 * British Columbia Cancer Agency
 */


#include <string>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <stdint.h>
#include <unistd.h>
#include <getopt.h>
#include <cassert>
#include <zlib.h>
#include "kseq.h"
#include "RollingHashIterator.h"
KSEQ_INIT(gzFile, gzread)

#ifdef _OPENMP
#include <omp.h>
#endif

#define PROGRAM "ssnttest"

using namespace std;

static const char VERSION_MESSAGE[] =
    PROGRAM " Version 1.0.0 \n"
    "Written by Justin Chu & Hamid Mohammadi.\n"
    "Copyright 2018 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
    "Usage: " PROGRAM " [OPTION]... QUERY\n"
    "Report bugs to cjustin@bcgsc.ca.\n";

namespace opt {
vector<string> sseeds;
bool fastq = false;
}

using namespace std;

static const char shortopts[] = "s:";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "sseed",	required_argument, NULL, 's' },
    { "help",	no_argument, NULL, OPT_HELP },
    { "version",	no_argument, NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

/*
 * Parses input string into separate strings, returning a vector.
 */
static vector<string> convertInputString(const string &inputString) {
	vector<string> currentString;
	string temp;
	stringstream converter(inputString);
	while (converter >> temp) {
		currentString.push_back(temp);
	}
	assert(currentString.size() > 0);
	return currentString;
}

void computeHashByMasking(const char * file,
		const vector<string> &sseeds) {
#ifdef _OPENMP
	double sTime = omp_get_wtime();
#else
	clock_t start = clock();
#endif
//populate sdsl bitvector (bloomFilter)
	gzFile fp;
	fp = gzopen(file, "r");
	kseq_t *seq = kseq_init(fp);
	int l;
//#pragma omp parallel private(l)
	for (;;) {
		string sequence;
//#pragma omp critical(seq)
		{
			l = kseq_read(seq);
			if (l >= 0) {
				sequence = string(seq->seq.s, seq->seq.l);
			}
		}
		if (l >= 0) {
			for (RollingHashIterator itr(sequence, sseeds[0].length(),
					RollingHash::parseSeedString(sseeds)); itr != itr.end();
					++itr) {
			}
		} else {
			break;
		}
	}
	kseq_destroy(seq);
	gzclose(fp);

#ifdef _OPENMP
	std::cerr << "hash_time=" <<setprecision(4) << fixed << omp_get_wtime() - sTime << "\n";
#else
	std::cerr << "hash_time=" << setprecision(4)
			<< (double) (clock() - start) / CLOCKS_PER_SEC << "\n";
#endif
}

int main(int argc, char** argv) {

    bool die = false;
    for (int c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
        case '?':
            die = true;
            break;
        case 's':
            opt::sseeds = convertInputString(optarg);
            break;
        case OPT_HELP:
            std::cerr << USAGE_MESSAGE;
            exit(EXIT_SUCCESS);
        case OPT_VERSION:
            std::cerr << VERSION_MESSAGE;
            exit(EXIT_SUCCESS);
        }
		//causing errors with -S command
//		if (optarg != NULL && !arg.eof()) {
//			cerr << PROGRAM ": invalid option: `-" << (char) c << optarg
//					<< "'\n";
//			exit(EXIT_FAILURE);
//		}
    }
    if (argc - optind != 1 && argc - optind != 2) {
        std::cerr << PROGRAM ": missing arguments\n";
        die = true;
    }

    if (opt::sseeds.empty()){
    	cerr << "need seed pattern" << endl;
    	die = true;
    }

    if (die) {
        std::cerr << "Try `" << PROGRAM
                  << " --help' for more information.\n";
        exit(EXIT_FAILURE);
    }

	//Stores fasta input file names
    const char *readName(argv[argc-1]);

	computeHashByMasking(readName, opt::sseeds);

    return 0;
}
