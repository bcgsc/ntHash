#include <string>
#include <cstring>
#include <cassert>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <getopt.h>
#include "rolling.h"
#include "city.h"

/* as per instructions at top of stb_image_write.h */
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

using namespace std;

#define PROGRAM "collision"

static const char USAGE_MESSAGE[] =
"Usage:" PROGRAM " [OPTIONS] [FASTA]...\n"
"Measure number of hash collisions using CityHash64 or rolling hash.\n"
"K-mers from [FASTA] are hashed to positions within a fixed-size\n"
"window (using the modulus operator).\n"
"\n"
"Options:\n"
"     -d, --max-density N     stop hashing k-mers when\n"
"                             DISTINCT_KMERS / WINDOW_SIZE > N [0.05]\n"
"     -h, --help              show this message\n"
"     -H, --hash HASH         hash function to use ('city' or 'rolling')\n"
"                             [rolling]\n"
"     -k, --kmer-size         size of k-mer [20]\n"
"     -p, --png FILE          write bitmap of hash positions\n"
"                             to FILE; dimensions of output PNG are\n"
"                             approximately sqrt(WINDOW_SIZE) by \n"
"                             sqrt(WINDOW_SIZE)\n"
"     -P, --progress-step N   show progress message after hashing\n"
"                             every N k-mers [10000]\n"
"     -v, --verbose           show progress messages\n"
"     -w, --window-size N     window size; this should normally\n"
"                             be a prime number [1048573]\n";

namespace opt {
	static float maxDensity = 0.05f;
	static int help = 0;
	static string hashFunc("rolling");
	static unsigned k = 20;
	static string pngPath;
	static size_t progressStep = 10000;
	static int verbose = 0;
	static size_t windowSize = 1048573;
}

static const char shortopts[] = "d:hH:k:p:P:vw:";

static const struct option longopts[] = {
	{ "max-density", required_argument, NULL, 'd' },
	{ "help", no_argument, NULL, 'h' },
	{ "hash", required_argument, NULL, 'H' },
	{ "kmer-size", required_argument, NULL, 'k' },
	{ "png", required_argument, NULL, 'p' },
	{ "progress-step", required_argument, NULL, 'P' },
	{ "verbose", no_argument, NULL, 'v' },
	{ "window-size", required_argument, NULL, 'w' },
	{ NULL, 0, NULL, 0 }
};

#define RGB_FORMAT 3
#define BYTES_PER_PIXEL 3

/** Dimensions and pixel data for output PNG file */
struct PNG {

	int width;
	int height;
	/* 3 bytes per pixel: R,G,B */
	void* data;

	void setPixel(unsigned x, unsigned y, uint8_t r, uint8_t g, uint8_t b) {
		size_t offset = (y * width + x) * BYTES_PER_PIXEL;
		uint8_t* bytes = (uint8_t*)data;
		bytes[offset] = r;
		bytes[offset + 1] = g;
		bytes[offset + 2] = b;
	}
};

typedef unordered_map<uint64_t, uint8_t> HashCountMap;

static inline bool hashSeq(string seq, unordered_set<string>& kmers,
	HashCountMap& hashCounts)
{
	transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
	uint64_t fwdHash, rcHash, hash;
	string prevKmer;
	for (unsigned i = 0; i < seq.length() - opt::k + 1; ++i) {
		string kmer = seq.substr(i, opt::k);
		size_t pos = kmer.find_first_not_of("ACGT");
		if (pos != string::npos) {
			i += pos;
			prevKmer.clear();
			continue;
		}
		if (opt::hashFunc == "city") {
			hash = CityHash64(kmer.c_str(), kmer.length());
		} else {
			/* hashFunc == "rolling" */
			if (prevKmer.empty()) {
				hash = initHashes(kmer, fwdHash, rcHash);
			} else {
				hash = rollHashesRight(fwdHash, rcHash,
						prevKmer.front(), kmer.back(), opt::k);
			}
		}
		hash %= opt::windowSize;
		kmers.insert(kmer);
		hashCounts[hash]++;
		prevKmer = kmer;
		if ((float)kmers.size() / opt::windowSize > opt::maxDensity)
			return false;
		if (opt::verbose && kmers.size() % opt::progressStep == 0) {
			cerr << "hashed " << kmers.size() << " k-mers" << endl;
		}
	}
	return true;
}

static inline void drawPNG(HashCountMap& hashCounts, PNG& png)
{
	/* white background */
	memset(png.data, 255, png.width * png.height * BYTES_PER_PIXEL);

	for(HashCountMap::iterator it = hashCounts.begin();
		it != hashCounts.end(); ++it) {
		uint64_t hash = it->first;
		uint8_t count = it->second;
		unsigned x = hash % png.width;
		unsigned y = hash / png.width;
		if (y >= png.height) {
			cerr << "y >= png.height!" << endl;
			cerr << "hash: " << hash << endl;
			cerr << "y: " << y << endl;
			cerr << "png.width: " << png.width << endl;
			cerr << "png.height: " << png.height << endl;
		}
		assert(y < png.height);
		assert(count > 0);
		if (count > 1) {
			/* red to indicate collision */
			png.setPixel(x, y, 255, 0, 0);
		} else {
			/* black to indicate non-collision */
			png.setPixel(x, y, 0, 0, 0);
		}
	}
}

int main(int argc, char** argv)
{
	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
		shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case 'd':
				arg >> opt::maxDensity; break;
			case 'h':
				cout << USAGE_MESSAGE;
				exit(EXIT_SUCCESS);
			case 'H':
				arg >> opt::hashFunc; break;
			case 'k':
				arg >> opt::k; break;
			case 'p':
				arg >> opt::pngPath; break;
			case 'P':
				arg >> opt::progressStep; break;
			case 'v':
				opt::verbose++; break;
			case 'w':
				arg >> opt::windowSize; break;
			case '?':
				die = true; break;
		}
		if (optarg != NULL && !arg.eof()) {
			cerr << PROGRAM ": invalid option: `-"
				<< (char)c << optarg << "'\n";
			exit(EXIT_FAILURE);
		}
	}

	if (opt::maxDensity < 0.0f || opt::maxDensity > 1.0f) {
		cerr << "error: value for -d (--max-density) must be in "
			" range [0,1]" << endl;
		die = true;
	}

	if (opt::hashFunc != "rolling" && opt::hashFunc != "city") {
		cerr << "error: unrecognized argument for -h (--hash); "
			" legal values are: 'rolling', 'city'" << endl;
	}

	if (die) {
		cerr << USAGE_MESSAGE;
		exit(EXIT_FAILURE);
	}


	ifstream fin;
	if (argc > optind)
		fin.open(argv[optind]);
	istream& in = (argc > optind) ? fin : cin;
	assert(in);

	/* k-mers inserted so far */
	unordered_set<string> kmers;
	/* number of occurrences of each hash value (usually 1) */
	HashCountMap hashCounts;


	/* read and hash FASTA sequences */
	string line;
	while(getline(in, line)) {
		if (line.empty() || line.at(0) == '>')
			continue;
		/* return false when window exceeds opt::maxDensity */
		if(!hashSeq(line, kmers, hashCounts))
			break;
	}

	/* count collisions */
	uint64_t collisions = kmers.size() - hashCounts.size();
	cout << "distinct k-mers hashed: " << kmers.size() << endl;
	cout << "hash collisions: " << collisions << endl;

	/* create PNG image */
	if (!opt::pngPath.empty()) {
		PNG png;
		png.width = (int)ceil(sqrt(opt::windowSize));
		png.height = (int)ceil((double)opt::windowSize / png.width);
		png.data = malloc(png.width * png.height * BYTES_PER_PIXEL);
		drawPNG(hashCounts, png);
		stbi_write_png(opt::pngPath.c_str(), png.width, png.height,
			RGB_FORMAT, png.data, png.width*BYTES_PER_PIXEL);
	}

	return 0;
}
