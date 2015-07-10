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
"                             be a prime number [1048573]\n"
"\n"
"Sample Prime Numbers (for window size):\n"
"\n"
"     1048573\n"
"     2047541\n"
"     3055471\n"
"     4051051\n"
"     5051143\n";

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

static char complement[256] = {
	0,  0,  0,  0,  0,  0,  0,  0, // 0..7
	0,  0,  0,  0,  0,  0,  0,  0, // 8..15
	0,  0,  0,  0,  0,  0,  0,  0, // 16..23
	0,  0,  0,  0,  0,  0,  0,  0, // 24..31
	0,  0,  0,  0,  0,  0,  0,  0, // 32..39
	0,  0,  0,  0,  0,  0,  0,  0, // 40..47
	0,  0,  0,  0,  0,  0,  0,  0, // 48..55
	0,  0,  0,  0,  0,  0,  0,  0, // 56..63
	0, 84,  0, 71,  0,  0,  0, 67, // 64..71 (A,C,G)
	0,  0,  0,  0,  0,  0,  0,  0, // 72..79
	0,  0,  0,  0, 65,  0,  0,  0, // 80..87 (T)
	0,  0,  0,  0,  0,  0,  0,  0, // 88..95
	0,116,  0,103,  0,  0,  0, 99, // 96..103 (a,c,g)
	0,  0,  0,  0,  0,  0,  0,  0, // 104..111
	0,  0,  0,  0, 97,  0,  0,  0, // 112..119 (t)
	0,  0,  0,  0,  0,  0,  0,  0, // 120..127
	0,  0,  0,  0,  0,  0,  0,  0, // 128..135
	0,  0,  0,  0,  0,  0,  0,  0, // 136..143
	0,  0,  0,  0,  0,  0,  0,  0, // 144..151
	0,  0,  0,  0,  0,  0,  0,  0, // 152..159
	0,  0,  0,  0,  0,  0,  0,  0, // 160..167
	0,  0,  0,  0,  0,  0,  0,  0, // 168..175
	0,  0,  0,  0,  0,  0,  0,  0, // 176..183
	0,  0,  0,  0,  0,  0,  0,  0, // 184..191
	0,  0,  0,  0,  0,  0,  0,  0, // 192..199
	0,  0,  0,  0,  0,  0,  0,  0, // 200..207
	0,  0,  0,  0,  0,  0,  0,  0, // 208..215
	0,  0,  0,  0,  0,  0,  0,  0, // 216..223
	0,  0,  0,  0,  0,  0,  0,  0, // 224..231
	0,  0,  0,  0,  0,  0,  0,  0, // 232..239
	0,  0,  0,  0,  0,  0,  0,  0, // 240..247
	0,  0,  0,  0,  0,  0,  0,  0  // 248..255
};

static inline void canonicalize(string& seq)
{
	unsigned k = seq.length();
	string rc(k, 'N');
	for (unsigned i = 0; i < k; ++i) {
		unsigned char rcChar = complement[seq.at(i)];
		rc.at(k-i-1) = rcChar;
		if (seq.at(i) != rcChar) {
			if (seq.at(i) < rcChar) {
				return;
			} else {
				/* finish constructing reverse complement */
				for (unsigned j = i + 1; j < k; ++j)
					rc.at(k-j-1) = complement[seq.at(j)];
				seq = rc;
				return;
			}
		}
	}
}

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
					prevKmer.at(0), kmer.at(opt::k-1), opt::k);
			}
			/*
			 * The rolling hash returns the same hash value for
			 * both orientations of a k-mer. In order to count
			 * the number of collisions correctly, we must store
			 * only the canonical version of each k-mer.
			 */
			canonicalize(kmer);
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
	cout << "distinct_kmers\tcollisions\n";
	cout << kmers.size() << "\t" << collisions << "\n";

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
