#define PROGRAM_NAME "ntHash"
#define VERSION "2.0"

#include "argparse/argparse.hpp"
#include "btllib/seq_reader.hpp"
#include "nthash/nthash.hpp"
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

enum OutputOrganization
{
  FILES,
  RECORDS,
  COLLECT
};

struct CommandLineArguments
{

  unsigned k;
  unsigned h;
  std::string out_path;
  std::vector<std::string> file_paths;
  std::vector<std::string> seeds;
  bool long_mode, binary_output, verbose;
  OutputOrganization out_org;

  CommandLineArguments(int argc, char** argv)
  {
    auto default_args = argparse::default_arguments::version;
    argparse::ArgumentParser parser(PROGRAM_NAME, VERSION, default_args);

    parser.add_argument("-k")
      .help("k-mer size")
      .required()
      .scan<'u', unsigned>();
    parser.add_argument("-o")
      .help("Output file (for -f collect) or directory path")
      .required();
    parser.add_argument("-f")
      .help("Output file organization (store hashes for each "
            "'file', 'record', or 'collect' all hashes into a single file")
      .default_value(std::string("file"));
    parser.add_argument("-h")
      .help("Number of hashes per k-mer/seed")
      .default_value((unsigned)1)
      .scan<'u', unsigned>();
    parser.add_argument("-s").help("Input spaced seed patterns separated by "
                                   "commas (e.g. 1110111,11011011). Performs "
                                   "k-mer hashing if no value provided.");
    parser.add_argument("files")
      .help("Input sequence files")
      .required()
      .remaining();
    parser.add_argument("--long")
      .help("Optimize file reader for long sequences (>5kbp)")
      .default_value(false)
      .implicit_value(true);
    parser.add_argument("--binary")
      .help("Output hashes in binary files (otherwise plain text)")
      .default_value(false)
      .implicit_value(true);
    parser.add_argument("--verbose")
      .help("Print progress to stdout")
      .default_value(false)
      .implicit_value(true);

    try {
      parser.parse_args(argc, argv);
    } catch (const std::runtime_error& err) {
      std::cerr << err.what() << std::endl;
      std::cerr << parser;
      std::exit(1); // NOLINT
    }

    this->k = parser.get<unsigned>("-k");
    this->h = parser.get<unsigned>("-h");
    this->long_mode = parser.get<bool>("--long");
    this->binary_output = parser.get<bool>("--binary");
    this->verbose = parser.get<bool>("--verbose");

    this->out_path = parser.get("-o");
    std::string out_org_str = parser.get("-f");
    if (out_org_str == "record") {
      this->out_org = OutputOrganization::RECORDS;
    } else if (out_org_str == "file") {
      this->out_org = OutputOrganization::FILES;
    } else if (out_org_str == "collect") {
      this->out_org = OutputOrganization::COLLECT;
    } else {
      std::cerr << "Invalid output organization: " << out_org_str << std::endl;
      std::cerr << parser;
      std::exit(1); // NOLINT
    }

    if (parser.is_used("-s")) {
      std::istringstream seed_parser(parser.get("-s"));
      std::string seed;
      while (std::getline(seed_parser, seed, ',')) {
        this->seeds.emplace_back(seed);
        btllib::check_error(!this->seeds.empty() &&
                              seed.size() != this->seeds[0].size(),
                            "All input seeds must have the same length");
      }
    }

    try {
      this->file_paths = parser.get<std::vector<std::string>>("files");
    } catch (std::logic_error& e) {
      std::cerr << "No files provided" << std::endl;
      std::exit(1); // NOLINT
    }
  }
};

class Logger
{
private:
  bool verbose;

public:
  Logger(bool verbose)
    : verbose(verbose)
  {}

  void print_parameter_report(const CommandLineArguments& args) const
  {
    if (verbose) {
      std::cout << PROGRAM_NAME << " " << VERSION << std::endl << std::endl;
      if (args.seeds.empty()) {
        std::cout << "Using k = " << args.k << ", " << args.h
                  << " hash(es) per k-mer" << std::endl
                  << std::endl;
      } else {
        std::cout << "Spaced seed patterns to be used:" << std::endl;
        for (const auto& seed : args.seeds) {
          std::cout << seed << std::endl;
        }
        std::cout << args.h << " hashes per spaced seed" << std::endl
                  << std::endl;
      }
    }
  }

  void print(const std::string& message, const bool endl = true) const
  {
    if (verbose) {
      std::cout << message;
      if (endl) {
        std::cout << std::endl;
      } else {
        std::cout << std::flush;
      }
    }
  }
};

class HashContainer
{
private:
  std::vector<const uint64_t*> hashes;

  void write_to_binary_file(const unsigned num_hashes, const std::string& path)
  {
    std::fstream out;
    out.open(path, std::ios::binary);
    uint64_t hash_value;
    for (const auto& hash_array : hashes) {
      for (unsigned i = 0; i < num_hashes; i++) {
        hash_value = hash_array[num_hashes];
        out.write(reinterpret_cast<char*>(&hash_value), sizeof(hash_value));
      }
    }
    out.close();
  }

  void write_to_text_file(const unsigned num_hashes, const std::string& path)
  {
    std::ofstream out;
    out.open(path);
    for (const auto& hash_array : hashes) {
      for (unsigned i = 0; i < num_hashes; i++) {
        out << hash_array[i] << '\t';
      }
      out << std::endl;
    }
    out.close();
  }

public:
  void reserve(const unsigned long seq_len, const unsigned long num_hashes)
  {
    hashes.reserve(hashes.size() + seq_len * num_hashes);
  }

  void add(const uint64_t* hash_array, const unsigned array_length)
  {
    auto* copy_array = new uint64_t[array_length];
    std::copy(hash_array, hash_array + array_length, copy_array);
    hashes.push_back(copy_array);
  }

  void write_to_file(const unsigned num_hashes,
                     const std::string& path,
                     const bool binary)
  {
    if (binary) {
      write_to_binary_file(num_hashes, path);
    } else {
      write_to_text_file(num_hashes, path);
    }
    hashes.clear();
  }
};

inline void
hash_kmers(const std::string& seq,
           unsigned num_hashes,
           unsigned kmer_size,
           HashContainer& hashes)
{
  nthash::NtHash nth(seq, num_hashes, kmer_size);
  while (nth.roll()) {
    hashes.add(nth.hashes(), nth.get_hash_num());
  }
}

inline void
hash_spaced_seed(const std::string& seq,
                 const std::vector<std::string>& seeds,
                 unsigned num_hashes,
                 HashContainer& hashes)
{
  nthash::SeedNtHash nth(seq, seeds, num_hashes, seeds[0].size());
  while (nth.roll()) {
    hashes.add(nth.hashes(), nth.get_hash_num());
  }
}

int
main(int argc, char** argv)
{
  CommandLineArguments args(argc, argv);
  Logger logger(args.verbose);

  unsigned reader_flag;
  if (args.long_mode) {
    reader_flag = btllib::SeqReader::Flag::LONG_MODE;
  } else {
    reader_flag = btllib::SeqReader::Flag::SHORT_MODE;
  }

  logger.print_parameter_report(args);

  HashContainer hashes;
  std::chrono::duration<long double> tms{};
  std::string file_type = args.binary_output ? "binary" : "text";
  for (const auto& file_path : args.file_paths) {
    logger.print("Opening " + file_path);
    btllib::SeqReader reader(file_path, reader_flag);
    for (const auto& record : reader) {
      logger.print("\tHashing sequence " + record.id + "...", false);
      hashes.reserve(record.seq.size(), args.h);
      if (args.seeds.empty()) {
        auto tp0 = std::chrono::steady_clock::now();
        hash_kmers(record.seq, args.h, args.k, hashes);
        auto tp1 = std::chrono::steady_clock::now();
        tms = std::chrono::duration_cast<std::chrono::seconds>(tp1 - tp0);
      } else {
        auto tp0 = std::chrono::steady_clock::now();
        hash_spaced_seed(record.seq, args.seeds, args.h, hashes);
        auto tp1 = std::chrono::steady_clock::now();
        tms = std::chrono::duration_cast<std::chrono::seconds>(tp1 - tp0);
      }
      logger.print(" DONE (" + std::to_string(tms.count()) + "s)");
      if (args.out_org == OutputOrganization::RECORDS) {
        std::filesystem::path dir(args.out_path);
        std::filesystem::path file(record.id);
        std::string out_path = dir / file;
        hashes.write_to_file(args.h, out_path, args.binary_output);
        logger.print("\tOutput " + file_type + " file: " + out_path);
      }
    }
    if (args.out_org == OutputOrganization::FILES) {
      std::filesystem::path dir(args.out_path);
      std::filesystem::path file(file_path);
      std::string out_path = dir / file;
      hashes.write_to_file(args.h, out_path, args.binary_output);
      logger.print("\tOutput " + file_type + " file: " + out_path);
    }
  }
  if (args.out_org == OutputOrganization::COLLECT) {
    hashes.write_to_file(args.h, args.out_path, args.binary_output);
    logger.print("Output " + file_type + " file: " + args.out_path);
  }

  return 0;
}