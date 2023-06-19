#include <nthash/nthash.hpp>

#include <chrono>
#include <fstream>
#include <string>

int
main(int argc, char* argv[])
{
  std::ifstream file(argv[1]);
  std::string seq;
  auto t1 = std::chrono::system_clock::now();
  while (file >> seq) {
    uint64_t sum = 0;
    nthash::NtHash h(seq, 3, 64);
    while (h.roll()) {
      sum += h.hashes()[2];
    }
  }
  auto t2 = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed = (t2 - t1);
  std::cout << elapsed.count() << std::endl;
  return 0;
}