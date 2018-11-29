bin_PROGRAMS = nttest nttest_avx

nttest_CPPFLAGS = -I$(top_srcdir)/lib

nttest_avx_CPPFLAGS = --std=c++11 -march=sandybridge -mavx512f -mavx512bw -mavx512vl -I$(top_srcdir)/lib

nttest_SOURCES = \
	lib/BloomFilter.hpp \
	lib/city.cc \
	lib/city.h \
	lib/murmur.hpp \
	lib/seedgen.cpp \
	lib/seqgen.hpp \
	lib/xxhash.c \
	lib/xxhash.h \
	nthash.hpp \
	ntHashIterator.hpp \
	nttest.cpp

nttest_avx_SOURCES = \
	nthash.hpp \
	nthash_avx.hpp \
	nttest_avx.cpp

dist_doc_DATA = \
	ChangeLog \
	CITATION.bib \
	LICENSE \
	README.md

EXTRA_DIST=autogen.sh
