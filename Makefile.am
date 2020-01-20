bin_PROGRAMS = nttest

nttest_CPPFLAGS = -I$(top_srcdir)/lib

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

dist_doc_DATA = \
	ChangeLog \
	CITATION.bib \
	LICENSE \
	README.md

EXTRA_DIST=autogen.sh

# Check the C++ source code for white-space errors with clang-format.
clang-format:
	for i in *pp; do clang-format -style=file $$i >$$i.fixed; done
	for i in *pp; do diff -su $$i $$i.fixed && rm -f $$i.fixed; done
	if ls *.fixed; then exit 1; fi
