bin_PROGRAMS = sttest sstest

sttest_CPPFLAGS = -I$(top_srcdir)/lib

sttest_SOURCES = \
	nthash.hpp \
	stHashIterator.hpp \
	sttest.cpp

sstest_CPPFLAGS = -I$(top_srcdir)/lib

sstest_SOURCES = \
	nthash.hpp \
	ssHashIterator.hpp \
	sstest.cpp

dist_doc_DATA = \
	ChangeLog \
	CITATION.bib \
	LICENSE \
	README.md
	
EXTRA_DIST=autogen.sh
