SHELL=/bin/sh

# turn off suffix-based rules
.SUFFIXES:

TESTS:=revcomp
revcomp_SOURCES:=revcomp.cpp
bloom_HEADERS:=../rolling.h

.PHONY: clean run_tests
default: $(TESTS) run_tests
clean:
	rm -f $(TESTS)

revcomp: $(revcomp_SOURCES) $(revcomp_HEADERS)
	$(CXX) -I.. $(CPPFLAGS) $(CXXFLAGS) $(revcomp_SOURCES) -o $@

run_tests:
	./revcomp
