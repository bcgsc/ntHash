SHELL=/bin/sh

# turn off suffix-based rules
.SUFFIXES:

CXXFLAGS:=-O3 -fopenmp $(CXXFLAGS)

PROGRAMS:=bloom collision
bloom_SOURCES:=bloom.cpp
bloom_HEADERS:=../rolling.h
collision_SOURCES:=lib/city.cc collision.cpp
collision_HEADERS:=../rolling.h lib/city.h

.PHONY: clean
default: $(PROGRAMS)
clean:
	rm -f $(PROGRAMS)

bloom: $(bloom_SOURCES) $(bloom_HEADERS)
	$(CXX) -I.. $(CPPFLAGS) $(CXXFLAGS) $(bloom_SOURCES) -o $@

collision: $(collision_SOURCES) $(collision_HEADERS)
	$(CXX) -I.. -Ilib $(CPPFLAGS) $(CXXFLAGS) $(collision_SOURCES) -o $@
