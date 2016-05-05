CXX=g++
CPPFLAGS=-c -fopenmp
OPTFLAGS=-O3
LIBPATH=-I. -Ilib
LDFLAGS=-fopenmp

EXEC=nttest

all: $(EXEC)

SRCS=nttest.cpp lib/city.cc lib/xxhash.c

$(EXEC):$(SRCS)
	$(CXX) $(OPTFLAGS) $(LDFLAGS) $(LIBPATH) -o $@ $^

clean:
	rm $(EXEC)

