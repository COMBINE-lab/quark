CXX=g++
CXXFLAGS=-g -std=c++11 -Wall -pedantic
all:kseq.h kseq_test.cpp randomtrav.cpp
		$(CXX) $(CXXFLAGS) kseq_test.cpp -o kseq_test -lz
		$(CXX) $(CXXFLAGS) randomtrav.cpp -o randomtrav

clean:
		rm -f *.o kseq_test
		rm -f *.o randomtrav
