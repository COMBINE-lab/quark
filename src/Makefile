RobSeriator: ReadSeriatorRob.cpp
	g++ -I./ -O3 -w -std=c++11 Kmer.cpp hash.cpp ReadSeriatorRob.cpp cigargen.cc -o RobSeriator -lz

ReadSeriatorChunk: ReadSeriatorChunk.cpp
	g++ -I./ -Ofast -march=native -funroll-loops -std=c++11 Kmer.cpp hash.cpp ReadSeriatorChunk.cpp cigargen.cc xxhash.c -o ReadSeriatorChunk -lz

