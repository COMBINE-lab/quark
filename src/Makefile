RobSeriator: ReadSeriatorRob.cpp
	g++ -I./ -O3 -w -std=c++11 Kmer.cpp hash.cpp ReadSeriatorRob.cpp cigargen.cc -o RobSeriator -lz

ReadSeriatorChunk: ReadSeriatorChunk.cpp
	g++ -I./ -Ofast -march=native -funroll-loops -std=c++11 Kmer.cpp hash.cpp ReadSeriatorChunk.cpp cigargen.cc xxhash.c -o ReadSeriatorChunk -lz

ReadSeriatorBothpair: ReadSeriatorBothpair.cpp
	g++ -I./ -Ofast -march=native -funroll-loops -std=c++11 Kmer.cpp hash.cpp ReadSeriatorBothpair.cpp cigargen.cc xxhash.c -o ReadSeriatorBothpair -lz

graphtrav: randomtrav.cpp
	g++ -g -std=c++11 -I../boost_1_60_0/ randomtrav.cpp -o randomtrav -L../boost_1_60_0/lib -lboost_filesystem -lboost_system

decoder: decoder.cpp
	g++ -g --std=c++11 -I./ -I../boost_1_60_0/ decoder.cpp -o decoder -L=../boost_1_60_0/lib -lboost_filesystem -lboost_system
