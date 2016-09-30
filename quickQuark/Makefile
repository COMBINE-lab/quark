#g++ -Iinclude -Lfmt/build/fmt/ -g -O2 -pthread -std=c++11 BinaryDecoder.cpp -o decoder -lfmt -lboost_system
CC=g++
IDIR=include
LDIR=fmt/build/fmt/
CFLAG= -g -pthread -std=c++11
LFLAG= -lfmt -lboost_system

all: BinaryDecoder.cpp
	$(CC) -I$(IDIR) -L$(LDIR) $(CFLAG) BinaryDecoder.cpp -o decoder $(LFLAG)

clean:
	$(RM) decoder
