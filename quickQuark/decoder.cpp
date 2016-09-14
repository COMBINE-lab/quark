#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <climits>
#include <cctype>
#include <thread>
#include <mutex>

#include "concurrentqueue.h"
//#include <boost/filesystem.hpp>

using namespace std;
#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

void printProgress (double percentage)
{
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf ("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush (stdout);
}
//using namespace boost::filesystem;

char complement(char& c){
    switch(c){
        case 'A': c = 'T';
                  return c;
        case 'T': c = 'A';
                  return c;
        case 'C': c = 'G';
                  return c;
        case 'G': c = 'C';
                  return c;
        case 'a': c = 't';
                  return c;
        case 't': c = 'a';
                  return c;
        case 'c': c = 'g';
                  return c;
        case 'g': c = 'c';
                  return c;
        default : c = 'N';
                  return c;

    }
}

std::string revComp(std::string s){
    int n = s.size();
    int halfLength = s.size() / 2;
    for (int i=0; i<halfLength; i++)
    {
        char temp = complement(s[i]);
        s[i] = complement(s[n-1-i]);
        s[n-1-i] = temp;
    }
    return s;
}
/*
void split(const string &s, char delim, vector<string> &elems) {
    stringstream ss;
    ss.str(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
}


vector<string> split(const string &s, char delim) {
    vector<string> elems;
    split(s, delim, elems);
    return elems;
}
*/
std::vector<std::string> split(const std::string &text, char sep) {
    std::vector<std::string> tokens;
    std::size_t start = 0, end = 0;
    while ((end = text.find(sep, start)) != std::string::npos) {
        tokens.push_back(text.substr(start, end - start));
        start = end + 1;
    }
    tokens.push_back(text.substr(start));
    return tokens;
}

std::string decode(std::string& island, std::string& encread,int pos){
	int ind = 0;
	std::string real = "";
	std::string read = "";
	if(encread[0] == '-'){
		std::vector<std::string> vec ;
		vec = split(encread,':');
		int overhang;
		read = vec[1];
		overhang = stoi(vec[0]);
		real.append(read.substr(0,abs(overhang)));

		ind = abs(overhang) ;
		//debug
		//std::cout<<read<<"\n";
		//std::cout<<real<<"\n";
		//std::cout<<ind<<"\n";
	}else{
		read = encread;
	}
	while(ind < read.size()){
		if(read[ind] == 'M'){
			std::string digit = "";
			int digitInd = ind+1;
			char c;
			c = read[digitInd];
			while(std::isdigit(c) != 0){
				digit.append(std::string(1,read[digitInd]));
				digitInd += 1;
				if (digitInd == read.size())
					break;
				c = read[digitInd];
			}
			int matches = std::stoi(digit);
			std::string matchChars = island.substr(pos,matches);
			real.append(matchChars);
			pos += matches;
			ind = digitInd;
		}else{
			real.append(std::string(1,read[ind]));
			ind += 1;
			pos += 1;
		}
	}

	return real;
}

struct ReadPair {
    std::string left;
    std::string right;
};

struct EncodedChunk {
    std::vector<std::string> islands;
    std::vector<ReadPair> encodedReads;
};
    
struct OutputWriter{
    std::mutex iomutex;
    std::ofstream lMap;
    std::ofstream rMap;
};

void decodeWorker(moodycamel::ConcurrentQueue<EncodedChunk*>* structQueuePtr, 
                  moodycamel::ConcurrentQueue<EncodedChunk*>* workQueuePtr, 
                  OutputWriter* owriterPtr, 
                  std::atomic<bool>& done) {
    auto& workQueue = *workQueuePtr;
    auto& structQueue = *structQueuePtr;
    auto& owriter = *owriterPtr;
    
    EncodedChunk* ec{nullptr};
    while (workQueue.try_dequeue(ec) or !done) {
        if (ec != nullptr) {
            auto& islands = ec->islands;
            auto nreads = ec->encodedReads.size();
            std::vector<std::string> lDecodedVec(nreads);
            std::vector<std::string> rDecodedVec(nreads);
            for (size_t i = 0; i < nreads; ++i) {
                auto& ep = ec->encodedReads[i];
                std::string& lDecoded = lDecodedVec[i]; 
                std::string& rDecoded = rDecodedVec[i];

                std::vector<std::string> meta ;
                std::vector<std::string> parts ;
			
                std::string left = "";
                std::string right = "";
                std::string leftSeq = "" ;
                std::string rightSeq = "" ;

                std::string oreleft ;
                std::string oreright ;
                int lIsland,lPos,rIsland,rPos;
                auto& code_0 = ep.left;
                auto& code_1 = ep.right;
                parts = split(code_0,'|');

                oreleft = parts[0][parts[0].size()-1];
                left = parts[0].substr(0,parts[0].size()-1);

                oreright = parts[1][parts[1].size()-1];
                right = parts[1].substr(0,parts[1].size()-1);

                meta = split(code_1,',');
                lIsland = std::stoi(meta[0]);
                lPos = std::stoi(meta[1]);
                rIsland = std::stoi(meta[2]);
                rPos = std::stoi(meta[3]);

                if(islands[lIsland] == "$"){
                    lDecoded = parts[0] ;
                }else{
                    lDecoded = decode(islands[lIsland],left,lPos);
                    lDecoded = (oreleft == "0") ? revComp(lDecoded) : lDecoded ;
                }
                if(islands[rIsland] == "$"){
                    rDecoded = parts[1] ;
                }else{
                    rDecoded = decode(islands[rIsland],right,rPos);
                    rDecoded = (oreright == "0") ? revComp(rDecoded) : rDecoded ;
                }
                //exit(0);

                if (lDecoded == "GCGCCTCTCACG"){
                    cout << "\n" << code_0 << "\n";
                    exit(0);
                }
            }

			owriter.iomutex.lock();
            for (size_t i = 0; i < nreads; ++i) {
                owriter.lMap << lDecodedVec[i] << '\n';
                owriter.rMap << rDecodedVec[i] << '\n';
            }
			owriter.iomutex.unlock();
            structQueue.enqueue(ec);
            ec = nullptr;
        }
    }
}



int readCompressed(moodycamel::ConcurrentQueue<EncodedChunk*>& structQueue, 
                   moodycamel::ConcurrentQueue<EncodedChunk*>& workQueue, 
                   std::string islandFile,
                   std::string quarkFile, 
                   std::atomic<bool>& done){
	std::ifstream iFile(islandFile);
	std::ifstream qFile(quarkFile);

	int numEqClasses ;
	iFile >> numEqClasses ;
	int eqClass ;
	for(eqClass = 0; eqClass < numEqClasses; eqClass++){
		printProgress(double(eqClass)/double(numEqClasses));
        EncodedChunk* encChunk{nullptr};
        while(!structQueue.try_dequeue(encChunk)) {
        }
		//read the new class
		//make an vector of islands
		int numOfIslands = 0;
		iFile >> numOfIslands ;
		encChunk->islands.resize(numOfIslands);
		for(int i = 0;i < numOfIslands;i++){
			iFile >> encChunk->islands[i];
		}
		//now read the encoded reads
		//and use the islands to decode them
		int numReads = 0;
		qFile>>numReads ;
        encChunk->encodedReads.resize(numReads);
		for(int i = 0; i < numReads; i++){
            qFile >> encChunk->encodedReads[i].left >> encChunk->encodedReads[i].right;
        }
        workQueue.enqueue(encChunk);
        encChunk = nullptr;
	}
    done = true;
}

int main(int argc, char* argv[]){
	std::vector<std::string> args(argv+1,argv+argc);
	std::cout << "\n Starting Quark Decoder Module ...\n";
	std::cout << "\n Island file : " << args[0];
	std::cout << "\n Quark read file : " << args[1];
	std::cout << "\n Output directory : " << args[2] << "\n\n";

        // must be a power-of-two
    /*
        size_t max_q_size = 2097152;
        spdlog::set_async_mode(max_q_size);

	auto fileSink = std::make_shared<spdlog::sinks::ostream_sink_mt>(logFile);
        auto consoleSink = std::make_shared<spdlog::sinks::stderr_sink_mt>();
        auto consoleLog = spdlog::create("stderrLog", {consoleSink});
        auto fileLog = spdlog::create("fileLog", {fileSink});
        auto jointLog = spdlog::create("jointLog", {fileSink, consoleSink});
    */
        
    std::atomic<bool> done{false};
        size_t nstruct{10000};
        moodycamel::ConcurrentQueue<EncodedChunk*> structQueue(nstruct);
        moodycamel::ConcurrentQueue<EncodedChunk*> workQueue(nstruct);
        for (size_t i = 0; i < nstruct; ++i) {
            structQueue.enqueue(new EncodedChunk);
        }
        std::thread producer(readCompressed, std::ref(structQueue), std::ref(workQueue), args[0],args[1], std::ref(done));

        std::string outDir(args[2]);
        OutputWriter owriter;
        owriter.lMap.open(outDir+"/mapped.1");
        owriter.rMap.open(outDir+"/mapped.2");

        std::vector<std::thread> decoders;
        size_t nthread = 8;
        for (size_t i = 0; i < nthread; ++i) {
            decoders.emplace_back(decodeWorker, &structQueue, &workQueue, &owriter, std::ref(done));
        }

        producer.join();
        for(auto& d : decoders) { d.join(); }
        owriter.lMap.close();
        owriter.rMap.close();
        std::cout << "\n";
        return 0;
}
