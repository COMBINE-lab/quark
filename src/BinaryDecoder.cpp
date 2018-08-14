#include <thread>
#include <string>
#include <atomic>
#include <vector>
#include <fstream>
#include <iostream>
#include <map>
#include <iterator>
#include <functional>
#include <cmath>
#include <mutex>

#include <boost/dynamic_bitset.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "pstream.h"
#include "spdlog/details/format.h"

namespace bfs = boost::filesystem ;

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
    if(s.size()%2 != 0)
        s[halfLength] = complement(s[halfLength]);
    return s;
}

std::string decode(std::string& island, std::string& encread,uint32_t pos){

	//debug

	int ind = 0;
	std::string real = "";
	std::string read = "";
	if(encread[0] == '0'){
		//std::vector<std::string> vec ;
		//vec = split(encread,':');
		std::string digit = "";
		int digitInd = 1;
		char c;
		c = encread[digitInd];
		while(std::isdigit(c) != 0){
			digit.append(std::string(1,encread[digitInd]));
			digitInd += 1;
			c = encread[digitInd];
		}
		int overhang;
		//read = vec[1];
		read = encread.substr(digitInd,encread.size()-digitInd+1);
		overhang = std::stoi(digit);
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
		if(std::isdigit(read[ind])){
			std::string digit = "";
			int digitInd = ind;
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
			//debug

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



char codes[16] = {
		'|','A','C','G',
		'T','N','0','1',
		'2','3','4','5',
		'6','7','8','9'
};


void readCompressedSingle(std::string &ifnamestr, std::string& ofname, bool qualityScore){

		bfs::path ifname(ifnamestr) ;

	    bfs::path seqPathLeft = ifname / "reads.quark.lz";
		bfs::path offsetPathLeft = ifname / "offset.quark.lz";
		bfs::path islandPath = ifname / "islands.txt";
		bfs::path qualityPath = ifname / "quality.quark";

		std::cout << "Sequence File : { " << seqPathLeft << " }\n";
		std::cout << "Offset File: { " << offsetPathLeft << " }\n";
		std::cout << "Island file : { " << islandPath << " }\n";
		if(qualityScore)
			std::cout << "Quality file : { " << qualityPath << " }\n";

		//set up output directory
		bfs::path outdir(ofname);
		if(!(bfs::exists(outdir))){
			if(bfs::create_directory(ofname))
				std::cout << "Output would be written in "<<outdir<<"\n";
		}

		std::ifstream iFile(islandPath.string());


		std::ifstream qualityFile;
		if(qualityScore){
			qualityFile.open(qualityPath.string(),std::ofstream::in);
		}


		std::ofstream seqLeftOut;
		seqLeftOut.open(ofname + "/1.fastq");



		bfs::path unmapped_1 = ifname / "um_.fa";
		std::ifstream um1(unmapped_1.string());

		fmt::MemoryWriter w;
		w.write("plzip -d -c -n {} {}", 2, seqPathLeft);
		redi::ipstream seqLeft(w.str());
		w.clear();

		w.write("plzip -d -c -n {} {}", 2, offsetPathLeft);
		redi::ipstream offLeft(w.str());
		w.clear();

		int numOfEquivClasses = 0;
		iFile >> numOfEquivClasses;
		int numOfReads = 0;



		std::string qScore ;

		for(int eqNum = 0; eqNum < numOfEquivClasses; eqNum++){

			printProgress(double(eqNum)/double(numOfEquivClasses));

			int numOfIslands = 0;
			iFile >> numOfIslands;
			std::vector<std::string> islands;
			islands.resize(numOfIslands);

			for(int i=0 ; i < numOfIslands; i++){
				iFile >> islands[i];
			}

			bool isImpLeft = false;
			//bool isImpRight = false;

			char leftOverHang;
			//char rightOverHang;

			uint32_t numOfCodes{0};
			seqLeft.read(reinterpret_cast<char*>(&numOfCodes), sizeof(numOfCodes));
			uint32_t codeCount{0};

			while(codeCount < numOfCodes){
				std::string leftEnc{""};
				uint16_t leftIsland{0};
				uint32_t leftPos{0};


				offLeft.read(reinterpret_cast<char*>(&leftIsland), sizeof(leftIsland));
				offLeft.read(reinterpret_cast<char*>(&leftPos), sizeof(leftPos));


				bool endIt = false;

				if(isImpLeft){
					leftEnc.append(std::string(1,leftOverHang));
					isImpLeft = false;
				}

				while(!endIt){
					uint8_t temp{0};
					seqLeft.read(reinterpret_cast<char*>(&temp), sizeof(temp));
					uint8_t f1{0};
					uint8_t f2{0};
					f2 = temp >> 4;
					f1 = temp & 0x0f ;
					if(codes[f1] == '|' || codes[f2] == '|'){
						//we should end it here
						//because the line ends here
						//std::cout << "\n should not be here for first time\n" ;
						if(codes[f2] != '|'){
							isImpLeft = true;
							leftOverHang = codes[f2];
						}else{
							if(codes[f1] != '|')
								leftEnc.append(std::string(1,codes[f1]));
							isImpLeft = false;
						}



						endIt = true;
						codeCount++;
					}else{
						leftEnc.append(std::string(1,codes[f1]));
						leftEnc.append(std::string(1,codes[f2]));
					}
				}

				std::string ldecoded ;

				if(islands[leftIsland] == "$"){
					ldecoded = leftEnc ;
				}else{
					std::string left = "";
					left = leftEnc.substr(0,leftEnc.size()-1);
					ldecoded = decode(islands[leftIsland],left,leftPos);
					std::string ore = std::string(1,leftEnc[leftEnc.size()-1]);
					ldecoded = (ore == "0") ? revComp(ldecoded) : ldecoded ;
				}
				//seqLeftOut << ldecoded << "\n" ;

				//read the quality file and then write one by one


				seqLeftOut << "@" << numOfReads << "\n";
				seqLeftOut << ldecoded << "\n" ;
				seqLeftOut << "+" << "\n";
				if(qualityScore){
					qualityFile >> qScore ;
					seqLeftOut << qScore << "\n";
				}
				else{
					std::string quality(ldecoded.size(),'I');
					seqLeftOut << quality << "\n" ;
				}
				numOfReads++;
			}

		}

		std::string content;

		std::string head;
		while(um1 >> head){
			um1 >> content;
			qualityFile >> qScore ;
			seqLeftOut << "@" << numOfReads << "\n";
			seqLeftOut << content << "\n" ;
			seqLeftOut << "+" << "\n";
			if(qualityScore){
				qualityFile >> qScore ;
				seqLeftOut << qScore << "\n";
			}
			else{
				std::string quality(content.size(),'I');
				seqLeftOut << quality << "\n" ;
			}
			numOfReads++;
		}
}

void readCompressed(std::string &ifname, std::string &ofname, bool qualityScore){

	    bfs::path seqPathLeft = ifname + "/reads_1.quark.lz";
		bfs::path seqPathRight = ifname + "/reads_2.quark.lz";
		bfs::path offsetPathLeft = ifname + "/offset_1.quark.lz";
		bfs::path offsetPathRight = ifname + "/offset_2.quark.lz";
		bfs::path qualityPathLeft = ifname + "quality_1.quark";
		bfs::path qualityPathRight = ifname + "quality_2.quark";

		std::string islandPath = ifname + "/islands.txt";



		std::string unmapped_1 = ifname + "/um_1.fa";
		std::string unmapped_2 = ifname + "/um_2.fa";
		std::ifstream um1(unmapped_1);
		std::ifstream um2(unmapped_2);


		std::cout << "Left end : { " << seqPathLeft << " }\n";
		std::cout << "Right end : { " << seqPathRight << " }\n";
		std::cout << "Left Offset : { " << offsetPathLeft << " }\n";
		std::cout << "Right Offset : { " << offsetPathRight << " }\n";
		std::cout << "Island file : { " << islandPath << " }\n";
		if(qualityScore){
					std::cout << "Left quality file : { " << qualityPathLeft << " }\n";
					std::cout << "Right quality file : { " << qualityPathRight << " }\n";
		}

		std::ifstream iFile(islandPath);

		//set up output directory
		bfs::path outdir(ofname);
		if(!(bfs::exists(outdir))){
			if(bfs::create_directory(ofname))
				std::cout << "Output would be written in "<<outdir<<"\n";
		}
		std::ofstream seqLeftOut;
		seqLeftOut.open(ofname + "1.fastq");
		std::ofstream seqRightOut;
		seqRightOut.open(ofname + "2.fastq");

		std::ifstream qualityFile1 ;
		std::ifstream qualityFile2 ;
		if(qualityScore){
			qualityFile1.open(qualityPathLeft.string(), std::ofstream::in);
			qualityFile2.open(qualityPathRight.string(), std::ofstream::in);
		}




		fmt::MemoryWriter w;
		w.write("plzip -d -c -n {} {}", 2, seqPathLeft.string());
		redi::ipstream seqLeft(w.str());
		w.clear();


		w.write("plzip -d -c -n {} {}", 2, seqPathRight.string());
		redi::ipstream seqRight(w.str());
		w.clear();

		w.write("plzip -d -c -n {} {}", 2, offsetPathLeft.string());
		redi::ipstream offLeft(w.str());
		w.clear();


		w.write("plzip -d -c -n {} {}", 2, offsetPathRight.string());
		redi::ipstream offRight(w.str());
		w.clear();

		int numOfEquivClasses = 0;
		iFile >> numOfEquivClasses;
		int numOfReads = 0;

		std::string qScore;



		for(int eqNum = 0; eqNum < numOfEquivClasses; eqNum++){

			printProgress(double(eqNum)/double(numOfEquivClasses));

		    //read the islands file
			int numOfIslands = 0;
			iFile >> numOfIslands;
			std::vector<std::string> islands;
			islands.resize(numOfIslands);

			for(int i=0 ; i < numOfIslands; i++){
				iFile >> islands[i];
			}
		    //read the file from the left end and right end and corresponding
		    //offset files
			//std::vector<std::string> leftEnd ;
			//std::vector<std::string> rightEnd ;


			bool isImpLeft = false;
			bool isImpRight = false;

			char leftOverHang;
			char rightOverHang;

			uint32_t numOfCodes{0};
			seqLeft.read(reinterpret_cast<char*>(&numOfCodes), sizeof(numOfCodes));
			uint32_t codeCount{0};

			while(codeCount < numOfCodes){

				numOfReads++;

				std::string leftEnc{""};
				std::string rightEnc{""};
				uint16_t leftIsland{0};
				uint32_t leftPos{0};
				uint16_t rightIsland{0};
				uint32_t rightPos{0};

				// read the position and island ids
				// we would just plug things in

				offLeft.read(reinterpret_cast<char*>(&leftIsland), sizeof(leftIsland));
				offLeft.read(reinterpret_cast<char*>(&leftPos), sizeof(leftPos));
				offRight.read(reinterpret_cast<char*>(&rightIsland), sizeof(rightIsland));
				offRight.read(reinterpret_cast<char*>(&rightPos), sizeof(rightPos));


				bool endIt = false; // end delimiter
				//check if something is left from the previous call
				if(isImpLeft){
					leftEnc.append(std::string(1,leftOverHang));
					isImpLeft = false;
				}

				while(!endIt){
					uint8_t temp{0};
					seqLeft.read(reinterpret_cast<char*>(&temp), sizeof(temp));
					uint8_t f1{0};
					uint8_t f2{0};
					f2 = temp >> 4;
					f1 = temp & 0x0f ;
					if(codes[f1] == '|' || codes[f2] == '|'){
						//we should end it here
						//because the line ends here
						//std::cout << "\n should not be here for first time\n" ;
						if(codes[f2] != '|'){
							isImpLeft = true;
							leftOverHang = codes[f2];
						}else{
							if(codes[f1] != '|')
								leftEnc.append(std::string(1,codes[f1]));
							isImpLeft = false;
						}




						endIt = true;
						codeCount++;
					}else{
						leftEnc.append(std::string(1,codes[f1]));
						leftEnc.append(std::string(1,codes[f2]));
					}
				}

				if(isImpRight){
					rightEnc.append(std::string(1,rightOverHang));
					isImpRight = false;
				}

				endIt = false;
				while(!endIt){
					uint8_t temp{0};
					seqRight.read(reinterpret_cast<char*>(&temp), sizeof(temp));
					uint8_t f1{0};
					uint8_t f2{0};
					f2 = temp >> 4;
					f1 = temp & 0x0f ;
					if(codes[f1] == '|' || codes[f2] == '|'){
						//we should end it here
						//because the line ends here
						if(codes[f2] != '|'){
							isImpRight = true;
							rightOverHang = codes[f2];
						}else{
							if(codes[f1] != '|')
								rightEnc.append(std::string(1,codes[f1]));
							isImpRight = false;
						}

						endIt = true;
					}else{
						rightEnc.append(std::string(1,codes[f1]));
						rightEnc.append(std::string(1,codes[f2]));
					}
				}

				std::string ldecoded;
				std::string rdecoded;


				if(islands[leftIsland] == "$"){
					ldecoded = leftEnc ;
				}else{
					std::string left = "";
					left = leftEnc.substr(0,leftEnc.size()-1);
					ldecoded = decode(islands[leftIsland],left,leftPos);
					std::string ore = std::string(1,leftEnc[leftEnc.size()-1]);
					ldecoded = (ore == "0") ? revComp(ldecoded) : ldecoded ;
				}

				if(islands[rightIsland] == "$"){
					rdecoded = rightEnc ;
				}else{
					std::string right = "";
					right = rightEnc.substr(0,rightEnc.size()-1);
					rdecoded = decode(islands[rightIsland],right,rightPos);
					std::string ore = std::string(1,rightEnc[rightEnc.size()-1]);
					rdecoded = (ore == "0") ? revComp(rdecoded) : rdecoded ;
				}


				seqLeftOut << "@" << numOfReads << "\n";
				seqLeftOut << ldecoded << "\n" ;
				seqLeftOut << "+" << "\n";
				if(qualityScore){
					qualityFile1 >> qScore ;
					seqLeftOut << qScore << "\n" ;
				}else{
					std::string quality(ldecoded.size(),'I');
					seqLeftOut << quality << "\n";
				}

				seqRightOut << "@" << numOfReads << "\n";
				seqRightOut << rdecoded << "\n" ;
				seqRightOut << "+" << "\n";
				if(qualityScore){
					qualityFile2 >> qScore ;
					seqRightOut << qScore << "\n" ;
				}else{
					std::string quality(rdecoded.size(),'I');
					seqRightOut << quality << "\n";
				}
				//seqRightOut << quality << "\n";

			}
		}

		int numOfReads2 = numOfReads;
		std::string content;
		std::string head;
		while(um1 >> head){
			um1 >> content;
			std::string quality(content.size(),'I');
			seqLeftOut << "@" << numOfReads << "\n";
			seqLeftOut << content << "\n" ;
			seqLeftOut << "+" << "\n";
			if(qualityScore){
				qualityFile1 >> qScore ;
				seqLeftOut << qScore << "\n" ;
			}else{
				//std::string quality(ldecoded.size(),'I');
				seqLeftOut << quality << "\n";
			}
			//seqLeftOut << quality << "\n";
			numOfReads++;
		}
		//numOfReads = 0;
		while(um2 >> head){
			um2 >> content;
			std::string quality(content.size(),'I');
			seqRightOut << "@" << numOfReads2 << "\n";
			seqRightOut << content << "\n" ;
			seqRightOut << "+" << "\n";
			if(qualityScore){
				qualityFile2 >> qScore ;
				seqRightOut << qScore << "\n" ;
			}else{
				//std::string quality(ldecoded.size(),'I');
				seqRightOut << quality << "\n";
			}
			//seqRightOut << quality << "\n";
			numOfReads2++;
		}



}


int main(int argc, char* argv[]){
	std::vector<std::string> args(argv+1,argv+argc);
	bool qualityScore{false} ;

	if(argc < 4){
		if(args[0] == "--help"){
			std::cout << "\n Correct syntax is \n"
					  << "./decoder <input dir> <output dir> <P/S> <Q/N>" ;
			return 0;
		}else{
			std::cout << "\n Invoked incorrectly \n" ;
			return 0;
		}
	}
	std::cout << "\n Starting Quark Decoder Module ...\n";
	std::cout << "\n Input directory : " << args[0];
	std::cout << "\n Output directory : " << args[1];


	if(args[3] == "Q"){
		std::cout << "\n Quality will be added ";
		qualityScore = true ;
	}else{
		std::cout << "\n No quality score file found ";
		qualityScore = false ;
	}
	//std::cout << "\n Library type : " << args[2];

	if(args[2] == "S"){
		std::cout << "\n Library type : SINGLE END\n" ;
	}else{
		std::cout << "\n Library type : PAIRED END\n" ;
	}
	//std::cout << "\n Number of theads : " << args[4] << "\n\n";

	if(args[2] != "S"){
		readCompressed(args[0],args[1],qualityScore);

	}else{
		readCompressedSingle(args[0],args[1],qualityScore);
	}

	std::cout<< "\n";
	return 0;

}
