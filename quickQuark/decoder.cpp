#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <climits>
#include <cctype>
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

int readCompressed(std::string islandFile,std::string quarkFile,std::string outDir){
	std::ifstream iFile(islandFile);
	std::ifstream qFile(quarkFile);
	std::ofstream lMap(outDir+"/mapped.1");
	std::ofstream rMap(outDir+"/mapped.2");

	int numEqClasses ;
	iFile >> numEqClasses ;
	int eqClass ;
	for(eqClass = 0; eqClass < numEqClasses; eqClass++){
		printProgress(double(eqClass)/double(numEqClasses));
		//read the new class
		//make an vector of islands
		std::vector<std::string> islands ;
		int numOfIslands = 0;
		iFile >> numOfIslands ;
		for(int i = 0;i < numOfIslands;i++){
			std::string s;
			iFile>>s;
			islands.push_back(s);
		}
		//now read the encoded reads
		//and use the islands to decode them
		int numReads = 0;
		qFile>>numReads ;
		for(int i = 0; i < numReads; i++){

			//print for debugging purpose
			//std::cout <<"\n Number of reads " << numReads << "\n";



			std::string lDecoded = "";
			std::string rDecoded = "";
			std::string encodedLine ;

			//std::vector<std::string> code ;
			std::string code_0;
			std::string code_1;
			std::vector<std::string> meta ;
			std::vector<std::string> parts ;


			std::string left = "";
			std::string right = "";
			std::string leftSeq = "" ;
			std::string rightSeq = "" ;

			std::string oreleft ;
			std::string oreright ;
			int lIsland,lPos,rIsland,rPos;

			//qFile>>encodedLine ;
			//std::getline(qFile,encodedLine);
			qFile>>code_0>>code_1;
			//std::cout << code_0 << code_1 << "\n";

			//code = split(encodedLine, '\t');
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

			lMap << lDecoded << '\n';
			rMap << rDecoded << '\n';


		}
	}
}

int main(int argc, char* argv[]){
	std::vector<std::string> args(argv+1,argv+argc);
	std::cout << "\n Starting Quark Decoder Module ...\n";
	std::cout << "\n Island file : " << args[0];
	std::cout << "\n Quark read file : " << args[1];
	std::cout << "\n Output directory : " << args[2] << "\n\n";


	readCompressed(args[0],args[1],args[2]);
	std::cout << "\n";
	return 0;

}
