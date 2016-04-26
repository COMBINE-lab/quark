#include <zlib.h>
#include <sstream>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <fstream>
#include <sys/time.h>
#include <boost/range/irange.hpp>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)
struct read_record{
  std::string name;
  std::string seq;
  std::string quality;
  std::string comment = "+";
  public:
  read_record(std::string &name, std::string &seq, std::string &quality):
      name(name), seq(seq), quality(quality) {

  }
};

//read the new file
int main(int argc,char *argv[]){
    //Start reading from sequence reads from shuffled fastq files i
    //load the tmp file to get the information of number of eq classes
    //and then process them
    //check if the first read id really comed in the top of the shuffled file
    //other wise it needs modification -- fixed the bug

    // read the list containing all numberes
    // you have to read the list in reverse
    //
    std::vector<int> numEqClass ;
    std::ifstream ifs(argv[2]);
    if(ifs.is_open()){
        for(std::string line; std::getline(ifs,line);){
            std::istringstream iss(line);
            int n ;
            iss >> n ;
            numEqClass.push_back(n);
        }
    }
    //reverse the vector
    std::reverse(numEqClass.begin(),numEqClass.end());
    //open the fast q file and start reading in chunks given by numEqClass
    gzFile fp;
    kseq_t *seq;
    int l;
    fp = gzopen(argv[1],"r");
    seq = kseq_init(fp);
    // I don't need kseq for now
    // I can go with ordinary input stream
    std::ifstream fqstream;
    fqstream.open(argv[1], std::ifstream::in);


    std::cout << "start reading from fastq \n";
    std::string leonfile = argv[1] + std::string(".leon.fq");
    std::string mincfile = argv[1] + std::string(".H.fq");
    std::ofstream oleon ;
    std::ofstream ominc ;
    oleon.open(leonfile, std::ofstream::out);
    ominc.open(mincfile, std::ofstream::out);

    int numReads = 0;


    for(auto chunkSize : numEqClass){
        //create a local vector
        std::vector<read_record> chunkReads ;
        //feed it to leon

        if(chunkSize > 1000){
            for(int i=0;i < chunkSize; ++i){
                std::string header;
                std::getline(fqstream,header);
                std::string sequ;
                std::getline(fqstream,sequ);
                std::string com;
                std::getline(fqstream,com);
                std::string qual;
                std::getline(fqstream,qual);
                oleon << header << "\n" << sequ << "\n" << com << "\n" << qual << "\n";
                if (numReads%10000 == 0){
                    std::cout << numReads+1 << " reads processed \r" << std::flush ;
                }
                numReads++;

            }
        }else{
             //don't go for leon
            for(int i=0;i < chunkSize; ++i){
                std::string header;
                std::getline(fqstream,header);
                std::string sequ;
                std::getline(fqstream,sequ);
                std::string com;
                std::getline(fqstream,com);
                std::string qual;
                std::getline(fqstream,qual);
                ominc << header << "\n" << sequ << "\n" << com << "\n" << qual << "\n";
                if (numReads%10000 == 0){
                    std::cout << numReads+1 << " reads processed \r" << std::flush ;
                }
                numReads++;
            }

        }
    }

    ifs.close();
    oleon.close();
    ominc.close();
    return 0;
    //fp = gzopen(argv[1])
}
