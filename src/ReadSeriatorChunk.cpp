/****************************************************************************
Usage: Takes both ends of the fastq file and the order given by
randomtrav program which is a requirement for this code to run, files needed:
tmpeq.aux and the order file. It also takes the unmapped reads to create
their sequences.

Output:
.h file: file for plzip
.l file: file for leon
****************************************************************************/
#include <zlib.h>
#include <stdio.h>
#include <chrono>
#include <set>
#include <cstring>
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <tuple>
#include <set>
#include <algorithm>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <sys/time.h>
#include "tclap/CmdLine.h"
//#include <Kmer.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/functional/hash.hpp>
//#include <Eigen/Dense>
//#include <hash.hpp>
#include "cigargen.h"
#include "kseq.h"
#include "edlib.h"
#include "xxhash.h"

typedef std::unordered_map<int,std::set<std::string> > EqReadID ;
//using idpos = std::pair<std::string, int> ;
typedef std::unordered_map<std::string, int> RevReadIDEq ;
typedef std::unordered_map<std::string, int> RevReadIDPos ;

struct StringHasher {
    size_t operator()(const std::string& s) const {
        return XXH64(s.c_str(), s.length(), 0);
    }
};

typedef std::pair<std::string,std::string> matepair ;
struct idpos{
    matepair seqlr;
     int pos;
};

//typedef std::unordered_map<int,std::vector<matepairPos> > EqSeq ;
int orderId = 0;

KSEQ_INIT(gzFile, gzread)
struct read_record{
    std::string name;
    std::string seq;
    std::string mateseq;
    std::string quality;
    std::string comment = "+";
    public:
    read_record(std::string &name, std::string &seq, std::string &mateseq, std::string &quality):
        name(name), seq(seq), mateseq(mateseq), quality(quality) {

    }
};


void encodeAsShift(std::vector<idpos> &l,
                   auto& ofs_mapped_l,
                   auto& ofs_mapped_r){

    if (orderId%100 == 0){
                std::cout << orderId << " equivalence class processed \r" << std::flush ;
            }
            ++orderId;

    std::sort(l.begin(), l.end(),
                    [](const idpos &p1, const idpos &p2) -> bool {
                    return p1.pos < p2.pos ;
                    });

    bool refFlag = true ;
    std::string prev;
    std::string curr;
    int prevpos ;
    int currpos ;

      for(auto& p : l){
                    if(refFlag){
                        ofs_mapped_l<<p.seqlr.first<< "\n";
                        prev = p.seqlr.first;
                        prevpos = p.pos ;
                        refFlag = false;
                    }else{
                        curr = p.seqlr.first;
                        currpos = p.pos ;
                        int match = 0;
                        int shift = currpos - prevpos;
                        if(shift < curr.size()){
                            if(shift == 0){
                                ofs_mapped_l<<"M"<<curr.size()<<"\t";
                            }else{
                                ofs_mapped_l<<"S"<<shift<<"\t";
                                while(prev[shift+match] == curr[match])
                                    match++;
                                ofs_mapped_l<<"M"<<match<<"\t";
                                if(match < curr.size())
                                    ofs_mapped_l<<curr.substr(match);
                            }
                            ofs_mapped_l<<"\n";
                        }else{
                            ofs_mapped_l<<curr<<"\n";
                        }
                        prev = curr;
                        prevpos = currpos;
                    }
                    ofs_mapped_r<<p.seqlr.second<<"\n";
        }




}

int main(int argc, char *argv[])
{

    std::cout << "Quark seriator:\n";

    std::string versionString = "0.0.1";
    TCLAP::CmdLine cmd(
            "Quark seriator",
            ' ',
            versionString);
    cmd.getProgramName() = "quarkSeriate";

    TCLAP::ValueArg<std::string> read1("1", "leftMates", "The location of the left paired-end reads", false, "", "path");
    TCLAP::ValueArg<std::string> read2("2", "rightMates", "The location of the right paired-end reads", false, "", "path");
    TCLAP::ValueArg<std::string> orderFile("p", "order", "The location of equivalence class path (i.e. order) file", false, "", "path");
    TCLAP::ValueArg<std::string> sizeFile("s", "sizes", "The location of the equivalence class size file", false, "", "path");
    TCLAP::ValueArg<std::string> unmappedFile("u", "unmapped", "The location of the unmapped read names", false, "", "path");
    TCLAP::ValueArg<uint32_t> numThreads("t", "numThreads", "Number of threads to use", false, 1, "positive integer");
    TCLAP::ValueArg<std::string> outname("o", "output", "The output file (default: stdout)", false, "", "path");

    cmd.add(read1);
    cmd.add(read2);
    cmd.add(unmappedFile);
    cmd.add(outname);
    cmd.add(numThreads);
    cmd.add(orderFile);
    cmd.add(sizeFile);
    /*
	if (argc == 1) {
        fprintf(stderr, "Usage: %s <left.fq> <right.fq> <order> <chunksize> <unmapped file> <outdir>\n", argv[0]);
		return 1;
	}*/


    //read the fastq file because it will be needed for everything
    // TODO: unordered_map is not a smart idea, have to do better
    //
    try {
        cmd.parse(argc, argv);
        std::cout << "\nGathering eqivalence classes \n";
        RevReadIDEq req ;

        std::string eqClassOrderFile(orderFile.getValue());
        std::string eqClassSizeFile(sizeFile.getValue());

        std::ifstream orderFile(eqClassOrderFile);
      if (!orderFile.good()) {
          std::cerr << "Could not open the order file " << eqClassOrderFile << '\n';
          std::exit(1);
      }


        std::ifstream sizeFile(eqClassSizeFile);
        if (!sizeFile.good()) {
          std::cerr << "Could not open the size file " << eqClassSizeFile << '\n';
          std::exit(1);
        }


        gzFile fp1, fp2;
        kseq_t *seq1, *seq2;
        int l1, l2;

        std::string outDir(outname.getValue());
        std::string outFilenameLeft = outDir + "/r1.enc";
        std::string outFilenameRight = outDir + "/r2.enc";
        std::ofstream outFileLeft(outFilenameLeft.c_str());
        std::ofstream outFileRight(outFilenameRight.c_str());

        // process the equivalence class
        // read the size of the next equivalence class
        const size_t chunkSize{12000000};
        size_t classSize{0};
        size_t eqID{0};
        bool done{false};

        //std::cout << "\n read this far\n" ;
        while(!done){
        size_t numToProcess{0};
          struct EqClassInfo {
              std::vector<idpos> readSeqs;
              //std::vector<idpos> readSeqsRight;
              size_t chunkSize;
          };

        std::vector<EqClassInfo> chunkSizes;
          while ((sizeFile >> classSize) and (numToProcess < chunkSize)) {
              chunkSizes.push_back({std::vector<idpos>(), classSize});
              chunkSizes.back().readSeqs.reserve(classSize);
              numToProcess += classSize;
          }
        if (numToProcess < chunkSize) { done = true; }

        std::unordered_map<std::string, size_t, StringHasher> readNames;
        std::unordered_map<std::string, size_t, StringHasher> readNamesPos;

        size_t classID{0};
        size_t nextTarget = chunkSizes.front().chunkSize;

        // the next eq class has classSize reads, whose names are given by the
          // next classSize lines in the order file.
          for (size_t i = 0; i < numToProcess; ++i) {
              if (i >= nextTarget) { classID++; nextTarget += chunkSizes[classID].chunkSize; }
                std::vector<std::string> vec3 ;
                std::string line ;
                std::getline(orderFile,line);
                //orderFile >> line ;
                char *token = std::strtok((char *)line.c_str()," ");
                while(token != NULL){
                    vec3.push_back(token);
                    token = std::strtok(NULL," ");
                }
                //std::cout << "\nvec size "<<line <<"\n";
                int pos =  std::atoi(vec3[2].c_str());

                //readPos[std::string(vec3[0])] = std::atoi(vec3[2].c_str());

                std::string rn = std::string(vec3[0]);
              readNames[rn] = classID;
              readNamesPos[rn] = pos;
          }

          fp1 = gzopen(read1.getValue().c_str(), "r"); // STEP 2: open the file handler
          seq1 = kseq_init(fp1); // STEP 3: initialize seq

          fp2 = gzopen(read2.getValue().c_str(), "r"); // STEP 2: open the file handler
          seq2 = kseq_init(fp2); // STEP 3: initialize seq

          std::chrono::time_point<std::chrono::system_clock> start, end;
          start = std::chrono::system_clock::now();

          uint64_t readID{0};
          while ((l1 = kseq_read(seq1)) >= 0 and (l2 = kseq_read(seq2)) >= 0) { // STEP 4: read sequence
              if (readID % 10000 == 1) {
                  std::cout << "\r\rparsed " << readID << " reads";
              }
              auto it = readNames.find(seq1->name.s);
              //std::cout << "r: " << seq1->name.s << '\n';
              if (it != readNames.end()) {
                  auto& eq = chunkSizes[it->second];
                  int matepos = 0;
                  int pos = readNamesPos[seq1->name.s];
                  //std::cout << pos <<"\n";
                  eq.readSeqs.push_back({std::make_pair(seq1->seq.s,seq2->seq.s),pos});
                  //eq.readSeqsRight.push_back(std::make_pair(seq2->seq.s,matepos));
              } else if (readNames.find(seq2->name.s) != readNames.end()) {
                  auto& eq = chunkSizes[it->second];
                  int matepos = 0;
                  int pos = readNamesPos[seq1->name.s];

                  eq.readSeqs.push_back({std::make_pair(seq1->seq.s,seq2->seq.s),pos});
                  //eq.readSeqsLeft.push_back(std::make_pair(seq1->seq.s,pos));
                  //eq.readSeqsRight.push_back(std::make_pair(seq2->seq.s,matepos));
              }
              ++readID;
          }
          gzclose(fp1);
          gzclose(fp2);

          end = std::chrono::system_clock::now();
          std::chrono::duration<double> elapsed_seconds = end-start;
          std::cout << "tool " << elapsed_seconds.count() << ", seconds\n";
          std::cerr << "gathered " << readNames.size() << " reads ... encoding\n";

          for (auto& eq : chunkSizes) {
              // If the number of reads is too small, just encode them "simply"
              if (eq.chunkSize <= 10) {
                  if (eq.chunkSize != eq.readSeqs.size()) {
                      std::cout << "classSize = " << eq.chunkSize << ", but only found " << eq.readSeqs.size() << " reads!\n";
                  }
                  for (size_t i = 0; i < eq.chunkSize; ++i) {
                      outFileLeft << eq.readSeqs[i].seqlr.first << '\n';
                      outFileRight << eq.readSeqs[i].seqlr.second << '\n';
                  }
              } else {
                  encodeAsShift(eq.readSeqs, outFileLeft, outFileRight);
              }
              ++eqID;
          }
          std::cerr << "done\n";
        }
        

        std::cout << "\n equivalence classes gathered\n" ;


        std::string lFile_l = outDir + std::string("/r1.l.fq");
        std::string lFile_r = outDir + std::string("/r2.l.fq");
        std::ofstream ofs_unmapped_l(lFile_l, std::ofstream::out);
        std::ofstream ofs_unmapped_r(lFile_r, std::ofstream::out);

        std::ifstream unmappedInput(unmappedFile.getValue());
        std::unordered_set<std::string> unmappedNames;
        std::string rn;
        while (unmappedInput >> rn) { unmappedNames.insert(rn); }
        unmappedInput.close();

        fp1 = gzopen(read1.getValue().c_str(), "r"); // STEP 2: open the file handler
        seq1 = kseq_init(fp1); // STEP 3: initialize seq

        fp2 = gzopen(read2.getValue().c_str(), "r"); // STEP 2: open the file handler
        seq2 = kseq_init(fp2); // STEP 3: initialize seq

        std::vector<std::string> readSeqsLeft;
        std::vector<std::string> readSeqsRight;
        uint64_t readID{0};
        while ((l1 = kseq_read(seq1)) >= 0 and (l2 = kseq_read(seq2)) >= 0) { // STEP 4: read sequence
        bool found{false};
        auto it = unmappedNames.find(seq1->name.s);
        if (it != unmappedNames.end()) {
            found = true;
        } else if (unmappedNames.find(seq2->name.s) != unmappedNames.end()) {
            found = true;
        }
        if (found) {
            ofs_unmapped_l << "@" << seq1->name.s << '\n';
            ofs_unmapped_l << seq1->seq.s << '\n';
            ofs_unmapped_l << "+" << '\n';
            ofs_unmapped_l << seq1->qual.s << '\n';

            ofs_unmapped_r << "@" << seq2->name.s << '\n';
            ofs_unmapped_r << seq2->seq.s << '\n';
            ofs_unmapped_r << "+" << '\n';
            ofs_unmapped_r << seq2->qual.s << '\n';
        }
        }

        gzclose(fp1);
        gzclose(fp2);
        return 0;
    }catch (TCLAP::ArgException& e) {
        std::cerr << "Exception [" << e.error() << "] when parsing argument " << e.argId() << "\n";
      return 1;
    }

	//printf("return value: %d\n", l);
	return 0;
}
//How to run
// ./kseq_test <original seq file> <custom seq of read headers> <shuffled file>
