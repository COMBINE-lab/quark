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
#include <cstdlib>
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


int numReads = 0;
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
     bool status;
     int matepos;
    bool matestatus ;
     //bool matestatus;
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

std::string revComp(std::string &s){
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


bool allN(std::string& s){
    int i = 0;
    while(i < s.size()){
        if(s[i] != 'N')
            return false;
        i++;
    }
    return true;
}

void encodeAsShiftMatch(std::vector<idpos> &l,
                   auto& ofs_mapped_l,
                   auto& ofs_mapped_r,
                   auto& ofs_perm){

    if (orderId%100 == 0){
                std::cout << orderId << " equivalence class processed \r" << std::flush ;
            }
            ++orderId;

    std::sort(l.begin(), l.end(),
                    [](const idpos &p1, const idpos &p2) -> bool {
                    return p1.pos < p2.pos ;
                    });


    using idposidx = std::pair<idpos *, int>;
    std::vector<idposidx> right;
    int i = 0;
    for(auto it = l.begin(); it != l.end() ; ++it, ++i){
        right.emplace_back(std::make_pair(&(*it),i));
    }

    //sort the right end but remember the order
    std::sort(right.begin(), right.end(),
            [](const idposidx &p1, idposidx &p2) -> bool {
                return p1.first->matepos < p2.first->matepos ;
            });


    {
        bool refFlag = true ;
        std::string prev;
        std::string curr;
        int prevpos ;
        int currpos ;
        bool currfwd ;
        bool prevfwd ;
        for(auto & r : right){
            idpos * p = r.first ;
            ofs_perm<<r.second<<"\n";
            if(refFlag){

                // for the right end
                prev = (p->matestatus) ? p->seqlr.second : revComp(p->seqlr.second);
                ofs_mapped_r<<prev<<p->status<< "\n";

                prevpos = p->matepos ;
                prevfwd = p->matestatus ;
                refFlag = false;

            }else{

                // for the right end
                //
                if(allN(p->seqlr.second)){
                    ofs_mapped_r<<"N"<<p->seqlr.second.size()<<"\n";
                    continue ;
                }
                curr = p->matestatus ? p->seqlr.second : revComp(p->seqlr.second);
                currpos = p->matepos ;
                currfwd = p->matestatus;

                {
                    int shift = currpos - prevpos;
                    //ofs_mapped_l<<"S"<<shift<<"\t";
                    //match = 0;
                    int counter=0;
                    int match = 0;
                    if(shift > 0 && shift < prev.size())
                        ofs_mapped_r<<"S"<<shift;
                    else
                        shift=0;
                    //ofs_mapped_l << "shit theory: "<<shift<< " prev.size(): "<< prev.size() <<"\t"<<curr  <<"\t";
                    while(shift+counter <= prev.size()-1 && counter < curr.size()){
                        if(prev[shift+counter] == curr[counter] ){
                            match++;

                        }else{
                            if(match < 3 ) {
                                ofs_mapped_r << (match == 0 ? (std::string(1,curr[counter])) : curr.substr(counter - match,match+1));
                            } else {
                                ofs_mapped_r << "M" << match;
                                ofs_mapped_r << curr[counter];
                            }
                            match = 0;
                        }
                        counter++;
                    }

                    if(match > 0)
                       ofs_mapped_r <<"M"<<match ;
                    //ofs_mapped_l<<"M"<<match<<"\t";
                    //if(match < curr.size())
                    //    ofs_mapped_l<<curr.substr(match);

                    if (prev == "GTTTCAATGCCAGCTTCCTGCTCTGCCCTTCAGATTTTGTTTTTAAGATCAACAAAGCCTGTAG"){
                        std::cout <<"S"<<shift<<" "<<"M"<<match<<"\n"<< prev.size() << "\n" << curr.size() << "\n";
                        exit (EXIT_FAILURE);
                    }
                    if (counter < curr.size())
                        ofs_mapped_r<<curr.substr(counter);
                    numReads++;

                    ofs_mapped_r<<currfwd<<"\n";
                }

                prev = curr;
                prevpos = currpos;
                prevfwd = currfwd;

            }
        }
    }

    {
        bool refFlag = true ;
        std::string prev;
        std::string curr;
        int prevpos ;
        int currpos ;
        bool currfwd ;
        bool prevfwd ;

        for(auto& p : l){
                    if(refFlag){
                        numReads++;

                        //for the left end
                        prev = (p.status) ? p.seqlr.first : revComp(p.seqlr.first);
                        ofs_mapped_l<<prev<<p.status<< "\n";
                        prevpos = p.pos ;
                        prevfwd = p.status ;


                        refFlag = false ;
                    }else{

                        // for the left end
                        if(allN(p.seqlr.first)){
                            ofs_mapped_l<<"N"<<p.seqlr.first.size()<<"\n";
                            continue ;
                        }
                        curr = p.status ? p.seqlr.first : revComp(p.seqlr.first);
                        currpos = p.pos ;
                        currfwd = p.status;

                        {
                            int shift = currpos - prevpos;
                            //ofs_mapped_l<<"S"<<shift<<"\t";
                            //match = 0;
                            int counter=0;
                            int match = 0;
                            if(shift > 0 && shift < prev.size())
                                ofs_mapped_l<<"S"<<shift;
                            else
                                shift=0;
                            //ofs_mapped_l << "shit theory: "<<shift<< " prev.size(): "<< prev.size() <<"\t"<<curr  <<"\t";
                            while(shift+counter <= prev.size()-1 && counter < curr.size()){
                                if(prev[shift+counter] == curr[counter] ){
                                    match++;

                                }else{
                                    if(match < 3 ) {
                                        ofs_mapped_l << (match == 0 ? (std::string(1,curr[counter])) : curr.substr(counter - match,match+1));
                                    } else {
                                        ofs_mapped_l << "M" << match;
                                        ofs_mapped_l << curr[counter];
                                    }
                                    match = 0;
                                }
                                counter++;
                            }

                            if(match > 0)
                               ofs_mapped_l <<"M"<<match ;
                            //ofs_mapped_l<<"M"<<match<<"\t";
                            //if(match < curr.size())
                            //    ofs_mapped_l<<curr.substr(match);

                            if (prev == "GTTTCAATGCCAGCTTCCTGCTCTGCCCTTCAGATTTTGTTTTTAAGATCAACAAAGCCTGTAG"){
                                std::cout <<"S"<<shift<<" "<<"M"<<match<<"\n"<< prev.size() << "\n" << curr.size() << "\n";
                                exit (EXIT_FAILURE);
                            }
                            if (counter < curr.size())
                                ofs_mapped_l<<curr.substr(counter);
                            numReads++;

                            ofs_mapped_l <<currfwd<<"\n";
                        }

                        prev = curr;
                        prevpos = currpos;
                        prevfwd = currfwd;
                    }

                    //ofs_mapped_r<<p.seqlr.second<<"\n";
        }

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
        std::string permFileName = outDir + "/perm.list";

        std::ofstream permFile(permFileName.c_str());
        std::ofstream outFileLeft(outFilenameLeft.c_str());
        std::ofstream outFileRight(outFilenameRight.c_str());

        // process the equivalence class
        // read the size of the next equivalence class
        const size_t chunkSize{12000000};
        size_t classSize{0};
        size_t eqID{0};
        bool done{false};

        //std::cout << "\n read this far\n" ;
        int cc = 0;
        //
        while(!done){
        size_t numToProcess{0};
          struct EqClassInfo {
              std::vector<idpos> readSeqs;
              //std::vector<idpos> readSeqsRight;
              size_t chunkSize;
          };

        std::vector<EqClassInfo> chunkSizes;
        while ((numToProcess < chunkSize) and (sizeFile >> classSize))  {
              chunkSizes.push_back({std::vector<idpos>(), classSize});
              chunkSizes.back().readSeqs.reserve(classSize);
              numToProcess += classSize;
          }
          std::cout << "\n I think so far " << numToProcess << " processed \n";
        if (numToProcess < chunkSize) { done = true; }

        std::unordered_map<std::string, size_t, StringHasher> readNames;
        std::unordered_map<std::string, size_t, StringHasher> readNamesPos;
        std::unordered_map<std::string, size_t, StringHasher> readNameStatus;
        std::unordered_map<std::string, size_t, StringHasher> readNamesMatePos;
        std::unordered_map<std::string, size_t, StringHasher> readNameMateStatus;

        size_t classID{0};
        size_t nextTarget = chunkSizes.front().chunkSize;

        // the next eq class has classSize reads, whose names are given by the
          // next classSize lines in the order file.
          for (size_t i = 0; i < numToProcess; ++i) {
              if (i >= nextTarget) { classID++; nextTarget += chunkSizes[classID].chunkSize; }
                std::vector<std::string> vec3 ;
                std::string line ;
                std::getline(orderFile,line);
                //std::cout << line << "\t";
                //orderFile >> line ;
                char *token = std::strtok((char *)line.c_str()," ");
                while(token != NULL){
                    vec3.push_back(token);
                    token = std::strtok(NULL," ");
                }
                //std::cout << "\nvec size "<<line <<"\n";
                //vec has the following way of parsing
                //<read name> <matepos> <pos> <txp> <mate flag> <flag>
                int matepos =  std::atoi(vec3[1].c_str());
                int pos =  std::atoi(vec3[2].c_str());

                bool status,matestatus ;

                status = (vec3[5] == std::string(1,'1')) ? true : false ;
                matestatus = (vec3[4] == std::string(1,'1')) ? true : false ;

                //std::cout <<"pos: "<<vec3[1]<<" "<<vec3[2]<< pos << "\n";
                //if(pos == 0){
                  //  exit (EXIT_FAILURE);
                //}
                //readPos[std::string(vec3[0])] = std::atoi(vec3[2].c_str());

                std::string rn = std::string(vec3[0]);
              readNames[rn] = classID;
              readNamesPos[rn] = pos;
              readNameStatus[rn] = status;
              readNamesMatePos[rn] = matepos;
              readNameMateStatus[rn] = matestatus;


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
                  //int matepos = 0;
                  int pos = readNamesPos[seq1->name.s];
                  bool status = readNameStatus[seq1->name.s];
                  int matepos = readNamesMatePos[seq1->name.s];
                  bool matestatus = readNameMateStatus[seq1->name.s];
                  //std::cout << pos <<"\n";
                  eq.readSeqs.push_back({std::make_pair(seq1->seq.s,seq2->seq.s),pos,status,matepos,matestatus});
                  //eq.readSeqsRight.push_back(std::make_pair(seq2->seq.s,matepos));
              } else if (readNames.find(seq2->name.s) != readNames.end()) {
                  auto& eq = chunkSizes[it->second];
                  //int matepos = 0;
                  int pos = readNamesPos[seq1->name.s];

                  bool status = readNameStatus[seq1->name.s];
                  int matepos = readNamesMatePos[seq1->name.s];
                  bool matestatus = readNameMateStatus[seq1->name.s];
                  eq.readSeqs.push_back({std::make_pair(seq1->seq.s,seq2->seq.s),pos,status,matepos,matestatus});
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
              if (eq.chunkSize <= 2) {
                  if (eq.chunkSize != eq.readSeqs.size()) {
                      std::cout << "classSize = " << eq.chunkSize << ", but only found " << eq.readSeqs.size() << " reads!\n";
                  }
                  for (size_t i = 0; i < eq.chunkSize; ++i) {
                      numReads++;
                      outFileLeft << eq.readSeqs[i].seqlr.first << '\n';
                      outFileRight << eq.readSeqs[i].seqlr.second << '\n';
                  }
              } else {
                  //encodeAsShift(eq.readSeqs, outFileLeft, outFileRight);
                  encodeAsShiftMatch(eq.readSeqs, outFileLeft, outFileRight, permFile);
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

        std::cout << "Number of Reasds "<< numReads<<"\n";

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
