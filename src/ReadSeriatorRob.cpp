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
#include <set>
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
//#include <Kmer.hpp>
#include <boost/functional/hash.hpp>
//#include <hash.hpp>
#include <cigargen.h>
#include "kseq.h"
#include "edlib.h"
#include "tclap/CmdLine.h"

typedef std::unordered_map<int,std::set<std::string> > EqReadID ;
typedef std::unordered_map<std::string, int> RevReadIDEq ;
typedef std::pair<std::string,std::string> matepair ;
typedef std::unordered_map<int,std::vector<matepair> > EqSeq ;


KSEQ_INIT(gzFile, gzread)

void refineKmers2(std::vector<std::string>& lefts, 
                  std::unordered_map<std::string, std::vector<int>>& final_kmers_left, 
                  std::vector<int>& unmappedIDs) {
  // maps kmers to read ids
  int k = 31;
  std::cout << "we are in new refineKmer2\n";
  std::unordered_set<std::string> invalidKmers;
  std::vector<int> coveredReads(lefts.size(), 0);

  int cutoff = 20; // ** arbitrary **
  int rid = 0;//
  //create a freq vector for each document
  //First do it for a read
  // populate the kmer set
  // that is get all the terms
  for(auto &val : lefts){
      //std::cout << val.length() << "\n";
      for(int i = 0; i < val.length()-k+1;++i){
          std::string mer = (val.substr(i,k));
          auto it = final_kmers_left.find(mer);
          if (it != final_kmers_left.end()) {
              if (it->second.size() > cutoff) {
                  invalidKmers.insert(mer);
                  continue;
              }
              // if the last read is the same as the current read, then skip it
              if (it->second.size() == 0) { std::cerr << "WTF!\n"; }
              if (it->second.back() != rid) {
                  it->second.push_back(rid);
                  coveredReads[rid]++;
              }
          } else { // otherwise
              final_kmers_left[mer].push_back(rid);
          }
      }
      rid++;
  }

  for (auto& m : invalidKmers) {
    std::vector<int>& readIds = final_kmers_left[m];
    for (int r : readIds) {
      coveredReads[r]--;
      if (coveredReads[r] == 0) { unmappedIDs.push_back(r); }
    }
    final_kmers_left.erase(m);
  }
}

#include <type_traits>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/undirected_dfs.hpp>

template <typename GraphT>
class MyVisitor : public boost::default_dfs_visitor
{
 public:
    MyVisitor(std::ofstream* left, std::ofstream* right, std::vector<std::string>* leftSeqs, std::vector<std::string>* rightSeqs) :
        _left(left), _right(right), _leftSeqs(leftSeqs), _rightSeqs(rightSeqs) {}

    void discover_vertex(typename boost::graph_traits<GraphT>::vertex_descriptor& v, const GraphT& g) const
    {
        (*_left) << (*_leftSeqs)[g[v].readID] << '\n';
        (*_right) << (*_rightSeqs)[g[v].readID] << '\n';
        return;
    }

private:
    std::vector<std::string>* _leftSeqs;
    std::vector<std::string>* _rightSeqs;
    std::ofstream* _left;
    std::ofstream* _right;
};

struct VertexInfo {
    int readID;
    bool discovered;
};

struct EdgeInfo {
    float weight;
};

template <typename GraphT>
void buildGraph(std::unordered_map<std::string, 
                std::vector<int>>& merMap,
                std::vector<int>& isolatedIDs,
                GraphT& G) {
    using VD = typename boost::graph_traits<GraphT>::vertex_descriptor;
    
    //Graph G;
    std::unordered_map<int, VD> vertMap;

    for (auto& kv : merMap) {
        auto& readIDs = kv.second;
        for (size_t i = 0; i < readIDs.size(); ++i) {
            auto vi = readIDs[i];
            auto itI = vertMap.find(vi);

            VD descI; 
            if (itI == vertMap.end()) {
                descI = boost::add_vertex(G);
                G[descI].readID = vi;
                G[descI].discovered = false;
                vertMap[vi] = descI;
            } else {
                descI = itI->second;
            }

            for (size_t j = i+1; j < readIDs.size(); ++j) {
                auto vj = readIDs[j];
                auto itJ = vertMap.find(vj);

                VD descJ; 
                if (itJ == vertMap.end()) {
                    descJ = boost::add_vertex(G);
                    G[descJ].readID = vj;
                    G[descJ].discovered = false;
                    vertMap[vj] = descJ;
                } else {
                    descJ = itJ->second;
                }
               
                auto edgeDesc = boost::edge(descI, descJ, G);
                if (edgeDesc.second) { // edge already exists
                    G[edgeDesc.first].weight += 1.0;
                } else {
                    auto e = boost::add_edge(descI, descJ, G);
                    G[e.first].weight = 1.0;
                }
            }
        }
    }
    // Add the disconnected nodes
   
    for (auto x : isolatedIDs) {
        auto vd = boost::add_vertex(G);
        G[vd].readID = x;
        G[vd].discovered = false;
    }
}

#include <boost/utility.hpp>

template <typename GraphT>
void visitViaDFS(GraphT& G, 
                 std::vector<std::string>& readSeqsLeft, 
                 std::vector<std::string>& readSeqsRight, 
                 std::ofstream& outFileLeft, 
                 std::ofstream& outFileRight) {

    using VD = typename boost::graph_traits<GraphT>::vertex_descriptor;
    using EdgeDesc = typename boost::graph_traits<GraphT>::edge_descriptor;

    using VertexIt = typename boost::graph_traits<GraphT>::vertex_iterator;
    using EdgeIt = typename boost::graph_traits<GraphT>::out_edge_iterator;

    std::unordered_set<VD> unvisited;
    // Attempt to remove all the vertices. Wrong!
    VertexIt vi, vi_end;
    for (boost::tie(vi, vi_end) = boost::vertices(G); vi != vi_end; ++vi) {
        unvisited.insert(*vi);
    }
    
    size_t visited{0};
    while (!unvisited.empty()) {
        // start from a random vertex
        auto vert = *(unvisited.begin());
        std::queue<VD> vertQueue;
        vertQueue.push(vert);
        while (!vertQueue.empty()) {
            auto v = vertQueue.front();
            vertQueue.pop();
            if (!G[v].discovered) {
                ++visited;
                G[v].discovered  = true;
                unvisited.erase(v);

                outFileLeft << readSeqsLeft[G[v].readID] << '\n';
                outFileRight<< readSeqsRight[G[v].readID] << '\n';

                EdgeIt ei, ei_end;
                std::vector<EdgeDesc> edges;
                for (boost::tie(ei, ei_end) = boost::out_edges(v, G); ei != ei_end; ++ei) {
                    edges.push_back(*ei);
                }
                std::sort(edges.begin(), edges.end(), 
                          [&G](const EdgeDesc& e1, const EdgeDesc& e2) -> bool {
                              return G[e1].weight > G[e2].weight;
                          });

                for (auto& e : edges) {
                    auto v1 = boost::source(e, G);
                    auto v2 = boost::target(e, G);
                    vertQueue.push((v1 == v) ? v2 : v1);
                }
            }
        }
    }
    

}

void encodeAsPath(std::vector<std::string>& readSeqsLeft, 
                  std::vector<std::string>& readSeqsRight, 
                  std::ofstream& outFileLeft, 
                  std::ofstream& outFileRight) {

    /** **/
    using EdgeWeightProp = boost::property<boost::edge_weight_t, int>;
    using Graph = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, VertexInfo, EdgeInfo>;
    using Edge = boost::graph_traits<Graph>::edge_descriptor;
    using VD = boost::graph_traits<Graph>::vertex_descriptor;
    
    Graph G;
    std::vector<int> isolatedLeft;
    std::vector<int> isolatedRight;
    std::unordered_map<std::string, std::vector<int>> merMapLeft;
    std::unordered_map<std::string, std::vector<int>> merMapRight;

    refineKmers2(readSeqsLeft, merMapLeft, isolatedLeft);
    buildGraph(merMapLeft, isolatedLeft, G);
    
    visitViaDFS(G, readSeqsLeft, readSeqsRight, outFileLeft, outFileRight);
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
      }
    */
    try {

      cmd.parse(argc, argv);
      std::cout << "\nGathering eqiovalence classes \n";
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
      size_t classSize{0};
      size_t eqID{0};
      while (sizeFile >> classSize) {
          // the next eq class has classSize reads, whose names are given by the
          // next classSize lines in the order file.
          std::unordered_set<std::string> readNames;
          for (size_t i = 0; i < classSize; ++i) {
              std::string rn;
              orderFile >> rn;
              readNames.insert(rn);
          }

          fp1 = gzopen(read1.getValue().c_str(), "r"); // STEP 2: open the file handler
          seq1 = kseq_init(fp1); // STEP 3: initialize seq
      
          fp2 = gzopen(read2.getValue().c_str(), "r"); // STEP 2: open the file handler
          seq2 = kseq_init(fp2); // STEP 3: initialize seq

          std::vector<std::string> readSeqsLeft;
          std::vector<std::string> readSeqsRight;
          uint64_t readID{0};
          while ((l1 = kseq_read(seq1)) >= 0 and (l2 = kseq_read(seq2)) >= 0) { // STEP 4: read sequence
              auto it = readNames.find(seq1->name.s);
              //std::cout << "r: " << seq1->name.s << '\n';
              if (it != readNames.end()) {
                  readSeqsLeft.push_back(seq1->seq.s); 
                  readSeqsRight.push_back(seq2->seq.s); 
              } else if (readNames.find(seq2->name.s) != readNames.end()) {
                  readSeqsLeft.push_back(seq1->seq.s); 
                  readSeqsRight.push_back(seq2->seq.s); 
              }
              ++readID;
          } 
          
          gzclose(fp1);
          gzclose(fp2);

          // If the number of reads is too small, just encode them "simply"
          if (classSize <= 10) {
              if (classSize != readSeqsLeft.size()) {
                  std::cout << "classSize = " << classSize << ", but only found " << readSeqsLeft.size() << " reads!\n";
              }
              for (size_t i = 0; i < classSize; ++i) {
                  outFileLeft << readSeqsLeft[i] << '\n';
                  outFileRight << readSeqsRight[i] << '\n';
              }
          } else {
              encodeAsPath(readSeqsLeft, readSeqsRight, outFileLeft, outFileRight);
          }

          std::cout << "equivalence class " << eqID << " had " << classSize << " reads\n";
          ++eqID;
      }
      //buildEqHeaders(eqClassOrderFile, eqClassSizeFile, req);
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
    } catch (TCLAP::ArgException& e) {
        std::cerr << "Exception [" << e.error() << "] when parsing argument " << e.argId() << "\n";
      return 1;
    }
}
//How to run
// ./kseq_test <original seq file> <custom seq of read headers> <shuffled file>
