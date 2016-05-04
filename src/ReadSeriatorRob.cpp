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

void buildEqHeaders(const std::string& orderFile, const std::string& chunkFile, RevReadIDEq &revidreq){
    // read the chunk sizes in reverse
    std::vector<int> chunkSizes ;
    std::ifstream ifs(chunkFile);
    if(ifs.is_open()){
        for(std::string line; std::getline(ifs,line);){
            std::istringstream iss(line);
            int n;
            iss >> n;
            chunkSizes.push_back(n);
        }
    }
    // Reverse the vector chunkSizes
    //std::reverse(chunkSizes.begin(),chunkSizes.end());
    ifs.close();
    std::ifstream ifsorder(orderFile);
    int orderId = 1;
    if(ifsorder.is_open()){
        for(int chunk : chunkSizes){
            for(int i=0;i < chunk; ++i){
                std::string line;
                std::getline(ifsorder,line);
                revidreq[line] = orderId;
            }
            if (orderId%100 == 0){
                std::cout << orderId << " equivalence processed \r" << std::flush ;
            }
            ++orderId;
        }
    }
    std::cout <<"\n" <<orderId << " total equivalence classes "<<"\n";
    ifs.close();

}

void populateEqSeq(auto& fqFile,
                   auto& mateFile,
                   EqSeq &eseq,
                   RevReadIDEq &req,
                   auto& unmappedHeaders,
                   auto& u_records){

    std::cout << "\n Start reading read from fastq "<<fqFile << " and " <<mateFile << "\n";
	gzFile fpleft, fpright;
	kseq_t *seql,*seqr;
	int ll,lr;
	fpleft = gzopen(fqFile, "r");
	fpright = gzopen(mateFile, "r");
    seql = kseq_init(fpleft);
    seqr = kseq_init(fpright);

    int numReads = 0;
	while ((ll = kseq_read(seql) >= 0) && (lr = kseq_read(seqr) >= 0)) {
        std::string header = seql->name.s;
        if(unmappedHeaders.find(header) != unmappedHeaders.end()){
            std::string name(header);
            std::string sequl(seql->seq.s);
            std::string quali(seql->qual.s);
            std::string sequr(seqr->seq.s);
            read_record v(name,sequl,sequr,quali);
            u_records.push_back(v);
            continue ;
        }
        std::string sequl = seql->seq.s;
        std::string sequr = seqr->seq.s;
        matepair p = std::make_pair(sequl,sequr);
        int classid = req[header];
        if(eseq.find(classid) != eseq.end()){
            eseq[classid].push_back(p);
        }
        else{
             eseq[classid] = {p};
        }
        if (numReads%1000000 == 0){
            std::cout << numReads/1000000 << " million reads processed \r" << std::flush ;
        }
        numReads++;
	}
	kseq_destroy(seqr);
	kseq_destroy(seql);
	gzclose(fpleft);
	gzclose(fpright);
}

/*
namespace std{
    template<typename T> struct myhash;
    template <> struct myhash<Kmer>
    {
        uint64_t operator()(const Kmer& o) const{
            return o.hash();
        }
    };
}
*/

typedef std::string Kmer;
typedef std::pair<int,int> Edge;
typedef std::vector<int> Path ;
typedef std::unordered_map< Edge, double, boost::hash< Edge > > EdgeList ;
typedef std::unordered_map< int, std::vector<int> > AdjList;
typedef std::unordered_set<int> Vertices;
typedef struct Graph{
    EdgeList e;
    Vertices v;
    AdjList adj;
};

// Adopted from https://github.com/sbebo/graph/blob/master/mst.cpp
struct UnionFind{
    std::vector<int> parent; //parent for each node
    std::vector<int> rank;  //rank for each node

    UnionFind(int n): parent(n), rank(n)
    {
        for(int i=0;i<n;i++)
            make_set(i);

        }

    void make_set(int i)
    {
        parent[i]=i;
        rank[i]=0;
    }


    int find_set(int i)
    {
        if(i!=parent[i])
            parent[i] = find_set(parent[i]); //path compression
        return parent[i];
    }

    void join(int i, int j)
    {
        i=find_set(i);
        j=find_set(j);
        if (i==j)
            return;
        if (rank[i]>rank[j])  //i longer. attach j to i
            parent[j]=i;
        else{
            parent[i]=j; //i not longer. attach i to j
            if(rank[i]==rank[j]) //rank grows if equal length
                rank[j]+=1;
        }

    }

};


void createGraph(auto& svec, Graph &lg){

    //std::cout << "\n in graph "<<svec.size() << "\n";
    for(int i=0;i<svec.size()-1;++i){
        for(int j=i+1;j<svec.size();++j){
            //std::cout << svec[i]<< " "<<svec[j]<<"\n";
            if(svec[i] < svec[j]){
                Edge epair = std::make_pair(svec[i],svec[j]);
                lg.v.insert(svec[i]);
                lg.v.insert(svec[j]);
                if(lg.e.find(epair) == lg.e.end()){
                    lg.e[epair] = 1;
                }else{
                    lg.e[epair] += 1;
                }
                if(lg.adj.find(svec[i]) == lg.adj.end() && lg.adj.find(svec[j]) == lg.adj.end()){
                    //seeing this edge for the first time
                    lg.adj[svec[i]] = {svec[j]};
                    lg.adj[svec[j]] = {svec[i]};
                }else{
                    // add
                    lg.adj[svec[i]].push_back(svec[j]) ;
                    lg.adj[svec[j]].push_back(svec[i]) ;
                }

            }
        }
    }
}

void kruskal(Graph &g){

}

int maxDegree(AdjList &adj){
    //std::cout <<"\n" <<adj.size() << "\n";
    int m = 0;
    int ind = 0;
    for(auto u : adj){
        if(u.second.size() > m){
            m = u.second.size();
            ind = u.first;
        }
    }
    return ind;
}

void createMST(Graph &g, Graph &tr){
    int n = g.e.size();
    //std::cout << "\n Forest size "<<n<<"\n";
    UnionFind forest(n);
    // sort the edges
    std::vector<Edge> eorder;
    for(auto u: g.e)
        eorder.push_back(u.first);
    std::sort(eorder.begin(),eorder.end(),
            [&g](const Edge& e1, const Edge& e2) -> bool {
                return g.e[e1] > g.e[e2];
            });
    int size = 0;
    std::unordered_set<int> seenV ;
    seenV.reserve(g.v.size());
    //std::cout << "\n eorder size "<<eorder.size()<<"\n";
    for(int i = 0; seenV.size() < g.v.size(); i++){
        Edge currE = eorder.at(i);
        int a = currE.first ;
        int b = currE.second ;
        //std::cout << "\n a "<<a<<" f(a) "<<forest.find_set(a)<<" b "<<b<<"\n";
        if(forest.find_set(a) != forest.find_set(b)){
            //std::cout << "\na "<<a<<"b "<<b<<"\n";
            tr.e[currE] = g.e[currE];
            forest.join(a,b);
            seenV.insert(a);
            seenV.insert(b);
            tr.v.insert(a);
            tr.v.insert(b);
            if(tr.adj.find(a) == tr.adj.end() && tr.adj.find(b) == tr.adj.end()){
                tr.adj[a] = {b};
                tr.adj[b] = {a};
            }else{
                tr.adj[a].push_back(b);
                tr.adj[b].push_back(a);
            }
        }
    }

}
typedef struct Rline{
    int a;
    int b;
    std::string strval;
    Rline(int a,int b, std::string val):a(a),b(b),strval(val){

    }
};

typedef std::vector<Rline> Robj;


void createObj(auto &tr,int start_node ,auto& seqs, Robj &obj){
    //std::cout << "\n Start Node "<< start_node<< "\n" ;
    std::set<int> seenNodes;
    std::queue<int> L;
    std::string refString = seqs.at(start_node);
    L.push(start_node);
    std::unordered_map<int, std::string > cigarMap ;
    //seenNodes.insert(start_node);
    //Robj obj;
    bool atfirst = true;
    while(L.size() > 0){
        //see the start node
        //std::cout << "\n L sz "<<L.size()<<"\n";
        int curr_start = L.front();
        L.pop();
        auto qvec = tr.adj[curr_start] ;
        seenNodes.insert(curr_start);
        for(auto v_i: qvec){
            if(seenNodes.find(v_i) == seenNodes.end()){
                std::string reference = seqs.at(curr_start);
                std::string query = seqs.at(v_i);
                L.push(v_i);
                if(atfirst){
                    Rline line(curr_start,v_i,refString);
                    obj.push_back(line);
                    atfirst = false ;
                }else{
                    std::string cigar;
                    std::string alignment = "";
                    uint32_t alignment_l = 0;
                    int32_t edit_distance = GenerateCigar((char *) query.c_str(),
                            query.size(), (char *) reference.c_str(),
                            reference.size(), &cigar, &alignment_l, &alignment);
                    Rline line(curr_start,v_i,cigar);
                    obj.push_back(line);
                }
            }
        }
    }
}

void buildReadGraph(auto &pvec, auto& lobj, auto& robj){
    //std::cout << "\nIn build graph \n";
    int k = 31;
    std::vector<std::string> lefts;
    std::vector<std::string> rights;

    //std::cout << "\nSetting kmer size to  "<<k<<"\n";
    //Kmer::set_k(k);
    for(auto p : pvec){
        lefts.push_back(p.first);
        rights.push_back(p.second);
    }
    std::unordered_map<Kmer, std::vector<int>, boost::hash<Kmer> > leftkmers;
    std::unordered_map<Kmer, std::vector<int>, boost::hash<Kmer> > rightkmers;
    int rid = 0;
    for(auto &val : lefts){
        //std::cout << val.length() << "\n";
        for(int i = 0; i < val.length()-k+1;++i){
            Kmer o(val.substr(i,k));
            if(leftkmers.find(o) == leftkmers.end()){
                leftkmers[o] = {rid};
            }else{
                leftkmers[o].push_back(rid);
            }
        }
        rid++;
    }
    rid = 0;

    for(auto &val : rights){
        //std::cout << val.length() << "\n";
        for(int i = 0; i < val.length()-k+1;++i){
            Kmer o((const char *)val.substr(i,k).c_str());
            if(rightkmers.find(o) == rightkmers.end()){
                rightkmers[o] = {rid};
            }else{
                rightkmers[o].push_back(rid);
            }
        }
        rid++;
    }

    Graph lg,rg,ltr,rtr;

    std::cout << "\n Number of kmers = " << leftkmers.size() << '\n';
    for(auto &kmere : leftkmers){
        auto svec = kmere.second;
        std::sort(svec.begin(),svec.end());
        createGraph(svec,lg);
        //std::cout << "\n Graph created with V:"<<lg.v.size()<<" E:"<<lg.e.size()<<"\n";
    }

    //std::cout << "\n Graph created with V:"<<lg.v.size()<<" E:"<<lg.e.size()<<"\n";
    //Done with making edge list and adj
    //
    //print graph : testing
    for(auto ed: lg.e){
        auto e1 = ed.first;
        double val = ed.second;
        //std::cout << "\n "<<e1.first<<" "<<e1.second<<" "<<val<<"\n";
    }
    createMST(lg,ltr);
    //std::cout << "\n MST created with V:"<<ltr.v.size()<<" E:"<<ltr.e.size()<<" adj sz: "<<ltr.adj.size()<<"\n";
    int start_node = maxDegree(ltr.adj);
    createObj(ltr,start_node,lefts,lobj);


    for(auto &kmere : rightkmers){
        auto svec = kmere.second;
        std::sort(svec.begin(),svec.end());
        createGraph(svec,rg);
    }

    //Done with making edge list and adj
    createMST(rg,rtr);
    start_node = maxDegree(rtr.adj);
    createObj(rtr,start_node,rights,robj);


    //buildKmerHash
}


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

    
    //std::vector<read_record > allreads ;

      /*
	gzFile fpleft, fpright;
	kseq_t *seql,*seqr;
	int ll,lr;
	fpleft = gzopen(argv[1], "r");
	fpright = gzopen(argv[2], "r");
	seql = kseq_init(fpleft);
	seqr = kseq_init(fpright);

	int numReads = 0;
	while ((ll = kseq_read(seql) >= 0) && (lr = kseq_read(seqr) >= 0)) {
	std::string name(seql->name.s);
	std::string sequl(seql->seq.s);
	std::string quali(seql->qual.s);
	std::string sequr(seqr->seq.s);
	read_record v(name,sequl,sequr,quali);
	allreads.push_back(v);
	if (numReads%10000 == 0){
	std::cout << numReads+1 << " reads processed \r" << std::flush ;
	}
	numReads++;
	}
      */

      //get done with unmapped file
      //You just have to write the read seq
      //read the unmapped headers
      //EqSeq eseq; // equivalence class contain <class id>: vector of pair sequence
      /*
      {
	std::vector<read_record> u_records ;

	std::cout << "\n Start reading unmapped read headers \n";
	std::set<std::string> unmappedHeaders;
	std::ifstream ifsunmapped;
	ifsunmapped.open(argv[5],std::ifstream::in);
	for(std::string line; std::getline(ifsunmapped,line) ;){
	  unmappedHeaders.insert(line);
	}
	//output to std::vector<matepair>
	//populateEqSeq(argv[1],argv[2],eseq,req,unmappedHeaders,u_records);// Populate eseq
    */
	/*
	  {
	  std::set<std::string> unmappedHeaders;
	  std::ifstream ifsunmapped;
	  ifsunmapped.open(argv[5],std::ifstream::in);
	  for(std::string line; std::getline(ifsunmapped,line) ;){
	  unmappedHeaders.insert(line);
	  }
	  //filter out these sequences
	  std::vector<read_record> u_records ;
	  std::copy_if(allreads.begin(),allreads.end(),
	  std::back_inserter(u_records),
	  [&unmappedHeaders](const read_record& r) -> bool {
	  return unmappedHeaders.find(r.name) != unmappedHeaders.end(); } );

	  allreads.erase(std::remove_if(allreads.begin(),
	  allreads.end(),
	  [&unmappedHeaders](const read_record& r) -> bool {
	  return unmappedHeaders.find(r.name) != unmappedHeaders.end();
	  }),
	  allreads.end() );
	  }*/

	//Write them
      /*
	std::ofstream ofs_unmapped_l;
	std::ofstream ofs_unmapped_r;

	std::string lFile_l = outDir + std::string("r1.l.fq");
	std::string lFile_r = outDir + std::string("r2.l.fq");
	std::cout << "\n Start writing unmapped to "<<lFile_l << " "<<lFile_r << "\n";
	ofs_unmapped_l.open(lFile_l, std::ofstream::out);
	ofs_unmapped_r.open(lFile_r, std::ofstream::out);

	for(auto &urid : u_records){
	  //Written the left end
	  ofs_unmapped_l << "@" << urid.name << "\n";
	  ofs_unmapped_l << urid.seq << "\n";
	  ofs_unmapped_l<< urid.comment << "\n";
	  ofs_unmapped_l<< urid.quality << "\n";

	  //Written the right end
	  ofs_unmapped_r << "@" << urid.name << "\n";
	  ofs_unmapped_r << urid.mateseq << "\n";
	  ofs_unmapped_r << urid.comment << "\n";
	  ofs_unmapped_r << urid.quality << "\n";
	}
	std::cout << "\n Done writing unmapped reads \n";

	ofs_unmapped_l.close();
	ofs_unmapped_r.close();
	ifsunmapped.close();
      }
    */
      // Done with unmapped reads hence .l files



      //create vector of all read sequence
      //std::vector<kseq_t *> allreads ;
      //std::vector<std::string > readids ;
      /*
	std::unordered_map<std::string,int > nameToRank ;
	std::ifstream ifs;
	ifs.open(argv[3],std::ifstream::in) ;
	int rank = 1;
	std::cout << "Start reading read id orders from "<<argv[3] << "\n";


	for(std::string line; std::getline(ifs, line) ;){
	//if (line.empty() || (line[0] == "-")) continue ;
	nameToRank[line] = rank ;
	if (rank%10000 == 0){
	std::cout << rank+1 << " reads processed \r" << std::flush ;
	}
	rank++;
	}


	std::cout << "\n Start sorting in reverse\n" ;
	timespec startt,endt,result;
	clock_gettime(CLOCK_REALTIME, &startt);
	std::sort(allreads.begin(), allreads.end(), [&nameToRank](const read_record& r1,
	const read_record& r2) -> bool {
	//std::cout << r1.name << "\n";
	//std::cout << r2.name << "\n";
	//std::cout<< nameToRank[r1.name] << nameToRank[r2.name] << "\n" ;
	return nameToRank[r1.name] > nameToRank[r2.name]; }) ;

	clock_gettime(CLOCK_REALTIME, &endt);
	//inttimeval_subtract(&result,&startt,&endt);
	std::cout << "\n Sorting done, took "<< (endt.tv_sec - startt.tv_sec) << " sec to finish \n";
      */

      /*
      std::cout << "Number of equivalence Classes "<<eseq.size()<<"\n";
      //At this point we have all the reads and the make a read GRaph and do mst
      {
	std::cout << "\n Start writing quark files \n";
	std::vector<Robj> finalLeftObjs,finalRightObjs;
	int eqCount = 0;
	for(auto eqClass : eseq){
	  int id = eqClass.first;
	  auto pvec = eqClass.second ;
	  if(pvec.size() > 4){
	    eqCount++;
	    Robj lobj,robj;
	    buildReadGraph(pvec,lobj,robj);
	    finalLeftObjs.push_back(lobj);
	    finalRightObjs.push_back(robj);
	    //eseq.erase(id);

	    std::cout << "\n Done for eq: "<<eqCount << "\n";
	  }
	}

	//For now write them in pure text in .quark files
	//
	std::cout << "Number of equivalence Classes remain"<<eseq.size()<<"\n";

	std::cout << "\n Still writing quark files \n";
	std::ofstream ofs_mapped_l;
	std::ofstream ofs_mapped_r;
	std::string File_l = outdir + std::string("r1.quark");
	std::string File_r = outdir + std::string("r2.quark");
	ofs_mapped_l.open(File_l, std::ofstream::out);
	ofs_mapped_r.open(File_r, std::ofstream::out);
	int objlines = 0;
	for(auto obj: finalLeftObjs){
	  ofs_mapped_l<<obj.size()<<"\n";
	  for(Rline l:obj){
	    ofs_mapped_l<<l.a<<" "<<l.b<<" "<<l.strval<<"\t" ;
	  }
	  ofs_mapped_l<<"\n";
	}
	for(auto obj: finalRightObjs){
	  ofs_mapped_r<<obj.size()<<"\n";
	  for(Rline l:obj){
	    ofs_mapped_r<<l.a<<" "<<l.b<<" "<<l.strval<<"\t" ;
	  }
	  ofs_mapped_r<<"\n";
	}
      }
      //Open two output streams
      //from eseq
      //for .h file
      {
	std::ofstream ofs_mapped_l;
	std::ofstream ofs_mapped_r;
	std::string File_l = outdir + std::string("r1.h");
	std::string File_r = outdir + std::string("r2.h");
	ofs_mapped_l.open(File_l, std::ofstream::out);
	ofs_mapped_r.open(File_r, std::ofstream::out);
	std::ofstream ofsl;

	std::cout << "\n Start writing to shuffled seq to "<< File_l << " " << File_r << "\n";

	int numReads = 0;
	for(auto &eqClass : eseq){
	  auto pvec = eqClass.second;
	  for(auto& p : pvec){
	    ofs_mapped_l<<p.first << "\n";
	    ofs_mapped_r<<p.second << "\n";
	    if (numReads%1000000 == 0){
	      std::cout << numReads / 1000000 << " million reads written \r" << std::flush ;
	    }
	    numReads++;
	  }
	}
	ofs_mapped_l.close();
	ofs_mapped_r.close();
      }
      */
      //printf("return value: %d\n", l);
      return 0;
    } catch (TCLAP::ArgException& e) {
        std::cerr << "Exception [" << e.error() << "] when parsing argument " << e.argId() << "\n";
      return 1;
    }
}
//How to run
// ./kseq_test <original seq file> <custom seq of read headers> <shuffled file>
