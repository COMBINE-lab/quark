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
#include <stack>
//#include <Kmer.hpp>
#include <ctime>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/functional/hash.hpp>
//#include <Eigen/Dense>
//#include <hash.hpp>
#include "cigargen.h"
#include "kseq.h"
#include "edlib.h"

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

void buildEqHeaders(auto& orderFile,auto& chunkFile, RevReadIDEq &revidreq){
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

typedef struct Rline{
    int a;
    int b;
    std::string strval;
    Rline(int a,int b, std::string val):a(a),b(b),strval(val){

    }
};

typedef std::vector<Rline> Robj;


int RandomSample(Vertices& vs){
        std::srand(std::time(0));
            int steps = std::rand() % vs.size();
                auto it = vs.begin();
                for(int i =0;i < steps;i++){
                            ++it;

                }
                    int elem = *it;
                        //vs.erase(it);
                        //return elem;
                        //
}

void createObj(auto& order,auto& seq, Robj &obj){
    //std::cout << "\n Start Node "<< start_node<< "\n" ;
    //
    //
    //
    //
    std::cout << "\n In createObj \n" ;
    std::set<int> seenNodes;
    for(auto p : order){
        for(int i = 0; i < p.size()-1;++i){
            auto u = p[i];
            auto v = p[i+1];
            /* edlib module
           auto query =(const unsigned char *)seq[u].c_str();
            auto target =(const unsigned char *)seq[v].c_str();
            int alphabetLength = 4;
            int queryLength = seq[u].size();
            int targetLength = seq[v].size();
            int score, numLocations, alignmentLength;
            int* startLocations;
            int* endLocations;
            unsigned char* alignment;
            edlibCalcEditDistance(query, queryLength, target, targetLength,
                                          alphabetLength, -1, EDLIB_MODE_NW, true, true,
                                                                &score, &endLocations, &startLocations, &numLocations,
                                                                                      &alignment, &alignmentLength);
            char* cigar;
            edlibAlignmentToCigar(alignment, alignmentLength, EDLIB_CIGAR_EXTENDED, &cigar);
            */
                    std::string reference = seq[u];
                    std::string query = seq[v];
                    std::string cigar;
                    std::string alignment = "";
                    uint32_t alignment_l = 0;
                    int32_t edit_distance = GenerateCigar((char *) query.c_str(),
                            query.size(), (char *) reference.c_str(),
                            reference.size(), &cigar, &alignment_l, &alignment);
            std::cout <<reference.size()<<" "<<query.size()<<" " <<  cigar <<std::endl;
        }
    }
    /*
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
    }*/
}

void return_tf (auto &val, auto& kmerlist, auto &v){
    int k = 31;
    for(int i = 0; i < val.length()-k+1;++i){
        std::string mer = (val.substr(i,k));
        v(kmerlist[mer]) += 1;
    }
    //v = v / v.sum();
}

void refineKmers2(auto& lefts, auto& final_kmers_left, auto& unmappedIDs) {
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

  //print rids
  /*
  for(auto & u: final_kmers_left){
      for(auto r : u.second){
          std::cout << r << "\t";
      }
      std::cout << "\n";
  }*/

  for (auto& m : invalidKmers) {
    std::vector<int>& readIds = final_kmers_left[m];
    for (int r : readIds) {
      coveredReads[r]--;
    }
    final_kmers_left.erase(m);
  }

  for (size_t i = 0; i < coveredReads.size(); ++i) {
    if (coveredReads[i] == 0) { unmappedIDs.push_back(i); }
  }
}

//void return idf()
void createGraph(auto& kmerhash, auto& G){
    std::cout << "\n Creating Graph \n";
    for(auto& u: kmerhash){
        std::vector<int> vind = u.second;
        std::sort(vind.begin(),vind.end());
        //create pairs of nodes
        if(vind.size() == 1)
            G.v.insert(vind[0]);
        for(size_t i = 0;i < vind.size()-1;i++){
            for(size_t j = i+1; j < vind.size();j++){
                if(vind[i] != vind[j]){
                    auto r = std::make_pair(vind[i],vind[j]);
                    G.e[r] += 1;
                    G.v.insert(vind[i]);
                    G.v.insert(vind[j]);
                    G.adj[vind[i]].push_back(vind[j]);
                    G.adj[vind[j]].push_back(vind[i]);
                }
            }
        }
    }
    std::cout<< "\n Graph creation complete \n" ;
}

void DFS(auto&g, auto& path, auto& lone){
    //Take the highest node;
    std::cout << "\n A Graph with v: "<<g.v.size()<<" e: "<<g.adj.size()<<"\n";
    auto vcopy = g.v ;
    //take the maximum edge
    //
    int maxval = 0;
    Edge maxe ;
   for(auto u:g.e){
       if(u.second > maxval){
           maxval = u.second;
           maxe = u.first;
       }
   }

   auto u = maxe.first ;
   auto v = maxe.second ;
   //std::cout << "seed node "<< u << "\n";
   //start doing dfs unless you see all nodes ;
   bool start = true;
   std::vector<int> path_t;
   while(!vcopy.empty()){
        int seed ;
        if(start){
            seed = u;
            start = false;
        }else{
            seed = RandomSample(vcopy);
        }
        std::unordered_set<int> seenV ;
        std::stack<int> dfsstack;
        dfsstack.push(seed);
        while(!dfsstack.empty()){
            int v = dfsstack.top();
            dfsstack.pop();
            if(seenV.find(v) == seenV.end()){

                seenV.insert(v);
                //std::cout<< "\n visited "<<v<<"\n";
                //std::cout<< "\n size of vcopy:  "<<vcopy.size()<<"\n";
                //std::cout<< "\n size of path_t:  "<<path_t.size()<<"\n";
                vcopy.erase(v);
                if(g.adj.find(v) == g.adj.end()){
                    std::cout << "\n no edge with "<<v << "\n";
                    lone.push_back(v);
                    continue;
                }
                path_t.push_back(v);
                std::vector<int> neighbors = g.adj[v];
                //std::cout << " \n WTF \n";
                std::vector<Edge> incEdges ;
                for(auto& i: neighbors){
                    if(seenV.find(i) == seenV.end()){
                        Edge e = std::make_pair(v,i);
                        incEdges.push_back(e);
                    }
                }
                //sort edges with weight
                std::sort(incEdges.begin(),incEdges.end(),
                        [&g](const Edge& e1, const Edge& e2) -> bool {
                                            return g.e[e1] < g.e[e2];
                        });
                for(auto &ei : incEdges){
                    dfsstack.push(ei.second);
                }
            }

        }
        if(!path_t.empty()){
            std::cout << "\n path_t size "<<path_t.size()<<std::endl;
            path.push_back(path_t);
            path_t.clear();
        }

    }
}

void buildReadGraph_tfidf(auto &pvec, auto& lobj, auto& robj){
    //use boost matrix and vector
    std::cout << "\nIn tf-idf calculation block\n" ;
    int k = 31;
    std::vector<std::string> lefts;
    std::vector<std::string> rights;

    //std::cout << "\nSetting kmer size to  "<<k<<"\n";
    //Kmer::set_k(k);
    for(auto p : pvec){
        lefts.push_back(p.first);
        rights.push_back(p.second);
    }
    //using kmersfreq = std::unordered_map<std::string ,long long int>;
    using kmerHash = std::unordered_map<std::string,std::vector<int> >;
    //Clear the vector
    //vec.clear();

    Graph lg,rg,ltr,rtr;
    kmerHash leftkmers;
    kmerHash rightkmers;
    std::vector<int> ul ;
    std::vector<int> ur ;
    refineKmers2(lefts,leftkmers,ul);
    refineKmers2(rights,rightkmers,ur);

    std::cout << "\n Number of kmers = " << leftkmers.size() << '\n';
    createGraph(leftkmers,lg);
    std::vector<std::vector<int> > order_l;
    std::vector<std::vector<int> > order_r;
    std::vector<int> lone_l;
    std::vector<int> lone_r;
    DFS(lg,order_l,lone_l);
    std::cout << "\n path size is "<<order_l.size() << "\n";
    DFS(lg,order_r,lone_r);

    //for order use quark objects

    if(!order_l.empty()){
        createObj(order_l,lefts,lobj);
    }
    if(!order_r.empty()){
        createObj(order_r,rights,robj);
    }

}



int main(int argc, char *argv[])
{


	if (argc == 1) {
        fprintf(stderr, "Usage: %s <left.fq> <right.fq> <order> <chunksize> <unmapped file> <outdir>\n", argv[0]);
		return 1;
	}
    //read the fastq file because it will be needed for everything
    // TODO: unordered_map is not a smart idea, have to do better
    //

    std::cout << "\nGathering eqiovalence classes \n";
    RevReadIDEq req ;

    buildEqHeaders(argv[3],argv[4],req);
    std::cout << "\n equivalence classes gathered\n" ;



    std::vector<read_record > allreads ;
    std::string outdir = argv[6];

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
    EqSeq eseq; // equivalence class contain <class id>: vector of pair sequence

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
        populateEqSeq(argv[1],argv[2],eseq,req,unmappedHeaders,u_records);// Populate eseq

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

        std::ofstream ofs_unmapped_l;
        std::ofstream ofs_unmapped_r;

        std::string lFile_l = outdir + std::string("r1.l.fq");
        std::string lFile_r = outdir + std::string("r2.l.fq");
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

        /*
        for(auto &urid : u_records){
            ofs_unmapped_l<<urid.seq << "\n";
            ofs_unmapped_r<<urid.mateseq << "\n";

        }
        */
        ofs_unmapped_l.close();
        ofs_unmapped_r.close();
        ifsunmapped.close();
    }
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
                //buildReadGraph(pvec,lobj,robj);
                buildReadGraph_tfidf(pvec,lobj,robj) ;
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

	//printf("return value: %d\n", l);
	return 0;
}
//How to run
// ./kseq_test <original seq file> <custom seq of read headers> <shuffled file>
