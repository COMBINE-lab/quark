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
#include <fstream>
#include <sys/time.h>
#include <Kmer.hpp>
#include <boost/functional/hash.hpp>
#include <hash.hpp>
#include <cigargen.h>
#include "kseq.h"

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

namespace std{
    template<typename T> struct myhash;
    template <> struct myhash<Kmer>
    {
        uint64_t operator()(const Kmer& o) const{
            return o.hash();
        }
    };
}


typedef std::pair<int,int> Edge;
typedef std::vector<int> Path ;
typedef std::unordered_map< Edge, double, boost::hash< Edge > > EdgeList ;
typedef std::unordered_map< int, std::list<int> > AdjList;
typedef std::set<int> Vertices;
typedef struct Graph{
    EdgeList e;
    Vertices v;
    AdjList adj;
};

// Adopted from https://github.com/sbebo/graph/blob/master/mst.cpp
typedef struct UnionFind{
    std::vector<int> parent;
    std::vector<int> rank;
    UnionFind(int n):parent(n),rank(n)
    {
        for(int i = 0; i < n;i++)
            make_set(n);
    }
    void make_set(int i){
        parent[i] = i;
        rank[i] = 0;
    }
    int find_set(int i){
        if(i != parent[i])
            parent[i] = find_set(parent[i]);
    }
    void join(int i,int j){
        i = find_set(i);
        j = find_set(j);
        if(i == j)
            return;
        if(rank[i]>rank[j]){
             parent[j] = i;
        }else{
            parent[i] = j;
            if(rank[i] == rank[j]){
                rank[j] += 1;
            }
        }
    }
};


void createGraph(auto& svec, Graph &lg){

    std::cout << "\n in graph "<<svec.size() << "\n";
    for(int i=0;i<svec.size()-1;++i){
        for(int j=i;j<svec.size();++j){
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

void kruskal(Graph &g){

}

int maxDegree(AdjList &adj){
    std::cout <<"\n" <<adj.size() << "\n";
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
    for(int i = 0; i < n-1; i++){
        Edge currE = eorder.at(i);
        int a = currE.first ;
        int b = currE.second ;
        if(forest.find_set(a) != forest.find_set(b)){
            tr.e[currE] = g.e[currE];
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
    std::cout << "\n "<< start_node<< "\n" ;
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
        int curr_start = L.front();
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
    std::cout << "\nIn build graph \n";
    int k = 31;
    std::vector<std::string> lefts;
    std::vector<std::string> rights;

    std::cout << "\nSetting kmer size to  "<<k<<"\n";
    Kmer::set_k(k);
    for(auto p : pvec){
        lefts.push_back(p.first);
        rights.push_back(p.second);
    }
    std::unordered_map<Kmer, std::vector<int>, std::myhash<Kmer> > leftkmers;
    std::unordered_map<Kmer, std::vector<int>, std::myhash<Kmer> > rightkmers;
    int rid = 0;
    for(auto &val : lefts){
        std::cout << val.length() << "\n";
        for(int i = 0; i < val.length()-k+1;++i){
            Kmer o((const char *)val.substr(i,k).c_str());
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
        std::cout << val.length() << "\n";
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

    for(auto &kmere : leftkmers){
        auto svec = kmere.second;
        std::sort(svec.begin(),svec.end());
        createGraph(svec,lg);
    }

    //Done with making edge list and adj
    createMST(lg,ltr);
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

    //At this point we have all the reads and the make a read GRaph and do mst
    {
        std::cout << "\n Start writing quark files \n";
        std::vector<Robj> finalLeftObjs,finalRightObjs;
        for(auto eqClass : eseq){
            int id = eqClass.first;
            auto pvec = eqClass.second ;
            if(pvec.size() > 4){
                Robj lobj,robj;
                buildReadGraph(pvec,lobj,robj);
                finalLeftObjs.push_back(lobj);
                finalRightObjs.push_back(robj);
                eseq.erase(id);
            }
        }

    //For now write them in pure text in .quark files
    //

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