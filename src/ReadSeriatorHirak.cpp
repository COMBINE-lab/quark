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

typedef std::unordered_map<int,std::set<std::string> > EqReadID ;
using idpos = std::pair<std::string, int> ;
typedef std::unordered_map<std::string, int> RevReadIDEq ;
typedef std::unordered_map<std::string, int> RevReadIDPos ;

typedef std::pair<std::string,std::string> matepair ;
struct matepairPos{
    matepair seqlr;
     int pos;
};
typedef std::unordered_map<int,std::vector<matepairPos> > EqSeq ;


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



void buildEqHeaders(auto& orderFile,
                    auto& chunkFile,
                    RevReadIDEq &revidreq,
                    RevReadIDPos &revidpos){
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
    //create a vector of eqClassHeader;
    //
    std::ifstream ifsorder(orderFile);
    int orderId = 1; // can be taken as eqID
    if(ifsorder.is_open()){
        for(int chunk : chunkSizes){
            //read the headers
            for(int i = 0; i < chunk; ++i){
                std::vector<std::string> vec3 ;
                std::string line ;
                std::getline(ifsorder,line);
                char *token = std::strtok((char *)line.c_str()," ");
                while(token != NULL){
                    vec3.push_back(token);
                    token = std::strtok(NULL," ");
                }
                int pos =  std::atoi(vec3[2].c_str());
                //readPos[std::string(vec3[0])] = std::atoi(vec3[2].c_str());
                revidreq[std::string(vec3[0])] = orderId ;
                revidpos[std::string(vec3[0])] = pos;
            }

            if (orderId%100 == 0){
                std::cout << orderId << " equivalence class processed \r" << std::flush ;
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
                   RevReadIDPos &rpos,
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
        int pos = rpos[header];
        eseq[classid].push_back({p,pos});
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
    RevReadIDPos rpos;

    buildEqHeaders(argv[3],argv[4],req,rpos);
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
        populateEqSeq(argv[1],argv[2],eseq,req,rpos,unmappedHeaders,u_records);// Populate eseq


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

    // Done with unmapped reads hence .l files

    std::cout << "Number of equivalence Classes "<<eseq.size()<<"\n";
    //

    //At this point we have all the reads and the make a read GRaph and do mst
    int numReads = 0;


        std::ofstream ofs_mapped_l;
        std::ofstream ofs_mapped_r;
        std::string File_l = outdir + std::string("r1.quark");
        std::string File_r = outdir + std::string("r2.quark");
        ofs_mapped_l.open(File_l, std::ofstream::out);
        ofs_mapped_r.open(File_r, std::ofstream::out);

        for(auto& eqClass : eseq){

            if (numReads%1000000 == 0){
                std::cout << numReads / 1000000 << " million reads written \r" << std::flush ;
            }
            auto pvec = eqClass.second; // a vector of structs sort it
            // first
            std::sort(pvec.begin(), pvec.end(),
                    [&pvec](const matepairPos &p1, const matepairPos &p2) -> bool {
                        return p1.pos < p2.pos ;
                    });

            /*
            std::cout << "\n >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
            for(auto& m : pvec){
                std::cout << m.seqlr.first << "\t" << m.seqlr.second << "\t" << m.pos << "\n";
            }*/
            if(pvec.size() > 10){
                bool refFlag = true ;
                std::string prev;
                std::string curr;
                int prevpos ;
                int currpos ;
                for(auto& p : pvec){
                    numReads++;
                    if(refFlag){
                        ofs_mapped_l<<p.seqlr.first<< "\n";
                        prev = p.seqlr.first;
                        prevpos = p.pos ;
                        refFlag = false;
                    }else{
                        curr = p.seqlr.first;
                        currpos = p.pos;
                        int match = 0;
                        int shift = currpos - prevpos;
                        if(shift < curr.size()){
                            if(shift == 0){
                                ofs_mapped_l<<"M"<<curr.size()<<"\t";
                            }else{
                                ofs_mapped_l<<"S"<<shift;
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
            }else{
                for(auto& p : pvec){
                    numReads++;
                    ofs_mapped_l<<p.seqlr.first<<"\n";
                    ofs_mapped_r<<p.seqlr.second<<"\n";
                }
            }

        }


        ofs_unmapped_l.close();
        ofs_unmapped_r.close();
        ifsunmapped.close();

        ofs_mapped_l.close();
        ofs_mapped_r.close();
    }
    //Open two output streams
    //from eseq
    //for .h file
    /*
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
                ofs_mapped_l<<p.seqlr.first << "\n";
                ofs_mapped_r<<p.seqlr.second << "\n";
                if (numReads%1000000 == 0){
                    std::cout << numReads / 1000000 << " million reads written \r" << std::flush ;
                }
                numReads++;
            }
        }
        ofs_mapped_l.close();
        ofs_mapped_r.close();
    }*/

	//printf("return value: %d\n", l);
	return 0;
}
//How to run
// ./kseq_test <original seq file> <custom seq of read headers> <shuffled file>
