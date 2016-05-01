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
#include <unordered_map>
#include <fstream>
#include <sys/time.h>

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
                if (numReads%10000 == 0){
                    std::cout << numReads+1 << " reads written \r" << std::flush ;
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
