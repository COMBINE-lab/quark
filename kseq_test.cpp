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
#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <unordered_map>
#include <fstream>
#include <sys/time.h>

#include "kseq.h"



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

int main(int argc, char *argv[])
{


	if (argc == 1) {
        fprintf(stderr, "Usage: %s <left.fq> <right.fq> <order> <chunksize> <unmapped file> <outdir>\n", argv[0]);
		return 1;
	}
    //read the fastq file because it will be needed for everything
    // TODO: unordered_map is not a smart idea, have to do better
    std::vector<read_record > allreads ;
    std::string outdir = argv[6];

    std::cout << "\n Start reading read from fastq "<<argv[1] << " and " <<argv[2] << "\n";
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


    //get done with unmapped file
    //You just have to write the read seq
    //read the unmapped headers
    std::cout << "\n Start writing unmapped reads \n";
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

        //Write them
        std::ofstream ofs_unmapped_l;
        std::ofstream ofs_unmapped_r;
        std::string lFile_l = outdir + std::string("r1.l");
        std::string lFile_r = outdir + std::string("r2.l");
        std::cout << "\n Start writing unmapped to "<<lFile_l << " "<<lFile_r << "\n";
        ofs_unmapped_l.open(lFile_l, std::ofstream::out);
        ofs_unmapped_r.open(lFile_r, std::ofstream::out);
        for(auto &urid : u_records){
            ofs_unmapped_l<<urid.seq << "\n";
            ofs_unmapped_r<<urid.mateseq << "\n";

        }
        ofs_unmapped_l.close();
        ofs_unmapped_r.close();
        ifsunmapped.close();
    }
    // Done with unmapped reads hence .l files



    //create vector of all read sequence
    //std::vector<kseq_t *> allreads ;
    //std::vector<std::string > readids ;
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



    //Open two output streams
    //for .h file
    {
        std::ofstream ofs_mapped_l;
        std::ofstream ofs_mapped_r;
        std::string File_l = outdir + std::string("r1.h");
        std::string File_r = outdir + std::string("r2.h");
        ofs_mapped_l.open(File_l, std::ofstream::out);
        ofs_mapped_r.open(File_r, std::ofstream::out);
        std::ofstream ofsl;

        std::cout << "\n Start writing to shuffled fastq "<<argv[3] << "\n";

        numReads = 0;
        for(auto &readid: allreads){
            ofs_mapped_l<<readid.seq << "\n";
            ofs_mapped_r<<readid.mateseq << "\n";
            if (numReads%10000 == 0){
                std::cout << numReads+1 << " reads written \r" << std::flush ;
            }
            numReads++;
        }
        ofs_mapped_l.close();
        ofs_mapped_r.close();
    }

	//printf("return value: %d\n", l);
	kseq_destroy(seqr);
	kseq_destroy(seql);
    ifs.close() ;
	gzclose(fpleft);
	gzclose(fpright);
	return 0;
}
//How to run
// ./kseq_test <original seq file> <custom seq of read headers> <shuffled file>
