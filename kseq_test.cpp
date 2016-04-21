#include <zlib.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <fstream>
#include <sys/time.h>

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

int main(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *seq;
	int l;

    //create vector of all read sequence
    //std::vector<kseq_t *> allreads ;
    std::vector<read_record > allreads ;
    //std::vector<std::string > readids ;
    std::unordered_map<std::string,int > nameToRank ;

	if (argc == 1) {
		fprintf(stderr, "Usage: %s <in.seq>\n", argv[0]);
		return 1;
	}
	fp = gzopen(argv[1], "r");
    std::ifstream ifs;
    ifs.open(argv[2],std::ifstream::in) ;
    int rank = 0;
    std::cout << "Start reading read ids from "<<argv[2] << "\n";


    for(std::string line; std::getline(ifs, line) ;){
        nameToRank[line] = rank ;
        if (rank%10000 == 0){
            std::cout << rank+1 << " reads processed \r" << std::flush ;
        }
        rank++;
    }


    std::cout << "\n Start reading read from fastq "<<argv[1] << "\n";
	seq = kseq_init(fp);
    int numReads = 0;
	while ((l = kseq_read(seq)) >= 0) {


		//printf("name: %s\n", seq->name.s);
        //if (seq->comment.l)
        //printf("comment: %s\n", seq->comment.s);
		//printf("seq: %s\n", seq->seq.s);
		//if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
        std::string name(seq->name.s);
        std::string sequ(seq->seq.s);
        std::string quali(seq->qual.s);
        read_record v(name,sequ,quali);
        allreads.push_back(v);
        if (numReads%10000 == 0){
            std::cout << numReads+1 << " reads processed \r" << std::flush ;
        }
        numReads++;
	}

    std::cout << "\n Start sorting \n" ;
    timespec startt,endt,result;
    clock_gettime(CLOCK_REALTIME, &startt);
    std::sort(allreads.begin(), allreads.end(), [&nameToRank](const read_record& r1,
                const read_record& r2) -> bool {
                //std::cout << r1.name << "\n";
                //std::cout << r2.name << "\n";
                //std::cout<< nameToRank[r1.name] << nameToRank[r2.name] << "\n" ;
                return nameToRank[r1.name] < nameToRank[r2.name]; }) ;

    clock_gettime(CLOCK_REALTIME, &endt);
    //inttimeval_subtract(&result,&startt,&endt);
    std::cout << "\n Sorting done, took "<< (startt.tv_sec - endt.tv_sec) << " ns to finish \n";

    std::ofstream ofs;
    ofs.open(argv[3], std::ofstream::out) ;

    std::cout << "\n Start writing to shuffled fastq "<<argv[3] << "\n";

    numReads = 0;
    for(auto &readid: allreads){
        ofs << "@" << readid.name << "\n" ;
        ofs << readid.seq<< "\n" ;
        ofs << readid.comment << "\n" ;
        ofs << readid.quality << "\n" ;
        if (numReads%10000 == 0){
            std::cout << numReads+1 << " reads written \r" << std::flush ;
        }
        numReads++;
    }

	//printf("return value: %d\n", l);
	kseq_destroy(seq);
    ifs.close() ;
    ofs.close() ;
	gzclose(fp);
	return 0;
}
