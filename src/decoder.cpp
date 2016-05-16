#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include "tclap/CmdLine.h"
#include <boost/filesystem.hpp>
#include <string>

void split(std::string& line, std::vector<std::string>& vec){
    char *token = std::strtok((char *)line.c_str()," ");
    while(token != NULL){
        vec.push_back(token);
        token = std::strtok(NULL,"\t");
    }
}

int main(int argc, char *argv[])
{
    std::cout << "Quark decoder module:\n";

    std::string versionString = "0.0.1";
    TCLAP::CmdLine cmd(
            "Quark seriator",
            ' ',
            versionString);
    cmd.getProgramName() = "quarkSeriate";

    TCLAP::ValueArg<std::string> read1("1", "leftMates", "The location of the left paired-end reads", false, "", "path");
    TCLAP::ValueArg<std::string> read2("2", "rightMates", "The location of the right paired-end reads", false, "", "path");

    cmd.add(read1);
    cmd.add(read2);

     try {

        cmd.parse(argc, argv);
        std::string left(read1.getValue());
        std::string right(read2.getValue());

        // create out dir from left
        boost::filesystem::path p(left);
        boost::filesystem::path outdir = p.parent_path();


        std::ifstream readL(left);
        std::ifstream readR(right);

        // create output file
        std::string outFilenameLeft = outdir.string() + "/r1.seq";
        std::string outFilenameRight = outdir.string() + "/r2.seq";
        std::ofstream outReadL(outFilenameLeft.c_str());
        std::ofstream outReadR(outFilenameRight.c_str());

        std::string ref ;

        for(std::string line; std::getline(readL,line);){
            if(line.empty()) continue ;
            std::vector<std::string> vec ;
            split(line,vec);
            if(vec.size() == 1){
                if(vec[0][0] != 'M'){
                    outReadL<<vec[0]<<"\n";
                    ref = vec[0];
                }else{
                    outReadL<<ref<<"\n";
                }
            }else{
                int shifts = std::atoi(vec[0].substr(1).c_str());
                int matches = std::atoi(vec[1].substr(1).c_str());
                std::string tmp = ref.substr(shifts+matches)+vec[2];
                ref = tmp;
                std::cout << ref << "\n";
                outReadL<<ref<<"\n";

            }
        }
    }catch (TCLAP::ArgException& e) {
        std::cerr << "Exception [" << e.error() << "] when parsing argument " << e.argId() << "\n";
      return 1;
    }
    return 0;
}
