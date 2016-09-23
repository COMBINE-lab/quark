#include <ctime>
#include <fstream>
#include <memory>

#include "cereal/archives/json.hpp"
#include "spdlog/details/format.h"

#include <boost/dynamic_bitset.hpp>
#include "pstream.h"

#include "GZipWriter.hpp"
#include "SailfishOpts.hpp"
#include "Transcript.hpp"
#include "LibraryFormat.hpp"

GZipWriter::GZipWriter(const boost::filesystem::path path, std::shared_ptr<spdlog::logger> logger) :
  path_(path), logger_(logger) {
}

GZipWriter::~GZipWriter() {
  if (bsStream_) {
    bsStream_->reset();
  }
}

/**
 * Creates a new gzipped file (path) and writes the contents
 * of the vector (vec) to the file in binary.
 */
template <typename T>
bool writeVectorToFile(boost::filesystem::path path,
                       const std::vector<T>& vec) {

    {
        bool binary = std::is_same<T, std::string>::value;
        auto flags = std::ios_base::out | std::ios_base::binary;

        boost::iostreams::filtering_ostream out;
        out.push(boost::iostreams::gzip_compressor(6));
        out.push(boost::iostreams::file_sink(path.string(), flags));

        size_t num = vec.size();
        size_t elemSize = sizeof(typename std::vector<T>::value_type);
        // We have to get rid of constness below, but this should be OK
        out.write(reinterpret_cast<char*>(const_cast<T*>(vec.data())),
                  num * elemSize);
        out.reset();
    }
    return true;
}

/**
 * Write the equivalence class information to file.
 * The header will contain the transcript / target ids in
 * a fixed order, then each equivalence class will consist
 * of a line / row.
 */

std::vector<std::string> split(const std::string &text, char sep) {
    std::vector<std::string> tokens;
    std::size_t start = 0, end = 0;
    while ((end = text.find(sep, start)) != std::string::npos) {
        tokens.push_back(text.substr(start, end - start));
        start = end + 1;
    }
    tokens.push_back(text.substr(start));
    return tokens;
}

void makeIslands(std::vector<std::pair<int32_t,int32_t>>& intervals, std::vector<int>& mapid, std::vector<int32_t>& relpos){

	//do a linear search on index
	typedef struct{
		std::pair<int32_t,int32_t> content ;
		int index ;
	}island;

	std::vector<island> fIslands ;
	int i = 0;
	for(auto c : intervals){
		island ci = {c,i} ;
		fIslands.push_back(ci) ;
		i++;
	}

	std::sort(fIslands.begin(),fIslands.end(),
				[](const island& p1, const island& p2) -> bool {
					return p1.content.first < p2.content.first ;
				});


	//std::vector<int> mapid ;
	std::vector<int> remapid ;
	std::vector<int> retmapid ;
	std::vector<int> retrelpos ;
	//std::vector<int32_t> relpos ;

	mapid.resize(intervals.size());
	retmapid.resize(intervals.size());
	remapid.resize(intervals.size());
	relpos.resize(intervals.size());
	retrelpos.resize(intervals.size());

	i = 0;
	for(auto ci :fIslands){
		intervals[i] = ci.content;
		mapid[ci.index] = i ;
		remapid[i] = i ;
		i++;
	}

	std::vector<std::pair<int32_t,int32_t>> correctedIslands;

	int currentindex;
	for(i = 0; i < intervals.size(); i++){
		if(i == 0){
			auto temp = intervals[i];
			correctedIslands.push_back(temp);
			remapid[i] = i;
			relpos[i] = 0 ;
		}else{
			auto temp = intervals[i];
			if(temp.first <= correctedIslands[correctedIslands.size()-1].second){
				// not new island relative position is offset
				relpos[i] = temp.first - correctedIslands[correctedIslands.size()-1].first ;
				//case 2
				if(temp.second > correctedIslands[correctedIslands.size()-1].second){
					correctedIslands[correctedIslands.size()-1].second = temp.second ;
				}
				remapid[i] = correctedIslands.size()-1;
			}else{
				correctedIslands.push_back(temp);
				remapid[i] = correctedIslands.size()-1;
				//new island
				relpos[i] = 0;
			}
		}
	}

	//std::cout << "du du";
	/*
	for(i = 0; i < intervals.size() ; i++){
		if(i == 0){
			auto temp = intervals[0];
			correctedIslands.push_back(temp);
			currentindex = 0;
			relpos[0] = 0;
		}else{
			//they are disjoint
			//time to create new island
			if(intervals[i].first > correctedIslands[currentindex].second){
				auto temp = intervals[i];
				correctedIslands.push_back(temp);
				currentindex++ ;
				remapid[i] = currentindex ;
				relpos[i] = 0;
			}else{
				//Now we have a tricky situation
				//where current island has an intersection with
				//the corrected current island so we extend our
				// old island to include this
				if(intervals[i].second > correctedIslands[currentindex].second){
					correctedIslands[currentindex].second = intervals[i].second;

				}
				//so we can rearrange whatever we have and see the relative positions
				remapid[i] = currentindex ;
				relpos[i] = intervals[i].first - correctedIslands[currentindex].first ;
				//otherwise we can swallow this small island

			}
		}

	}*/

	intervals = correctedIslands ;
	//get the remapped value
	for(i = 0; i < mapid.size(); i++){
		retmapid[i] = remapid[mapid[i]];
		retrelpos[i] = relpos[mapid[i]];
	}

	mapid = retmapid;
	relpos = retrelpos;
	//done
	//intervals = {{0,0}};
	// correctedIslands ;
}

std::vector<uint8_t> threeBitEncode(const std::string& str) {
	boost::dynamic_bitset<uint8_t> bitmap;
	for (auto c : str) {
		switch(c) {
			case 'A':
				bitmap.push_back(0);
				bitmap.push_back(0);
				bitmap.push_back(0);
				break;
			case 'C':
				bitmap.push_back(0);
				bitmap.push_back(0);
				bitmap.push_back(1);
				break;
			case 'G':
				bitmap.push_back(0);
				bitmap.push_back(1);
				bitmap.push_back(0);
				break;
			case 'T':
				bitmap.push_back(0);
				bitmap.push_back(1);
				bitmap.push_back(1);
				break;
			case 'N':
				bitmap.push_back(1);
				bitmap.push_back(0);
				bitmap.push_back(0);
				break;
			case '$':
				bitmap.push_back(1);
				bitmap.push_back(0);
				bitmap.push_back(1);
				break;
		}
	}
	std::vector<uint8_t> bytes;
	bitmap.push_back(1);
	bitmap.push_back(1);
	bitmap.push_back(0);
	boost::to_block_range(bitmap, std::back_inserter(bytes));
	return bytes;
}

std::vector<uint8_t> fourBitEncode(const std::string& str) {
	boost::dynamic_bitset<uint8_t> bitmap;
	for (auto c : str) {
		switch(c) {
			case 'A':
				bitmap.push_back(0);
				bitmap.push_back(0);
				bitmap.push_back(0);
				bitmap.push_back(1);
				break;
			case 'C':
				bitmap.push_back(0);
				bitmap.push_back(0);
				bitmap.push_back(1);
				bitmap.push_back(0);
				break;
			case 'G':
				bitmap.push_back(0);
				bitmap.push_back(0);
				bitmap.push_back(1);
				bitmap.push_back(1);
				break;
			case 'T':
				bitmap.push_back(0);
				bitmap.push_back(1);
				bitmap.push_back(0);
				bitmap.push_back(0);
				break;
			case 'N':
				bitmap.push_back(0);
				bitmap.push_back(1);
				bitmap.push_back(0);
				bitmap.push_back(1);
				break;
			case '0':
				bitmap.push_back(0);
				bitmap.push_back(1);
				bitmap.push_back(1);
				bitmap.push_back(0);
				break;
			case '1':
				bitmap.push_back(0);
				bitmap.push_back(1);
				bitmap.push_back(1);
				bitmap.push_back(1);
				break;
			case '2':
				bitmap.push_back(1);
				bitmap.push_back(0);
				bitmap.push_back(0);
				bitmap.push_back(0);
				break;
			case '3':
				bitmap.push_back(1);
				bitmap.push_back(0);
				bitmap.push_back(0);
				bitmap.push_back(1);
				break;
			case '4':
				bitmap.push_back(1);
				bitmap.push_back(0);
				bitmap.push_back(1);
				bitmap.push_back(0);
				break;
			case '5':
				bitmap.push_back(1);
				bitmap.push_back(0);
				bitmap.push_back(1);
				bitmap.push_back(1);
				break;
			case '6':
				bitmap.push_back(1);
				bitmap.push_back(1);
				bitmap.push_back(0);
				bitmap.push_back(0);
				break;
			case '7':
				bitmap.push_back(1);
				bitmap.push_back(1);
				bitmap.push_back(0);
				bitmap.push_back(1);
				break;
			case '8':
				bitmap.push_back(1);
				bitmap.push_back(1);
				bitmap.push_back(1);
				bitmap.push_back(0);
				break;
			case '9':
				bitmap.push_back(1);
				bitmap.push_back(1);
				bitmap.push_back(1);
				bitmap.push_back(1);
				break;
		}
	}
	std::vector<uint8_t> bytes;
	bitmap.push_back(0);
	bitmap.push_back(0);
	bitmap.push_back(0);
	bitmap.push_back(0);
	boost::to_block_range(bitmap, std::back_inserter(bytes));
	return bytes;
}

bool GZipWriter::writeEquivCounts(
    const SailfishOpts& opts,
    ReadExperiment& experiment,
	std::vector<std::vector<std::string>>& unmapped) {

  namespace bfs = boost::filesystem;

  bfs::path auxDir = path_ / opts.auxDir;
  bool auxSuccess = boost::filesystem::create_directories(auxDir);
  bfs::path eqFilePath = auxDir / "eq_classes.txt";

  //for paired end
  bfs::path quarkFilePath_1  = auxDir/ "reads_1.quark";
  bfs::path quarkFilePath_2  = auxDir/ "reads_2.quark";
  bfs::path offsetFilePath_1  = auxDir/ "offset_1.quark";
  bfs::path offsetFilePath_2  = auxDir/ "offset_2.quark";

  std::unique_ptr<redi::opstream> leftSeqPtr(nullptr);
  std::unique_ptr<redi::opstream> rightSeqPtr(nullptr);
  std::unique_ptr<redi::opstream> leftOffsetPtr(nullptr);
  std::unique_ptr<redi::opstream> rightOffsetPtr(nullptr);

  //for single end
  bfs::path quarkFilePath = auxDir/ "reads.quark";
  bfs::path offsetFilePath = auxDir/ "offset.quark";


  std::unique_ptr<redi::opstream> seqPtr(nullptr);
  std::unique_ptr<redi::opstream> offsetPtr(nullptr);

  bfs::path islandFile = auxDir/"islands.quark";
  fmt::MemoryWriter w;
  w.write("plzip -o {} -f -n {} -", islandFile.string(), opts.numThreads - 3);
  redi::opstream iFile(w.str());
  w.clear();

  if(experiment.readLibraries().front().format().type != ReadType::SINGLE_END){
	fmt::MemoryWriter w1;
	w1.write("plzip -o {} -f -n {} -", quarkFilePath_1.string(), opts.numThreads - 3);
	leftSeqPtr.reset(new redi::opstream(w1.str()));
	w1.clear();

	fmt::MemoryWriter w2;
	w2.write("plzip -o {} -f -n {} -", quarkFilePath_2.string(), opts.numThreads - 3);
	rightSeqPtr.reset(new redi::opstream(w2.str()));
	w2.clear();

	fmt::MemoryWriter w3;
	w3.write("plzip -o {} -f -n {} -", offsetFilePath_1.string(), opts.numThreads - 3);
	leftOffsetPtr.reset(new redi::opstream(w3.str()));
	w3.clear();

	fmt::MemoryWriter w4;
	w4.write("plzip -o {} -f -n {} -", offsetFilePath_2.string(), opts.numThreads - 3);
	rightOffsetPtr.reset(new redi::opstream(w4.str()));
	w4.clear();
  }else{


	fmt::MemoryWriter w2;
	w2.write("plzip -o {} -f -n {} -", quarkFilePath.string(), opts.numThreads - 3);
	seqPtr.reset(new redi::opstream(w2.str()));
	w2.clear();

	fmt::MemoryWriter w3;
	w3.write("plzip -o {} -f -n {} -", offsetFilePath.string(), opts.numThreads - 3);
	offsetPtr.reset(new redi::opstream(w3.str()));
	w3.clear();

  }

  //bfs::path quarkFilePath  = auxDir/ "read_names.quark";



  std::ofstream equivFile(eqFilePath.string());

  //**************previous way of writing*********
  //std::ofstream qFile(quarkFilePath.string());
  //std::ofstream iFile(islandFile.string());
  //**********************************************


  auto& transcripts = experiment.transcripts();
  std::vector<std::pair<const TranscriptGroup, TGValue>>& eqVec =
        experiment.equivalenceClassBuilder().eqVec();

  std::vector<std::pair<const TranscriptGroup, QStrings>>& qVec =
		experiment.quarkEqClassBuilder().eqVec();

  // Number of transcripts
  equivFile << transcripts.size() << '\n';

  // Number of equivalence classes
  equivFile << eqVec.size() << '\n';


  for (auto& t : transcripts) {
    equivFile << t.RefName << '\n';
  }

  for (auto& eq : eqVec) {
    uint64_t count = eq.second.count;
    // for each transcript in this class
    const TranscriptGroup& tgroup = eq.first;
    const std::vector<uint32_t>& txps = tgroup.txps;
    // group size
    equivFile << txps.size() << '\t';
    // each group member
    for (auto tid : txps) { equivFile << tid << '\t'; }
    // count for this class
    equivFile << count << '\n' ;
  }

  // Number of equivalence classes
  //iFile << eqVec.size() << '\n';

  for (auto& q : qVec){
	  const TranscriptGroup& tgroup = q.first;
	  const std::vector<uint32_t>& txps = tgroup.txps;
	  uint64_t count= q.second.count;
	  auto& qcodes = q.second.qcodes;
	  auto intervals = q.second.intervals;

	  std::vector<int> mapid;
	  std::vector<int32_t> relpos;
	  makeIslands(intervals,mapid,relpos);

	  struct quarkStruct{
	  		  std::string qcode;
	  		  int lIslandId;
	  		  int32_t lpos;
	  		  int rIslandId;
	  		  int32_t rpos;
	  	  };

	  int i = 0;
	  	  std::vector<quarkStruct> quarkStructVec;
	  	  for(auto qcode : qcodes) {
	  		  quarkStructVec.push_back({qcode,mapid[i],relpos[i],mapid[i+1],relpos[i+1]});
	  		  //qFile << qcode << "\t" << mapid[i] << ","<< relpos[i] << "," << mapid[i+1] << ","<< relpos[i+1]  << "\n";
	  		  i = i + 2;
	  	  }

	  	std::sort(quarkStructVec.begin(),quarkStructVec.end(),
	  				  [](const quarkStruct& q1, const quarkStruct& q2) -> bool {
	  			  	  	  if(q1.lIslandId != q2.lIslandId)
	  			  	  		  return q1.lIslandId < q2.lIslandId;
	  			  	  	  return q1.lpos < q2.lpos ;

	  		  });






	  if(experiment.readLibraries().front().format().type != ReadType::SINGLE_END){
		  leftSeqPtr->write(reinterpret_cast<char*>(&count), sizeof(count));
		  rightSeqPtr->write(reinterpret_cast<char*>(&count), sizeof(count));
		  int i=0;

		  int oldIsland = -1;


		  for(auto& qS : quarkStructVec) {
			  std::vector<std::string> seqs;
			  auto qcode = qS.qcode;
			  seqs = split(qcode,'|');
			  //std::string leftSeq = seqs[0];
			  //std::string rightSeq = seqs[1];
			  //encode left sequence
			  auto leftBytes = fourBitEncode(seqs[0]);
			  auto rightBytes = fourBitEncode(seqs[1]);

			  leftSeqPtr->write(reinterpret_cast<char*>(&leftBytes[0]), sizeof(leftBytes[0])* leftBytes.size());
			  rightSeqPtr->write(reinterpret_cast<char*>(&rightBytes[0]), sizeof(rightBytes[0])* rightBytes.size());

			  //qFile << qcode << "\t" << mapid[i] << ","<< relpos[i] << "," << mapid[i+1] << ","<< relpos[i+1]  << "\n";
			  uint8_t leftIsland = qS.lIslandId;
			  uint32_t leftPos = qS.lpos;
			  uint8_t rightIsland = qS.rIslandId;
			  uint32_t rightPos = qS.rpos;

			  leftOffsetPtr->write(reinterpret_cast<char*>(&leftIsland), sizeof(leftIsland));
			  leftOffsetPtr->write(reinterpret_cast<char*>(&leftPos), sizeof(leftPos));
			  rightOffsetPtr->write(reinterpret_cast<char*>(&rightIsland), sizeof(rightIsland));
			  rightOffsetPtr->write(reinterpret_cast<char*>(&rightPos), sizeof(rightPos));

			  i = i + 2;
		  }
	  }else{
		  seqPtr->write(reinterpret_cast<char*>(&count), sizeof(count));
		  int i = 0;
		  for(auto qS : quarkStructVec) {
			  auto qcode = qS.qcode;
			  auto bytes = fourBitEncode(qcode);

			  seqPtr->write(reinterpret_cast<char*>(&bytes[0]), sizeof(bytes[0])* bytes.size());
			  //qFile << qcode << "\t" << mapid[i] << ","<< relpos[i] << "\n";
			  uint8_t leftIsland = qS.lIslandId;
			  uint32_t leftPos = qS.lpos;

			  offsetPtr->write(reinterpret_cast<char*>(&leftIsland), sizeof(leftIsland));
			  offsetPtr->write(reinterpret_cast<char*>(&leftPos), sizeof(leftPos));

			  i = i + 2;
		  }
	  }




	  /*
	   * Old style encoding
	  if(experiment.readLibraries().front().format().type != ReadType::SINGLE_END){
		  for(auto qcode : qcodes) {
			  qFile << qcode << "\t" << mapid[i] << ","<< relpos[i] << "," << mapid[i+1] << ","<< relpos[i+1]  << "\n";
			  i = i + 2;
		  }
	  }else{
		  for(auto qcode : qcodes) {
			  qFile << qcode << "\t" << mapid[i] << ","<< relpos[i] << "\n";
			  i = i + 2;
		  }
	  }
	  */


	  /*
	  for(auto qcode : qcodes){
		  qFile << qcode << "\n";

	  }*/

	  //iFile << intervals.size()  << "\n";
	  uint32_t intervalSize = intervals.size();
	  //iFile.write(reinterpret_cast<char*>(&intervalSize), sizeof(intervalSize));

	  //iFile << intervals.size() << "," << transcripts[txps[0]].id << ","<< transcripts[txps[0]].RefLength << "\n";
	  for(auto interval : intervals) {
		  const char *txpSeq = transcripts[txps[0]].Sequence();
		  //block for $ symbols

		  /*
		  if(interval.first == interval.second){
			  std::string txpStr("$");
			  auto bytes = threeBitEncode(txpStr);
			  iFile.write(reinterpret_cast<char*>(&bytes[0]), sizeof(bytes[0])* bytes.size());
		  }else{
			  std::string txpStr(txpSeq+interval.first,txpSeq+interval.second+1);
			  auto bytes = threeBitEncode(txpStr);
			  iFile.write(reinterpret_cast<char*>(&bytes[0]), sizeof(bytes[0])* bytes.size());
		  }*/


		  boost::dynamic_bitset<uint8_t> bitmap;
		  for(int ind=0 ; ind <= interval.second ; ind++){
			switch(txpSeq[ind]) {
				case 'A':
					bitmap.push_back(0);
					bitmap.push_back(0);
					bitmap.push_back(0);
					break;
				case 'C':
					bitmap.push_back(0);
					bitmap.push_back(0);
					bitmap.push_back(1);
					break;
				case 'G':
					bitmap.push_back(0);
					bitmap.push_back(1);
					bitmap.push_back(0);
					break;
				case 'T':
					bitmap.push_back(0);
					bitmap.push_back(1);
					bitmap.push_back(1);
					break;
				case 'N':
					bitmap.push_back(1);
					bitmap.push_back(0);
					bitmap.push_back(0);
					break;
				case '$':
					bitmap.push_back(1);
					bitmap.push_back(0);
					bitmap.push_back(1);
					break;
			}
		  }


		std::vector<uint8_t> bytes;
		bitmap.push_back(1);
		bitmap.push_back(1);
		bitmap.push_back(0);
		boost::to_block_range(bitmap, std::back_inserter(bytes));
		iFile.write(reinterpret_cast<char*>(&bytes[0]), sizeof(bytes[0])* bytes.size());
		  /*
		  for(int ind = interval.first; ind <= interval.second;ind++)
			  iFile << txpSeq[ind];
		  iFile << "\n";
		  */

		  //iFile << interval.first << "\t" << interval.second << "\n";
	  }
	  boost::dynamic_bitset<uint8_t> bitmapD;
	  std::vector<uint8_t> delimBytes;
	  bitmapD.push_back(1);
	  bitmapD.push_back(1);
	  bitmapD.push_back(1);
	  boost::to_block_range(bitmapD, std::back_inserter(delimBytes));
	  iFile.write(reinterpret_cast<char*>(&delimBytes[0]), sizeof(delimBytes[0])* delimBytes.size());



  }

  //iFile.close();

  //write unmapped sequences
  if(unmapped.size() > 0){
	  int uid = 0;

	  if(experiment.readLibraries().front().format().type != ReadType::SINGLE_END){
		  bfs::path unMappedFile_l = auxDir/"unmapped_1.fastq";
		  bfs::path unMappedFile_r = auxDir/"unmapped_2.fastq";

		  /*
		  bfs::path unMappedFile = auxDir/"unmapped.fa";
		  std::ofstream uFile(unMappedFile.string());
		  for(auto& seqvec : unmapped){
			  for(auto& seq : seqvec){
				  uFile << seq << "\n";
			  }
		  }
		  */

		  std::ofstream uFile_l(unMappedFile_l.string());
		  std::ofstream uFile_r(unMappedFile_r.string());
		  for(auto seqvec : unmapped)
			  for(auto seq : seqvec){
				  int il = 0;
				  int len = (seq.size()-1)/2;
				  uFile_l << "@"<<uid<< "\n";
				  for(il = 0; il < len ; il++){
					  uFile_l << seq[il];
				  }
				  uFile_l << "\n";
				  uFile_l << "+"<< "\n";
				  for(il = 0; il < len ; il++){
					  uFile_l << "I";
				  }
				  uFile_l << "\n";



				  uFile_r << "@"<<uid<< "\n";
				  for(il = len+1; il < seq.size() ; il++){
					  uFile_r << seq[il] ;
				  }
				  uFile_r << "\n";
				  uFile_r << "+"<< "\n";
				  for(il = len+1; il < seq.size() ; il++){
					  uFile_r << "I" ;
				  }
				  uFile_r << "\n";
				  uid++;

		  }
	  }else{
		  bfs::path unMappedFile_l = auxDir/"unmapped.fastq";
		  std::ofstream uFile_l(unMappedFile_l.string());
		  for(auto seqvec : unmapped){
			  for(auto seq : seqvec){
				  uFile_l << "@"<<uid<< "\n";
				  for(int il = 0; il < seq.size() ; il++){
					  uFile_l << seq[il];
				  }
				  uFile_l << "\n";
				  uFile_l << "+"<< "\n";
				  for(int il = 0; il < seq.size() ; il++){
					  uFile_l << "I";
				  }

				  uFile_l << "\n";
				  uid++;
			  }
		  }

	  }
  }

  equivFile.close();
  return true;
}

/**
 * Write the ``main'' metadata to file.  Currently this includes:
 *   -- Names of the target id's if bootstrapping / gibbs is performed
 *   -- The fragment length distribution
 *   -- The expected and observed bias values
 *   -- A json file with information about the run
 */
bool GZipWriter::writeMeta(
    const SailfishOpts& opts,
    const ReadExperiment& experiment,
    const std::string& tstring // the start time of the run
    ) {

  namespace bfs = boost::filesystem;

  bfs::path auxDir = path_ / opts.auxDir;
  bool auxSuccess = boost::filesystem::create_directories(auxDir);

  auto numBootstraps = opts.numBootstraps;
  auto numSamples = (numBootstraps > 0) ? numBootstraps : opts.numGibbsSamples;
  if (numSamples > 0) {
      bsPath_ = auxDir / "bootstrap";
      bool bsSuccess = boost::filesystem::create_directories(bsPath_);
      {

          boost::iostreams::filtering_ostream nameOut;
          nameOut.push(boost::iostreams::gzip_compressor(6));
          auto bsFilename = bsPath_ / "names.tsv.gz";
          nameOut.push(boost::iostreams::file_sink(bsFilename.string(), std::ios_base::out));

          auto& transcripts = experiment.transcripts();
          size_t numTxps = transcripts.size();
          if (numTxps == 0) { return false; }
          for (size_t tn = 0; tn < numTxps; ++tn) {
              auto& t  = transcripts[tn];
              nameOut << t.RefName;
              if (tn < numTxps - 1) {
                  nameOut << '\t';
              }
          }
          nameOut << '\n';
          nameOut.reset();
      }

  }

  bfs::path fldPath = auxDir / "fld.gz";
  auto* fld = experiment.fragLengthDist();
  auto fragLengthSamples = fld->realize();
  writeVectorToFile(fldPath, fragLengthSamples);

  bfs::path normBiasPath = auxDir / "expected_bias.gz";
  writeVectorToFile(normBiasPath, experiment.expectedSeqBias());

  bfs::path obsBiasPath = auxDir / "observed_bias.gz";
  const auto& bcounts = experiment.readBias().counts;
  std::vector<int32_t> observedBias(bcounts.size(), 0);
  std::copy(bcounts.begin(), bcounts.end(), observedBias.begin());
  writeVectorToFile(obsBiasPath, observedBias);

  bfs::path normGCPath = auxDir / "expected_gc.gz";
  writeVectorToFile(normGCPath, experiment.expectedGCBias());

  bfs::path obsGCPath = auxDir / "observed_gc.gz";
  const auto& gcCounts = experiment.observedGC();
  std::vector<int32_t> observedGC(gcCounts.size(), 0);
  std::copy(gcCounts.begin(), gcCounts.end(), observedGC.begin());
  writeVectorToFile(obsGCPath, observedGC);

  bfs::path info = auxDir / "meta_info.json";

  {
      std::ofstream os(info.string());
      cereal::JSONOutputArchive oa(os);

      std::string sampType = "none";
      if (numBootstraps == 0 and numSamples > 0) {
          sampType = "gibbs";
      }
      if (numBootstraps > 0) {
          sampType = "bootstrap";
      }

      auto& transcripts = experiment.transcripts();
      oa(cereal::make_nvp("sf_version", std::string(sailfish::version)));
      oa(cereal::make_nvp("samp_type", sampType));
      oa(cereal::make_nvp("frag_dist_length", fld->maxValue()));
      oa(cereal::make_nvp("bias_correct", opts.biasCorrect));
      oa(cereal::make_nvp("num_bias_bins", bcounts.size()));
      oa(cereal::make_nvp("num_targets", transcripts.size()));
      oa(cereal::make_nvp("num_bootstraps", numBootstraps));
      oa(cereal::make_nvp("num_processed", experiment.numObservedFragments()));
      oa(cereal::make_nvp("num_mapped", experiment.numMappedFragments()));
      oa(cereal::make_nvp("percent_mapped", experiment.mappingRate() * 100.0));
      oa(cereal::make_nvp("call", std::string("quant")));
      oa(cereal::make_nvp("start_time", tstring));
  }
  return true;
}

bool GZipWriter::writeAbundances(
    const SailfishOpts& sopt,
    ReadExperiment& readExp) {

  using sailfish::math::LOG_0;
  using sailfish::math::LOG_1;
  namespace bfs = boost::filesystem;

  bool useScaledCounts = false;//(sopt.allowOrphans == false);
  bfs::path fname = path_ / "quant.sf";

  std::unique_ptr<std::FILE, int (*)(std::FILE *)> output(std::fopen(fname.c_str(), "w"), std::fclose);

  /* No comments for now
  if (headerComments.length() > 0) {
    fmt::print(output.get(), "{}", headerComments);
  }
  */

  // The header
  fmt::print(output.get(), "Name\tLength\tEffectiveLength\tTPM\tNumReads\n");

  double numMappedFrags = readExp.numMappedFragments();

  std::vector<Transcript>& transcripts_ = readExp.transcripts();
  for (auto& transcript : transcripts_) {
    transcript.projectedCounts = useScaledCounts ?
      (transcript.mass() * numMappedFrags) : transcript.estCount();
  }

  double tfracDenom{0.0};
  for (auto& transcript : transcripts_) {
    double refLength = sopt.noEffectiveLengthCorrection ?
      transcript.RefLength :
      transcript.EffectiveLength;
    tfracDenom += (transcript.projectedCounts / numMappedFrags) / refLength;
  }

  double million = 1000000.0;
  // Now posterior has the transcript fraction
  for (auto& transcript : transcripts_) {
    auto effLen = sopt.noEffectiveLengthCorrection ?
      transcript.RefLength :
      transcript.EffectiveLength;
    double count = transcript.projectedCounts;
    double npm = (transcript.projectedCounts / numMappedFrags);
    double tfrac = (npm / effLen) / tfracDenom;
    double tpm = tfrac * million;
    fmt::print(output.get(), "{}\t{}\t{}\t{}\t{}\n",
	transcript.RefName, transcript.RefLength, effLen,
	tpm, count);
  }

  return true;
}

template <typename T>
bool GZipWriter::writeBootstrap(const std::vector<T>& abund) {
#if defined __APPLE__
            spin_lock::scoped_lock sl(writeMutex_);
#else
            std::lock_guard<std::mutex> lock(writeMutex_);
#endif
	    if (!bsStream_) {
	      bsStream_.reset(new boost::iostreams::filtering_ostream);
	      bsStream_->push(boost::iostreams::gzip_compressor(6));
	      auto bsFilename = bsPath_ / "bootstraps.gz";
	      bsStream_->push(
                  boost::iostreams::file_sink(bsFilename.string(),
                                              std::ios_base::out | std::ios_base::binary));
	    }

	    boost::iostreams::filtering_ostream& ofile = *bsStream_;
	    size_t num = abund.size();
        size_t elSize = sizeof(typename std::vector<T>::value_type);
        ofile.write(reinterpret_cast<char*>(const_cast<T*>(abund.data())),
                    elSize * num);
        /*
        for (size_t tn = 0; tn < num; ++tn) {
            auto& a  = abund[tn];
            ofile << a;
            if (tn < num - 1) {
                ofile << '\t';
            }
        }
        ofile << '\n';
        */
        logger_->info("wrote {} bootstraps", numBootstrapsWritten_.load()+1);
        ++numBootstrapsWritten_;
        return true;
}

template
bool GZipWriter::writeBootstrap<double>(const std::vector<double>& abund);

template
bool GZipWriter::writeBootstrap<int>(const std::vector<int>& abund);

