#include <ctime>
#include <fstream>
#include <chrono>
#include <thread>
#include <memory>
#include <functional>
#include <stdio.h>
#include <unistd.h>

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
				bitmap.push_back(1);
				bitmap.push_back(0);
				bitmap.push_back(0);
				bitmap.push_back(0);
				break;
			case 'C':
				bitmap.push_back(0);
				bitmap.push_back(1);
				bitmap.push_back(0);
				bitmap.push_back(0);
				break;
			case 'G':
				bitmap.push_back(1);
				bitmap.push_back(1);
				bitmap.push_back(0);
				bitmap.push_back(0);
				break;
			case 'T':
				bitmap.push_back(0);
				bitmap.push_back(0);
				bitmap.push_back(1);
				bitmap.push_back(0);
				break;
			case 'N':
				bitmap.push_back(1);
				bitmap.push_back(0);
				bitmap.push_back(1);
				bitmap.push_back(0);
				break;
			case '0':
				bitmap.push_back(0);
				bitmap.push_back(1);
				bitmap.push_back(1);
				bitmap.push_back(0);
				break;
			case '1':
				bitmap.push_back(1);
				bitmap.push_back(1);
				bitmap.push_back(1);
				bitmap.push_back(0);
				break;
			case '2':
				bitmap.push_back(0);
				bitmap.push_back(0);
				bitmap.push_back(0);
				bitmap.push_back(1);
				break;
			case '3':
				bitmap.push_back(1);
				bitmap.push_back(0);
				bitmap.push_back(0);
				bitmap.push_back(1);
				break;
			case '4':
				bitmap.push_back(0);
				bitmap.push_back(1);
				bitmap.push_back(0);
				bitmap.push_back(1);
				break;
			case '5':
				bitmap.push_back(1);
				bitmap.push_back(1);
				bitmap.push_back(0);
				bitmap.push_back(1);
				break;
			case '6':
				bitmap.push_back(0);
				bitmap.push_back(0);
				bitmap.push_back(1);
				bitmap.push_back(1);
				break;
			case '7':
				bitmap.push_back(1);
				bitmap.push_back(0);
				bitmap.push_back(1);
				bitmap.push_back(1);
				break;
			case '8':
				bitmap.push_back(0);
				bitmap.push_back(1);
				bitmap.push_back(1);
				bitmap.push_back(1);
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

	/*
	std::cout << str << "\n";
	for(auto b : bytes){
		std::cout << unsigned(b) << "\t";
	}
	exit(0);
	*/
	return bytes;
}

template <typename T, typename S, typename P>
void pushIslandBitmap(T &intervals,S &txpSeq ,P &islandPtr){

	  for(auto interval : intervals) {
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
		  auto intervalSize = interval.second - interval.first + 1;
		  bitmap.resize((intervalSize+1) * 3);
		  auto offset = 0;
		  for(int ind=interval.first ; ind <= interval.second ; ind++){
			//iFile<<txpSeq[ind];
			auto c = txpSeq[ind];
			switch(c) {
				case 'A':
					bitmap[offset] =   0;
					bitmap[offset+1] = 0;
					bitmap[offset+2] = 0;
					break;
				case 'C':
					bitmap[offset] =   0;
					bitmap[offset+1] = 0;
					bitmap[offset+2] = 1;
					break;
				case 'G':
					bitmap[offset] =   0;
					bitmap[offset+1] = 1;
					bitmap[offset+2] = 0;
					break;
				case 'T':
					bitmap[offset] =   0;
					bitmap[offset+1] = 1;
					bitmap[offset+2] = 1;
					break;
				case 'N':
					bitmap[offset] =   1;
					bitmap[offset+1] = 0;
					bitmap[offset+2] = 0;
					break;
				case '$':
					bitmap[offset] =   1;
					bitmap[offset+1] = 0;
					bitmap[offset+2] = 1;
					break;
			}
			offset += 3;
		  }


			std::vector<uint8_t> bytes;
			bitmap[offset]=1;
			bitmap[offset+1]=1;
			bitmap[offset+2]=0;
			boost::to_block_range(bitmap, std::back_inserter(bytes));
			islandPtr->write(reinterpret_cast<char*>(&bytes[0]), sizeof(bytes[0])* bytes.size());

			//numOfEquivClasses++;

		  //iFile << "\n";


		  //iFile << interval.first << "\t" << interval.second << "\n";
	  }
}

std::function<void(redi::opstream*)> closeStreamDeleter(std::string s)  {
	  std::function<void(redi::opstream*)> deleter =  [s](redi::opstream* ptr) -> void {
	  ptr->rdbuf()->peof();
	  //fmt::print(stderr, "\n start writing unmapped\n");
	  fmt::print(stderr,"\nPtr deleted {}\n", s);
	  fmt::print(stderr,"status = {}\n", ptr->rdbuf()->status());
	  std::chrono::seconds wt(2);
	  std::this_thread::sleep_for(wt);
	  fmt::print(stderr,"status = {}\n", ptr->rdbuf()->status());
	  ptr->close();
	  fmt::print(stderr,"status = {}\n", ptr->rdbuf()->status());
	  delete ptr;
	  };
	  return deleter;
}

bool GZipWriter::writeEncoding(
    const SailfishOpts& opts,
    ReadExperiment& experiment,
	std::vector<std::vector<std::string>>& unmapped,
	bool& qualityScore) {


  namespace bfs = boost::filesystem;

  using pstream_ptr = std::unique_ptr<redi::opstream, std::function<void(redi::opstream*)>>;

  bfs::path auxDir = path_ / opts.auxDir;
  bool auxSuccess = boost::filesystem::create_directories(auxDir);

  bfs::path islandFile = auxDir/"islands.quark";
  bfs::path islandTxtFile = auxDir/"islands.txt";
  std::ofstream iFile(islandTxtFile.string());

  /*
  pstream_ptr islandPtr(nullptr, closeStreamDeleter("islandPtr"));
  {
	  fmt::MemoryWriter w;
	  w.write("plzip -o {} -f -n {} -", islandFile.string(), 1);
	  //w.write("gzip -1 - > {}", islandFile.string()+".gz");
	  islandPtr.reset(new redi::opstream(w.str()));
	  w.clear();
  }
  */
  bfs::path eqFilePath = auxDir / "eq_classes.txt";

  //for debugging
  bfs::path txtFilePath = auxDir/ "reads.txt";


  //for paired end

  bfs::path quarkFilePath_1  = auxDir/ "reads_1.quark";
  bfs::path quarkFilePath_2  = auxDir/ "reads_2.quark";
  bfs::path qualityFilePath_1  = auxDir/ "quality_1.quark";
  bfs::path qualityFilePath_2  = auxDir/ "quality_2.quark";
  bfs::path offsetFilePath_1  = auxDir/ "offset_1.quark";
  bfs::path offsetFilePath_2  = auxDir/ "offset_2.quark";
  bfs::path chunkFilePath_1  = auxDir/ "chunk_1.quark";
  bfs::path chunkFilePath_2  = auxDir/ "chunk_2.quark";
  //bfs::path mince  = auxDir/ "chunk_2.quark";





  pstream_ptr leftSeqPtr(nullptr, closeStreamDeleter("leftSeqPtr"));
  pstream_ptr rightSeqPtr(nullptr, closeStreamDeleter("rightSeqPtr"));
  pstream_ptr leftOffsetPtr(nullptr, closeStreamDeleter("leftOffsetPtr"));
  pstream_ptr rightOffsetPtr(nullptr, closeStreamDeleter("rightOffsetPtr"));
  //pstream_ptr leftChunkPtr(nullptr, closeStreamDeleter("leftChunkPtr"));
  //pstream_ptr rightChunkPtr(nullptr, closeStreamDeleter("rightChunkPtr"));

  //for single end
  bfs::path quarkFilePath = auxDir/ "reads.quark";
  bfs::path qualityFilePath = auxDir/ "quality.quark";
  bfs::path offsetFilePath = auxDir/ "offset.quark";

  //debug
  //bfs::path normalSeq = auxDir/"mapped.reads";


  pstream_ptr seqPtr(nullptr, closeStreamDeleter("seqPtr"));
  pstream_ptr offsetPtr(nullptr, closeStreamDeleter("offsetPtr"));





  if(experiment.readLibraries().front().format().type != ReadType::SINGLE_END){
	fmt::MemoryWriter w1;
	w1.write("plzip -o {} -f -n {} -", quarkFilePath_1.string(), 1);
	leftSeqPtr.reset(new redi::opstream(w1.str()));
	w1.clear();

	fmt::MemoryWriter w2;
	w2.write("plzip -o {} -f -n {} -", quarkFilePath_2.string(), 1);
	rightSeqPtr.reset(new redi::opstream(w2.str()));
	w2.clear();

	fmt::MemoryWriter w3;
	w3.write("plzip -o {} -f -n {} -", offsetFilePath_1.string(), 1);
	leftOffsetPtr.reset(new redi::opstream(w3.str()));
	w3.clear();

	fmt::MemoryWriter w4;
	w4.write("plzip -o {} -f -n {} -", offsetFilePath_2.string(), 1);
	rightOffsetPtr.reset(new redi::opstream(w4.str()));
	w4.clear();

	/*
	w4.write("plzip -o {} -f -n {} -", chunkFilePath_1.string(), 1);
	leftChunkPtr.reset(new redi::opstream(w4.str()));
	w4.clear();


	w4.write("plzip -o {} -f -n {} -", chunkFilePath_2.string(), 1);
	rightChunkPtr.reset(new redi::opstream(w4.str()));
	w4.clear();
	*/


  }else{


	fmt::MemoryWriter w2;
	w2.write("plzip -o {} -f -n {} -", quarkFilePath.string(),1);
	seqPtr.reset(new redi::opstream(w2.str()));
	w2.clear();

	fmt::MemoryWriter w3;
	w3.write("plzip -o {} -f -n {} -", offsetFilePath.string(),1);
	offsetPtr.reset(new redi::opstream(w3.str()));
	w3.clear();

  }

  //bfs::path quarkFilePath  = auxDir/ "read_names.quark";



  std::ofstream equivFile(eqFilePath.string());
  std::ofstream quality_1;
  std::ofstream quality_2;
  std::ofstream quality_single;

  //write quality score
  if(qualityScore){
	  if(experiment.readLibraries().front().format().type != ReadType::SINGLE_END){
		  quality_1.open(qualityFilePath_1.string(), std::ofstream::out | std::ofstream::app);
		  quality_2.open(qualityFilePath_2.string(), std::ofstream::out | std::ofstream::app);

	  }else{
		  quality_single.open(qualityFilePath.string(), std::ofstream::out | std::ofstream::app);

	  }
  }

  //for debugging

  //**************previous way of writing*********
  std::ofstream qFile(txtFilePath.string());
  //std::ofstream seqFile(normalSeq.string());
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
  iFile << eqVec.size() << '\n';

  int numOfEquivClasses = 0;
  for (auto& q : qVec){
	  const TranscriptGroup& tgroup = q.first;
	  const std::vector<uint32_t>& txps = tgroup.txps;
	  uint32_t count= q.second.count;
	  auto& qcodes = q.second.qcodes;
	  auto& qualityscores = q.second.qualityscores ;
	  auto intervals = q.second.intervals;

	  std::vector<int> mapid;
	  std::vector<int32_t> relpos;
	  makeIslands(intervals,mapid,relpos);

	  struct quarkStruct{
	  		  std::string qcode;
	  		  std::string quality;
	  		  int lIslandId;
	  		  int32_t lpos;
	  		  int rIslandId;
	  		  int32_t rpos;
	  	  };

	  int i = 0;
	  int quality_index = 0;
	  	  std::vector<quarkStruct> quarkStructVec;
	  	  for(auto qcode : qcodes) {
	  		  if(!qualityScore){
	  			  quarkStructVec.push_back({qcode,"",mapid[i],relpos[i],mapid[i+1],relpos[i+1]});
	  		  }else{
	  			  quarkStructVec.push_back({qcode,qualityscores[quality_index],mapid[i],relpos[i],mapid[i+1],relpos[i+1]});
	  			  quality_index++;
	  		  }
	  		  //qFile << qcode << "\t" << mapid[i] << ","<< relpos[i] << "," << mapid[i+1] << ","<< relpos[i+1]  << "\n";
	  		  i = i + 2;
	  	  }

	  	std::sort(quarkStructVec.begin(),quarkStructVec.end(),
	  				  [](const quarkStruct& q1, const quarkStruct& q2) -> bool {
	  			  	  	  if(q1.lIslandId != q2.lIslandId)
	  			  	  		  return q1.lIslandId < q2.lIslandId;
	  			  	  	  return q1.lpos < q2.lpos ;

	  		  });





      qFile << count << "\n";
	  if(experiment.readLibraries().front().format().type != ReadType::SINGLE_END){
		  leftSeqPtr->write(reinterpret_cast<char*>(&count), sizeof(count));

		  //rightSeqPtr->write(reinterpret_cast<char*>(&count), sizeof(count));
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

			  uint8_t leftChunk = leftBytes.size();
			  uint8_t rightChunk = rightBytes.size();

			  leftSeqPtr->write(reinterpret_cast<char*>(&leftBytes[0]), sizeof(leftBytes[0])* leftBytes.size());
			  rightSeqPtr->write(reinterpret_cast<char*>(&rightBytes[0]), sizeof(rightBytes[0])* rightBytes.size());

			  //leftChunkPtr->write(reinterpret_cast<char*>(&leftChunk), sizeof(leftChunk));
			  //rightChunkPtr->write(reinterpret_cast<char*>(&rightChunk), sizeof(rightChunk));

			  qFile << qS.qcode << "\t" << qS.lIslandId << ","<< qS.lpos << "," << qS.rIslandId << ","<< qS.rpos  << "\n";

			  uint8_t leftIsland = qS.lIslandId;
			  uint32_t leftPos = qS.lpos;
			  uint8_t rightIsland = qS.rIslandId;
			  uint32_t rightPos = qS.rpos;

			  leftOffsetPtr->write(reinterpret_cast<char*>(&leftIsland), sizeof(leftIsland));
			  leftOffsetPtr->write(reinterpret_cast<char*>(&leftPos), sizeof(leftPos));
			  rightOffsetPtr->write(reinterpret_cast<char*>(&rightIsland), sizeof(rightIsland));
			  rightOffsetPtr->write(reinterpret_cast<char*>(&rightPos), sizeof(rightPos));

			  i = i + 2;

			  if(qualityScore){
				  std::vector<std::string> qualities;
				  auto qualityscores = qS.quality ;
				  qualities = split(qualityscores,'|');
				  quality_1 << qualities[0] << "\n";
				  quality_2 << qualities[1] << "\n";

			  }
		  }
	  }else{
		  seqPtr->write(reinterpret_cast<char*>(&count), sizeof(count));
		  int i = 0;
		  for(auto qS : quarkStructVec) {
			  auto qcode = qS.qcode;
			  auto bytes = fourBitEncode(qcode);

			  seqPtr->write(reinterpret_cast<char*>(&bytes[0]), sizeof(bytes[0])* bytes.size());

			  qFile << qcode << "\t" << qS.lIslandId << "," << qS.lpos << "\n" ;
			  //qFile << qcode << "\t" << mapid[i] << ","<< relpos[i] << "\n";
			  uint8_t leftIsland = qS.lIslandId;
			  uint32_t leftPos = qS.lpos;

			  offsetPtr->write(reinterpret_cast<char*>(&leftIsland), sizeof(leftIsland));
			  offsetPtr->write(reinterpret_cast<char*>(&leftPos), sizeof(leftPos));

			  i = i + 2;

			  if(qualityScore){
				  auto qualityscores = qS.quality ;
				  quality_single << qualityscores << "\n";

			  }
		  }
	  }




	  /*
	   * Old style encoding
	   */
	  /*
	  int j = 0;
	  qFile << count << "\n";
	  if(experiment.readLibraries().front().format().type != ReadType::SINGLE_END){
		  for(auto qcode : qcodes) {
			  qFile << qcode << "\t" << mapid[j] << ","<< relpos[j] << "," << mapid[j+1] << ","<< relpos[j+1]  << "\n";
			  j = j + 2;
		  }
	  }else{
		  for(auto qcode : qcodes) {
			  qFile << qcode << "\t" << mapid[j] << ","<< relpos[j] << "\n";
			  j = j + 2;
		  }
	  }
	  */
	 /* */


	  /*
	  for(auto qcode : qcodes){
		  qFile << qcode << "\n";

	  }*/


	  //************** Write islands *********************
	  //iFile << intervals.size()  << "\n";
	  //uint32_t intervalSize = intervals.size();
	  //iFile.write(reinterpret_cast<char*>(&intervalSize), sizeof(intervalSize));
	  //iFile << intervals.size() << "," << transcripts[txps[0]].id << ","<< transcripts[txps[0]].RefLength << "\n";
	  /*
	  if (numOfEquivClasses % 10000 == 0) {

		  fmt::print(stderr, "\n processed {} equivalence classes (cumulative)", numOfEquivClasses);

	  }

	  numOfEquivClasses++;
	  */


	  const char *txpSeq = transcripts[txps[0]].Sequence();
	  //write islands in bitmap fashion
	  //pushIslandBitmap(intervals,txpSeq,islandPtr);

	  iFile << intervals.size()  << "\n";
	  //iFile << intervals.size() << "," << transcripts[txps[0]].id << ","<< transcripts[txps[0]].RefLength << "\n";
	  for(auto interval:intervals){
		  for(int ind = interval.first; ind <= interval.second;ind++)
			  iFile << txpSeq[ind];
		  iFile <<"\n";
		  //iFile << interval.first << "\t" << interval.second <<"\n";
	  }


	  boost::dynamic_bitset<uint8_t> bitmapD;
	  std::vector<uint8_t> delimBytes;
	  bitmapD.push_back(1);
	  bitmapD.push_back(1);
	  bitmapD.push_back(1);
	  boost::to_block_range(bitmapD, std::back_inserter(delimBytes));
	  //islandPtr->write(reinterpret_cast<char*>(&delimBytes[0]), sizeof(delimBytes[0])* delimBytes.size());

	  //*****************************************************




  }


  iFile.close();
  fmt::MemoryWriter w5;
  w5.write("plzip -k -f -n {} {}",opts.numThreads,islandTxtFile.string());
  std::system(w5.c_str());
  w5.clear();




	  //boost::dynamic_bitset<uint8_t> bitmapC;
	  //std::vector<uint8_t> delimBytes;
	  //bitmapC.push_back(1);
	  //bitmapC.push_back(1);
	  //bitmapC.push_back(1);
	  //boost::to_block_range(bitmapC, std::back_inserter(delimBytes));
	  //islandPtr->write(reinterpret_cast<char*>(&delimBytes[0]), sizeof(delimBytes[0])* delimBytes.size());
  //iFile.close();
  //(*islandPtr) << redi::peof ;

  fmt::print(stderr, "\n start writing unmapped\n");

  bfs::path minceScriptFile = auxDir/"run_mince.sh";
  std::ofstream minceScript(minceScriptFile.string());
  //write unmapped sequences
  if(unmapped.size() > 0){
	  int uid = 0;

	  if(experiment.readLibraries().front().format().type != ReadType::SINGLE_END){
		  bfs::path unMappedFile_l = auxDir/"unmapped_1.fastq";
		  bfs::path unMappedFile_r = auxDir/"unmapped_2.fastq";


		  std::string mincePrefix = "unmin" ;

		  //pstream_ptr unlPtr(nullptr, closeStreamDeleter("leftUnmappedPtr"));
		  //pstream_ptr unrPtr(nullptr, closeStreamDeleter("rightUnmappedPtr"));

		  fmt::MemoryWriter w1;
		  w1.write("/home/rob/mince/bin/mince -e -l IU -1 {} -2 {} -p {} -o {}",unMappedFile_l.string(),unMappedFile_r.string(),opts.numThreads,mincePrefix);


		  //fmt::MemoryWriter w2;
		  //w2.write("cd {}",auxDir.string());

		  //std::cout << "\n " << w1.str() << "\n";


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

				  std::vector<std::string> qualities;
				  if(qualityScore){
					  qualities = split(seq,'|');
				  }


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

				 if(qualityScore){
					quality_1 << qualities[2] << "\n";
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
				  if(qualityScore){
					  quality_2 << qualities[3] <<"\n";
				  }
				  uFile_r << "\n";
				  uid++;

		  }


		  char cwd[1024];
		  if (getcwd(cwd, sizeof(cwd)) != NULL)
		         fprintf(stdout, "Current working dir: %s\n", cwd);

		  auto currDir = cwd ;


		  //std::system(w2.c_str());
		  //chdir(auxDir.string().c_str());
		  minceScript << w1.str() << "\n";

		  if (getcwd(cwd, sizeof(cwd)) != NULL)
		         fprintf(stdout, "Changed working dir: %s\n", cwd);

		  w1.clear();
		  w1.write("chmod +x {}",minceScriptFile.string());
		  std::system(w1.c_str());

		  //chdir(currDir);


		  w1.clear();
		  //w2.clear();
	  }else{
		  bfs::path unMappedFile_l = auxDir/"unmapped.fastq";

		  std::string mincePrefix = "unmin" ;
		  std::ofstream uFile_l(unMappedFile_l.string());

		  fmt::MemoryWriter w1;
		  w1.write("/home/rob/mince/bin/mince -e -l U -r {} -p {} -o {}",unMappedFile_l.string(),opts.numThreads,mincePrefix);



		  for(auto seqvec : unmapped){
			  for(auto seq : seqvec){

				  std::vector<std::string> qualities ;
				  if(qualityScore){
					  qualities = split(seq,'|');
				  }

				  uFile_l << "@"<<uid<< "\n";
				  for(int il = 0; il < seq.size() ; il++){
					  uFile_l << seq[il];
				  }
				  uFile_l << "\n";
				  uFile_l << "+"<< "\n";

				  for(int il = 0; il < seq.size() ; il++){
					  uFile_l << "I";
				  }
				  if(qualityScore){
					  quality_single << qualities[1] << "\n";
				  }

				  uFile_l << "\n";
				  uid++;
			  }
		  }


		  minceScript << w1.str() << "\n";
		  w1.clear();
		  w1.write("chmod +x {}",minceScriptFile.string());
		  std::system(w1.c_str());
		  w1.clear();
	  }

  }


  if(qualityScore){
	  if(experiment.readLibraries().front().format().type != ReadType::SINGLE_END){
		  quality_1.close();
		  fmt::MemoryWriter w6;
		  w6.write("plzip -k -f -n {} {}",opts.numThreads,qualityFilePath_1.string());
		  std::system(w6.c_str());
		  w6.clear();
		  w6.write("rm {}",qualityFilePath_1.string());
		  std::system(w6.c_str());
		  w6.clear();


		  quality_2.close();
		  fmt::MemoryWriter w7;
		  w7.write("plzip -k -f -n {} {}",opts.numThreads,qualityFilePath_2.string());
		  std::system(w7.c_str());
		  w7.clear();
		  w7.write("rm {}",qualityFilePath_2.string());
		  std::system(w7.c_str());
		  w7.clear();

	  }else{

		  quality_single.close();
		  fmt::MemoryWriter w6;
		  w6.write("plzip -k -f -n {} {}",opts.numThreads,qualityFilePath.string());
		  std::system(w6.c_str());
		  w6.clear();
		  w6.write("rm {}",qualityFilePath.string());
		  std::system(w6.c_str());
		  w6.clear();
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
      oa(cereal::make_nvp("quark_version", std::string(sailfish::version)));
      oa(cereal::make_nvp("samp_type", sampType));
      oa(cereal::make_nvp("frag_dist_length", fld->maxValue()));
      oa(cereal::make_nvp("num_targets", transcripts.size()));
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

