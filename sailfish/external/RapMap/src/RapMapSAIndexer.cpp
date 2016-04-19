#include <iostream>
#include <mutex>
#include <vector>
#include <random>
#include <unordered_map>
#include <type_traits>
#include <fstream>
#include <cctype>
#include <cstdio>

#include "tclap/CmdLine.h"

#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>

#include "xxhash.h"
//#include "btree/btree_map.h"

#include "spdlog/spdlog.h"

// Jellyfish 2 include
#include "jellyfish/mer_dna.hpp"
#include "jellyfish/stream_manager.hpp"
#include "jellyfish/whole_sequence_parser.hpp"

#include "divsufsort.h"
#include "divsufsort64.h"
//#include "sais.h"
//#include "RSDic.hpp"
//#include "RSDicBuilder.hpp"
#include "bit_array.h"
#include "RapMapUtils.hpp"
#include "RapMapFileSystem.hpp"
#include "ScopedTimer.hpp"

#include "jellyfish/file_header.hpp"
#include "jellyfish/binary_dumper.hpp"
#include "jellyfish/thread_exec.hpp"
#include "jellyfish/hash_counter.hpp"
#include "jellyfish/mer_overlap_sequence_parser.hpp"
#include "jellyfish/mer_iterator.hpp"
#include "JFRaw.hpp"
//#include "bitmap.h"
#include "rank9b.h"

#include "sparsehash/dense_hash_map"

#include "IndexHeader.hpp"

#include <chrono>

using stream_manager = jellyfish::stream_manager<std::vector<std::string>::const_iterator>;
using single_parser = jellyfish::whole_sequence_parser<stream_manager>;
using TranscriptID = uint32_t;
using TranscriptIDVector = std::vector<TranscriptID>;
using KmerIDMap = std::vector<TranscriptIDVector>;
using MerMapT = jellyfish::cooperative::hash_counter<rapmap::utils::my_mer>;

bool buildSA(const std::string& outputDir,
             std::string& concatText,
             size_t tlen,
             std::vector<int64_t>& SA) {
// IndexT is the signed index type
// UIndexT is the unsigned index type
  using IndexT = int64_t;
  using UIndexT = uint64_t;
  bool success{false};

  std::ofstream saStream(outputDir + "sa.bin", std::ios::binary);
  {
    ScopedTimer timer;
    SA.resize(tlen, 0);
    IndexT textLen = static_cast<IndexT>(tlen);
    std::cerr << "Building suffix array . . . ";
    auto ret = divsufsort64(reinterpret_cast<unsigned char*>(
      const_cast<char*>(concatText.data())),
      SA.data(), tlen);

      success = (ret == 0);
      if (success) {
        std::cerr << "success\n";
        {
          ScopedTimer timer2;
          std::cerr << "saving to disk . . . ";
          cereal::BinaryOutputArchive saArchive(saStream);
          saArchive(SA);
          std::cerr << "done\n";
        }
      } else {
        std::cerr << "FAILURE: return code from libdivsufsort64() was " << ret << "\n";
        saStream.close();
        std::exit(1);
      }
      std::cerr << "done\n";
    }
    saStream.close();
    return success;
}

bool buildSA(const std::string& outputDir,
             std::string& concatText,
             size_t tlen,
             std::vector<int32_t>& SA) {
// IndexT is the signed index type
// UIndexT is the unsigned index type
  using IndexT = int32_t;
  using UIndexT = uint32_t;
  bool success{false};

  std::ofstream saStream(outputDir + "sa.bin", std::ios::binary);
  {
    ScopedTimer timer;
    SA.resize(tlen, 0);
    IndexT textLen = static_cast<IndexT>(tlen);
    std::cerr << "Building suffix array . . . ";
    auto ret = divsufsort(reinterpret_cast<unsigned char*>(
      const_cast<char*>(concatText.data())),
      SA.data(), tlen);

      success = (ret == 0);
      if (success) {
        std::cerr << "success\n";
        {
          ScopedTimer timer2;
          std::cerr << "saving to disk . . . ";
          cereal::BinaryOutputArchive saArchive(saStream);
          saArchive(SA);
          std::cerr << "done\n";
        }
      } else {
        std::cerr << "FAILURE: return code from libdivsufsort() was " << ret << "\n";
        saStream.close();
        std::exit(1);
      }
      std::cerr << "done\n";
    }
    saStream.close();
    return success;
}

// IndexT is the index type.
// int32_t for "small" suffix arrays
// int64_t for "large" ones
template <typename IndexT>
bool buildHash(const std::string& outputDir,
               std::string& concatText,
               size_t tlen,
               uint32_t k,
               std::vector<IndexT>& SA
              ) {
  // Now, build the k-mer lookup table
  google::dense_hash_map<uint64_t,
                      rapmap::utils::SAInterval<IndexT>,
                      rapmap::utils::KmerKeyHasher> khash;
  khash.set_empty_key(std::numeric_limits<uint64_t>::max());

  // The start and stop of the current interval
  IndexT start = 0, stop = 0;
  // An iterator to the beginning of the text
  auto textB = concatText.begin();
  auto textE = concatText.end();
  // The current k-mer as a string
  rapmap::utils::my_mer mer;
  bool currentValid{false};
  std::string currentKmer;
  std::string nextKmer;
  while (stop < tlen) {
      // Check if the string starting at the
      // current position is valid (i.e. doesn't contain $)
      // and is <= k bases from the end of the string
      nextKmer = concatText.substr(SA[stop], k);
      if (nextKmer.length() == k and
          nextKmer.find_first_of('$') == std::string::npos) {
          // If this is a new k-mer, then hash the current k-mer
          if (nextKmer != currentKmer) {
              if (currentKmer.length() == k and
                  currentKmer.find_first_of('$') == std::string::npos) {
                  mer = rapmap::utils::my_mer(currentKmer);
                  auto bits = mer.get_bits(0, 2*k);
                  auto hashIt = khash.find(bits);
                  if (hashIt == khash.end()) {
                      if (start > 1) {
                          if (concatText.substr(SA[start-1], k) ==
                              concatText.substr(SA[start], k)) {
                              std::cerr << "T[SA["
                                        << start-1 << "]:" << k << "] = "
                                        << concatText.substr(SA[start-1], k)
                                        << " = T[SA[" << start << "]:" << k << "]\n";
                              std::cerr << "start = " << start << ", stop = " << stop << "\n";
                              std::cerr << "[fatal (1)] THIS SHOULD NOT HAPPEN\n";
                              std::exit(1);
                          }
                      }
                      if (start == stop) {
                          std::cerr << "[fatal (2)] Interval is empty! (start = " << start
                                    << ") = (stop =  " << stop << ")\n";
                      }
                      if (start == stop) {
                          std::cerr << "[fatal (3)] Interval is empty! (start = " << start
                              << ") = (stop =  " << stop << ")\n";
                      }

                      khash[bits] = {start, stop};
                  } else {
                      std::cerr << "\nERROR (1): trying to add same suffix "
                                << currentKmer << " (len = "
                                << currentKmer.length() << ") multiple times!\n";
                      auto prevInt = hashIt->second;
                      std::cerr << "existing interval is ["
                                << prevInt.begin << ", " << prevInt.end << ")\n";
                      for (auto x = prevInt.begin; x < prevInt.end; ++x) {
                          auto suff = concatText.substr(SA[x], k);
                          for (auto c : suff) {
                              std::cerr << "*" << c << "*";
                          }
                          std::cerr << " (len = " << suff.length() <<")\n";
                      }
                      std::cerr << "new interval is ["
                                << start << ", " << stop << ")\n";
                      for (auto x = start; x < stop; ++x) {
                          auto suff = concatText.substr(SA[x], k);
                          for (auto c : suff) {
                              std::cerr << "*" << c << "*";
                          }
                          std::cerr << "\n";
                      }
                  }
              }
              currentKmer = nextKmer;
              start = stop;
          }
      } else {
          // If this isn't a valid suffix (contains a $)

          // If the previous interval was valid, put it
          // in the hash.
          if (currentKmer.length() == k and
              currentKmer.find_first_of('$') == std::string::npos) {
              mer = rapmap::utils::my_mer(currentKmer);
              auto bits = mer.get_bits(0, 2*k);
              auto hashIt = khash.find(bits);
              if (hashIt == khash.end()) {
                  if (start > 2) {
                      if (concatText.substr(SA[start-1], k) ==
                          concatText.substr(SA[start], k)) {
                          std::cerr << "T[SA["
                              << start-1 << "]:" << k << "] = "
                              << concatText.substr(SA[start-1], k)
                              << " = T[SA[" << start << "]:" << k << "]\n";
                          std::cerr << "start = " << start << ", stop = " << stop << "\n";
                          std::cerr << "[fatal (4)] THIS SHOULD NOT HAPPEN\n";
                          std::exit(1);
                      }
                  }
                  khash[bits] = {start, stop};
              } else {
                  std::cerr << "\nERROR (2): trying to add same suffix "
                            << currentKmer << "multiple times!\n";
                  auto prevInt = hashIt->second;
                  std::cerr << "existing interval is ["
                      << prevInt.begin << ", " << prevInt.end << ")\n";
                  for (auto x = prevInt.begin; x < prevInt.end; ++x) {
                      std::cerr << concatText.substr(SA[x], k) << "\n";
                  }
                  std::cerr << "new interval is ["
                      << start << ", " << stop << ")\n";
                  for (auto x = start; x < stop; ++x) {
                      std::cerr << concatText.substr(SA[x], k) << "\n";
                  }
              }

          }
          // The current interval is invalid and empty
          currentKmer = nextKmer;
          start = stop;
      }
      if (stop % 1000000 == 0) {
          std::cerr << "\r\rprocessed " << stop << " positions";
      }
      // We always update the end position
      ++stop;
  }
  if (start < tlen) {
      if (currentKmer.length() == k and
          currentKmer.find_first_of('$') != std::string::npos) {
          mer = rapmap::utils::my_mer(currentKmer);
          khash[mer.get_bits(0, 2*k)] = {start, stop};
      }
  }
  std::cerr << "\nkhash had " << khash.size() << " keys\n";
  std::ofstream hashStream(outputDir + "hash.bin", std::ios::binary);
  {
      ScopedTimer timer;
      std::cerr << "saving hash to disk . . . ";
      cereal::BinaryOutputArchive hashArchive(hashStream);
      hashArchive(k);
      khash.serialize(typename google::dense_hash_map<uint64_t,
                      rapmap::utils::SAInterval<IndexT>,
                      rapmap::utils::KmerKeyHasher>::NopointerSerializer(), &hashStream);
      //hashArchive(khash);
      std::cerr << "done\n";
  }
  hashStream.close();
  return true;
}

// To use the parser in the following, we get "jobs" until none is
// available. A job behaves like a pointer to the type
// jellyfish::sequence_list (see whole_sequence_parser.hpp).
template <typename ParserT>//, typename CoverageCalculator>
void indexTranscriptsSA(ParserT* parser,
            		std::string& outputDir,
        			bool noClipPolyA,
                    std::mutex& iomutex) {
    // Seed with a real random value, if available
    std::random_device rd;

    // Create a random uniform distribution
    std::default_random_engine eng(rd());

    std::uniform_int_distribution<> dis(0, 3);

    uint32_t n{0};
    uint32_t k = rapmap::utils::my_mer::k();
    std::vector<std::string> transcriptNames;
    std::vector<int64_t> transcriptStarts;
    //std::vector<uint32_t> positionIDs;
    constexpr char bases[] = {'A', 'C', 'G', 'T'};
    uint32_t polyAClipLength{10};
    uint32_t numPolyAsClipped{0};
    uint32_t numNucleotidesReplaced{0};
    std::string polyA(polyAClipLength, 'A');

    using TranscriptList = std::vector<uint32_t>;
    using eager_iterator = MerMapT::array::eager_iterator;
    using KmerBinT = uint64_t;

    bool clipPolyA = !noClipPolyA;
    size_t numDistinctKmers{0};
    size_t numKmers{0};
    size_t currIndex{0};
    std::cerr << "\n[Step 1 of 4] : counting k-mers\n";

    //rsdic::RSDicBuilder rsdb;
    std::vector<uint64_t> onePos; // Positions in the bit array where we should write a '1'
    fmt::MemoryWriter txpSeqStream;
    {
        ScopedTimer timer;
        while(true) {
            typename ParserT::job j(*parser);
            if(j.is_empty()) break;
            for(size_t i = 0; i < j->nb_filled; ++i) { // For each sequence
                std::string& readStr = j->data[i].seq;
                readStr.erase(std::remove_if(readStr.begin(), readStr.end(),
                               [](const char a) -> bool {
                                    return !(isprint(a));
                                }), readStr.end());

                uint32_t readLen = readStr.size();
                // First, replace non ATCG nucleotides
                for (size_t b = 0; b < readLen; ++b) {
                    readStr[b] = ::toupper(readStr[b]);
                    int c = jellyfish::mer_dna::code(readStr[b]);
                    // Replace non-ACGT bases with pseudo-random bases
                    if (jellyfish::mer_dna::not_dna(c)) {
                        char rbase = bases[dis(eng)];
                        c = jellyfish::mer_dna::code(rbase);
                        readStr[b] = rbase;
                        ++numNucleotidesReplaced;
                    }
                }

                // Now, do Kallisto-esque clipping of polyA tails
                if (clipPolyA) {
                    if (readStr.size() > polyAClipLength and
                            readStr.substr(readStr.length() - polyAClipLength) == polyA) {

                        auto newEndPos = readStr.find_last_not_of("Aa");
                        // If it was all As
                        if (newEndPos == std::string::npos) {
                            readStr.resize(0);
                        } else {
                            readStr.resize(newEndPos + 1);
                        }
                        ++numPolyAsClipped;
                    }
                }

                readLen  = readStr.size();
                uint32_t txpIndex = n++;

                // The name of the current transcript
                auto& recHeader = j->data[i].header;
                transcriptNames.emplace_back(recHeader.substr(0, recHeader.find_first_of(" \t")));

                // The position at which this transcript starts
                transcriptStarts.push_back(currIndex);

                txpSeqStream << readStr;
                txpSeqStream << '$';
                currIndex += readLen + 1;
                onePos.push_back(currIndex - 1);
            }
            if (n % 10000 == 0) {
                std::cerr << "\r\rcounted k-mers for " << n << " transcripts";
            }
        }
    }
    std::cerr << "\n";

    std::cerr << "Replaced " << numNucleotidesReplaced << " non-ATCG nucleotides\n";
    std::cerr << "Clipped poly-A tails from " << numPolyAsClipped << " transcripts\n";

    // Put the concatenated text in a string
    std::string concatText = txpSeqStream.str();
    // And clear the stream
    txpSeqStream.clear();

    // Build the suffix array
    size_t tlen = concatText.length();
    size_t maxInt = std::numeric_limits<int32_t>::max();
    bool largeIndex = (tlen + 1 > maxInt);

    // Make our dense bit arrray
    BIT_ARRAY* bitArray = bit_array_create(concatText.length());
    for (auto p : onePos) {
	    bit_array_set_bit(bitArray, p);
    }

    onePos.clear();
    onePos.shrink_to_fit();

    std::string rsFileName = outputDir + "rsd.bin";
    FILE* rsFile = fopen(rsFileName.c_str(), "w");
    {
      ScopedTimer timer;
      std::cerr << "Building rank-select dictionary and saving to disk ";
      bit_array_save(bitArray, rsFile);
      std::cerr << "done\n";
    }
    fclose(rsFile);
    bit_array_free(bitArray);


    std::ofstream seqStream(outputDir + "txpInfo.bin", std::ios::binary);
    {
        ScopedTimer timer;
        std::cerr << "Writing sequence data to file . . . ";
        cereal::BinaryOutputArchive seqArchive(seqStream);
        seqArchive(transcriptNames);
        if (largeIndex) {
          seqArchive(transcriptStarts);
        } else {
          std::vector<int32_t> txpStarts(transcriptStarts.size(), 0);
          size_t numTranscriptStarts = transcriptStarts.size();
          for (size_t i = 0; i < numTranscriptStarts; ++i) {
            txpStarts[i] = static_cast<int32_t>(transcriptStarts[i]);
          }
          transcriptStarts.clear();
          transcriptStarts.shrink_to_fit();
          {
            seqArchive(txpStarts);
          }
        }
        //seqArchive(positionIDs);
        seqArchive(concatText);
        std::cerr << "done\n";
    }
    seqStream.close();

    // clear stuff we no longer need
    //positionIDs.clear();
    //positionIDs.shrink_to_fit();
    transcriptStarts.clear();
    transcriptStarts.shrink_to_fit();
    transcriptNames.clear();
    transcriptNames.shrink_to_fit();
    // done clearing


    if (largeIndex) {
        largeIndex = true;
        std::cerr << "[info] Building 64-bit suffix array "
                     "(length of generalized text is " << tlen << " )\n";
        using IndexT = int64_t;
        std::vector<IndexT> SA;
        bool success = buildSA(outputDir, concatText, tlen, SA);
        if (!success) {
          std::cerr << "[fatal] Could not build the suffix array!\n";
          std::exit(1);
        }

        success = buildHash<IndexT>(outputDir, concatText, tlen, k, SA);
        if (!success) {
          std::cerr << "[fatal] Could not build the suffix interval hash!\n";
          std::exit(1);
        }
    } else {
        std::cerr << "[info] Building 32-bit suffix array "
                     "(length of generalized text is " << tlen << ")\n";
        using IndexT = int32_t;
        std::vector<IndexT> SA;
        bool success = buildSA(outputDir, concatText, tlen, SA);
        if (!success) {
          std::cerr << "[fatal] Could not build the suffix array!\n";
          std::exit(1);
        }

        success = buildHash<IndexT>(outputDir, concatText, tlen, k, SA);
        if (!success) {
          std::cerr << "[fatal] Could not build the suffix interval hash!\n";
          std::exit(1);
        }
    }

/*
    std::vector<int> SA(tlen, 0);
    std::ofstream saStream(outputDir + "sa.bin", std::ios::binary);
    {
        ScopedTimer timer;
        std::cerr << "Building suffix array . . . ";
        auto ret = sais(reinterpret_cast<unsigned char*>(
                        const_cast<char*>(concatText.c_str())),
                        SA.data(), tlen + 1);
        if (ret == 0) {
            std::cerr << "success\n";
            {
                ScopedTimer timer2;
                std::cerr << "saving to disk . . . ";
                cereal::BinaryOutputArchive saArchive(saStream);
                saArchive(SA);
		// don't actually need the LCP right now
                // saArchive(LCP);
                std::cerr << "done\n";
            }
        } else {
            std::cerr << "FAILURE: return code from sais() was " << ret << "\n";
	    std::exit(1);
        }
        std::cerr << "done\n";
    }
    saStream.close();
    */


    std::string indexVersion = "q1";
    IndexHeader header(IndexType::QUASI, indexVersion, true, k, largeIndex);
    // Finally (since everything presumably succeeded) write the header
    std::ofstream headerStream(outputDir + "header.json");
    {
	cereal::JSONOutputArchive archive(headerStream);
	archive(header);
    }
    headerStream.close();

}



int rapMapSAIndex(int argc, char* argv[]) {
    std::cerr << "RapMap Indexer\n";

    TCLAP::CmdLine cmd("RapMap Indexer");
    TCLAP::ValueArg<std::string> transcripts("t", "transcripts", "The transcript file to be indexed", true, "", "path");
    TCLAP::ValueArg<std::string> index("i", "index", "The location where the index should be written", true, "", "path");
    TCLAP::ValueArg<uint32_t> kval("k", "klen", "The length of k-mer to index", false, 31, "positive integer less than 32");
    TCLAP::SwitchArg  	      noClip("n", "noClip", "Don't clip poly-A tails from the ends of target sequences", false);
    cmd.add(transcripts);
    cmd.add(index);
    cmd.add(kval);
    cmd.add(noClip);

    cmd.parse(argc, argv);

    // stupid parsing for now
    std::string transcriptFile(transcripts.getValue());
    std::vector<std::string> transcriptFiles({ transcriptFile });

    uint32_t k = kval.getValue();
    if (k % 2 == 0) {
        std::cerr << "Error: k must be an odd value, you chose " << k << '\n';
        std::exit(1);
    } else
    if (k > 31) {
        std::cerr << "Error: k must not be larger than 31, you chose " << k << '\n';
        std::exit(1);
    }
    rapmap::utils::my_mer::k(k);

    std::string indexDir = index.getValue();
    if (indexDir.back() != '/') {
	indexDir += '/';
    }
    bool dirExists = rapmap::fs::DirExists(indexDir.c_str());
    bool dirIsFile = rapmap::fs::FileExists(indexDir.c_str());
    if (dirIsFile) {
        std::cerr << "The requested index directory already exists as a file.";
        std::exit(1);
    }
    if (!dirExists) {
        rapmap::fs::MakeDir(indexDir.c_str());
    }

    size_t maxReadGroup{1000}; // Number of reads in each "job"
    size_t concurrentFile{2}; // Number of files to read simultaneously
    size_t numThreads{2};
    stream_manager streams(transcriptFiles.begin(), transcriptFiles.end(), concurrentFile);
    std::unique_ptr<single_parser> transcriptParserPtr{nullptr};
    transcriptParserPtr.reset(new single_parser(4 * numThreads, maxReadGroup,
                              concurrentFile, streams));

    bool noClipPolyA = noClip.getValue();
    std::mutex iomutex;
    indexTranscriptsSA(transcriptParserPtr.get(), indexDir, noClipPolyA, iomutex);
    return 0;
}


  /** SANITY CHECKS RELATED TO THE RANK structure **/
    /*
    uint64_t nextSetBit{0};
    uint64_t offset{0};
    auto numBits = bit_array_length(bitArray);
    while (offset < numBits and bit_array_find_next_set_bit(bitArray, offset, &nextSetBit)) {
	if (concatText[nextSetBit] != '$') {
		std::cerr << "Bit # " << nextSetBit << " is set to 1, but the "
			  << "corresponding character in the text is " << concatText[nextSetBit] << "\n";
	}
	offset = nextSetBit + 1;
    }

    if (bit_array_num_bits_set(bitArray) != onePos.size()) {
	std::cerr << "ERROR: Bit array has " << bit_array_num_bits_set(bitArray)
		  << " bits set, but this should be " << onePos.size() << "!\n";
	std::exit(1);
    }

    rank9b bitmap(bitArray->words, bitArray->num_of_bits);
    for (size_t i = 0; i < onePos.size() - 1; ++i) {
	    auto pos = onePos[i];
	    auto r = bitmap.rank(pos+1);

	    if (r != i+1) {
		std::cerr << "rank should be " << i+1 << " but it's " << r << "\n";
		std::cerr << "text is " << concatText[pos] < "\n\n";
		std::cerr << "bit vector says " << (bit_array_get_bit(bitArray, pos) ? '1' : '0') << "\n";
	    }
    }

    std::ofstream rsStream(outputDir + "rsdSafe.bin", std::ios::binary);
    {
        ScopedTimer timer;
        rsdic::RSDic rsd;
        rsdb.Build(rsd);
        rsd.Save(rsStream);
        std::cerr << "done\n";
    }
    rsStream.close();
    */
    /** END OF SANITY CHECK **/
