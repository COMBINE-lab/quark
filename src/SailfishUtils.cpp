/**
>HEADER
    Copyright (c) 2013 Rob Patro robp@cs.cmu.edu

    This file is part of Sailfish.

    Sailfish is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Sailfish is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Sailfish.  If not, see <http://www.gnu.org/licenses/>.
<HEADER
**/

#include <boost/thread/thread.hpp>
#include <boost/filesystem.hpp>
#include <algorithm>
#include <iostream>
#include <tuple>
#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/range/join.hpp>

#include "spdlog/spdlog.h"

#include "gff.h"

#include "jellyfish/stream_manager.hpp"
#include "jellyfish/whole_sequence_parser.hpp"
#include "jellyfish/mer_dna.hpp"

#include "SailfishUtils.hpp"
#include "ReadExperiment.hpp"

//S_AYUSH_CODE
#include "UtilityFunctions.hpp"
//T_AYUSH_CODE

namespace sailfish {
    namespace utils {

        using std::string;
        using NameVector = std::vector<string>;
        using IndexVector = std::vector<size_t>;
        using KmerVector = std::vector<uint64_t>;

        /**
         * This function parses the library format string that specifies the format in which
         * the reads are to be expected.
         */
        LibraryFormat parseLibraryFormatString(std::string& fmt) {
            using std::vector;
            using std::string;
            using std::map;
            using std::stringstream;

            map<string, LibraryFormat> formatMap = {
                {"IU", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::U)},
                {"ISF", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::SA)},
                {"ISR", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::AS)},
                {"OU", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::U)},
                {"OSF", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::SA)},
                {"OSR", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::AS)},
                {"MU", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::U)},
                {"MSF", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::S)},
                {"MSR", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::A)},
                {"U", LibraryFormat(ReadType::SINGLE_END, ReadOrientation::NONE, ReadStrandedness::U)},
                {"SF", LibraryFormat(ReadType::SINGLE_END, ReadOrientation::NONE, ReadStrandedness::S)},
                {"SR", LibraryFormat(ReadType::SINGLE_END, ReadOrientation::NONE, ReadStrandedness::A)}};

            // inspired by http://stackoverflow.com/questions/236129/how-to-split-a-string-in-c
            // first convert the string to upper-case
            for (auto& c : fmt) { c = std::toupper(c); }


            auto libFmtIt = formatMap.find(fmt);

            if (libFmtIt == formatMap.end()) {
                stringstream errstr;
                errstr << "unknown library format string : " << fmt;
                throw std::invalid_argument(errstr.str());
            }

            return libFmtIt->second;
        }

        /**
         * Parses a set of __ordered__ command line options and extracts the relevant
         * read libraries from them.
         */
        std::vector<ReadLibrary> extractReadLibraries(boost::program_options::parsed_options& orderedOptions) {
            // The current (default) format for paired end data
            LibraryFormat peFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::U);
            // The current (default) format for single end data
            LibraryFormat seFormat(ReadType::SINGLE_END, ReadOrientation::NONE, ReadStrandedness::U);

            std::vector<ReadLibrary> peLibs{ReadLibrary(peFormat)};
            std::vector<ReadLibrary> seLibs{ReadLibrary(seFormat)};
            for (auto& opt : orderedOptions.options) {
                // Update the library type
                if (opt.string_key == "libType") {
                    auto libFmt = parseLibraryFormatString(opt.value[0]);
                    if (libFmt.type == ReadType::PAIRED_END) {
                        peFormat = libFmt;
                        peLibs.emplace_back(libFmt);
                    } else {
                        seFormat = libFmt;
                        seLibs.emplace_back(libFmt);
                    }
                }
                if (opt.string_key == "mates1") {
                    peLibs.back().addMates1(opt.value);
                }
                if (opt.string_key == "mates2") {
                    peLibs.back().addMates2(opt.value);
                }
                if (opt.string_key == "unmatedReads") {
                    seLibs.back().addUnmated(opt.value);
                }
            }

            std::vector<ReadLibrary> libs;
            libs.reserve(peLibs.size() + seLibs.size());
            for (auto& lib : boost::range::join(seLibs, peLibs)) {
                if (lib.format().type == ReadType::SINGLE_END) {
                    if (lib.unmated().size() == 0) {
                        // Didn't use default single end library type
                        continue;
                    }
                } else if (lib.format().type == ReadType::PAIRED_END) {
                    if (lib.mates1().size() == 0 or lib.mates2().size() == 0) {
                        // Didn't use default paired-end library type
                        continue;
                    }
                }
                libs.push_back(lib);
            }
            size_t numLibs = libs.size();
            std::cerr << "there " << ((numLibs > 1) ? "are " : "is ") << libs.size() << ((numLibs > 1) ? " libs\n" : " lib\n");
            return libs;
        }


        // for single end reads or orphans
        bool compatibleHit(LibraryFormat expected,
                           int32_t start, bool isForward, MateStatus ms) {
            auto expectedStrand = expected.strandedness;
            switch (ms) {
                case MateStatus::SINGLE_END:
                    if (isForward) { // U, SF
                        return (expectedStrand == ReadStrandedness::U or
                                expectedStrand == ReadStrandedness::S);
                    } else { // U, SR
                        return (expectedStrand == ReadStrandedness::U or
                                expectedStrand == ReadStrandedness::A);
                    }
                    break;
                case MateStatus::PAIRED_END_LEFT:
                    // "M"atching or same orientation is a special case
                    if (expected.orientation == ReadOrientation::SAME) {
                        return (expectedStrand == ReadStrandedness::U
                                or
                                (expectedStrand == ReadStrandedness::S and isForward)
                                or
                                (expectedStrand == ReadStrandedness::A and !isForward));
                    } else if (isForward) { // IU, ISF, OU, OSF, MU, MSF
                        return (expectedStrand == ReadStrandedness::U or
                                expectedStrand == ReadStrandedness::S);
                    } else { // IU, ISR, OU, OSR, MU, MSR
                        return (expectedStrand == ReadStrandedness::U or
                                expectedStrand == ReadStrandedness::A);
                    }
                    break;
                case MateStatus::PAIRED_END_RIGHT:
                    // "M"atching or same orientation is a special case
                    if (expected.orientation == ReadOrientation::SAME) {
                        return (expectedStrand == ReadStrandedness::U
                                or
                                (expectedStrand == ReadStrandedness::S and isForward)
                                or
                                (expectedStrand == ReadStrandedness::A and !isForward));
                    } else if (isForward) { // IU, ISR, OU, OSR, MU, MSR
                        return (expectedStrand == ReadStrandedness::U or
                                expectedStrand == ReadStrandedness::A);
                    } else { // IU, ISF, OU, OSF, MU, MSF
                        return (expectedStrand == ReadStrandedness::U or
                                expectedStrand == ReadStrandedness::S);
                    }
                    break;
                default:
                    // SHOULD NOT GET HERE
                    fmt::print(stderr, "WARNING: Could not associate known library type with read!\n");
                    return false;
                    break;
            }
            // SHOULD NOT GET HERE
            fmt::print(stderr, "WARNING: Could not associate known library type with read!\n");
            return false;
        }


        // for paired-end reads
        bool compatibleHit(LibraryFormat expected, LibraryFormat observed) {
            if (observed.type != ReadType::PAIRED_END) {
                // SHOULD NOT GET HERE
                fmt::print(stderr, "WARNING: PE compatibility function called with SE read!\n");
                return false;
            }

            auto es = expected.strandedness;
            auto eo = expected.orientation;

            auto os = observed.strandedness;
            auto oo = observed.orientation;

            // If the orientations are different, they are incompatible
            if (eo != oo) {
                return false;
            } else { // In this branch, the orientations are always compatible
                return (es == ReadStrandedness::U or
                        es == os);
            }
            // SHOULD NOT GET HERE
            fmt::print(stderr, "WARNING: Could not determine strand compatibility!");
            fmt::print(stderr, "please report this.\n");
            return false;
        }


        // Determine the library type of paired-end reads
        LibraryFormat hitType(int32_t end1Start, bool end1Fwd, uint32_t len1,
                              int32_t end2Start, bool end2Fwd, uint32_t len2, bool canDovetail) {

            // If the reads come from opposite strands
            if (end1Fwd != end2Fwd) {
                // and if read 1 comes from the forward strand
                if (end1Fwd) {
                    // then if read 1 start < read 2 start ==> ISF
                    // NOTE: We can't really delineate between inward facing reads that stretch
                    // past each other and outward facing reads --- the purpose of stretch is to help
                    // make this determinateion.
                    int32_t stretch = canDovetail ? len2 : 0;
                    if (end1Start <= end2Start + stretch) {
                        return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::SA);
                    } // otherwise read 2 start < read 1 start ==> OSF
                    else {
                        return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::SA);
                    }
                }
                // and if read 2 comes from the forward strand
                if (end2Fwd) {
                    // then if read 2 start <= read 1 start ==> ISR
                    // NOTE: We can't really delineate between inward facing reads that stretch
                    // past each other and outward facing reads --- the purpose of stretch is to help
                    // make this determinateion.
                    int32_t stretch = canDovetail ? len1 : 0;
                    if (end2Start <= end1Start + stretch) {
                        return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::AS);
                    } // otherwise, read 2 start > read 1 start ==> OSR
                    else {
                        return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::AS);
                    }
                }
            } else { // Otherwise, the reads come from the same strand
                if (end1Fwd) { // if it's the forward strand ==> MSF
                    return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::S);
                } else { // if it's the reverse strand ==> MSR
                    return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::A);
                }
            }
            // SHOULD NOT GET HERE
            spdlog::get("jointLog")->error("ERROR: Could not associate any known library type with read! "
                                           "Please report this bug!\n");
            std::this_thread::sleep_for(std::chrono::seconds(1));
            std::exit(-1);
            return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::NONE, ReadStrandedness::U);
        }

        uint64_t encode(uint64_t tid, uint64_t offset) {
            uint64_t res = (((tid & 0xFFFFFFFF) << 32) | (offset & 0xFFFFFFFF));
            return res;
        }

        uint32_t transcript(uint64_t enc) {
            uint32_t t = (enc & 0xFFFFFFFF00000000) >> 32;
            return t;
        }

        uint32_t offset(uint64_t enc) {
            uint32_t o = enc & 0xFFFFFFFF;
            return o;
        }

        size_t numberOfReadsInFastaFile(const std::string& fname) {
            constexpr size_t bufferSize = 16184;
            char buffer[bufferSize];
            std::ifstream ifile(fname, std::ifstream::in);
            ifile.rdbuf()->pubsetbuf(buffer, bufferSize);

            size_t numReads = 0;
            std::string s;
            while (ifile >> s) { if (s.front() == '>') { ++numReads; } }

            ifile.close();

            return numReads;
        }

        class ExpressionRecord {
            public:
                ExpressionRecord(const std::string& targetIn, uint32_t lengthIn, double effLengthIn,
                        std::vector<double>& expValsIn) :
                    target(targetIn), length(lengthIn), effLength(effLengthIn), expVals(expValsIn) {}

                ExpressionRecord( ExpressionRecord&& other ) {
                    std::swap(target, other.target);
                    length = other.length;
		    effLength = other.effLength;
                    std::swap(expVals, other.expVals);
                }

                ExpressionRecord(std::vector<std::string>& inputLine) {
                    if (inputLine.size() < 3) {
                        std::string err ("Any expression line must contain at least 3 tokens");
                        throw std::invalid_argument(err);
                    } else {
                        auto it = inputLine.begin();
                        target = *it; ++it;
                        length = std::stoi(*it); ++it;
			effLength = std::stod(*it); ++it;
                        for (; it != inputLine.end(); ++it) {
                            expVals.push_back(std::stod(*it));
                        }
                    }
                }

                std::string target;
                uint32_t length;
		double effLength;
                std::vector<double> expVals;
        };

        // From : http://stackoverflow.com/questions/9435385/split-a-string-using-c11
        std::vector<std::string> split(const std::string& str, int delimiter(int) = ::isspace){
            using namespace std;
            vector<string> result;
            auto e=str.end();
            auto i=str.begin();
            while (i != e) {
                i = find_if_not(i,e, delimiter);
                if (i == e) break;
                auto j = find_if(i,e, delimiter);
                result.push_back(string(i,j));
                i = j;
            }
            return result;
        }
    }
}

//=======
