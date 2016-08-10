#ifndef QUARK_EQUIVALENCE_CLASS_HPP
#define QUARK_EQUIVALENCE_CLASS_HPP

#include <unordered_map>
#include <vector>
#include <thread>
#include <memory>
#include <mutex>

// Logger includes
#include "spdlog/spdlog.h"

#include "cuckoohash_map.hh"
#include "concurrentqueue.h"
#include "TranscriptGroup.hpp"

//structure that contains the encoded strings
//and counts
//

//typedef struct islands{
//	std::pair<int32_t,int32_t> content;
//	int id;
//};

struct QStrings {
	QStrings(const QStrings& o){
		qcodes = o.qcodes;
		count.store(o.count.load());
		intervals = o.intervals;

	}

	QStrings(std::string qcode,std::pair<int32_t,int32_t> interval,uint64_t countIn) :
		qcodes({qcode}),intervals({interval}){
		count.store(countIn);
	}

	std::vector<std::pair<int32_t,int32_t>> makeIslands() const{
		std::sort(intervals.begin(),intervals.end(),
					[](const std::pair<int32_t,int32_t>& p1, const std::pair<int32_t,int32_t>& p2) -> bool {
						return p1.first < p2.first ;
					});

		//do a linear search on index
		std::vector<std::pair<int32_t,int32_t>> correctedIslands;
		int currentindex;
		//std::cout << "du du";
		for(int i = 0; i < intervals.size(); i++){
			if(i == 0){
				auto temp = intervals[0];
				correctedIslands.push_back(temp);
				currentindex = 0;
			}else{
				//they are disjoint
				//time to create new island
				if(intervals[i].first > correctedIslands[currentindex].second){
					auto temp = intervals[i];
					correctedIslands.push_back(temp);
					currentindex++ ;
				}else{
					//Now we have a tricky situation
					//where current island has an intersection with
					//the corrected current island so we extend our
					// old island to include this
					if(intervals[i].second > correctedIslands[currentindex].second){
						correctedIslands[currentindex].second = intervals[i].second;
					}
					//otherwise we can swallow this small island
				}
			}

		}
		//intervals = correctedIslands ;
		//intervals = {{0,0}};
		return correctedIslands ;
	}

	mutable std::vector<std::string> qcodes ;
	std::atomic<uint64_t> count{0};
	mutable std::vector<std::pair<int32_t,int32_t>> intervals;

#if defined __APPLE__
        spin_lock writeMutex_;
#else
        std::mutex writeMutex_;
#endif

};

class QuarkEquivalenceClassBuilder{
	public:
		QuarkEquivalenceClassBuilder(std::shared_ptr<spdlog::logger> loggerIn) :
			logger_(loggerIn) {
	            countMap_.reserve(1000000);
	        }
		~QuarkEquivalenceClassBuilder() {}

		  void start() { active_ = true; }

		        bool finish() {
		            active_ = false;
		            size_t totalCount{0};
		            auto lt = countMap_.lock_table();
		            for (auto& kv : lt) {
		            //for (auto kv = countMap_.begin(); !kv.is_end(); ++kv) {
		                //kv.second.normalizeAux();
		                totalCount += kv.second.count;
		                countVec_.push_back(kv);
		                //code to merge intervals and create a map
		                // that can also be done while writing
		                // whatever
		                kv.second.intervals = kv.second.makeIslands();
		            }

		    	    logger_->info("Computed {} rich equivalence classes "
					  "for further processing", countVec_.size());
		            logger_->info("Counted {} total reads in the equivalence classes ",
		                    totalCount);
		            return true;
		        }


		        inline void addGroup(TranscriptGroup&& g,
		                                     std::string qcode,
											 std::pair<int32_t,int32_t> interval) {

		                    auto upfn = [&qcode,&interval](QStrings& x) -> void {
		                     // update the count
		                        // x.scoped_lock {
#if defined __APPLE__
		                    spin_lock::scoped_lock sl(x.writeMutex_);
#else
		                    std::lock_guard<std::mutex> lock(x.writeMutex_);
#endif
		                        x.count++;
		                        x.qcodes.push_back(qcode);
		                        x.intervals.push_back(interval);
		                    };
		                    QStrings v(qcode,interval,1);
		                    countMap_.upsert(g, upfn, v);
		                }

		        		// return the equivalence class
		        		// I assume having the same name for the
		        		// function won't bring distress
		                // to the world
		                std::vector<std::pair<const TranscriptGroup, QStrings>>& eqVec() {
		                    return countVec_;
		                }




	private:
        std::atomic<bool> active_;
	    cuckoohash_map<TranscriptGroup, QStrings, TranscriptGroupHasher> countMap_;
        std::vector<std::pair<const TranscriptGroup, QStrings> > countVec_;
    	std::shared_ptr<spdlog::logger> logger_;

};
//class that actually hold all the equivalence classes
#endif
