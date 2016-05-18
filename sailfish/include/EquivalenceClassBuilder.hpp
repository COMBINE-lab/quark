#ifndef EQUIVALENCE_CLASS_BUILDER_HPP
#define EQUIVALENCE_CLASS_BUILDER_HPP

#include <unordered_map>
#include <vector>
#include <thread>
#include <memory>
#include <mutex>

// Logger includes
#include "spdlog/spdlog.h"

#include "cuckoohash_map.hh"
#include "concurrentqueue.h"
#include "Transcript.hpp"
#include "TranscriptGroup.hpp"
#include "SailfishSpinLock.hpp"


struct TGValue {
    TGValue(const TGValue& o) {
        weights = o.weights;
        count.store(o.count.load());
        readNames = o.readNames ;
        txpNames = o.txpNames ;
        positions = o.positions;
        matePositions = o.matePositions;
        status = o.status;
        matestatus = o.matestatus;
    }

    TGValue(std::vector<double>& weightIn, uint64_t countIn, std::string& readNamesIn, uint32_t txpNameIn, int32_t& pos, int32_t& matepos, bool& stat, bool& matestat) :
        weights(weightIn), readNames({readNamesIn}), txpNames({txpNameIn}), positions({pos}), matePositions({matepos}), status({stat}), matestatus({matestat}) {
    count.store(countIn);
    }
    // const is a lie
    void normalizeAux() const {
        double sumOfAux{0.0};
        for (size_t i = 0; i < weights.size(); ++i) {
            sumOfAux += weights[i];
        }
        double norm = 1.0 / sumOfAux;
        for (size_t i = 0; i < weights.size(); ++i) {
            weights[i] *= norm;
        }
        /* LOG SPACE
        double sumOfAux = salmon::math::LOG_0;
        for (size_t i = 0; i < weights.size(); ++i) {
            sumOfAux = salmon::math::logAdd(sumOfAux, weights[i]);
        }
        for (size_t i = 0; i < weights.size(); ++i) {
            weights[i] = std::exp(weights[i] - sumOfAux);
        }
        */
    }

    // forget synchronizing this for the time being
    mutable std::vector<double> weights;
    std::atomic<uint64_t> count{0};
    mutable std::vector<std::string> readNames ;
    mutable std::vector<uint32_t> txpNames ;
    mutable std::vector<int32_t> positions ;
    mutable std::vector<int32_t> matePositions ;
    mutable std::vector<bool> status;
    mutable std::vector<bool> matestatus;
    //spin_lock l;
    // only one writer thread at a time
    //
#if defined __APPLE__
        spin_lock writeMutex_;
#else
        std::mutex writeMutex_;
#endif
};

class EquivalenceClassBuilder {
    public:
        EquivalenceClassBuilder(std::shared_ptr<spdlog::logger> loggerIn) :
		logger_(loggerIn) {
            countMap_.reserve(1000000);
        }

        ~EquivalenceClassBuilder() {}

        void start() { active_ = true; }

        bool finish() {
            active_ = false;
            size_t totalCount{0};
            auto lt = countMap_.lock_table();
            for (auto& kv : lt) {
            //for (auto kv = countMap_.begin(); !kv.is_end(); ++kv) {
                kv.second.normalizeAux();
                totalCount += kv.second.count;
                countVec_.push_back(kv);
            }

    	    logger_->info("Computed {} rich equivalence classes "
			  "for further processing", countVec_.size());
            logger_->info("Counted {} total reads in the equivalence classes ",
                    totalCount);
            return true;
        }

        inline void insertGroup(TranscriptGroup g, uint32_t count) {
            std::string dummy ;
            uint32_t dummy1;
            int32_t dummy2;
            bool dummy3;
            std::vector<double> weights(g.txps.size(), 0.0);
            //auto updatefn = [count](TGValue& x) { x.count = count; };
            TGValue v(weights, count,dummy,dummy1,dummy2,dummy2,dummy3,dummy3);
            countVec_.push_back(std::make_pair(g, v));
            //countMap_.upsert(g, updatefn, v);
        }

        inline void addGroup(TranscriptGroup&& g,
                             std::vector<double>& weights,
                             std::string &readName,
                             uint32_t &txpName,
                             int32_t &position,
                             int32_t &matepos,
                             bool &stat,
                             bool &matestat) {

            auto upfn = [&weights, &readName, &txpName, &position, &matepos, &stat, &matestat](TGValue& x) -> void {
             // update the count
                // x.scoped_lock {
#if defined __APPLE__
            spin_lock::scoped_lock sl(x.writeMutex_);
#else
            std::lock_guard<std::mutex> lock(x.writeMutex_);
#endif
                x.count++;
                x.readNames.push_back(readName);
                x.txpNames.push_back(txpName);
                x.positions.push_back(position);
                x.matePositions.push_back(matepos);
                x.status.push_back(stat);
                x.matestatus.push_back(matestat);
                // update the weights
                for (size_t i = 0; i < x.weights.size(); ++i) {
                    // Possibly atomicized in the future
                    weights[i] += x.weights[i];
                    /* LOG SPACE
                    x.weights[i] =
                        salmon::math::logAdd(x.weights[i], weights[i]);
                    */
                }
            };
            TGValue v(weights, 1,readName,txpName,position,matepos,stat,matestat);
            countMap_.upsert(g, upfn, v);
        }

        std::vector<std::pair<const TranscriptGroup, TGValue>>& eqVec() {
            return countVec_;
        }

    private:
        std::atomic<bool> active_;
	    cuckoohash_map<TranscriptGroup, TGValue, TranscriptGroupHasher> countMap_;
        std::vector<std::pair<const TranscriptGroup, TGValue>> countVec_;
    	std::shared_ptr<spdlog::logger> logger_;
};

#endif // EQUIVALENCE_CLASS_BUILDER_HPP
