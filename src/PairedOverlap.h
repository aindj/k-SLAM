/* Copyright 2014 David Ainsworth
 *
 * This file is part of SLAM
 *
 * SLAM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SLAM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.

 * You should have received a copy of the GNU Affero General Public License
 * along with SLAM.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef PAIREDOVERLAP_H_
#define PAIREDOVERLAP_H_
#include "Overlap.h"
//#include <unordered_map>
#include <parallel/algorithm>
#include <inttypes.h>
#include <sstream>
namespace SLAM {
/*
 * Represents an overlap between a paired end read and a genbank entry
 * If the overlap is not supported by both reads then it can be stored as
 * an R1 overlap (hasR1=true) or as an R2 overlap (hasR2=true).
 * The combined score and inferred insert size are stored.
 */
class PairedOverlap {
 public:
  PairedOverlap() {}
  ;
  PairedOverlap(uint16_t combinedScore, uint32_t entryPosInArray, int refStart,
                int refEnd, uint32_t insertSize, bool hasR1, bool hasR2,
                Overlap r1Overlap, Overlap r2Overlap)
      : combinedScore(combinedScore),
        entryPosInArray(entryPosInArray),
        refStart(refStart),
        refEnd(refEnd),
        insertSize(insertSize),
        hasR1(hasR1),
        hasR2(hasR2),
        r1Overlap(r1Overlap),
        r2Overlap(r2Overlap) {}
  ;
  uint32_t combinedScore = 0;
  uint32_t entryPosInArray = 0;
  int refStart = 0;
  int refEnd = 0;
  uint32_t insertSize = 0;
  bool hasR1 = false;
  bool hasR2 = false;
  Overlap r1Overlap;
  Overlap r2Overlap;
};
/*
 * Contains all overlaps for a read pair
 */
class ReadPairAndOverlaps {
 public:
  ReadPairAndOverlaps() {}
  ;
  ReadPairAndOverlaps(uint32_t r1PosInArray, uint32_t r2PosInArray,
                      std::vector<PairedOverlap> alignmentPairs)
      : r1PosInArray(r1PosInArray),
        r2PosInArray(r2PosInArray),
        alignmentPairs(alignmentPairs) {}
  ;
  uint32_t r1PosInArray = 0;
  uint32_t r2PosInArray = 0;
  std::vector<PairedOverlap> alignmentPairs;
};

template <typename FASTQType>
PairedOverlap makePair(const Overlap &R1, const Overlap &R2, bool orientation,
                       const std::vector<FASTQType> &reads);
template <typename ITType, typename FASTQType>
ITType getPairsFromRead(const ITType first, const ITType last,
                        const std::vector<FASTQType> &reads,
                        std::vector<PairedOverlap> &pairedOverlaps);
template <typename ITType, typename FASTQType>
std::vector<PairedOverlap> getPairedOverlaps(
    ITType first, const ITType last, const std::vector<FASTQType> &reads);
template <typename FASTQType>
std::vector<ReadPairAndOverlaps> getDummyAlignmentPairsFromSingleEndReads(
    const std::vector<ReadAndOverlaps> &alignments,
    const std::vector<FASTQType> &reads);
uint32_t getMaxAllowedInsertSize(
    const std::vector<ReadPairAndOverlaps> &readsAndAlignments);
void screenPairedAlignmentsByScore(
    std::vector<ReadPairAndOverlaps> &readsAndAlignments, double fraction);
template <typename ITType>
std::vector<ReadPairAndOverlaps> getPerReadOverlaps(ITType begin, ITType end,
                                                    const uint32_t midpoint);
template <typename FASTQType>
void pseudoAssembly(std::vector<ReadPairAndOverlaps> &pairedAlignments,
                    const std::vector<FASTQType> &reads,
                    const GenbankIndex &index);
/*
 * Creates a pairedOverlap from an R1 and R2 overlap
 * the bool 'orientation' refers to whether the pair maps to the reference
 * R1....R2 (true) or R2....R1 (false)
 */
template <typename FASTQType>
inline PairedOverlap makePair(const Overlap &R1, const Overlap &R2,
                              bool orientation,
                              const std::vector<FASTQType> &reads) {
  int32_t refStart = std::min(R1.alignment.ref_begin, R2.alignment.ref_begin);
  int32_t refEnd = std::max(R1.alignment.ref_end, R2.alignment.ref_end);
  uint32_t insertSize = 0;
  if (orientation)
    insertSize = R2.relativePosition - R1.relativePosition +
                 reads[R2.readPosInArray].bases.size();
  else
    insertSize = R1.relativePosition - R2.relativePosition +
                 reads[R1.readPosInArray].bases.size();
  return PairedOverlap(R1.alignment.sw_score + R2.alignment.sw_score,
                       R2.entryPosInArray, refStart, refEnd, insertSize, true,
                       true, R1, R2);
}
/*
 * Function for finding alignment pairs (two alignments, one from each read in a pair)
 * Overlaps are sorted by read ID then by position on the genome.
 * The overlaps are iterated over, finding adjacent R1 and R2 overlaps (that have opposite orientation) and pairing them
 * Note: Each overlap may be placed in more than one pair. This is to ensure that the correct pairs are found.
 * this also avoids pairing each R1 with all R2s (N^2 scaling). Overlaps that aren't paired are placed
 * in dummy pairs with a blank overlap.
 */
template <typename ITType, typename FASTQType>
inline ITType getPairsFromRead(const ITType first, const ITType last,
                               const std::vector<FASTQType> &reads,
                               std::vector<PairedOverlap> &pairedOverlaps) {
  const auto midpoint = reads.size() / 2;
  auto readID = first->readPosInArray % midpoint;
  auto entry = first->entryPosInArray;
  auto lastR1 = last;
  auto lastR2 = last;
  bool lastR1Used = false;
  bool lastR2Used = false;
  auto lastR1RC = last;
  auto lastR2RC = last;
  bool lastR1UsedRC = false;
  bool lastR2UsedRC = false;
  PairedOverlap overlap;
  if (first == last) return last;
  auto current = first;
  while (current != last && (current->readPosInArray % midpoint == readID) &&
         (current->entryPosInArray == entry)) {
    if (current->readPosInArray < midpoint) {
      //R1
      if (current->revComp) {
        if (!lastR1UsedRC && lastR1RC != last) {
          pairedOverlaps.push_back(PairedOverlap(
              lastR1RC->alignment.sw_score, lastR1RC->entryPosInArray,
              lastR1RC->alignment.ref_begin, lastR1RC->alignment.ref_end, 0,
              true, false, *lastR1RC, Overlap()));
        }
        lastR1RC = current;
        lastR1UsedRC = false;
        if (lastR2 != last) {
          pairedOverlaps.push_back(makePair(*current, *lastR2, false, reads));
          lastR1UsedRC = true;
          lastR2Used = true;
        }
      } else {
        if (!lastR1Used && lastR1 != last) {
          pairedOverlaps.push_back(PairedOverlap(
              lastR1->alignment.sw_score, lastR1->entryPosInArray,
              lastR1->alignment.ref_begin, lastR1->alignment.ref_end, 0, true,
              false, *lastR1, Overlap()));
        }
        lastR1 = current;
        lastR1Used = false;
        if (lastR2RC != last) {
          pairedOverlaps.push_back(makePair(*current, *lastR2RC, false, reads));
          lastR1Used = true;
          lastR2UsedRC = true;
        }
      }
    } else {
      //R2
      if (current->revComp) {
        if (!lastR2UsedRC && lastR2RC != last) {
          pairedOverlaps.push_back(PairedOverlap(
              lastR2RC->alignment.sw_score, lastR2RC->entryPosInArray,
              lastR2RC->alignment.ref_begin, lastR2RC->alignment.ref_end, 0,
              false, true, Overlap(), *lastR2RC));
        }
        lastR2RC = current;
        lastR2UsedRC = false;
        if (lastR1 != last) {
          pairedOverlaps.push_back(makePair(*lastR1, *current, true, reads));
          lastR1Used = true;
          lastR2UsedRC = true;
        }
      } else {
        if (!lastR2Used && lastR2 != last) {
          pairedOverlaps.push_back(PairedOverlap(
              lastR2->alignment.sw_score, lastR2->entryPosInArray,
              lastR2->alignment.ref_begin, lastR2->alignment.ref_end, 0, false,
              true, Overlap(), *lastR2));
        }
        lastR2 = current;
        lastR2Used = false;
        if (lastR1RC != last) {
          pairedOverlaps.push_back(makePair(*lastR1RC, *current, true, reads));
          lastR1UsedRC = true;
          lastR2Used = true;
        }
      }
    }
    current++;
  }
  if (!lastR2Used && lastR2 != last) {
    pairedOverlaps.push_back(
        PairedOverlap(lastR2->alignment.sw_score, lastR2->entryPosInArray,
                      lastR2->alignment.ref_begin, lastR2->alignment.ref_end, 0,
                      false, true, Overlap(), *lastR2));
  }
  if (!lastR2UsedRC && lastR2RC != last) {
    pairedOverlaps.push_back(PairedOverlap(
        lastR2RC->alignment.sw_score, lastR2RC->entryPosInArray,
        lastR2RC->alignment.ref_begin, lastR2RC->alignment.ref_end, 0, false,
        true, Overlap(), *lastR2RC));
  }
  if (!lastR1Used && lastR1 != last) {
    pairedOverlaps.push_back(
        PairedOverlap(lastR1->alignment.sw_score, lastR1->entryPosInArray,
                      lastR1->alignment.ref_begin, lastR1->alignment.ref_end, 0,
                      true, false, *lastR1, Overlap()));
  }
  if (!lastR1UsedRC && lastR1RC != last) {
    pairedOverlaps.push_back(PairedOverlap(
        lastR1RC->alignment.sw_score, lastR1RC->entryPosInArray,
        lastR1RC->alignment.ref_begin, lastR1RC->alignment.ref_end, 0, true,
        false, *lastR1RC, Overlap()));
  }
  return current;
}
template <typename ITType, typename FASTQType>
inline std::vector<PairedOverlap> getPairedOverlaps(
    ITType first, const ITType last, const std::vector<FASTQType> &reads) {
  log("Pairing alignments");
  const auto midpoint = reads.size() / 2;
  __gnu_parallel::sort(first, last, [&](const Overlap & i, const Overlap & j) {
    if ((i.readPosInArray % midpoint) == (j.readPosInArray % midpoint)) {
      if (i.entryPosInArray == j.entryPosInArray)
        return i.relativePosition < j.relativePosition;
      else
        return i.entryPosInArray < j.entryPosInArray;

    } else
      return (i.readPosInArray % midpoint) < (j.readPosInArray % midpoint);
  });
  auto tempFn = [&](ITType first, ITType last) {
    std::vector<PairedOverlap> pairedOverlaps;
    while (first != last)
      first = getPairsFromRead(first, last, reads, pairedOverlaps);
    return pairedOverlaps;
  }
  ;
  auto breakAt = [&](ITType i, ITType j) {
    return (i->readPosInArray % midpoint) != (j->readPosInArray % midpoint);
  }
  ;
  std::vector<PairedOverlap> pairedOverlaps;
  parallelize(first, last, pairedOverlaps, tempFn, breakAt);
  return pairedOverlaps;
}

/*
 * As the taxonomy step of SLAM requires ReadPairAndOverlaps, single end reads
 * are turned into a dummy pair with and their overlaps are set as the R1
 * overlaps
 * of the pair
 */
template <typename FASTQType>
inline std::vector<ReadPairAndOverlaps>
getDummyAlignmentPairsFromSingleEndReads(
    const std::vector<ReadAndOverlaps> &alignments,
    const std::vector<FASTQType> &reads) {
  std::vector<ReadPairAndOverlaps> dummyPairs(alignments.size());
  __gnu_parallel::transform(alignments.begin(), alignments.end(),
                            dummyPairs.begin(), [](const ReadAndOverlaps & i) {
    std::vector<PairedOverlap> dummy(i.overlaps.size());
    std::transform(i.overlaps.begin(), i.overlaps.end(), dummy.begin(),
                   [](const Overlap & i) {
      return PairedOverlap(i.alignment.sw_score, i.entryPosInArray,
                           i.alignment.ref_begin, i.alignment.ref_end, 0, true,
                           false, i, Overlap());
    });
    return ReadPairAndOverlaps(i.readPosInArray, 0, dummy);
  });
  return dummyPairs;
}

/*
 * Calculates the distribution of inferred insert sizes of read pairs.
 * At this stage, some of the inferred pair overlaps will be incorrect
 * (due to both reads in a pair aligning to the same reference but not near
 * to one another). These alignments can be filtered out if the insert size is
 * known.
 * Firstly the lower/upper quartile and median of the insert size distribution
 * is
 * calculated. Then sizes are picked from the range 0 to
 * UQ + 2 * (UQ - LQ) (unless there is a significant spike in insert size).
 * The mean and standard deviations of these insert sizes are used to calculate
 * an
 * upper limit of the correct insert size. The limit used is mean + 6*stddev.
 */
inline uint32_t getMaxAllowedInsertSize(
    const std::vector<ReadPairAndOverlaps> &readsAndAlignments) {
  log("Calculating insert size distribution");
  std::vector<int32_t> insertSizes;
  for (auto &read : readsAndAlignments) {
    for (auto &alignmentPair : read.alignmentPairs) {
      if (alignmentPair.insertSize != 0)
        insertSizes.push_back(alignmentPair.insertSize);
    }
  }
  if (insertSizes.size() == 0) return UINT32_MAX;
  __gnu_parallel::sort(insertSizes.begin(), insertSizes.end());
  int32_t limit = 0;
  for (int i = 0; i < 99; i++) {
    if ((insertSizes[floor(insertSizes.size() * (i + 1) / 100.0)] -
         insertSizes[floor(insertSizes.size() * (i) / 100.0)]) > 1000) {
      limit = insertSizes[floor(insertSizes.size() * (i) / 100)];
      break;
    }
  }
  int32_t LQ = insertSizes[floor(insertSizes.size() * 0.25)];
  int32_t Med = insertSizes[floor(insertSizes.size() * 0.5)];
  int32_t UQ = insertSizes[floor(insertSizes.size() * 0.75)];
  int32_t lowerLimit = 0;  // LQ - 2 * (UQ - LQ);
  int32_t upperLimit = UQ + 2 * (UQ - LQ);
  if (limit) upperLimit = limit;
  if (upperLimit == 0) upperLimit = INT32_MAX;
  log("Lower quartile = " + std::to_string(LQ) + ", median = " +
      std::to_string(Med) + ", upper quartile = " + std::to_string(UQ));
  log("Calculating mean and standard deviation from pairs " +
      std::to_string(lowerLimit) + " < insert size < " +
      std::to_string(upperLimit));
  auto endPos = std::remove_if(insertSizes.begin(), insertSizes.end(),
                               [&](const int32_t i) {
    return i < lowerLimit || i > upperLimit;
  });
  insertSizes.resize(std::distance(insertSizes.begin(), endPos));
  double sum = std::accumulate(insertSizes.begin(), insertSizes.end(), 0.0);
  double mean = sum / insertSizes.size();
  double sqSum = std::inner_product(insertSizes.begin(), insertSizes.end(),
                                    insertSizes.begin(), 0.0);
  double stdDev = std::sqrt(sqSum / insertSizes.size() - mean * mean);
  log("Mean = " + std::to_string(mean) + ", standard deviation = " +
      std::to_string(stdDev));
  auto result = floor(mean + 6 * stdDev);
  return std::isnan(result) ? UINT_MAX : result;
}
inline void screenPairedAlignmentsByScore(
    std::vector<ReadPairAndOverlaps> &readsAndAlignments, double fraction) {
  size_t orignalNumber=0;
  for(auto & read : readsAndAlignments){
    orignalNumber+=read.alignmentPairs.size();
  }
  log("Screening all " + std::to_string(orignalNumber) + " alignment pairs by score");
  __gnu_parallel::for_each(readsAndAlignments.begin(), readsAndAlignments.end(),
                           [&](ReadPairAndOverlaps & read) {
    if (read.alignmentPairs.size() == 0) return;
    std::sort(read.alignmentPairs.begin(), read.alignmentPairs.end(),
              [](const PairedOverlap & i, const PairedOverlap & j) {
      return i.combinedScore > j.combinedScore;
    });

    unsigned topScore = read.alignmentPairs[0].combinedScore;
    auto cutoff =
        std::find_if(read.alignmentPairs.begin(), read.alignmentPairs.end(),
                     [&](const PairedOverlap & i) {
      return i.combinedScore < topScore * fraction;
    });
    read.alignmentPairs.erase(cutoff, read.alignmentPairs.end());
  });
  size_t newNumber=0;
  for(auto & read : readsAndAlignments){
    newNumber+=read.alignmentPairs.size();
  }
  log("Screened " + std::to_string(orignalNumber-newNumber) +
      " overlaps");
}
/*
 * Removes alignment pairs which have an insert size > insertSize
 * If remove==true then the removed alignment pairs are split up
 * and added back to the vector
 */
inline void screenPairedAlignmentsByInsertSize(
    std::vector<ReadPairAndOverlaps> &readsAndAlignments,
    const uint32_t insertSize, bool replace) {
  log("Screening all alignment pairs with insert size >= " +
      std::to_string(insertSize));
  __gnu_parallel::for_each(readsAndAlignments.begin(), readsAndAlignments.end(),
                           [&](ReadPairAndOverlaps & read) {
    std::sort(read.alignmentPairs.begin(), read.alignmentPairs.end(),
              [](const PairedOverlap & i, const PairedOverlap & j) {
      return i.insertSize < j.insertSize;
    });
    auto cutoff =
        std::find_if(read.alignmentPairs.begin(), read.alignmentPairs.end(),
                     [&](const PairedOverlap & i) {
      return i.insertSize > insertSize;
    });
    auto cutoffPos = std::distance(read.alignmentPairs.begin(), cutoff);
    if (replace) {
      read.alignmentPairs
          .reserve(read.alignmentPairs.size() +
                   std::distance(cutoff, read.alignmentPairs.end()));
      cutoff = read.alignmentPairs.begin() + cutoffPos;
      auto oldEnd = read.alignmentPairs.end();
      for (auto current = cutoff; current < oldEnd; current++) {
        read.alignmentPairs.emplace_back(
            current->r1Overlap.alignment.sw_score, current->entryPosInArray,
            current->r1Overlap.alignment.ref_begin,
            current->r1Overlap.alignment.ref_end, 0, true, false,
            current->r1Overlap, Overlap());
        current->combinedScore = current->r2Overlap.alignment.sw_score;
        current->hasR1 = false;
        current->insertSize = 0;
        current->r1Overlap = Overlap();
        current->refStart = current->r2Overlap.alignment.ref_begin;
        current->refEnd = current->r2Overlap.alignment.ref_end;
      }
    } else {
      read.alignmentPairs.erase(cutoff, read.alignmentPairs.end());
    }
  });
}
template <typename ITType>
inline std::vector<ReadPairAndOverlaps> getPerReadOverlaps(
    ITType begin, ITType end, const uint32_t midpoint) {
  log("Getting per read overlaps");
  size_t numOverlaps = 0;
  std::vector<ReadPairAndOverlaps> readsAndOverlaps;
  readsAndOverlaps.reserve(midpoint * 2);
  ReadPairAndOverlaps readAndOverlaps;
  uint32_t readPos = 0;
  for (auto overlap = begin; overlap != end; overlap++) {
    uint32_t thisReadPos =
        overlap->hasR1 ? overlap->r1Overlap.readPosInArray
                       : (overlap->r2Overlap.readPosInArray - midpoint);
    if (thisReadPos == readPos) {
    } else {
      if (readAndOverlaps.alignmentPairs.size()) {
        numOverlaps += readAndOverlaps.alignmentPairs.size();
        readsAndOverlaps.push_back(std::move(readAndOverlaps));
        readAndOverlaps.alignmentPairs.clear();
      }
      readPos = thisReadPos;
    }
    readAndOverlaps.alignmentPairs.push_back(std::move(*overlap));
    readAndOverlaps.r1PosInArray = thisReadPos;
    readAndOverlaps.r2PosInArray = thisReadPos + midpoint;
  }
  if (readAndOverlaps.alignmentPairs.size()) {
    numOverlaps += readAndOverlaps.alignmentPairs.size();
    readsAndOverlaps.push_back(readAndOverlaps);
  }
  log(std::to_string(readsAndOverlaps.size()) +
      " entries have k-mer overlaps");
  return readsAndOverlaps;
}
/*
 * For each genome: find alignments that overlap (by at least 20 bases) along the genome and form them into a
 * pseudo-assembly chain.
 * This allows reads that are in conserved genome regions to be joined together into long
 * chains that extend into unique sequence. This allows all reads to be classified to species
 * level instead of genus.
 * The new score is calculated by taking the average score per base of each read in the chain and multiplying by
 * the total number of bases in the chain.
 */
template <typename FASTQType>
void pseudoAssembly(std::vector<ReadPairAndOverlaps> &pairedAlignments,
                    const std::vector<FASTQType> &reads,
                    const GenbankIndex &index) {
  log("Performing a pseudo-assembly");
  class coverage {
   public:
    int start = 0;
    int stop = 0;
  };
  class entryAndOverlaps {
   public:
    uint32_t entryPos = 0;
    std::vector<std::pair<coverage, PairedOverlap *>> reads;
  };
  std::unordered_map<uint32_t, entryAndOverlaps> entriesAndOverlaps;
  for (auto &read : pairedAlignments) {
    for (auto &overlap : read.alignmentPairs) {
      auto it = entriesAndOverlaps.find(overlap.entryPosInArray);
      if (it == entriesAndOverlaps.end()) {
        entryAndOverlaps temp;
        temp.entryPos = overlap.entryPosInArray;
        coverage temp2;
        temp2.start = overlap.refStart;
        temp2.stop = overlap.refEnd;
        temp.reads.push_back({
          temp2, &overlap
        });
        entriesAndOverlaps.insert({
          overlap.entryPosInArray, temp
        });
      } else {
        coverage temp2;
        temp2.start = overlap.refStart;
        temp2.stop = overlap.refEnd;
        it->second.reads.push_back({
          temp2, &overlap
        });
      }
    }
  }
  for (auto &entry : entriesAndOverlaps) {
    std::sort(entry.second.reads.begin(), entry.second.reads.end(),
              [](const std::pair<coverage, PairedOverlap *> & i,
                 const std::pair<coverage, PairedOverlap *> & j) {
      return i.first.start < j.first.start;
    });
    auto chainStart = entry.second.reads.begin();
    int highestPos = -1000000;
    uint32_t score = 0;
    uint32_t numBases = 0;
    double perBaseScore = 0;
    for (auto overlap = entry.second.reads.begin();
         overlap != entry.second.reads.end(); overlap++) {
      if (overlap->first.start > highestPos - 20) {
        auto chainLength = std::distance(chainStart, overlap);
        if (chainLength > 1) {
          double length = highestPos - chainStart->first.start;
          double coverage = numBases / length;
          double avgScorePerBase = perBaseScore / chainLength;
          double score = coverage * avgScorePerBase * length;
          for (auto overlap2 = chainStart; overlap2 != overlap; overlap2++) {
            overlap2->second->combinedScore = score;
          }
        }
        chainStart = overlap;
        highestPos = overlap->first.stop;
        score = overlap->second->combinedScore;
        perBaseScore = overlap->second->combinedScore * 1.0 /
                       abs(overlap->second->refEnd - overlap->second->refStart);
        numBases = abs(overlap->second->refEnd - overlap->second->refStart);

      } else {
        if (overlap->first.stop > highestPos) highestPos = overlap->first.stop;
        score += overlap->second->combinedScore;
        perBaseScore +=
            overlap->second->combinedScore * 1.0 /
            abs(overlap->second->refEnd - overlap->second->refStart);
        numBases += abs(overlap->second->refEnd - overlap->second->refStart);
      }
    }
    auto chainLength = std::distance(chainStart, entry.second.reads.end());
    if (chainLength > 1) {
      double length = highestPos - chainStart->first.start;
      double coverage = numBases / length;
      double avgScorePerBase = perBaseScore / chainLength;
      //todo check if this gives too great a score to chains (consider eliminating number of bases
      // and using only length.)
      double score = coverage * avgScorePerBase * length;
      for (auto overlap2 = chainStart; overlap2 != entry.second.reads.end();
           overlap2++) {
        overlap2->second->combinedScore = score;
      }
    }
  }
  return;
}
}

#endif /* PAIREDOVERLAP_H_ */
