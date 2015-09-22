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
#ifndef OVERLAP_H_
#define OVERLAP_H_
#include <vector>
#include "ssw_cpp.h"
#include "Globals.h"
#include <algorithm>
#include "ParallelTools.h"
#include <inttypes.h>
namespace SLAM {

/*
 * Class representing an overlap between a FASTQRead and a Genbank Entry
 * Stores the locations of the read and entry in their respective arrays,
 *  the position of the read relative to the entry and whether the read is
 *  reverse complemented. In addition to this the Smith-Waterman alignment
 *  corresponding to this overlap can be computed and stored
 */

class OverlapTemp {
 public:
  OverlapTemp(const uint32_t readPosInArray, const uint32_t entryPosInArray,
          const int32_t relativePosition, const bool revComp)
      : readPosInArray(readPosInArray),
        entryPosInArray(entryPosInArray),
        relativePosition(relativePosition),
        revComp(revComp) {}
  ;
  OverlapTemp() {}
  ;
  uint32_t readPosInArray = 0;
  uint32_t entryPosInArray = 0;
  int32_t relativePosition = 0;
  bool revComp = false;
//  StripedSmithWaterman::Alignment alignment;
};
class Overlap {
 public:
  Overlap(const uint32_t readPosInArray, const uint32_t entryPosInArray,
          const int32_t relativePosition, const bool revComp)
      : readPosInArray(readPosInArray),
        entryPosInArray(entryPosInArray),
        relativePosition(relativePosition),
        revComp(revComp) {}
  Overlap(const OverlapTemp & overlap)
      : readPosInArray(overlap.readPosInArray),
        entryPosInArray(overlap.entryPosInArray),
        relativePosition(overlap.relativePosition),
        revComp(overlap.revComp) {}
  ;
  Overlap() {}
  ;
  uint32_t readPosInArray = 0;
  uint32_t entryPosInArray = 0;
  int32_t relativePosition = 0;
  bool revComp = false;
  StripedSmithWaterman::Alignment alignment;
};
/*
 * Tests equality of overlaps
 */
//template<typename overlap>
struct overlapEqual {
  bool operator()(const OverlapTemp &i, const OverlapTemp &j) const {
    return (i.readPosInArray == j.readPosInArray) &&
           (i.entryPosInArray == j.entryPosInArray) &&
           (abs(i.relativePosition - j.relativePosition) < 3);
  }
};
//template<typename overlap>
struct overlapSort {
  bool operator()(const OverlapTemp &i, const OverlapTemp &j) const {
    if (i.readPosInArray == j.readPosInArray) {
      if (i.entryPosInArray == j.entryPosInArray)
        return i.relativePosition < j.relativePosition;
      else
        return i.entryPosInArray < j.entryPosInArray;

    } else
      return i.readPosInArray < j.readPosInArray;
  }
};
/*
 * Represents all of the overlaps for a read
 */
class ReadAndOverlaps {
 public:
  ReadAndOverlaps() {}
  ;
  ReadAndOverlaps(uint32_t readPosInArray, std::vector<Overlap> alignments)
      : readPosInArray(readPosInArray), overlaps(alignments) {}
  ;
  uint32_t readPosInArray = 0;
  std::vector<Overlap> overlaps;
};

template <typename FASTQType, typename ITType>
ITType processPileUp(ITType first, ITType last,
                     const std::vector<FASTQType> &reads,
                     std::vector<OverlapTemp> &overlaps);
template <typename ITType, typename FASTQType>
std::vector<OverlapTemp> findOverlaps(ITType first, ITType last,
                                  const std::vector<FASTQType> &reads);
template <typename ITType, typename FASTQType>
std::vector<OverlapTemp> findOverlaps_parallel(const ITType begin, const ITType end,
                                           const std::vector<FASTQType> &reads);
template <typename ITType>
std::vector<ReadAndOverlaps> getPerReadOverlaps(const ITType begin,
                                                const ITType end);
void screenAndSortOverlaps(std::vector<ReadAndOverlaps> &readsAndAlignments,
                           const double cutoff);
void screenOverlapsByScoreThreshold(std::vector<Overlap> &overlaps,
                                    const double cutoff);

/*
 * Takes a range [begin -> end) of kMerAndData whose "kMer" is identical.
 * If sequences share kMers we infer an overlap between those sequences.
 * This function finds overlaps between the FASTQ sequences and entries in
 * Genbank.
 * The function iterates over kMers and produces a number of "Overlaps", a class
 * which describes a kMer overlap between a read and a Genbank entry.
 * These overlaps are added to a vector.
 * The variable kMer->revComp represents whether the kMer was complemented from
 * the input sequence when they were added to the vector. The variable sameComp
 * represents whether both kMers in this particular overlap have the same
 * revComp.
 * The variable "offset" represents the distance between the start of the
 * genbank
 * entry and the start of the read.
 * If sameComp==false then the reverse complement of the read is what overlaps
 * with the genbank entry
 * This function relies on the kmers being sorted first alphabetically and
 * secondly by whether they are from a Genbank entry or a read (entry kmers at
 * the beginning).
 */

template <typename FASTQType, typename ITType>
inline ITType processPileUp(ITType first, ITType last,
                            const std::vector<FASTQType> &reads,
                            std::vector<OverlapTemp> &overlaps) {
  if (!first->kMerData.getisFromGB()) {
    for (auto kMer = first; kMer != last; kMer++) {
      if (kMer->kMerInt != first->kMerInt) return kMer;
    }
    return last;
  }
//  uint32_t numGB=0;
//  auto ret=last;
//  for (auto GBkMer = first; GBkMer != last; GBkMer++) {
//    if (GBkMer->kMerInt != first->kMerInt){
//      ret=GBkMer;
//      break;
//    }
//    if (!GBkMer->kMerData.getisFromGB()) continue;
//    numGB++;
//  }
//  if(numGB>5)
//    return ret;
  for (auto notGBkMer = first; notGBkMer != last; notGBkMer++) {
    if (notGBkMer->kMerInt != first->kMerInt) return notGBkMer;
    if (notGBkMer->kMerData.getisFromGB()) continue;
//    uint32_t numGB=0;
    for (auto GBkMer = first; GBkMer != notGBkMer; GBkMer++) {
      if (!GBkMer->kMerData.getisFromGB()) break;
//      numGB++;
      const bool sameComp =
          GBkMer->kMerData.getrevComp() == notGBkMer->kMerData.getrevComp();
      //todo encode read length for locality
      const uint32_t offset =
          !GBkMer->kMerData.getrevComp()
              ? notGBkMer->kMerData.offset
              : reads[notGBkMer->kMerData.getreadIDOrGenbankID()].bases.size() -
                    notGBkMer->kMerData.offset - k;
      OverlapTemp overlap(notGBkMer->kMerData.getreadIDOrGenbankID(),
                      GBkMer->kMerData.getreadIDOrGenbankID(),
                      GBkMer->kMerData.offset - offset, !sameComp);
      overlaps.push_back(overlap);
//      if(numGB>numGBKMers)
//        break;
    }
  }
  return last;
}
//template <typename FASTQType, typename ITType>
//inline ITType processPileUp(ITType first, ITType last,
//                            const std::vector<FASTQType> &reads,
//                            std::vector<Overlap> &overlaps) {
//  for (auto notGBkMer = first; notGBkMer != last; notGBkMer++) {
//    if (notGBkMer->kMerInt != first->kMerInt) return notGBkMer;
//    if (notGBkMer->kMerData.getisFromGB()) continue;
//    for (auto GBkMer = first; GBkMer != last; GBkMer++) {
//      if (GBkMer->kMerInt != first->kMerInt) break;
//      if (!GBkMer->kMerData.getisFromGB()) continue;
//      const bool sameComp =
//          GBkMer->kMerData.getrevComp() == notGBkMer->kMerData.getrevComp();
//      //todo encode read length for locality
//      const uint32_t offset =
//          !GBkMer->kMerData.getrevComp()
//              ? notGBkMer->kMerData.offset
//              : reads[notGBkMer->kMerData.getreadIDOrGenbankID()].bases.size()
// -
//                    notGBkMer->kMerData.offset - k;
//      Overlap overlap(notGBkMer->kMerData.getreadIDOrGenbankID(),
//                      GBkMer->kMerData.getreadIDOrGenbankID(),
//                      GBkMer->kMerData.offset - offset, !sameComp);
//      overlaps.push_back(overlap);
//    }
//  }
//  return last;
//}
/*
 * Finds kmers with occurence > 1 and applies processPileUp to get Overlaps
 */
template <typename ITType, typename FASTQType>
inline std::vector<OverlapTemp> findOverlaps(ITType first, ITType last,
                                         const std::vector<FASTQType> &reads) {
  typedef typename ITType::value_type type;
  std::vector<OverlapTemp> overlaps;
  while (first != last) {
    if (first->kMerInt == 0) {
      first++;
      continue;
    }
    first = std::adjacent_find(first, last, [](type i, type j) {
      return i.kMerInt == j.kMerInt;
    });
    first = processPileUp(first, last, reads, overlaps);
  }
  return overlaps;
}
inline void removeLowQualityOverlaps(std::vector<OverlapTemp> & overlaps){
  log("Removing low quality overlaps");
  uint32_t readPos = 0;
  uint32_t entryPos = 0;
  uint32_t count = 0;
  auto writePos = overlaps.begin();
  for (auto overlapIt = overlaps.begin(); overlapIt != overlaps.end();overlapIt++) {
    if ((overlapIt->readPosInArray == readPos) && (overlapIt->entryPosInArray == entryPos)) {
      count++;
      if(count>200){
      }
      else{
        *writePos=*overlapIt;
        writePos++;
      }
    }
    else {
      readPos = overlapIt->readPosInArray;
      entryPos = overlapIt->entryPosInArray;
      count=1;
      *writePos=*overlapIt;
      writePos++;
    }
  }
  overlaps.resize(std::distance(overlaps.begin(),writePos+1));
}
/*
 * Parallelises the function "findOverlaps" ensuring that the range is never
 * split at a point which separates equal kMers (as this would miss an overlap)
 */
template <typename ITType, typename FASTQType>
inline std::vector<OverlapTemp> findOverlaps_parallel(
    const ITType begin, const ITType end, const std::vector<FASTQType> &reads) {
  log("Finding overlaps");
  std::vector<OverlapTemp> overlaps;
  parallelize(begin, end, overlaps, [&](ITType beg, ITType end) {
    return findOverlaps(beg, end, reads);
  },
              [](const ITType i, const ITType j) {
    return i->kMerInt != j->kMerInt;
  });
  //todo sort using midpoint data here
  __gnu_parallel::sort(overlaps.begin(), overlaps.end(), overlapSort());
  auto it = std::unique(overlaps.begin(), overlaps.end(), overlapEqual());
  overlaps.resize(std::distance(overlaps.begin(), it));
//  removeLowQualityOverlaps(overlaps);
  log("Found " + std::to_string(overlaps.size()) + " k-mer overlaps");
  return overlaps;
}
/*
 * Iterates over a collection of "Overlap"s and creates a vector
 * of "ReadAndOverlaps", each element contains an index of a read and
 * a vector of all of its overlaps. This allows each read to be analysed
 * individually
 */
template <typename ITType>
inline std::vector<ReadAndOverlaps> getPerReadOverlaps(const ITType begin,
                                                       const ITType end) {
  log("Getting per read overlaps");
  std::vector<ReadAndOverlaps> readsAndOverlaps;
  ReadAndOverlaps readAndOverlap;
  uint32_t readPos = 0;
  for (auto overlap = begin; overlap != end; overlap++) {
    if (overlap->readPosInArray == readPos) {
    } else {
      if (readAndOverlap.overlaps.size()) {
        readsAndOverlaps.push_back(readAndOverlap);
        readAndOverlap.overlaps.clear();
      }
      readPos = overlap->readPosInArray;
    }
    readAndOverlap.overlaps.push_back(*overlap);
    readAndOverlap.readPosInArray = overlap->readPosInArray;
  }
  if (readAndOverlap.overlaps.size()) {
    readsAndOverlaps.push_back(readAndOverlap);
  }
  log(std::to_string(readsAndOverlaps.size()) +
      " segments have k-mer overlaps");
  return readsAndOverlaps;
}

inline void screenOverlapsByScoreThreshold(std::vector<Overlap> &overlaps,
                                           const double cutoff) {
  log("Screening all alignments with score < " + std::to_string(cutoff));
  auto originalSize = overlaps.size();
  auto end =
      std::remove_if(overlaps.begin(), overlaps.end(), [&](const Overlap & i) {
    return i.alignment.sw_score < scoreThreshold;
  });
  overlaps.resize(std::distance(overlaps.begin(), end));
  log("Screened " + std::to_string(originalSize - overlaps.size()) +
      " overlaps");
  return;
}
}

#endif /* OVERLAP_H_ */
