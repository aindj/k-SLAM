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
#ifndef SMITHWATERMAN_H_
#define SMITHWATERMAN_H_
#include "Overlap.h"
#include "GenbankTools.h"
#include <string>
#include "ssw_cpp.h"
#include <vector>
#include <inttypes.h>
#include "Globals.h"
namespace SLAM {
StripedSmithWaterman::Alignment smithWaterman(
    const std::string &ref, const std::string &query, const uint32_t offset,
    const uint32_t length, const StripedSmithWaterman::Aligner &aligner,
    const StripedSmithWaterman::Filter &filter);
// template <typename ITType, typename FASTQType>
// void performSmithWatermanOnRange_parallel(
//    ITType begin, ITType end, const std::vector<FASTQType> &reads,
//    const GenbankIndex &index);
//template <typename FASTQType>
//inline void pairwiseAlign(Overlap & overlap,
//                   const std::vector<FASTQType> &reads,
//                   const GenbankIndex &index){
//  const std::string &queryBases = reads[overlap.readPosInArray].bases;
//  const std::string &refBases = index.entries[overlap.entryPosInArray].bases;
//  int refStartPos = std::max(overlap.relativePosition, 0);
//  int length=queryBases.size();
//  int32_t score=0;
//  std::cout<<(overlap.revComp ? reverseComplement(queryBases) : queryBases)<<"\n";
//  std::cout<<refBases.substr(refStartPos,length)<<"\n";
//  if(overlap.revComp){
//    for(int i=0;i<queryBases.size();i++){
//        switch(refBases[refStartPos+i]){
//          case 'A':
//            if(queryBases[length-1-i]=='T')
//              score+=2;
//            break;
//          case 'C':
//            if(queryBases[length-1-i]=='G')
//              score+=2;
//            break;
//          case 'T':
//            if(queryBases[length-1-i]=='A')
//              score+=2;
//            break;
//          case 'G':
//            if(queryBases[length-1-i]=='C')
//              score+=2;
//            break;
//          default:
//            score-=3;
//            break;
//        }
//    }
//  }
//  else{
//    for(int i=0;i<queryBases.size();i++){
//      if(queryBases[i] == refBases[refStartPos+i])
//        score+=2;
//      else
//        score-=3;
//    }
//  }
//  std::cout<<score<<"\n";
//  overlap.alignment.cigar=nullptr;
//  overlap.alignment.cigarLen=0;
//  overlap.alignment.query_begin=0;
//  overlap.alignment.query_end=queryBases.size()-1;
//  overlap.alignment.ref_begin=refStartPos;
//  overlap.alignment.ref_end=refStartPos+queryBases.size()-1;
//  overlap.alignment.sw_score=score;
//}
template <typename RAIter, typename FASTQType>
void performSmithWatermanOnRange(RAIter begin, RAIter end,
                                 const std::vector<FASTQType> &reads,
                                 const GenbankIndex &index);
/*
 * Takes a range (specified by random access iterators) of "ReadAndOverlaps"
 * and performs a Smith-Waterman alignment on each.
 *
 * The iterators point to a "ReadAndOverlaps" which is a class representing a
 * single read's overlaps.
 * The function performs a Smith-Waterman alignment of each of these overlaps.
 * The function is optimised in that it sorts each alignment by reference
 *sequence
 * so that multiple overlaps with the same reference sequence only have to be
 * aligned once
 * The functions performs the alignment and stores it in each Overlap's
 *alignment
 * field
 */
template <typename RAIter, typename FASTQType>
inline void performSmithWatermanOnRange(RAIter begin, RAIter end,
                                        const std::vector<FASTQType> &reads,
                                        const GenbankIndex &index) {
  //  log("Performing pairwise Smith-Waterman");
  if (begin == end) return;
  typedef typename RAIter::value_type type;
  class overlapAndRefBases {
   public:
    overlapAndRefBases(Overlap *overlap, std::string &&bases)
        : overlap(overlap), bases(bases) {}
    ;
    Overlap *overlap = nullptr;
    std::string bases;
  };
  const StripedSmithWaterman::Aligner aligner(match, misMatch, gapOpen,
                                              gapExtend);
  StripedSmithWaterman::Filter filter;
  filter.report_begin_position = true;
  //todo this doesnt stop cigar
  filter.report_cigar = reportCigar;
  filter.score_filter = scoreThreshold;
  std::vector<overlapAndRefBases> overlapsAndBases;
  auto current = begin;
  while (current != end) {
    auto initial = current;
    const std::string &queryBases = reads[current->readPosInArray].bases;
    overlapsAndBases.clear();
    for (; current != end && current->readPosInArray == initial->readPosInArray;
         current++) {
      auto &overlap = *current;
      int refStartPos = std::max(overlap.relativePosition, 0);
      if (refStartPos >= 2) refStartPos -= 2;
      auto refBases = index.entries[overlap.entryPosInArray].bases
          .substr(refStartPos, queryBases.size() + 4);
      if (overlap.revComp) inPlaceReverseComplement(refBases);
      overlapsAndBases.emplace_back(&overlap, std::move(refBases));
    }
    std::sort(overlapsAndBases.begin(), overlapsAndBases.end(),
              [](const overlapAndRefBases & i, const overlapAndRefBases & j) {
      return i.bases < j.bases;
    });
    std::string tempBases;
    StripedSmithWaterman::Alignment alignment;
    for (auto overlap = overlapsAndBases.begin();
         overlap != overlapsAndBases.end(); overlap++) {
      int refStartPos = std::max(overlap->overlap->relativePosition, 0);
      if (refStartPos >= 2) refStartPos -= 2;
      if (overlap->bases != tempBases) {
        alignment = smithWaterman(overlap->bases, queryBases, 0,
                                  queryBases.size() + 4, aligner, filter);
        tempBases = overlap->bases;
      }
      overlap->overlap->alignment = alignment;
      if (overlap->overlap->revComp) {
        if (reportCigar && overlap->overlap->alignment.cigar) {
          std::reverse(overlap->overlap->alignment.cigar,
                       overlap->overlap->alignment.cigar +
                           overlap->overlap->alignment.cigarLen);
        }
        auto temp = overlap->overlap->alignment.ref_begin;
        overlap->overlap->alignment.ref_begin =
            overlap->bases.size() - (overlap->overlap->alignment.ref_end + 1);
        overlap->overlap->alignment.ref_end =
            overlap->bases.size() - (temp + 1);
        temp = overlap->overlap->alignment.query_begin;
        overlap->overlap->alignment.query_begin =
            queryBases.size() - (overlap->overlap->alignment.query_end + 1);
        overlap->overlap->alignment.query_end = queryBases.size() - (temp + 1);
      }
      overlap->overlap->alignment.ref_begin += refStartPos;
      overlap->overlap->alignment.ref_end += refStartPos;
    }
  }
  return;
}
template <typename RAIter, typename FASTQType>
inline void performSmithWatermanOnRange2(RAIter begin, RAIter end,
                                        const std::vector<FASTQType> &reads,
                                        const GenbankIndex &index) {
  //  log("Performing pairwise Smith-Waterman");
  if (begin == end) return;
  typedef typename RAIter::value_type type;
  const StripedSmithWaterman::Aligner aligner(match, misMatch, gapOpen,
                                              gapExtend);
  StripedSmithWaterman::Filter filter;
  filter.report_begin_position = true;
  //todo this doesnt stop cigar
  filter.report_cigar = reportCigar;
  filter.score_filter = scoreThreshold;
  auto current = begin;
  while (current != end) {
    auto initial = current;
    const std::string &queryBases = reads[current->readPosInArray].bases;
    for (; current != end && current->readPosInArray == initial->readPosInArray;
         current++) {
      auto &overlap = *current;
      int refStartPos = std::max(overlap.relativePosition, 0);
      auto refBases = index.entries[overlap.entryPosInArray].bases
          .substr(refStartPos, queryBases.size());
      if (overlap.revComp) inPlaceReverseComplement(refBases);
      aligner.Align(queryBases.c_str(), refBases.c_str(),
                    std::min(int(queryBases.size()), int(refBases.size())), filter,
                    &overlap.alignment);
      if (overlap.revComp) {
        if (reportCigar && overlap.alignment.cigar) {
          std::reverse(overlap.alignment.cigar,
                       overlap.alignment.cigar +
                       overlap.alignment.cigarLen);
        }
        auto temp = overlap.alignment.ref_begin;
        overlap.alignment.ref_begin =
            refBases.size() - (overlap.alignment.ref_end + 1);
        overlap.alignment.ref_end =
            refBases.size() - (temp + 1);
        temp = overlap.alignment.query_begin;
        overlap.alignment.query_begin =
            queryBases.size() - (overlap.alignment.query_end + 1);
        overlap.alignment.query_end = queryBases.size() - (temp + 1);
      }
      overlap.alignment.ref_begin += refStartPos;
      overlap.alignment.ref_end += refStartPos;
    }
  }
  return;
}
template <typename RAIter, typename FASTQType>
inline void performSmithWatermanOnRange_parallel(
    RAIter begin, RAIter end, const std::vector<FASTQType> &reads,
    const GenbankIndex &index) {
  log("Performing pairwise Smith-Waterman");
  ;
  parallelForEachWithSplit(begin, end,
                           [&](const std::vector<Overlap>::iterator begin,
                               const std::vector<Overlap>::iterator end) {
    return performSmithWatermanOnRange2(begin, end, reads, index);
  },
                           [](const std::vector<Overlap>::iterator begin,
                              const std::vector<Overlap>::iterator end) {
    return begin->readPosInArray != end->readPosInArray;
  });
}
/*
 * Perform a pairwise Smith Waterman alignment of "query" against "ref"
 * "ref" refers to the reference substring, "offset" refers to the
 * position along the reference substring (usually 0) that the alignment starts
 * at
 * and "length" refers to the length of the reference substring that
 * the query will be aligned to
 */
inline StripedSmithWaterman::Alignment smithWaterman(
    const std::string &ref, const std::string &query, const uint32_t offset,
    const uint32_t length, const StripedSmithWaterman::Aligner &aligner,
    const StripedSmithWaterman::Filter &filter) {
  StripedSmithWaterman::Alignment alignment;
  aligner.Align(query.c_str(), ref.c_str() + offset,
                std::min(int(length), int(ref.size()) - int(offset)), filter,
                &alignment);
  //  std::cout<<alignment.sw_score<<"\n";
  return alignment;
}
/*
 * Parallelises performSmithWaterman on range.
 * The reason that this is not done with a standard for_each is that
 * each SmithWaterman alignment needs a "filter", This is expensive to
 * create for each alignment so this method creates a filter for each
 * thread. todo look at whether it is quicker to do a simple for_each
 */
// template <typename ITType, typename FASTQType>
// inline void performSmithWatermanOnRange_parallel(
//    ITType begin, ITType end, const std::vector<FASTQType> &reads,
//    const GenbankIndex &index) {
//  log("Performing pairwise Smith-Waterman");
//  performSmithWatermanOnRange(
//          begin, end, reads, index);
//  return;
//
//  auto breakAt = [](ITType i, ITType j) { return true; };
//  auto functor = [&](ITType beg, ITType end) {
//    performSmithWatermanOnRange(
//        beg, end, reads, index);
//  };
//  parallelForEachWithSplit(begin, end, functor, breakAt);
//}
//template <typename RAIter, typename FASTQType>
//inline void performSmithWatermanOnRange(RAIter begin, RAIter end,
//                                        const std::vector<FASTQType> &reads,
//                                        const GenbankIndex &index) {
//  log("Performing pairwise Smith-Waterman");
//  typedef typename RAIter::value_type type;
//  class overlapAndRefBases {
//   public:
//    overlapAndRefBases(Overlap *overlap, std::string &&bases)
//        : overlap(overlap), bases(bases) {};
//    Overlap *overlap = nullptr;
//    std::string bases;
//  };
//  const StripedSmithWaterman::Aligner aligner(match, misMatch, gapOpen,
//                                              gapExtend);
//  StripedSmithWaterman::Filter filter;
//  filter.report_begin_position = true;
//  filter.report_cigar = cigar;
//  filter.score_filter = scoreThreshold;
//  __gnu_parallel::for_each(begin, end, [&, filter](type &readAndAlignments) {
//    //  std::for_each(begin,end,[&](type & readAndAlignments){
//    const std::string &queryBases =
//        reads[readAndAlignments.readPosInArray].bases;
//    std::vector<overlapAndRefBases> overlapsAndBases;
//    for (auto &overlap : readAndAlignments.overlaps) {
//      int refStartPos = std::max(overlap.relativePosition, 0);
//      if (refStartPos >= 2) refStartPos -= 2;
//      auto refBases = index.entries[overlap.entryPosInArray]
//                          .bases.substr(refStartPos, queryBases.size() + 4);
//      if (overlap.revComp) inPlaceReverseComplement(refBases);
//      overlapsAndBases.emplace_back(&overlap, std::move(refBases));
//    }
//    std::sort(overlapsAndBases.begin(), overlapsAndBases.end(),
//              [](const overlapAndRefBases &i,
//                 const overlapAndRefBases &j) { return i.bases < j.bases; });
//    std::string tempBases;
//    StripedSmithWaterman::Alignment alignment;
//    for (auto overlap = overlapsAndBases.begin();
//         overlap != overlapsAndBases.end(); overlap++) {
//      int refStartPos = std::max(overlap->overlap->relativePosition, 0);
//      if (refStartPos >= 2) refStartPos -= 2;
//      if (overlap->bases != tempBases) {
//        alignment = smithWaterman(overlap->bases, queryBases, 0,
//                                  queryBases.size() + 4, aligner, filter);
//        tempBases = overlap->bases;
//      }
//      overlap->overlap->alignment = alignment;
//      if (overlap->overlap->revComp) {
//        std::reverse(overlap->overlap->alignment.cigar.begin(),
//                     overlap->overlap->alignment.cigar.end());
//        auto temp = overlap->overlap->alignment.ref_begin;
//        overlap->overlap->alignment.ref_begin =
//            overlap->bases.size() - (overlap->overlap->alignment.ref_end + 1);
//        overlap->overlap->alignment.ref_end =
//            overlap->bases.size() - (temp + 1);
//        temp = overlap->overlap->alignment.query_begin;
//        overlap->overlap->alignment.query_begin =
//            queryBases.size() - (overlap->overlap->alignment.query_end + 1);
//        overlap->overlap->alignment.query_end = queryBases.size() - (temp +
// 1);
//      }
//      overlap->overlap->alignment.ref_begin += refStartPos;
//      overlap->overlap->alignment.ref_end += refStartPos;
//    }
//  });
//}
}

#endif /* SMITHWATERMAN_H_ */
