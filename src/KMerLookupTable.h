/* Copyright 2013 David Ainsworth
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
#ifndef KMERLOOKUPTABLE_H_
#define KMERLOOKUPTABLE_H_

#include "KMer.h"
#include "LookupTable.h"
#include "GenbankTools.h"
#include <atomic>
#include <limits>
/*
 * Implements a specialisation of the LookupTable
 * class for kMers
 * Contains functions for adding/lookup of kMers and
 * for the adding/lookup of strings which can be decomposed
 * into kMers
 *
 * This class is designed to be used to filter a set of FASTQ reads
 * by taxonomy ID. The lookup table would be initialised with a set of
 * genbank entries from a single taxID using the "addToTable" function.
 * This table would then be used to lookup each FASTQ read using the
 * "isFromTable"
 * function to test whether they map to any entries that are in the table.
 */
namespace SLAM {
template <typename KMerInt, unsigned K>
class KMerLookupTable : public LookupTable<KMerInt> {
 private:
  using LookupTable<KMerInt>::lookup;
  using LookupTable<KMerInt>::set;

 public:
  KMerLookupTable();
  KMerLookupTable(const char* fileName);
  bool lookup(const KMer<KMerInt, K>& kMer) const;
  void set(const KMer<KMerInt, K>& kMer);
  bool isHost(const std::string& bases, const unsigned gap,
              const unsigned cutoff, const unsigned secondarycutoff,
              const unsigned complexityCutoff) const;
  void addToTable(const std::string& bases, const unsigned gap);
  std::vector<int> getVectorOfTableQueryResults(const std::string& bases) const;
  std::vector<std::pair<unsigned, unsigned>> getChainStartPositionsAndLengths(
      std::vector<int>& foundPositions, const unsigned gap) const;
  bool areLongestTwoChainsValid(
      const std::vector<std::pair<unsigned, unsigned>>& startPosAndLengths,
      const unsigned cutoff) const;
  bool isLowComplexity(const std::string& bases,
                       const unsigned complexityCutoff) const;
};
template <typename KMerInt, unsigned K>
KMerLookupTable<KMerInt, K>::KMerLookupTable()
    : LookupTable<KMerInt>(getMask<KMerInt, K>()) {}

/*
 * See LookupTable.h for guidelines here
 */
template <typename KMerInt, unsigned K>
KMerLookupTable<KMerInt, K>::KMerLookupTable(const char* fileName)
    : LookupTable<KMerInt>(getMask<KMerInt, K>(), fileName) {}

/*
 * Inserts a kMer into the table
 */
template <typename KMerInt, unsigned K>
void KMerLookupTable<KMerInt, K>::set(const KMer<KMerInt, K>& kMer) {
  set(kMer.kMerInt);
}

/*
 * Lookup kMer in table
 */
template <typename KMerInt, unsigned K>
bool KMerLookupTable<KMerInt, K>::lookup(const KMer<KMerInt, K>& kMer) const {
  return lookup(kMer.kMerInt);
}

/*
 * Takes a string of bases and splits into K-mers.
 * Each kmer is looked up in the table and these lookup results are stored in a
 * boolean vector "foundPositions"
 */
template <typename KMerInt, unsigned K>
std::vector<int> KMerLookupTable<KMerInt, K>::getVectorOfTableQueryResults(
    const std::string& bases) const {
  KMerInt kMer = 0, rckMer = 0;
  std::vector<int> foundPositions(bases.size() - K + 1, 0);
  for (unsigned i = 0; i < bases.size(); i++) {
    addBaseToKMers<KMerInt, K>(bases[i], kMer, rckMer);
    if (i < (K - 1)) continue;
    bool found = lookup(kMer < rckMer ? kMer : rckMer);
    foundPositions[i - (K - 1)] = found;
  }
  return foundPositions;
}

/*
 * Goes through a vector of boolean lookup results and finds all chains longer
 * than 3 kmers.
 * Returns a vector of chain start positions and lengths
 */
template <typename KMerInt, unsigned K>
std::vector<std::pair<unsigned, unsigned>>
KMerLookupTable<KMerInt, K>::getChainStartPositionsAndLengths(
    std::vector<int>& foundPositions, const unsigned gap) const {
  std::vector<std::pair<unsigned, unsigned>> startPosAndLengths;
  for (unsigned i = 0; i < foundPositions.size(); i++) {
    if (foundPositions[i]) {
      unsigned chainLength = 0;
      for (unsigned j = i;; j += gap) {
        if (j >= foundPositions.size()) break;
        if (foundPositions[j]) {
          chainLength++;
          foundPositions[j] = 0;
        } else
          break;
      }
      if (chainLength > 2)
        startPosAndLengths.push_back({ i, chainLength
        });
    }
  }
  std::sort(startPosAndLengths.begin(), startPosAndLengths.end(),
            [](const std::pair<unsigned, unsigned> & i,
               const std::pair<unsigned, unsigned> & j) {
    return i.second > j.second;
  });
  return startPosAndLengths;
}

/*
 * Iff the longest two chains do not overlap and have a combined length
 * >=cutoff, true is returned
 */
template <typename KMerInt, unsigned K>
bool KMerLookupTable<KMerInt, K>::areLongestTwoChainsValid(
    const std::vector<std::pair<unsigned, unsigned>>& startPosAndLengths,
    const unsigned cutoff) const {
  unsigned startPos1 = startPosAndLengths[0].first;
  unsigned startPos2 = startPosAndLengths[1].first;
  unsigned endPos1 = startPos1 + startPosAndLengths[0].second * 8 + 7;
  unsigned endPos2 = startPos2 + startPosAndLengths[1].second * 8 + 7;
  if (startPosAndLengths[0].second + startPosAndLengths[1].second >= cutoff) {
    if (!(startPos1 <= endPos2 && startPos2 <= endPos1)) return true;
  }
  return false;
}

/*
 * Operates on an initialised table.
 * Takes a string of bases and uses a kMer based method to test whether
 * the bases map to any of the entries in the table.
 * The function looks up kMers of length K (stored in type KMerInt) which
 * are separated by "gap" bases
 * eg: ACTGACTACGACTG with K=4 and gap=3 would lookup these kMers:
 *     ACTG
 *        GACT
 *           TACG
 *              GACT
 *
 * The function tries to find chains of overlapping K-mers which have hits in
 * the lookup table.
 *
 * 1: The function first tests whether the string of bases is low complexity, it
 * does this by calculating the number of distinct 3mers in the string, if this
 * is below a certain cutoff then the read is not screened. This step is done as
 * low complexity reads are more likely to have kmer hits and be unnecessarily
 * screened.
 *
 * 2: At each position along the read, the kmer is looked up in the table and
 * the result stored in a vector.
 *
 * 3: This vector is analysed to find all chains of K-mers (where the found
 * kmers
 * are separated by "gap" bases).
 *
 * 4: If the longest chain is longer than "cutoff" kmers then the read is
 * assumed to be from the host organism.
 *
 * 5: Otherwise, if the longest two chains do not overlap and have a combined
 * length >= "secondaryCutoff" then they are screened. Otherwise they are
 * assumed to not be from the host
 *
 * Note that the parameters for this function must match those used in the
 * "addToTable"
 * function
 */
template <typename KMerInt, unsigned K>
bool KMerLookupTable<KMerInt, K>::isHost(
    const std::string& bases, const unsigned gap, const unsigned cutoff,
    const unsigned secondarycutoff, const unsigned complexityCutoff) const {
  if (isLowComplexity(bases, complexityCutoff)) return false;
  if (bases.size() < K) return false;
  std::vector<int> foundPositions = getVectorOfTableQueryResults(bases);
  std::vector<std::pair<unsigned, unsigned>> startPosAndLengths =
      getChainStartPositionsAndLengths(foundPositions, gap);
  if (startPosAndLengths.size() && (startPosAndLengths[0].second >= cutoff))
    return true;
  if (startPosAndLengths.size() < 2) return false;
  return areLongestTwoChainsValid(startPosAndLengths, secondarycutoff);
}

template <typename KMerInt, unsigned K>
bool KMerLookupTable<KMerInt, K>::isLowComplexity(
    const std::string& bases, const unsigned complexityCutoff) const {
  unsigned numDistinct3Mers = getNumDistinct3Mers(bases);
  return numDistinct3Mers < complexityCutoff;
}

/*
 * Adds a string to the table.
 * kMers are added, each separated by "gap" bases
 */
template <typename KMerInt, unsigned K>
void KMerLookupTable<KMerInt, K>::addToTable(const std::string& bases,
                                             const unsigned gap) {
  KMerInt kMer = 0, rckMer = 0;
  if (bases.size() < K) return;
  for (unsigned i = 0; i < bases.size(); i++) {
    addBaseToKMers<KMerInt, K>(bases[i], kMer, rckMer);
    if (i < (K - 1)) continue;
    if (((i - (K - 1)) % gap) == 0) {
      set(kMer < rckMer ? kMer : rckMer);
    }
  }
}

template <typename KMerInt, unsigned K>
void writeLookupTable(const GenbankIndex& index, const std::string fileName);
template <typename FASTQType>
void labelHostReads(std::vector<FASTQType>& reads,
                    const KMerLookupTable<uint32_t, 16>& table);
template <typename FASTQType>
void labelLowComplexityReads(std::vector<FASTQType>& reads,
                             const unsigned cutoff,
                             const KMerLookupTable<uint32_t, 16>& table);
/*
 * Tries to identify reads which may be from a species with tax id = "hostTaxID"
 * Labels them with isHost=true
 * Uses a linear model to determine parameters for the isHost function from the
 * length of the read
 */
template <typename FASTQType>
inline void labelHostReads(std::vector<FASTQType>& reads,
                           const KMerLookupTable<uint32_t, 16>& table) {
  log("Host screening using k = 16");
  std::atomic<unsigned> numHost(0);
  __gnu_parallel::for_each(reads.begin(), reads.end(), [&](FASTQType & read) {
    size_t size = read.bases.size();
    bool isHost =
        table.isHost(read.bases, 8, floor(0.09 * size + 1.1),
                     floor(0.082 * size + 3.04), floor(0.1 * size + 26));
    numHost += isHost;
    read.isHost = isHost;
  });
  auto endPos =
      std::remove_if(reads.begin(), reads.end(), [&](const FASTQType & i) {
    return i.isHost;
  });
  reads.resize(std::distance(reads.begin(), endPos));
  log("Screened " + std::to_string(numHost) + " reads, now got " +
      std::to_string(reads.size()));
  return;
}

template <typename FASTQType>
inline void labelLowComplexityReads(
    std::vector<FASTQType>& reads, const unsigned cutoff,
    const KMerLookupTable<uint32_t, 16>& table) {
  log("Low complexity screening using num3Mers = " + std::to_string(cutoff));
    //  std::atomic<unsigned> numLowComplexity(0);
  __gnu_parallel::for_each(reads.begin(), reads.end(), [&](FASTQType & read) {
    bool isLowComplexity = table.isLowComplexity(read.bases, cutoff);
    //    numLowComplexity += isLowComplexity;
    read.isLowComplexity = isLowComplexity;
  });
  //  auto endPos =
  //      std::remove_if(reads.begin(), reads.end(),
  //                     [&](const FASTQType &i) { return i.isLowComplexity; });
  //  reads.resize(std::distance(reads.begin(), endPos));
  //  log("Screened " + std::to_string(numLowComplexity) + " reads, now got " +
  //      std::to_string(reads.size()));
  return;
}
/*
 * Produces a lookup table from a genbank index and writes it to file
 */
template <typename KMerInt, unsigned K>
inline void writeLookupTable(const GenbankIndex& index,
                             const std::string fileName) {
  KMerLookupTable<KMerInt, K> table;
  for (auto& entry : index.entries) {
    table.addToTable(entry.bases, K / 2);
  }
  table.writeToFile(fileName);
}
}

#endif /* KMERLOOKUPTABLE_H_ */
