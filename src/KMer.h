/* Copyright 2012 David Ainsworth
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

/*
 * Contains classes that represent compressed k-mers and tools
 * to work with these classes.
 *
 * Nucleotides are represented as follows: A=0b00, C=0b01, T=0b10, G=0b11
 * such that complement <-> flip of 2nd bit.
 * Bases are stored such that the final base of the k-mer is the LSB of the
 * integer
 * e.g: for k=3 with an 8 bit integer, TAG would be stored as 0b 00 10 00 11
 */

#ifndef KMER_H_
#define KMER_H_

#include <inttypes.h>
#include <string>
#include <iostream>
#include "sequenceTools.h"
#include <parallel/algorithm>
#include <thread>
namespace SLAM {
// Type used to store kMer offset
typedef uint32_t offsetType;
template <typename KMerInt, unsigned K> constexpr KMerInt getMask() {
  return std::numeric_limits<KMerInt>::max() >> (2 * (sizeof(KMerInt) * 4 - K));
}
/*
 * Class used to represent the metadata associated with kMers.
 * Offset: position of the k-mer in the sequence.
 * There are four pieces of metadata associated with a kMer:
 * readIDOrGenbankID: Represents the ID (position in array)
 * of the read or genbank entry that the k-mer came from
 * isFromGB: true if the k-mer came from a genbank entry, false if it came from
 * a read
 * revComp: true if the k-mer comes from the reverse-complement of the sequence.
 * Offset: the distance in bases from the start of the sequence (zero indexed)
 * The variable ID_isFromGB_RC stores readIDOrGenbankID in the least
 * significant 30 bits, isFromGB in the 31st bit and revComp in the 32nd
 */
class KMerData {
 public:
  uint32_t ID_isFromGB_RC;  // = 0;
  KMerData() {}
  ;
  KMerData(uint32_t readIDOrGenbankID, offsetType offset, bool isFromGenbank,
           bool revComp)
      : ID_isFromGB_RC(((readIDOrGenbankID &
                         0b00111111111111111111111111111111) |
                        (isFromGenbank << 31)) | (revComp << 30)),
        offset(offset) {}
  ;
  inline uint32_t getreadIDOrGenbankID() const {
    return ID_isFromGB_RC & 0b00111111111111111111111111111111;
  }
  ;
  inline bool getisFromGB() const {
    return ID_isFromGB_RC & 0b10000000000000000000000000000000;
  }
  ;
  inline bool getrevComp() const {
    return ID_isFromGB_RC & 0b01000000000000000000000000000000;
  }
  ;
  offsetType offset;  // = 0;
};

/*
 * Class used to define an arbitrary length k-mer, the k-mer itself
 * is stored in the KMerInt type
 */
template <typename KMerInt, unsigned K> class KMer {
 public:
  KMer(KMerInt kMer) : kMerInt(kMer) {}
  ;
  KMer() {}
  ;
  KMerInt kMerInt;  // = 0;
};

/*
 * Class used to represent a k-mer that has come from a FAST(A/Q) entry.
 * The kMer itself is stored in the KMerInt, the metadata stored in
 * a KMerData
 */
template <typename KMerInt, unsigned K>
class KMerAndData : public KMer<KMerInt, K> {
 public:
  KMerAndData(KMerInt kMerInt, uint32_t readIDOrGenbankID, offsetType offset,
              bool isFromGenbank, bool revComp)
      : KMer<KMerInt, K>(kMerInt),
        kMerData(readIDOrGenbankID, offset, isFromGenbank, revComp) {}
  ;
  KMerAndData(KMerInt kMerInt) : KMer<KMerInt, K>(kMerInt) {}
  ;
  KMerAndData() : KMer<KMerInt, K>() {}
  ;
  KMerData kMerData;
};

template <typename KMerInt, unsigned K>
void addBaseToKMers(const char base, KMerInt &kMer, KMerInt &rckMer);
template <typename KMerInt> inline KMerInt getTwoBits(const char base);
char convert(const char input);
template <typename KMerInt, unsigned K>
std::string decompress(const KMerAndData<KMerInt, K> &kMer);
template <typename KMerInt, unsigned K>
void printKMer(const KMerAndData<KMerInt, K> &kMer);
template <typename KMerInt, unsigned K>
void splitIntoKMersAndAddToVector(
    const std::string &bases, std::vector<KMerAndData<KMerInt, K> > &kMerVector,
    const bool isFromGenbank, const uint32_t ID, const unsigned gap,
    const size_t pos);
template <typename T, typename KMerInt, unsigned K>
void getKMers_parallel(const std::vector<T> &entries,
                       std::vector<KMerAndData<KMerInt, K> > &kMers,
                       bool isFromGenbank, unsigned gap);
unsigned getNumDistinct3Mers(const std::string &bases);
template <typename FASTQType, typename KMerType>
void getKMersFromReads(const std::vector<FASTQType> &reads,
                       std::vector<KMerType> &kMers);
template <typename KMerType> void sortKMers(std::vector<KMerType> &kMers);
/*
 * Splits a sequence into k-mers. The distance between each kMer is "gap" bases
 * For each k-mer, the k-mer and its reverse complement are generated
 * The smallest (numerically) of the k-mers is then added to the vector
 * This allows overlaps between reverse-complement strands to be found,
 * using only 1/2 of the memory required to store both the k-mer and
 * its revComp.
 * The variable "offset" represents the distance between the start of the
 * string and the start of the kMer.
 * In the case of the reverse complement kMer (rcKMer) being added,
 * for FASTQ reads (isFromGenbank==false), "offset"
 * refers to the offset of the rcKMer from the start of
 * reverseComplement(string).
 * For strings from Genbank (isFromGenbank==true), "offset"
 * refers to the offset of the kMer (not the rcKMer) from the start of the
 * string
 * Pos refers to the position in "kMerVector" that the kMers should be inserted
 * at, the vector should already be resized so that this is possible
 */

template <typename KMerInt, unsigned K>
inline void splitIntoKMersAndAddToVector(
    const std::string &bases, std::vector<KMerAndData<KMerInt, K> > &kMerVector,
    const bool isFromGenbank, const uint32_t ID, const unsigned gap,
    const size_t pos) {
  KMerInt kMer = 0, rckMer = 0;
  size_t offset = 0;
  if (bases.size() < K) return;
  // todo possibly a low complexity filter here
  for (unsigned i = 0; i < bases.size(); i++) {
    addBaseToKMers<KMerInt, K>(bases[i], kMer, rckMer);
    if (i < (K - 1)) continue;
    if (((i - (K - 1)) % gap) == 0) {
      kMer < rckMer ? kMerVector[pos + offset] =
          (KMerAndData<KMerInt, K>(kMer, ID, i - (K - 1), isFromGenbank, false))
: kMerVector[pos + offset] = (KMerAndData<KMerInt, K>(
              rckMer, ID, isFromGenbank ? i - (K - 1) : bases.size() - 1 - i,
              isFromGenbank, true));
      offset++;
    }
  }
}

/*
 * Applies "splitIntoKMersAndAddToVector" to a range of sequences in parallel
 * Load balancing is performed by calculating how many kmers will be produced
 * per entry and splitting entries amongst threads so that each thread has a
 * similar amount of kMers.
 * Each thread writes to a portion of a pre-resized vector
 */
template <typename T, typename KMerInt, unsigned K>
inline void getKMers_parallel(const std::vector<T> &entries,
                              std::vector<KMerAndData<KMerInt, K> > &kMers,
                              bool isFromGenbank, unsigned gap) {
  size_t oldSize = kMers.size();
  int numThreads = omp_get_max_threads();
  std::vector<size_t> numKMersInThread;
  std::vector<unsigned> startPositions { 0 }
  ;
  size_t numKMers = 0;
  for (int i = 0; i < entries.size(); i++) {
    auto &entry = entries[i];
    if (entry.bases.size() >= K) numKMers += (entry.bases.size() - K) / gap + 1;
  }
  size_t numKMersPerThread = numKMers / numThreads + 1;
  size_t runningNumKMers = 0;
  for (int i = 0; i < entries.size(); i++) {
    auto &entry = entries[i];
    if (entry.bases.size() >= K)
      runningNumKMers += (entry.bases.size() - K) / gap + 1;
    if (runningNumKMers > numKMersPerThread) {
      if (i + 1 < entries.size()) {
        startPositions.push_back(i + 1);
        numKMersInThread.push_back(runningNumKMers);
      }
      runningNumKMers = 0;
    }
  }
  kMers.resize(numKMers + oldSize);
  std::vector<std::thread *> threads;
  size_t kMerOffset = oldSize;
  for (int i = 0; i < startPositions.size(); i++) {
    unsigned startPos = startPositions[i];
    unsigned endPos =
        i == startPositions.size() - 1 ? entries.size() : startPositions[i + 1];
    threads.push_back(new std::thread([ =, &kMers, &entries]() {
      size_t kMerOffsetNew = kMerOffset;
      for (unsigned j = startPos; j < endPos; j++) {
        auto &entry = entries[j];
        if (entry.bases.size() < K) continue;

        splitIntoKMersAndAddToVector(entry.bases, kMers, isFromGenbank, j, gap,
                                     kMerOffsetNew);
        kMerOffsetNew += (entry.bases.size() - K) / gap + 1;
      }
    }));
    if (i < numKMersInThread.size()) kMerOffset += numKMersInThread[i];
  }
  for (int i = 0; i < threads.size(); i++) {
    threads[i]->join();
  }
}

/*
 * Converts an upper case base to its compressed representation
 */
template <typename KMerInt> inline KMerInt getTwoBits(const char base) {
  KMerInt twoBits;
  switch (base) {
    case 'A':
      twoBits = 0b00;
      break;
    case 'C':
      twoBits = 0b01;
      break;
    case 'T':
      twoBits = 0b10;
      break;
    case 'G':
      twoBits = 0b11;
      break;
    default:
      twoBits = 0;
      break;
  }
  return twoBits;
}

/*
 * Compresses a base and adds the it to a kmer and the base's RC to the
 * revComp kmer
 */
template <typename KMerInt, unsigned K>
inline void addBaseToKMers(const char base, KMerInt &kMer, KMerInt &rckMer) {
  kMer <<= 2;
  rckMer >>= 2;
  KMerInt twoBits = getTwoBits<KMerInt>(base);
  kMer |= twoBits;
  kMer &= getMask<KMerInt, K>();
  rckMer |= (twoBits ^ 0b10) << (2 * (K - 1));
}

/*
 * Iterates over a string to find the number of distinct 3mers
 * This is helpful when determining whether a read is low-complexity
 */
inline unsigned getNumDistinct3Mers(const std::string &bases) {
  std::array<char, 64> table;
  table.fill(0);
  uint32_t kMer = 0;
  uint32_t rckMer = 0;
  for (unsigned i = 0; i < bases.size(); i++) {
    addBaseToKMers<uint32_t, 3>(bases[i], kMer, rckMer);
    if (i < 2) continue;
    table[kMer]++;
  }
  unsigned numDistinct3Mers = 0;
  for (int i = 0; i < 64; i++) {
    if (table[i] > 0) numDistinct3Mers++;
  }
  return numDistinct3Mers;
}
template <typename KMerInt,unsigned K>
inline bool isLowComplexity(KMerInt kMerInt,const unsigned cutoff) {
//  KMerAndData<KMerInt, K> kMer;
//  kMer.kMerInt=kMerInt;
  std::array<char, 64> table;
  table.fill(0);
  for (int i = 0; i < K-2; i++) {
    table[kMerInt&0b111111]++;
    kMerInt>>=2;
  }
  unsigned numDistinct3Mers = 0;
  for (int i = 0; i < 64; i++) {
    if (table[i] > 0) numDistinct3Mers++;
  }
//  std::string temp=decompress(kMer);
//  if(getNumDistinct3Mers(temp)!=numDistinct3Mers)
//    std::cout<<"wrong"<<std::endl;
//  std::cout<<decompress(kMer)<<"\t"<<numDistinct3Mers<<std::endl;
  return numDistinct3Mers<cutoff;
}
/*
 *	Writes a decompressed kmer and metadata to stdout
 */
template <typename KMerInt, unsigned K>
inline void printKMer(const KMerAndData<KMerInt, K> &kMer) {
  std::cout << std::boolalpha << decompress<KMerInt, K>(kMer) << "\t"
            << kMer.kMerData.getreadIDOrGenbankID() << "\t"
            << kMer.kMerData.offset << "\t" << kMer.kMerData.getisFromGB()
            << "\t" << kMer.kMerData.getrevComp() << std::endl;
}

/*
 * Converts the compressed representation of a kmer to its string form
 */
template <typename KMerInt, unsigned K>
inline std::string decompress(const KMerAndData<KMerInt, K> &kMer) {
  std::string bigString;
  for (unsigned i = K; i > 0; i--) {
    bigString.push_back(convert((kMer.kMerInt >> 2 * (i - 1)) & 0b11));
  }
  return bigString;
}
template <typename KMerInt, unsigned K>
inline std::string decompress(const KMerInt &kMer) {
  std::string bigString;
  for (unsigned i = K; i > 0; i--) {
    bigString.push_back(convert((kMer >> 2 * (i - 1)) & 0b11));
  }
  return bigString;
}
/*
 * Converts the compressed representation of a base to its character
 * representation
 */
inline char convert(const char input) {
  switch (input) {
    case 0:
      return 'A';
    case 1:
      return 'C';
    case 2:
      return 'T';
    case 3:
      return 'G';
    default:
      return 'N';
  }
}
/*
 * Splits a vector of reads into kmers and appends to vector "kMers"
 */
template <typename FASTQType, typename KMerType>
inline void getKMersFromReads(const std::vector<FASTQType> &reads,
                              std::vector<KMerType> &kMers) {
  log("Getting k-mers from reads");
  size_t numKMers = kMers.size();
  getKMers_parallel(reads, kMers, false, 1);
  log("Obtained " + std::to_string(kMers.size() - numKMers) +
      " k-mers from reads");
}

/*
 * Sorts kMers based on value of kMerInt, this places identical kMers next
 * to each other.
 */

template <typename KMerType>
inline void sortKMers(std::vector<KMerType> &kMers) {
  log("Sorting k-mers");
  __gnu_parallel::sort(kMers.begin(), kMers.end(),
                       [](const KMerType & i, const KMerType & j) {
    if (i.kMerInt == j.kMerInt)
      return i.kMerData.ID_isFromGB_RC > j.kMerData.ID_isFromGB_RC;
    return i.kMerInt < j.kMerInt;
  });
  //  log("Done sorting k-mers");
}
//template <typename KMerType>
//inline void sortKMers(std::vector<KMerType> &kMers) {
//  log("Sorting k-mers");
//  __gnu_parallel::sort(kMers.begin(), kMers.end(),
//                       [](const KMerType & i, const KMerType & j) {
//    return i.kMerInt < j.kMerInt;
//  });
//  //  log("Done sorting k-mers");
//}
}
#endif /* KMER_H_ */
