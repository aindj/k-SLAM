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
 * This file contains classes which represents a single entry in a fastq file,
 * with nucleotide sequence and associated read quality data
 */
#ifndef FASTQSEQUENCE_H_
#define FASTQSEQUENCE_H_

#include "sequenceTools.h"
#include <string>
#include <fstream>
#include <utility>
#include <vector>
namespace SLAM {
/*
 * Base class representing a FASTQ sequence
 * Like all FASTQ entries it contains sequence identifier and
 *  a string of bases, each with associated quality values.
 *
 */
class FASTQSequence {
 public:
  FASTQSequence(std::string sequenceIdentifier, std::string bases,
                std::string quality);
  FASTQSequence() {}
  ;
  void resize(const size_t length);
  std::string sequenceIdentifier;
  std::string bases;
  std::string quality;
  std::string getFASTQEntry() const;
  std::string getFASTAEntry() const;
};
void FASTQSequence::resize(const size_t length) {
  if (length > bases.size()) return;
  bases.resize(length);
  quality.resize(length);
}

/*
 * Initialises the FASTQSequence, stripping the '@' sign from the
 * seqID line of the input
 */
inline FASTQSequence::FASTQSequence(std::string seqID, std::string inBases,
                                    std::string inQuality)
    : bases(std::move(inBases)), quality(std::move(inQuality)) {
  if (seqID.size() > 1) {
    auto spacePos = seqID.find(' ');
    if (spacePos > 0) spacePos--;
    sequenceIdentifier = seqID.substr(1, spacePos);
    auto slashPos = sequenceIdentifier.find('/');
    sequenceIdentifier = sequenceIdentifier.substr(0, slashPos);
  }
}
// Tested
/*
 * Returns an entry in FASTQ format for this read
 */
inline std::string FASTQSequence::getFASTQEntry() const {
  return "@" + sequenceIdentifier + "\n" + bases + "\n+\n" + quality + "\n";
}
inline std::string FASTQSequence::getFASTAEntry() const {
  return ">" + sequenceIdentifier + "\n" + bases + "\n";
}
template <typename FASTQSequenceType>
void getSequencesFromFASTQFile(std::ifstream &input,
                               std::vector<FASTQSequenceType> &reads,
                               const unsigned numReads);
template <typename FASTQType>
void writeReadsToFASTQ(const std::string FASTQFileName,
                       const std::vector<FASTQType> &reads);
template <typename FASTQType>
std::vector<FASTQType> getReadsFromFile(const std::string FastQFileName,
                                        const unsigned numReads);
template <typename FASTQType>
double getAverageQuality(const std::vector<FASTQType> &reads);

// Tested
/*
 * Writes all entries in a vector of FASTQ reads to a file in
 * FASTQ format
 */
template <typename FASTQType>
inline void writeReadsToFASTQ(const std::string FASTQFileName,
                              const std::vector<FASTQType> &reads) {
  std::ofstream outFile(FASTQFileName);
  for (auto &read : reads) {
    outFile << read.getFASTQEntry();
  }
}

// Tested
template <typename FASTQSequenceType>
inline void getPairedSequencesFromFASTQFiles(
    std::ifstream &R1File, std::ifstream &R2File,
    std::vector<FASTQSequenceType> &reads, const unsigned maxNumReads) {
  getSequencesFromFASTQFile(R1File, reads, maxNumReads);
  const uint32_t numReads = reads.size();
  if (numReads == 0) return;
  getSequencesFromFASTQFile(R2File, reads, maxNumReads);
  if (reads.size() / numReads != 2) {
    std::cout << reads.size() / numReads << "\t" << reads.size() << "\t"
              << numReads << "\n";
    throw std::runtime_error("mismatch in R1 and R2 size");
  }
}
// Tested
/*
 * Reads an ifstream and extracts the FASTQ reads (with seq ID, sequence and
 * quality data)
 */
template <typename FASTQSequenceType>
inline void getSequencesFromFASTQFile(std::ifstream &input,
                                      std::vector<FASTQSequenceType> &reads,
                                      const unsigned numReads) {
  if (!input.good()) return;
  std::string line;
  // 0th line=SeqID, 1st line=Bases, 2nd line="+", 3rd line=Quality;
  unsigned lineType = 0;
  std::string seqID, bases, quality;
  unsigned readsAdded = 0;
  while (safeGetline(input, line)) {
    switch (lineType) {
      case 0:
        seqID = std::move(line);
        lineType++;
        break;
      case 1:
        bases = std::move(line);
        lineType++;
        break;
      case 2:
        lineType++;
        break;
      case 3:
        quality = std::move(line);
        reads.emplace_back(seqID, bases, quality);
        readsAdded++;
        lineType = 0;
        break;
      default:
        break;
    }
    if (readsAdded >= numReads) break;
  }
  log(std::to_string(reads.size()) + " reads");
  return;
}
// Tested
/*
 * Returns a vector of reads from a FASTQ file. Throws exception if file doesn't
 * exist
 */
template <typename FASTQType>
inline std::vector<FASTQType> getReadsFromFile(const std::string FastQFileName,
                                               const unsigned numReads) {
  std::vector<FASTQType> reads;
  std::ifstream readFile(FastQFileName);
  if (!readFile.good()) {
    log("FASTQ file " + FastQFileName + " bad");
    return std::vector<FASTQType>();
  }
  log("Getting reads from FASTQ file " + std::string(FastQFileName));
  getSequencesFromFASTQFile<FASTQType>(readFile, reads, numReads);
  return reads;
}
// Tested
template <typename FASTQType>
inline std::vector<FASTQType> getPairedReadsFromFiles(
    const std::string R1FileName, const std::string R2FileName,
    const unsigned numReads) {
  std::vector<FASTQType> reads;
  std::ifstream R1File(R1FileName);
  if (!R1File.good()) {
    log("FASTQ file " + R1FileName + " bad");
    return std::vector<FASTQType>();
  }
  std::ifstream R2File(R2FileName);
  if (!R2File.good()) {
    log("FASTQ file " + R2FileName + " bad");
    return std::vector<FASTQType>();
  }
  log("Getting reads from FASTQ files " + R1FileName + " and " + R2FileName);
  getPairedSequencesFromFASTQFiles<FASTQType>(R1File, R2File, reads, numReads);
  return reads;
}
template <typename FASTQType>
inline double getAverageQuality(const std::vector<FASTQType> &reads) {
  log("Finding average base quality");
  long double qualityTotal = 0;
  unsigned __int128 numBases = 0;
  for (auto &read : reads) {
    for (auto &base : read.quality)
      qualityTotal += (base - 33);
    numBases += read.quality.size();
  }
  return qualityTotal / numBases;
}
}
#endif /* FASTQSEQUENCE_H_ */
