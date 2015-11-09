
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

#ifndef SLAM_H_
#define SLAM_H_
#include "KMerLookupTable.h"
#include "KMer.h"
#include "sequenceTools.h"
#include "FASTQsequence.h"
#include "GenbankTools.h"
#include "MetagenomicFASTQSequence.h"
#include "ParallelTools.h"
#include "TaxonomyDatabase.h"
#include "MetagenomicResults.h"
#include "Overlap.h"
#include "Globals.h"
#include "SmithWaterman.h"
//#include "comparison.h"
#include "PairedOverlap.h"
#include <inttypes.h>
#include <vector>
#include "SAM.h"
#include <fstream>
namespace SLAM {
void forceParallel();
inline void metagenomicAnalysis(const std::string R1FileName,
                                const std::string R2FileName,
                                const std::string databaseFileName,
                                const std::string outFileName,
                                const std::string samFileName,
                                const unsigned maxNumReads);
void metagenomicAnalysis_Low_Mem(const std::string R1FileName,
                                 const std::string R2FileName,
                                 const std::string databaseFile,
                                 const std::string outFileName,
                                 const unsigned readsPerGo,
                                 const unsigned maxNumReads);
template <typename FASTQType>
std::vector<Overlap> alignToDatabase(const std::vector<FASTQType> &reads,
                                     const GenbankIndex &genbankIndex);
//inline void serverMode(const std::string fileName, const unsigned maxNumReads);

template <typename FASTQType>
inline std::vector<Overlap> alignToDatabase(const std::vector<FASTQType> &reads,
                                            const GenbankIndex &genbankIndex) {
  log("Aligning reads to database using k = " + std::to_string(k));
  std::vector<KMerAndData<KMerInt, k>> kMers;
  getKMersFromReads(reads, kMers);
  genbankIndex.getKMers<KMerInt, k>(kMers, k / 2);
  sortKMers(kMers);
  std::vector<OverlapTemp> overlapTemp =
      findOverlaps_parallel(kMers.begin(), kMers.end(), reads);
  std::vector<Overlap> overlaps;
  overlaps.reserve(overlapTemp.size());
  for(auto & overlap : overlapTemp){
    overlaps.emplace_back(overlap);
  }
  kMers.resize(0);
  kMers.shrink_to_fit();
    performSmithWatermanOnRange_parallel(overlaps.begin(), overlaps.end(), reads,
                                       genbankIndex);
  return overlaps;
}

inline void metagenomicAnalysis(const std::string R1FileName,
                                const std::string R2FileName,
                                const std::string databaseFileName,
                                const std::string outFileName,
                                const std::string samFileName,
                                const unsigned maxNumReads) {
  forceParallel();
  log("Performing metagenomic analysis");
  if (samFileName.size()) reportCigar = true;
  bool isPaired = (R2FileName.size());
  pairedData=isPaired;
  std::vector<MetagenomicFASTQSequence> reads =
      isPaired
          ? getPairedReadsFromFiles<MetagenomicFASTQSequence>(
                R1FileName, R2FileName, maxNumReads)
          : getReadsFromFile<MetagenomicFASTQSequence>(R1FileName, maxNumReads);
  const size_t numReads = isPaired ? reads.size() / 2 : reads.size();
  if (reads.size() == 0) return;
  GenbankIndex genbankIndex =
      getIndexFromBoostSerial(databaseFileName + "/database");
  auto overlaps = alignToDatabase(reads, genbankIndex);
  screenOverlapsByScoreThreshold(overlaps, scoreThreshold);
  std::vector<ReadPairAndOverlaps> readPairsAndAlignments;
  if (isPaired) {
    auto pairedOverlaps =
        getPairedOverlaps(overlaps.begin(), overlaps.end(), reads);
    overlaps.resize(0);
    overlaps.shrink_to_fit();
    readPairsAndAlignments = getPerReadOverlaps(
        pairedOverlaps.begin(), pairedOverlaps.end(), reads.size() / 2);
    pairedOverlaps.resize(0);
    pairedOverlaps.shrink_to_fit();
    uint32_t maxInsertSize = getMaxAllowedInsertSize(readPairsAndAlignments);
    screenPairedAlignmentsByInsertSize(readPairsAndAlignments, maxInsertSize,
                                       true);
    screenPairedAlignmentsByScore(readPairsAndAlignments, scoreFractionThreshold);
  } else {
    auto perReadOverlaps = getPerReadOverlaps(overlaps.begin(), overlaps.end());
    readPairsAndAlignments =
        getDummyAlignmentPairsFromSingleEndReads(perReadOverlaps, reads);
    screenPairedAlignmentsByScore(readPairsAndAlignments, scoreFractionThreshold);
  }
  if (performPseudoAssembly) {
    pseudoAssembly(readPairsAndAlignments, reads, genbankIndex);
    screenPairedAlignmentsByScore(readPairsAndAlignments, scoreFractionThreshold);
  }
  if (samFileName.size()) {
    std::ofstream sam(samFileName);
    log("Writing SAM output");
    sam<<getHeader(genbankIndex);
    for (auto &read : readPairsAndAlignments)
      writeSAMOutputPairs(sam, read, reads, genbankIndex);
  }
  if(justAlign){
    log("Done");
    return;
  }
  TaxonomyDB taxDB(databaseFileName + "/taxDB");
  auto identifiedTaxonomies = convertAlignmentsToIdentifiedTaxonomies_parallel(
      readPairsAndAlignments.begin(), readPairsAndAlignments.end(), reads,
      genbankIndex, taxDB);
  std::ofstream perReadout(outFileName + "PerRead");
  writePerReadResults(identifiedTaxonomies, perReadout);
  identifiedTaxonomies = combineTaxonomies(identifiedTaxonomies);
  if (outFileName.size()) {
    std::ofstream outFile(outFileName);
    writeResults(identifiedTaxonomies, outFile, taxDB, numReads);
    //writePerReadResults(identifiedTaxonomies, outFileName + "_per_read");
    writeAbbreviatedResultsFile(identifiedTaxonomies,
                                outFileName + "_abbreviated", taxDB, numReads);
  } else
    writeResults(identifiedTaxonomies, std::cout, taxDB, numReads);
  log("Done");
  //  evaluateResults("metahitPerRead",
  //   "publication/data/reads/metahit/locusToTaxID",taxDB,"eval_32_16_meta");
  return;
}

inline void metagenomicAnalysis_Low_Mem(const std::string R1FileName,
                                        const std::string R2FileName,
                                        const std::string databaseFile,
                                        const std::string outFileName,
                                        const std::string samFileName,
                                        const unsigned readsPerGo,
                                        const unsigned maxNumReads) {
  unsigned numReads = 0;
  forceParallel();
  log("Performing metagenomic analysis");
  if (samFileName.size()) reportCigar = true;
  bool isPaired = (R2FileName.size());
  pairedData=isPaired;
  TaxonomyDB taxDB(databaseFile + "/taxDB");
  GenbankIndex genbankIndex =
      getIndexFromBoostSerial(databaseFile + "/database");
  std::vector<IdentifiedTaxonomy> identifiedTaxonomies;
  std::ifstream R1File(R1FileName);
  if (!R1File.good()) {
    log("FASTQ file " + R1FileName + " bad");
  }
  std::ifstream R2File;
  if (isPaired) {
    R2File.open(R2FileName);
    if (!R2File.good()) {
      log("FASTQ file " + R2FileName + " bad");
    }
  }
  std::ofstream sam;
  if (samFileName.size()) {
    sam.open(samFileName);
    sam<<getHeader(genbankIndex);
  }
  while (numReads < maxNumReads) {
    std::vector<MetagenomicFASTQSequence> reads;
    if (isPaired)
      log("Getting reads from FASTQ files " + R1FileName + " and " +
          R2FileName);
    else
      log("Getting reads from FASTQ file " + R1FileName);
    unsigned readsPerGoTemp = readsPerGo;
    if (numReads + readsPerGo > maxNumReads)
      readsPerGoTemp = maxNumReads - numReads;
    isPaired ? getPairedSequencesFromFASTQFiles(R1File, R2File, reads,
                                                readsPerGoTemp)
             : getSequencesFromFASTQFile(R1File, reads, readsPerGoTemp);
    if (reads.size() == 0) break;
    numReads += isPaired ? reads.size() / 2 : reads.size();
    auto overlaps = alignToDatabase(reads, genbankIndex);
    screenOverlapsByScoreThreshold(overlaps, scoreThreshold);
    std::vector<ReadPairAndOverlaps> readPairsAndAlignments;
    if (isPaired) {
      auto pairedOverlaps =
          getPairedOverlaps(overlaps.begin(), overlaps.end(), reads);
      overlaps.resize(0);
      overlaps.shrink_to_fit();
      readPairsAndAlignments = getPerReadOverlaps(
          pairedOverlaps.begin(), pairedOverlaps.end(), reads.size() / 2);
      pairedOverlaps.resize(0);
      pairedOverlaps.shrink_to_fit();
      uint32_t maxInsertSize = getMaxAllowedInsertSize(readPairsAndAlignments);
      screenPairedAlignmentsByInsertSize(readPairsAndAlignments, maxInsertSize,
                                             true);
      screenPairedAlignmentsByScore(readPairsAndAlignments, scoreFractionThreshold);
    } else {
      auto perReadOverlaps = getPerReadOverlaps(overlaps.begin(), overlaps.end());
      readPairsAndAlignments =
          getDummyAlignmentPairsFromSingleEndReads(perReadOverlaps, reads);
      screenPairedAlignmentsByScore(readPairsAndAlignments, scoreFractionThreshold);
    }
    if (performPseudoAssembly) {
      pseudoAssembly(readPairsAndAlignments, reads, genbankIndex);
      screenPairedAlignmentsByScore(readPairsAndAlignments, scoreFractionThreshold);
    }
    if (samFileName.size()) {
      log("Writing SAM output");
      for (auto &read : readPairsAndAlignments)
        writeSAMOutputPairs(sam, read, reads, genbankIndex);
    }
    if(justAlign){
      continue;
    }
    auto newIdentifiedTaxonomies =
        convertAlignmentsToIdentifiedTaxonomies_parallel(
            readPairsAndAlignments.begin(), readPairsAndAlignments.end(), reads,
            genbankIndex, taxDB);
    identifiedTaxonomies.insert(identifiedTaxonomies.end(),
                                newIdentifiedTaxonomies.begin(),
                                newIdentifiedTaxonomies.end());
    log("Processed\t" + std::to_string(numReads) + "\t reads");
  }
  std::ofstream perReadout(outFileName + "_PerRead");
  writePerReadResults(identifiedTaxonomies, perReadout);
  identifiedTaxonomies = combineTaxonomies(identifiedTaxonomies);
  if (outFileName.size()) {
    std::ofstream outFile(outFileName);
    writeResults(identifiedTaxonomies, outFile, taxDB, numReads);
    writeAbbreviatedResultsFile(identifiedTaxonomies,
                                outFileName + "_abbreviated", taxDB, numReads);
  } else
      writeResults(identifiedTaxonomies, std::cout, taxDB, numReads);
  log("Done");
  return;
}

// inline void serverMode(const std::string fileName, const unsigned
// maxNumReads) {
//  forceParallel();
//  std::vector<std::string> databaseFileNames{"PATRICGenomes"};
//  TaxonomyDB taxDB("taxonomyIndex");
//  GenbankIndex genbankIndex;
//  genbankIndex.buildIndex(databaseFileNames);
//  const KMerLookupTable<uint32_t, 16> table("humanHalfOverlap16");
//  std::ifstream pipe(fileName);
//  std::string line;
//  while (1) {
//    while (std::getline(pipe, line)) {
//      auto tokens = tokenise(line, "\t");
//      if (tokens.size() == 2) {
//        log("reset");
//        log("Performing metagenomic analysis");
//        std::vector<MetagenomicFASTQSequence> reads =
//            getReadsFromFile<MetagenomicFASTQSequence>(tokens[0].c_str(),
//                                                       maxNumReads);
//        const size_t numReads = reads.size();
//        if (reads.size() == 0) continue;
//        hostScreen(reads, table);
//        lowComplexityScreen(reads, 36, table);
//        auto alignments = alignToDatabase(reads, genbankIndex, taxDB);
//        screenAndSortOverlaps(alignments, 100);
//        auto identifiedTaxonomies =
//            convertAlignmentsToIdentifiedTaxonomies_parallel(
//                alignments.begin(), alignments.end(), reads, genbankIndex,
//                taxDB);
//        identifiedTaxonomies = concatenateTaxonomies(identifiedTaxonomies);
//        if(tokens[1].size()){
//          std::ofstream outFile(tokens[1]);
//          writeResults(identifiedTaxonomies, outFile, taxDB, numReads);
//        }
//        else
//          writeResults(identifiedTaxonomies, std::cout, taxDB, numReads);
//        //        auto results = convertAlignmentsToResults(
//        //            alignments.begin(), alignments.end(), reads,
// genbankIndex,
//        // taxDB);
//        //        writeResultsFile(results, tokens[1].c_str(), taxDB,
// numReads);
//      } else {
//        log("Quitting");
//        return;
//      }
//    }
//    if (pipe.eof()) {
//      sleep(1);
//      pipe.clear();
//    } else
//      return;
//  }
//
//  return;
//}

//  std::unordered_set<std::string> IDs;
//  std::ifstream idFile(databaseFileName);
//  std::string line;
//  while(getline(idFile,line)){
//    IDs.insert(line);
//  }
//  auto dotPos=R1FileName.find('.');
//  auto R1OutFileName=R1FileName;
//  R1OutFileName.insert(dotPos,"_filtered");
//  dotPos=R2FileName.find('.');
//  auto R2OutFileName=R2FileName;
//  R2OutFileName.insert(dotPos,"_filtered");
//  std::ofstream R1Out(R1OutFileName);
//  std::ofstream R2Out(R2OutFileName);
//  std::ifstream R1In(R1FileName);
//  std::ifstream R2In(R2FileName);
//  unsigned lineType = 0;
//  std::string seqID, bases, quality;
//  unsigned readsAdded = 0;
//  while (safeGetline(R1In, line)) {
//    switch (lineType) {
//      case 0:
//        seqID = std::move(line);
//        lineType++;
//        break;
//      case 1:
//        bases = std::move(line);
//        lineType++;
//        break;
//      case 2:
//        lineType++;
//        break;
//      case 3:{
//        quality = std::move(line);
//        FASTQSequence read(seqID, bases, quality);
//        if(IDs.find(read.sequenceIdentifier)==IDs.end()){
//          read.sequenceIdentifier=seqID.substr(1,std::string::npos);
//          R1Out<<read.getFASTQEntry();
//          readsAdded++;
//        }
//        lineType = 0;
//      }
//        break;
//      default:
//        break;
//    }
//  }
//  log("Wrote "+std::to_string(readsAdded) + " reads");
//  lineType = 0;
//  readsAdded=0;
//  while (safeGetline(R2In, line)) {
//    switch (lineType) {
//      case 0:
//        seqID = std::move(line);
//        lineType++;
//        break;
//      case 1:
//        bases = std::move(line);
//        lineType++;
//        break;
//      case 2:
//        lineType++;
//        break;
//      case 3:{
//        quality = std::move(line);
//        FASTQSequence read(seqID, bases, quality);
//        if(IDs.find(read.sequenceIdentifier)==IDs.end()){
//          read.sequenceIdentifier=seqID.substr(1,std::string::npos);
//          R2Out<<read.getFASTQEntry();
//          readsAdded++;
//        }
//      }
//        lineType = 0;
//        break;
//      default:
//        break;
//    }
//  }
//  log("Wrote "+std::to_string(readsAdded) + " reads");
//  return;
}

#endif /* SLAM_H_ */
