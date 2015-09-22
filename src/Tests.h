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
#ifndef TESTS_H_
#define TESTS_H_

#include "SLAM.h"
#include "KMer.h"
#include "FASTQsequence.h"
#include "MetagenomicFASTQSequence.h"
#include "GenbankTools.h"
#include "sequenceTools.h"
#include "LookupTable.h"
//#include "Globals.h"
#include "KMerLookupTable.h"

#include <stdlib.h>
#include <iostream>
#include <inttypes.h>
#include <vector>
#include <algorithm>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <array>
#include <random>

namespace SLAM {
class connection {
 public:
  connection(unsigned readID, uint32_t GID, int offset, bool revComp)
      : readID(readID), GID(GID), offset(offset), revComp(revComp) {}
  connection() {}
  ;
  unsigned readID = 0;
  uint32_t GID = 0;
  int offset = 0;
  bool revComp = false;
  int expectedSWScore = 0;
  std::vector<std::string> geneOverlaps;
};
bool getReadsFromFile_Test(const std::string R1FileName,
                           const std::string R2FileName);
bool findOverlaps_Test();
bool performSmithWatermanOnRange_parallel_test();
bool test(const std::string R1FileName, const std::string R2FileName,
          const std::string databaseFileName, const std::string outFileName) {
  //  return findOverlaps_Test();
  return performSmithWatermanOnRange_parallel_test();
  return getReadsFromFile_Test(R1FileName, R2FileName);
}
inline bool getReadsFromFile_Test(const std::string R1FileName,
                                  const std::string R2FileName) {
  bool isPaired = (R2FileName.size());
  std::vector<MetagenomicFASTQSequence> reads =
      isPaired
          ? getPairedReadsFromFiles<MetagenomicFASTQSequence>(
                R1FileName, R2FileName, UINT32_MAX)
          : getReadsFromFile<MetagenomicFASTQSequence>(R1FileName, UINT32_MAX);
  std::vector<MetagenomicFASTQSequence> R1, R2;
  for (int i = 0; i < reads.size(); i++) {
    if (i < reads.size() / 2) {
      R1.push_back(reads[i]);
    } else
      R2.push_back(reads[i]);
  }
  writeReadsToFASTQ("temp_test1.fq", R1);
  writeReadsToFASTQ("temp_test2.fq", R2);
  bool R1pass = system(std::string(
      "diff " + R1FileName + " temp_test1.fq  > /dev/null 2>&1").c_str()) == 0;
  bool R2pass = system(std::string(
      "diff " + R2FileName + " temp_test2.fq  > /dev/null 2>&1").c_str()) == 0;
  return R1pass && R2pass;
}
std::string getRandomString() {
  static std::random_device rd;
  static std::mt19937 gen(rd());
  static std::uniform_int_distribution<> dis(32, 126);
  std::string returnString;
  for (int i = 0; i < 10; i++) {
    returnString.push_back(dis(gen));
  }
  return returnString;
}
GenbankEntry generateGenbankEntry(const uint32_t ID, const unsigned length) {
  std::string bases = "ACTG";
  GenbankEntry genbankEntry;
  genbankEntry.genbankID = ID;
  for (unsigned i = 0; i < length; i++) {
    genbankEntry.bases.push_back(bases[rand() % 4]);
  }
  for (int i = 0; i < 10; i++) {
    CDS cds;
    cds.start = rand() % (length - 160);
    cds.stop = cds.start + 150;
    Gene gene(getRandomString(), getRandomString(), getRandomString(),
              getRandomString(), std::to_string(ID), cds);
    genbankEntry.genes.push_back(gene);
  }
  genbankEntry.taxonomyID = ID;
  return genbankEntry;
}

std::pair<FASTQSequence, connection> getOverlappingRead(
    const GenbankEntry &entry, unsigned readLength, unsigned k) {
  std::string bases = "ACTG";
  //  std::pair<int, bool> offsetAndRevComp = { 0, false };

  FASTQSequence read;
  connection con;
  std::vector<Gene> genes;
  con.offset = rand() % (entry.bases.length() + readLength - 2 * k + 1) -
               (readLength - k);
  //  std::cout<<offsetAndRevComp.first<<"\n";
  if (con.offset >= 0) {
    read.bases = entry.bases.substr(con.offset, readLength);
    genes = entry.getGenesInRange(con.offset, con.offset + readLength);
  } else {
    read.bases = entry.bases.substr(0, readLength - abs(con.offset));
    genes = entry.getGenesInRange(0, readLength - abs(con.offset));
  }
  con.expectedSWScore = 2 * read.bases.size();
  if (read.bases.size() != readLength) {
    if (con.offset < 0) {
      std::string tempBases;
      while (read.bases.size() + tempBases.size() < readLength) {
        tempBases.push_back(bases[rand() % 4]);
      }
      read.bases = tempBases + read.bases;
    } else {
      while (read.bases.size() < readLength) {
        read.bases.push_back(bases[rand() % 4]);
      }
    }
  }
  if (rand() % 2) {
    read.bases = reverseComplement(read.bases);
    con.revComp = true;
  }

  for (auto &gene : genes)
    con.geneOverlaps.push_back(gene.geneName);
  //  read.sequenceIdentifier = std::to_string(ID);
  con.GID = entry.genbankID;
  return { read, con };
}
bool findOverlaps_Test() {
  std::vector<connection> connections;
  std::vector<connection> connections2;
  std::vector<connection> connections3;
  std::vector<Overlap> overlaps;
  GenbankIndex genbankIndex;
  std::vector<FASTQSequence> reads;
  std::vector<KMerAndData<uint64_t, 32> > kMers;
  for (int i = 0; i < 1000; i++) {
    GenbankEntry genbankEntry = generateGenbankEntry(i, 3000);
    genbankIndex.entries.push_back(genbankEntry);
    for (int j = 0; j < 100; j++) {
      std::pair<FASTQSequence, connection> readAndData =
          getOverlappingRead(genbankEntry, 150, 32);
      readAndData.second.readID = i * 100 + j;
      readAndData.second.GID = i;
      connections.push_back(readAndData.second);
      reads.push_back(readAndData.first);
    }
  }
  genbankIndex.getKMers(kMers, 1);
  getKMersFromReads(reads, kMers);
  sortKMers(kMers);
  overlaps = findOverlaps_parallel(kMers.begin(), kMers.end(), reads);

  for (auto &overlap : overlaps) {
    connections2.push_back(
        connection(overlap.readPosInArray, overlap.entryPosInArray,
                   overlap.relativePosition, overlap.revComp));
  }
    //  for (unsigned i = 0; i < reads.size(); i++) {
    //    for (auto &entry : reads[i].overlappingGenbankEntries) {
    //      connections2.push_back(connection(i,
    // entry.genbankEntry->genbankID,
    //                                        entry.relativePosition,
    // entry.revComp));
    //    }
    //  }
  std::sort(connections.begin(), connections.end(),
            [](connection i, connection j) {
    if (i.readID == j.readID) {
      if (i.GID == j.GID) {
        return i.offset < j.offset;
      } else {
        return i.GID < j.GID;
      }
    } else {
      return i.readID < j.readID;
    }
  });
  std::sort(connections2.begin(), connections2.end(),
            [](connection i, connection j) {
    if (i.readID == j.readID) {
      if (i.GID == j.GID) {
        return i.offset < j.offset;
      } else {
        return i.GID < j.GID;
      }
    } else {
      return i.readID < j.readID;
    }
  });
  connections3 = connections2;
  auto it =
      std::unique_copy(connections2.begin(), connections2.end(),
                       connections3.begin(), [](connection i, connection j) {
    return i.readID == j.readID;
  });
  connections3.resize(std::distance(connections3.begin(), it));
  std::cout << connections.size() << "\t" << connections2.size() << "\t"
            << connections3.size() << std::endl;
  for (unsigned i = 0; i < connections.size(); i++) {
    connection &i2 = connections[i];
    connection &j2 = connections3[i];
    if ((i2.GID == j2.GID) && (i2.offset == j2.offset) &&
        (i2.readID == j2.readID) && (i2.revComp == j2.revComp)) {

    } else {
      std::string bases = reads[connections[i].readID].bases;
      std::string rc = reverseComplement(bases);

      std::cout << "Read ID\t" << connections[i].readID << "\n";
      std::cout << "GID\t" << connections[i].GID << "\n";
      std::cout << "Offset\t" << connections[i].offset << "\n";
      std::cout << "RevComp\t" << connections[i].revComp << "\n";
      std::cout << (connections[i].revComp ? rc : bases) << "\n";
      std::cout << genbankIndex.entries[connections[i].GID].bases << "\n";

      bases = reads[connections3[i].readID].bases;
      rc = reverseComplement(bases);
      std::cout << "Read ID\t" << connections3[i].readID << "\n";
      std::cout << "GID\t" << connections3[i].GID << "\n";
      std::cout << "Offset\t" << connections3[i].offset << "\n";
      std::cout << "RevComp\t" << connections3[i].revComp << "\n";
      std::cout << (connections3[i].revComp ? rc : bases) << "\n";
      std::cout << genbankIndex.entries[connections3[i].GID].bases << "\n\n";
    }
  }
  return std::equal(connections.begin(), connections.end(),
                    connections3.begin(), [](connection i, connection j) {
    return (i.GID == j.GID) && (i.offset == j.offset) &&
           (i.readID == j.readID) && (i.revComp == j.revComp);
  });
}
bool performSmithWatermanOnRange_parallel_test() {
  std::vector<connection> connections;
  std::vector<Overlap> overlaps;
  GenbankIndex genbankIndex;
  std::vector<FASTQSequence> reads;
  std::vector<KMerAndData<uint64_t, 32> > kMers;
  for (int i = 0; i < 100; i++) {
    GenbankEntry genbankEntry = generateGenbankEntry(i, 300);
    genbankIndex.entries.push_back(genbankEntry);
    for (int j = 0; j < 100; j++) {
      std::pair<FASTQSequence, connection> readAndData =
          getOverlappingRead(genbankEntry, 150, 32);
      readAndData.second.readID = i * 100 + j;
      readAndData.second.GID = i;
      connections.push_back(readAndData.second);
      reads.push_back(readAndData.first);
    }
  }
  genbankIndex.getKMers(kMers, 1);
  getKMersFromReads(reads, kMers);
  sortKMers(kMers);
  overlaps = findOverlaps_parallel(kMers.begin(), kMers.end(), reads);

  auto readsAndOverlaps = getPerReadOverlaps(overlaps.begin(), overlaps.end());
  performSmithWatermanOnRange(readsAndOverlaps.begin(), readsAndOverlaps.end(),
                              reads, genbankIndex);
  std::vector<Overlap> overlapsAndSW;
  for (auto &readAndAlignments : readsAndOverlaps) {
    for (auto &overlapAndSW : readAndAlignments.overlaps)
      overlapsAndSW.push_back(overlapAndSW);
  }
  std::sort(connections.begin(), connections.end(),
            [](connection i, connection j) {
    if (i.readID == j.readID) {
      if (i.GID == j.GID) {
        return i.offset < j.offset;
      } else {
        return i.GID < j.GID;
      }
    } else {
      return i.readID < j.readID;
    }
  });
  std::sort(overlapsAndSW.begin(), overlapsAndSW.end(),
            [](Overlap i, Overlap j) {
    if (i.readPosInArray == j.readPosInArray) {
      if (i.entryPosInArray == j.entryPosInArray) {
        return i.relativePosition < j.relativePosition;
      } else {
        return i.entryPosInArray < j.entryPosInArray;
      }
    } else {
      return i.readPosInArray < j.readPosInArray;
    }
  });
  std::cout << overlapsAndSW.size() << "\t" << connections.size() << std::endl;
  for (unsigned i = 0; i < overlapsAndSW.size(); i++) {
    auto &overlap = overlapsAndSW[i];
    if (overlap.aasdent->sw_score != connections[i].expectedSWScore) {
      std::cout << (reads[overlap.readPosInArray].bases) << "\n"
                << genbankIndex.entries[overlap.entryPosInArray].bases << "\n"
                << overlap.aasdent->sw_score << "\t"
                << connections[i].expectedSWScore << std::endl;
      std::cout << overlap.relativePosition << std::endl;
      return false;
    }
  }
  return true;
}
bool kMerCompare(const std::string &i, const std::string &j) {
  std::string i2, j2;
  for (auto &c : i) {
    switch (c) {
      case 'A':
        i2.push_back('a');
        break;
      case 'C':
        i2.push_back('b');
        break;
      case 'T':
        i2.push_back('c');
        break;
      case 'G':
        i2.push_back('d');
        break;
    }
  }
  for (auto &c : j) {
    switch (c) {
      case 'A':
        j2.push_back('a');
        break;
      case 'C':
        j2.push_back('b');
        break;
      case 'T':
        j2.push_back('c');
        break;
      case 'G':
        j2.push_back('d');
        break;
    }
  }
  return std::lexicographical_compare(i2.begin(), i2.end(), j2.begin(),
                                      j2.end());
}
std::vector<std::string> getKMersFromString(const std::string bases, unsigned k,
                                            unsigned gap) {

  std::vector<std::string> kMers;
  for (int i = 0; (i + k) <= bases.size(); i += gap) {
    std::string kMer = bases.substr(i, k);
    std::string rcKMer = reverseComplement(kMer);
    kMers.push_back(kMerCompare(kMer, rcKMer) ? kMer : rcKMer);
  }
  return kMers;
}
bool splitIntoKMersAndAddToVector_Test() {
  GenbankIndex genbankIndex;
  for (int i = 0; i < 300; i++) {
    GenbankEntry genbankEntry = generateGenbankEntry(i, 300);
    genbankIndex.entries.push_back(genbankEntry);
  }
  for (int gap = 1; gap < 15; gap++) {
    //    std::cout << gap << std::endl;
    std::vector<KMerAndData<uint64_t, 30> > kMers;
    getKMers_parallel(genbankIndex.entries, kMers, false, gap);
    std::vector<std::string> kMers2;
    for (auto &entry : genbankIndex.entries) {
      std::vector<std::string> kMerstemp =
          getKMersFromString(entry.bases, 30, gap);
      kMers2.insert(kMers2.end(), kMerstemp.begin(), kMerstemp.end());
    }
    std::sort(kMers2.begin(), kMers2.end(),
              [](const std::string & i, const std::string & j) {
      return kMerCompare(i, j);
    });
    sortKMers(kMers);
    for (unsigned i = 0; i < kMers2.size(); i++) {
      if (kMers2[i] != decompress(kMers[i])) {
        return false;
      }
    }
  }
  for (int gap = 1; gap < 15; gap++) {
    //    std::cout << gap << std::endl;
    std::vector<KMerAndData<uint64_t, 16> > kMers;
    getKMers_parallel(genbankIndex.entries, kMers, false, gap);
    std::vector<std::string> kMers2;
    for (auto &entry : genbankIndex.entries) {
      std::vector<std::string> kMerstemp =
          getKMersFromString(entry.bases, 16, gap);
      kMers2.insert(kMers2.end(), kMerstemp.begin(), kMerstemp.end());
    }
    std::sort(kMers2.begin(), kMers2.end(),
              [](const std::string & i, const std::string & j) {
      return kMerCompare(i, j);
    });
    sortKMers(kMers);
    for (unsigned i = 0; i < kMers2.size(); i++) {
      if (kMers2[i] != decompress(kMers[i])) {
        return false;
      }
    }
  }
  for (int gap = 1; gap < 15; gap++) {
    //    std::cout << gap << std::endl;
    std::vector<KMerAndData<uint64_t, 32> > kMers;
    getKMers_parallel(genbankIndex.entries, kMers, false, gap);
    std::vector<std::string> kMers2;
    for (auto &entry : genbankIndex.entries) {
      std::vector<std::string> kMerstemp =
          getKMersFromString(entry.bases, 32, gap);
      kMers2.insert(kMers2.end(), kMerstemp.begin(), kMerstemp.end());
    }
    std::sort(kMers2.begin(), kMers2.end(),
              [](const std::string & i, const std::string & j) {
      return kMerCompare(i, j);
    });
    sortKMers(kMers);
    for (unsigned i = 0; i < kMers2.size(); i++) {
      if (kMers2[i] != decompress(kMers[i])) {
        return false;
      }
    }
  }
  return true;
}
}

// bool kMerCompare(const std::string &i, const std::string &j);
// std::vector<std::string> getKMersFromString(const std::string bases, unsigned
// k,
//                                            unsigned gap);
// bool splitIntoKMersAndAddToVector_Test();
// void test();
// bool getReadsFromFile_Test();
// GenbankEntry generateGenbankEntry(const uint32_t ID, const unsigned length);
// std::pair<FASTQSequence, connection> getOverlappingRead(
//    const GenbankEntry &entry, unsigned readLength, unsigned k,
//    unsigned ID = 0);
// bool findOverlaps_Test();
// bool lookupTable_Test();
// bool kMerLookup_Test();
//
//// void test() {
////  std::cout << "getReadsFromFile\t"
////            << (getReadsFromFile_Test() ? "pass" : "fail") << std::endl;
////  std::cout << "splitIntoKMersAndAddToVector\t"
////            << (splitIntoKMersAndAddToVector_Test() ? "pass" : "fail")
////            << std::endl;
////}
//
// bool lookupTable_Test() {
//  for (unsigned tableSize = 100000; tableSize < 10000000; tableSize += 100000)
// {
//    LookupTable<uint64_t> table(tableSize - 1);
//    std::unordered_set<unsigned> uniqueInts;
//    for (unsigned i = 0; i < 10000; i++)
//      uniqueInts.insert(rand() % (tableSize - 1));
//    for (auto number : uniqueInts) {
//      if (table.lookup(number)) {
//        std::cout << "Found\t" << number << std::endl;
//        return false;
//      } else {
//        table.set(number);
//        if (!table.lookup(number)) {
//          std::cout << "Not inserted\t" << number << std::endl;
//          return false;
//        }
//      }
//    }
//    for (auto &number : uniqueInts) {
//      if (!table.lookup(number)) {
//        std::cout << "Not found\t" << number << std::endl;
//        return false;
//      }
//    }
//    table.writeToFile("tempTable");
//    LookupTable<uint64_t> table2(tableSize - 1, "tempTable");
//    if (!(table == table2)) {
//      std::cout << "Tables not equal" << std::endl;
//      return false;
//    }
//  }
//  return true;
//}
//

//

//
// template <typename KMerInt, unsigned k>
// bool kMerLookup_Test_Sub(unsigned gap, unsigned readLength) {
//  std::vector<GenbankEntry> entries;
//  std::vector<std::pair<FASTQSequence, bool> > basesAndIsHost;
//  KMerLookupTable<KMerInt, k> table;
//  std::string bases = "ACTG";
//  bool good = true;
//  for (int i = 0; i < 10; i++) {
//    entries.push_back(generateGenbankEntry(i, 10000));
//    table.addToTable(entries[entries.size() - 1].bases, gap);
//  }
//  for (auto &entry : entries) {
//    for (int i = 0; i < 10; i++)
//      basesAndIsHost.push_back({
//        getOverlappingRead(entry, readLength, readLength).first, true
//      });
//  }
//  unsigned predictedMaxHits = (readLength - k) / gap + 1;
//  unsigned predictedMinHits = (readLength - (k - 1)) / gap;
//  for (auto &read : basesAndIsHost) {
//    int i = 1;
//    for (; i < predictedMaxHits + 10; i++) {
//      if (!table.isHost(read.first.bases, gap, i, 100, 0)) break;
//    }
//    if ((i - 1) < predictedMinHits) {
//      std::cout << "Got " << i - 1 << " hits, expected " << predictedMinHits
//                << "\n";
//      good = false;
//    }
//    if ((i - 1) > predictedMaxHits) {
//      std::cout << "Got " << i - 1 << " hits, expected " << predictedMaxHits
//                << "\n";
//      good = false;
//    }
//  }
//  return good;
//}
// bool kMerLookup_Test() {
//  for (int i = 1; i < 13; i++)
//    kMerLookup_Test_Sub<uint32_t, 13>(i, 103);
//  return true;
//}
//

//

//
// bool predictedScore_test() {
//  std::random_device rd;
//  std::mt19937 gen(rd());
//  std::normal_distribution<> dis(15, 3);
//  std::uniform_real_distribution<> real(0, 1);
//  int totalNumErrors = 0;
//  int totalNumErrorSquareds = 0;
//  double expectedNumErrors = 0;
//  double averageStdDev = 0;
//  const int numReads = 100000;
//  for (auto i = 0; i < numReads; i++) {
//    std::string quality;
//    int numErrors = 0;
//    for (int j = 0; j < 150; j++) {
//      int qual = round(dis(gen));
//      double prob = pow(10, qual * (-0.1));
//      bool hasError = real(gen) < prob;
//      numErrors += hasError;
//      quality.push_back(qual + 33);
//    }
//    totalNumErrors += numErrors;
//    totalNumErrorSquareds += numErrors * numErrors;
//    std::pair<double, double> expectedNumErrorsAndStandardDeviation =
//        getExpectedNumErrorsAndStandardDeviation(quality);
//    //    std::cout << expectedNumErrorsAndStandardDeviation.first << "\t"
//    //              << expectedNumErrorsAndStandardDeviation.second << "\t"
//    //              <<numErrors<<std::endl;
//    expectedNumErrors += expectedNumErrorsAndStandardDeviation.first;
//    averageStdDev += expectedNumErrorsAndStandardDeviation.second;
//  }
//  double meanNoErrors = totalNumErrors * 1.0 / numReads;
//  double actualVariance =
//      totalNumErrorSquareds * 1.0 / numReads - meanNoErrors * meanNoErrors;
//  std::cout << totalNumErrors * 1.0 / numReads << "\t" << sqrt(actualVariance)
//            << "\t" << expectedNumErrors / numReads << "\t"
//            << averageStdDev / numReads << "\n";
//}
// bool getTaxonomyFromReadAndOverlaps_test() {
//  std::random_device rd;
//  std::mt19937 gen(rd());
//  int quality = 20;
//  std::normal_distribution<> dis(quality, 3);
//  std::uniform_real_distribution<> real(0, 1);
//  const int numReads = 100000;
//  for (auto i = 0; i < numReads; i++) {
//    std::string quality;
//    //    int numErrors = 0;
//    for (int j = 0; j < 150; j++) {
//      int qual = round(dis(gen));
//      quality.push_back(qual + 33);
//    }
//    std::pair<double, double> expectedScoreAndStandardDeviation =
//        getExpectedScoreAndStandardDeviation(quality);
//    ReadAndAlignments readAndAlignments, goodReadAndAlignments,
//        goodReadAndAlignments2;
//    int numGoodHits = rand() % 5;
//    double cutoff = expectedScoreAndStandardDeviation.first -
//                    2 * expectedScoreAndStandardDeviation.second;
//    for (int j = 0; j < numGoodHits; j++) {
//      OverlapAndSmithWaterman overlap;
//      overlap.alignment.sw_score = (real(gen) + 1.01) * cutoff;
//      readAndAlignments.overlaps.push_back(overlap);
//      goodReadAndAlignments.overlaps.push_back(overlap);
//    }
//    int numBadHits = rand() % 5;
//    for (int j = 0; j < numBadHits; j++) {
//      OverlapAndSmithWaterman overlap;
//      overlap.alignment.sw_score = real(gen) * cutoff;
//      readAndAlignments.overlaps.push_back(overlap);
//    }
//    goodReadAndAlignments2 = getGoodAlignments(readAndAlignments, cutoff);
//    std::sort(goodReadAndAlignments.overlaps.begin(),
//              goodReadAndAlignments.overlaps.end(),
//              [](OverlapAndSmithWaterman i, OverlapAndSmithWaterman j) {
//      return i.alignment.sw_score < j.alignment.sw_score;
//    });
//    std::sort(goodReadAndAlignments2.overlaps.begin(),
//              goodReadAndAlignments2.overlaps.end(),
//              [](OverlapAndSmithWaterman i, OverlapAndSmithWaterman j) {
//      return i.alignment.sw_score < j.alignment.sw_score;
//    });
//    //    std::cout << goodReadAndAlignments.overlaps.size() << "\t"
//    //              << goodReadAndAlignments2.overlaps.size() << std::endl;
//    for (unsigned j = 0; j < goodReadAndAlignments2.overlaps.size(); j++) {
//      if (goodReadAndAlignments2.overlaps[j].alignment.sw_score !=
//          goodReadAndAlignments.overlaps[j].alignment.sw_score)
//        return false;
//    }
//  }
//  return true;
//}
// bool convertAlignmentsToIdentifiedTaxonomy_test() {
//  TaxonomyDB taxDB("taxonomyIndex");
//  std::vector<connection> connections;
//  std::vector<Overlap> overlaps;
//  GenbankIndex genbankIndex;
//  std::vector<FASTQSequence> reads;
//  std::vector<KMerAndData<uint64_t, 32> > kMers;
//  for (int i = 0; i < 100; i++) {
//    GenbankEntry genbankEntry = generateGenbankEntry(i, 3000);
//    genbankIndex.entries.push_back(genbankEntry);
//    for (int j = 0; j < 100; j++) {
//      unsigned readID = i * 100 + j;
//      std::pair<FASTQSequence, connection> readAndData =
//          getOverlappingRead(genbankEntry, 150, 32, readID);
//      readAndData.second.readID = readID;
//      readAndData.second.GID = i;
//      connections.push_back(readAndData.second);
//      reads.push_back(readAndData.first);
//    }
//  }
//  genbankIndex.getKMers(kMers, 1);
//  getKMersFromReads(reads, kMers);
//  sortKMers(kMers);
//  overlaps = findOverlaps_parallel(kMers.begin(), kMers.end(), reads);
//  auto readsAndOverlaps = getPerReadOverlaps(overlaps.begin(),
// overlaps.end());
//  auto perReadAlignments = performSmithWatermanOnRange_parallel(
//      readsAndOverlaps.begin(), readsAndOverlaps.end(), reads, genbankIndex);
//  auto identifiedTaxonomies = convertAlignmentsToResults(
//      perReadAlignments.begin(), perReadAlignments.end(), reads, genbankIndex,
//      taxDB);
//
//  std::sort(connections.begin(), connections.end(),
//            [](connection i, connection j) {
//    return i.readID < j.readID;
//  });
//  std::sort(identifiedTaxonomies.begin(), identifiedTaxonomies.end(),
//            [](IdentifiedTaxonomy i, IdentifiedTaxonomy j) {
//    return stoi(i.reads[0]) < stoi(j.reads[0]);
//  });
//  std::cout << connections.size() << "\t" << identifiedTaxonomies.size()
//            << std::endl;
//  for (unsigned i = 0; i < connections.size(); i++) {
//    if (connections[i].readID != stoi(identifiedTaxonomies[i].reads[0]))
//      return false;
//    std::string gid = identifiedTaxonomies[i].genes.size() > 0
//                          ? identifiedTaxonomies[i].genes[0].referenceSequence
//                          : std::string();
//    if ((gid.size()) && (connections[i].GID != stoi(gid))) return false;
//    if (connections[i].geneOverlaps.size() !=
//        identifiedTaxonomies[i].genes.size())
//      return false;
//    for (unsigned j = 0; j < connections[i].geneOverlaps.size(); j++) {
//      if (connections[i].geneOverlaps[j] !=
//          identifiedTaxonomies[i].genes[j].geneName)
//        return false;
//    }
//
//    //
//    //
// std::cout<<connections[i].readID<<"\t"<<identifiedTaxonomies[i].reads[0]<<std::endl;
//    //	  std::string gid=identifiedTaxonomies[i].genes.size()>0 ?
//    // identifiedTaxonomies[i].genes[0].referenceSequence : std::string();
//    //	  std::cout<<connections[i].GID<<"\t"<<gid<<std::endl;
//    //	  for(auto & gene :connections[i].geneOverlaps)
//    //	  	std::cout<<gene<<"\t";
//    //	  if(connections[i].geneOverlaps.size()) std::cout<<std::endl;
//    //	  for(auto & gene : identifiedTaxonomies[i].genes)
//    //		  std::cout<<gene.geneName<<"\t";
//    //	  if(identifiedTaxonomies[i].genes.size())std::cout<<std::endl;
//    //	  std::cout<<std::endl;
//  }
//  return true;
//}
// bool pairingTest(const std::vector<ReadAndAlignments> &alignments,
//                 const std::vector<ReadAndAlignments> &pairedAlignments,
//                 const std::vector<PairedMetagenomicFASTQSequence> &reads) {
//  std::unordered_map<uint32_t, ReadAndAlignments> alignmentsMap;
//  std::unordered_map<uint32_t, ReadAndAlignments> pairedAlignmentsMap;
//  uint32_t divPt = reads.size() / 2;
//  for (auto &alignment : alignments) {
//    alignmentsMap.insert({alignment.readPosInArray, alignment});
//  }
//  for (auto &alignment : pairedAlignments) {
//    pairedAlignmentsMap.insert({alignment.readPosInArray, alignment});
//  }
//  //  std::cout<<"Checking single against paired\n";
//  for (auto &alignment : alignmentsMap) {
//    uint32_t readID = alignment.second.readPosInArray;
//    //    std::cout<<"Checking\t"<<readID<<"\n";
//    bool isR1 = true;
//    if (readID >= divPt) {
//      readID -= divPt;
//      isR1 = false;
//    }
//    auto it = pairedAlignmentsMap.find(readID);
//    if (it == pairedAlignmentsMap.end()) {
//      std::cout << "Paired alignment not found\n";
//      return false;
//    } else {
//      if (isR1) {
//        if (alignment.second.readPosInArray != it->second.readPosInArray) {
//          std::cout << "R1/2 mismatch1\t" << alignment.second.readPosInArray
//                    << "\t" << it->second.readPosInArray << "\n";
//          return false;
//        }
//      } else {
//        if (alignment.second.readPosInArray != it->second.r2PosInArray) {
//          std::cout << "R1/2 mismatch2\t" << alignment.second.readPosInArray
//                    << "\t" << it->second.r2PosInArray << "\n";
//          return false;
//        }
//      }
//      auto &read = reads[alignment.second.readPosInArray];
//      if (read.isR1 != isR1) {
//        std::cout << "R1 not equal\n";
//        return false;
//      }
//      auto &overlaps = read.isR1 ? it->second.overlaps :
// it->second.R2Overlaps;
//      if (alignment.second.overlaps.size() == overlaps.size()) {
//        std::sort(alignment.second.overlaps.begin(),
//                  alignment.second.overlaps.end(),
//                  [](OverlapAndSmithWaterman i, OverlapAndSmithWaterman j) {
//          if (i.overlap.entryPosInArray == j.overlap.entryPosInArray)
//            return i.overlap.readPosInArray < j.overlap.readPosInArray;
//          else
//            return i.overlap.entryPosInArray < j.overlap.entryPosInArray;
//        });
//        std::sort(overlaps.begin(), overlaps.end(),
//                  [](OverlapAndSmithWaterman i, OverlapAndSmithWaterman j) {
//          if (i.overlap.entryPosInArray == j.overlap.entryPosInArray)
//            return i.overlap.readPosInArray < j.overlap.readPosInArray;
//          else
//            return i.overlap.entryPosInArray < j.overlap.entryPosInArray;
//        });
//        if (!std::equal(alignment.second.overlaps.begin(),
//                        alignment.second.overlaps.end(), overlaps.begin(),
//                        [](OverlapAndSmithWaterman i,
//                           OverlapAndSmithWaterman j) {
//               return i.overlap.entryPosInArray == j.overlap.entryPosInArray
// &&
//                      i.overlap.readPosInArray == j.overlap.readPosInArray;
//             })) {
//          std::cout << "Overlaps not equal\n";
//          return false;
//        }
//      } else {
//        std::cout << "Overlap size not equal\n";
//        return false;
//      }
//    }
//  }
//  for (auto &paired : pairedAlignmentsMap) {
//    //    uint32_t readID=paired.second.readID;
//    uint32_t R1pos = paired.second.readPosInArray;
//    uint32_t R2pos = paired.second.r2PosInArray;
//    auto &R1 = reads[R1pos];
//    auto &R2 = reads[R2pos];
//    auto R1it = alignmentsMap.find(R1pos);
//    auto R2it = alignmentsMap.find(R2pos);
//    if (R1it == alignmentsMap.end()) {
//      if (paired.second.overlaps.size() != 0) {
//        std::cout << "Couldn't find R1 overlaps\n";
//        return false;
//      }
//    } else {
//      if (paired.second.overlaps.size() == 0 && R1it->second.overlaps.size())
// {
//        std::cout << "Found R1 but dont have any\n";
//        return false;
//      }
//      if (paired.second.overlaps.size() == R1it->second.overlaps.size()) {
//        std::sort(paired.second.overlaps.begin(),
// paired.second.overlaps.end(),
//                  [](OverlapAndSmithWaterman i, OverlapAndSmithWaterman j) {
//          if (i.overlap.entryPosInArray == j.overlap.entryPosInArray)
//            return i.overlap.readPosInArray < j.overlap.readPosInArray;
//          else
//            return i.overlap.entryPosInArray < j.overlap.entryPosInArray;
//        });
//        std::sort(R1it->second.overlaps.begin(), R1it->second.overlaps.end(),
//                  [](OverlapAndSmithWaterman i, OverlapAndSmithWaterman j) {
//          if (i.overlap.entryPosInArray == j.overlap.entryPosInArray)
//            return i.overlap.readPosInArray < j.overlap.readPosInArray;
//          else
//            return i.overlap.entryPosInArray < j.overlap.entryPosInArray;
//        });
//        if (!std::equal(
//                 R1it->second.overlaps.begin(), R1it->second.overlaps.end(),
//                 paired.second.overlaps.begin(),
//                 [](OverlapAndSmithWaterman i, OverlapAndSmithWaterman j) {
//                   return i.overlap.entryPosInArray ==
//                              j.overlap.entryPosInArray &&
//                          i.overlap.readPosInArray ==
// j.overlap.readPosInArray;
//                 })) {
//          std::cout << "R1 Overlaps not equal\n";
//          return false;
//        }
//      }
//    }
//    if (R2it == alignmentsMap.end()) {
//      if (paired.second.R2Overlaps.size() != 0) {
//        std::cout << "Couldn't find R2 overlaps\n";
//        return false;
//      }
//    } else {
//      if (paired.second.R2Overlaps.size() == 0 &&
//          R2it->second.overlaps.size()) {
//        std::cout << "Found R2 but dont have any\n";
//        return false;
//      }
//      if (paired.second.R2Overlaps.size() == R2it->second.overlaps.size()) {
//        std::sort(paired.second.R2Overlaps.begin(),
//                  paired.second.R2Overlaps.end(),
//                  [](OverlapAndSmithWaterman i, OverlapAndSmithWaterman j) {
//          if (i.overlap.entryPosInArray == j.overlap.entryPosInArray)
//            return i.overlap.readPosInArray < j.overlap.readPosInArray;
//          else
//            return i.overlap.entryPosInArray < j.overlap.entryPosInArray;
//        });
//        std::sort(R2it->second.overlaps.begin(), R2it->second.overlaps.end(),
//                  [](OverlapAndSmithWaterman i, OverlapAndSmithWaterman j) {
//          if (i.overlap.entryPosInArray == j.overlap.entryPosInArray)
//            return i.overlap.readPosInArray < j.overlap.readPosInArray;
//          else
//            return i.overlap.entryPosInArray < j.overlap.entryPosInArray;
//        });
//        if (!std::equal(
//                 R2it->second.overlaps.begin(), R2it->second.overlaps.end(),
//                 paired.second.R2Overlaps.begin(),
//                 [](OverlapAndSmithWaterman i, OverlapAndSmithWaterman j) {
//                   return i.overlap.entryPosInArray ==
//                              j.overlap.entryPosInArray &&
//                          i.overlap.readPosInArray ==
// j.overlap.readPosInArray;
//                 })) {
//          std::cout << "R2 Overlaps not equal\n";
//          return false;
//        }
//      }
//    }
//  }
//  return true;
//}
// bool getTopScoringEntriesFromPairs_test(
//    ReadAndOverlaps &singleAlignments,
//    std::vector<PairedOverlap> &pairedAlignments) {
//  ReadAndOverlaps convertedBack;
//  for (auto &el : pairedAlignments) {
//    if (el.hasR1) convertedBack.overlaps.push_back(el.r1Overlap);
//    if (el.hasR2) convertedBack.R2Overlaps.push_back(el.r2Overlap);
//  }
//  std::vector<decltype(convertedBack.overlaps) *> ovvec{
//      &convertedBack.overlaps,    &convertedBack.R2Overlaps,
//      &singleAlignments.overlaps, &singleAlignments.R2Overlaps};
//  for (auto singleReadPtr : ovvec) {
//    auto &singleReadOverlaps = *singleReadPtr;
//    std::sort(
//        singleReadOverlaps.begin(), singleReadOverlaps.end(),
//        [](const OverlapAndSmithWaterman &i, const OverlapAndSmithWaterman &j)
// {
//          if (i.overlap.entryPosInArray == j.overlap.entryPosInArray)
//            return i.overlap.overlapInfo.relativePosition <
//                   j.overlap.overlapInfo.relativePosition;
//          else
//            return i.overlap.entryPosInArray < j.overlap.entryPosInArray;
//        });
//  }
//  auto it = std::unique_copy(
//      convertedBack.overlaps.begin(), convertedBack.overlaps.end(),
//      convertedBack.overlaps.begin(),
//      [](const OverlapAndSmithWaterman &i, const OverlapAndSmithWaterman &j) {
//        return i.overlap.entryPosInArray == j.overlap.entryPosInArray &&
//               i.overlap.overlapInfo.relativePosition ==
//                   j.overlap.overlapInfo.relativePosition;
//      });
//  convertedBack.overlaps.resize(
//      std::distance(convertedBack.overlaps.begin(), it));
//  it = std::unique_copy(
//      convertedBack.R2Overlaps.begin(), convertedBack.R2Overlaps.end(),
//      convertedBack.R2Overlaps.begin(),
//      [](const OverlapAndSmithWaterman &i, const OverlapAndSmithWaterman &j) {
//        return i.overlap.entryPosInArray == j.overlap.entryPosInArray &&
//               i.overlap.overlapInfo.relativePosition ==
//                   j.overlap.overlapInfo.relativePosition;
//      });
//  convertedBack.R2Overlaps.resize(
//      std::distance(convertedBack.R2Overlaps.begin(), it));
//  if (convertedBack.overlaps.size() != singleAlignments.overlaps.size()) {
//    std::cout << "r1 size\n";
//    return false;
//  }
//  if (convertedBack.R2Overlaps.size() != singleAlignments.R2Overlaps.size()) {
//    std::cout << "r2 size\n";
//    return false;
//  }
//  return (std::equal(convertedBack.overlaps.begin(),
//                     convertedBack.overlaps.end(),
//                     singleAlignments.overlaps.begin(),
//                     [](OverlapAndSmithWaterman i, OverlapAndSmithWaterman j)
// {
//           return (i.overlap.entryPosInArray == j.overlap.entryPosInArray) &&
//                  (i.overlap.readPosInArray == j.overlap.readPosInArray) &&
//                  (i.overlap.overlapInfo.relativePosition ==
//                   j.overlap.overlapInfo.relativePosition) &&
//                  (i.alignment.sw_score == j.alignment.sw_score);
//         })) &&
//         (std::equal(convertedBack.R2Overlaps.begin(),
//                     convertedBack.R2Overlaps.end(),
//                     singleAlignments.R2Overlaps.begin(),
//                     [](OverlapAndSmithWaterman i, OverlapAndSmithWaterman j)
// {
//           return (i.overlap.entryPosInArray == j.overlap.entryPosInArray) &&
//                  (i.overlap.readPosInArray == j.overlap.readPosInArray) &&
//                  (i.overlap.overlapInfo.relativePosition ==
//                   j.overlap.overlapInfo.relativePosition) &&
//                  (i.alignment.sw_score == j.alignment.sw_score);
//         }));
//}
//}
#endif /* TESTS_H_ */
