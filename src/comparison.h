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
#ifndef COMPARISON_H_
#define COMPARISON_H_
#include "GenbankTools.h"
#include "TaxonomyDatabase.h"
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
namespace SLAM {
class mcnemar {
 public:
  bool krakenCorrectSp = false;
  bool krakenCorrectGe = false;
  bool phymmBLCorrectSp = false;
  bool phymmBLCorrectGe = false;
  bool SLAMCorrectSp = false;
  bool SLAMCorrectGe = false;
  bool NBCCorrectSp = false;
  bool NBCCorrectGe = false;
  bool RITACorrectSp = false;
  bool RITACorrectGe = false;
};
void mcnemar_eval(std::string locusToTaxIDFileName, const TaxonomyDB taxDB) {
  std::unordered_map<std::string, uint32_t> locusToTaxID;
  std::unordered_map<std::string, mcnemar> seqIDToResult;
  std::string line;
  std::ifstream inFile2("locusToTaxID");
  if (!inFile2.good()) throw std::runtime_error("unable to open map file");
  while (safeGetline(inFile2, line)) {
    std::vector<std::string> tokens = tokenise(line, "\t");
    if (tokens.size() < 2) break;
    uint32_t taxID = std::stoi(tokens[1]);
    locusToTaxID.insert({
      tokens[0], taxID
    });
  }
  std::vector<std::string> fileNames {
    "kraken", "phymmBL", "SLAM", "NBC", "RITA"
  }
  ;
  for (int i = 0; i < fileNames.size(); i++) {
    std::ifstream inFile(fileNames[i]);
    if (!inFile.good()) throw std::runtime_error("unable to open results file");
    while (safeGetline(inFile, line)) {
      std::vector<std::string> tokens = tokenise(line, "\t");
      if (tokens.size() < 2) continue;
      std::string seqID = tokens[0];
      if (seqID.size()) {
        if (seqID[seqID.size() - 1] == ' ') seqID.resize(seqID.size() - 1);
      }
      uint32_t predictedTaxID = std::stoi(tokens[1]);
      auto dotPos = seqID.find_first_of(". \t");
      std::string locus = seqID.substr(0, dotPos);
      uint32_t actualTaxID = 0;
      auto it = locusToTaxID.find(locus);
      if (it == locusToTaxID.end()) {
        std::cout << "Not found\t" << locus << "\n";
      } else {
        actualTaxID = it->second;
      }

      uint32_t actualSpeciesTaxID =
          taxDB.getTaxIDAtRank(actualTaxID, "species");
      uint32_t actualGenusTaxID = taxDB.getTaxIDAtRank(actualTaxID, "genus");
      uint32_t predictedSpeciesTaxID =
          taxDB.getTaxIDAtRank(predictedTaxID, "species");
      uint32_t predictedGenusTaxID =
          taxDB.getTaxIDAtRank(predictedTaxID, "genus");

      auto it2 = seqIDToResult.find(seqID);
      if (it2 == seqIDToResult.end()) {
        mcnemar temp;
        it2 = seqIDToResult.insert({
          seqID, temp
        }).first;
      }
      bool* genus = nullptr;
      bool* species = nullptr;
      switch (i) {
        case 0:
          genus = &it2->second.krakenCorrectGe;
          species = &it2->second.krakenCorrectSp;
          break;
        case 1:
          genus = &it2->second.phymmBLCorrectGe;
          species = &it2->second.phymmBLCorrectSp;
          break;
        case 2:
          genus = &it2->second.SLAMCorrectGe;
          species = &it2->second.SLAMCorrectSp;
          break;
        case 3:
          genus = &it2->second.NBCCorrectGe;
          species = &it2->second.NBCCorrectSp;
          break;
        case 4:
          genus = &it2->second.RITACorrectGe;
          species = &it2->second.RITACorrectSp;
          break;
        default:
          std::cout << "i error\t" << i << "\n";
          break;
      }
      bool isSubspecies = false;
      if (predictedSpeciesTaxID) {
        if (predictedSpeciesTaxID == actualSpeciesTaxID) *species = true;
      }
      if (predictedGenusTaxID) {
        if (predictedGenusTaxID == actualGenusTaxID) *genus = true;
      }
    }
  }
  //  std::cout<<"Acidobacterium_capsulatum_ATCC_51196_uid59127.000000000\tK_Sp\tK_Ge\tP_Sp\tP_ge\tS_Sp\tS_ge\tN_Sp\tN_Ge\tR_Sp\tR_Ge\n";
  double n01 = 0;
  double n10 = 0;
  int count = 0;
  for (auto& entry : seqIDToResult) {
    if (entry.second.SLAMCorrectSp && !entry.second.phymmBLCorrectSp)
      n01++;
    else if (!entry.second.SLAMCorrectSp && entry.second.phymmBLCorrectSp)
      n10++;
    //    std::cout<<entry.first<<"\t"<<entry.second.krakenCorrectSp<<"\t"<<entry.second.krakenCorrectGe<<"\t"<<
    //        entry.second.phymmBLCorrectSp<<"\t"<<entry.second.phymmBLCorrectGe<<"\t"<<
    //        entry.second.SLAMCorrectSp<<"\t"<<entry.second.SLAMCorrectGe<<"\t"<<
    //        entry.second.NBCCorrectSp<<"\t"<<entry.second.NBCCorrectGe<<"\t"<<
    //        entry.second.RITACorrectSp<<"\t"<<entry.second.RITACorrectGe<<"\n";
    //    count++;
    //    if(count>1000)
    //      break;
  }
  double chi = pow((abs(n01 - n10) - 1), 2) / (n01 + n10);
  std::cout << n01 << "\t" << n10 << "\t" << chi << "\n";
}
template <typename FASTQType>
void evaluateOverlaps(std::vector<ReadPairAndOverlaps>& readPairsAndAlignments,
                      const GenbankIndex& genbankIndex,
                      const std::vector<FASTQType>& reads,
                      std::string locusToTaxIDFileName,
                      const TaxonomyDB taxDB) {
  class Result {
   public:
    uint32_t taxID = 0;
    uint32_t speciesTaxID = 0;
    uint32_t genusTaxID = 0;
    uint32_t numSubSpecies = 0;
    uint32_t numSpecies = 0;
    uint32_t numGenus = 0;
  };
  std::unordered_map<uint32_t, Result> taxIDAndResult;
  std::unordered_map<std::string, uint32_t> locusToTaxID;
  unsigned numIncorrect = 0;
  std::string line;
  std::ifstream inFile2(locusToTaxIDFileName);
  if (!inFile2.good()) throw std::runtime_error("unable to open map file");
  while (safeGetline(inFile2, line)) {
    std::vector<std::string> tokens = tokenise(line, "\t");
    if (tokens.size() < 2) break;
    uint32_t taxID = std::stoi(tokens[1]);
    locusToTaxID.insert({
      tokens[0], taxID
    });
  }
  for (auto& readPair : readPairsAndAlignments) {
    auto& read = reads[readPair.r1PosInArray];
    const std::string& seqID = read.sequenceIdentifier;
    auto dotPos = seqID.find_first_of(". \t");
    std::string locus = seqID.substr(0, dotPos);
    uint32_t actualTaxID = 0;
    auto it = locusToTaxID.find(locus);
    if (it == locusToTaxID.end()) {
      std::cout << "Not found\t" << locus << "\n";
    } else {
      actualTaxID = it->second;
    }
    bool containsCorrectSpecies = false;
    std::sort(readPair.alignmentPairs.begin(), readPair.alignmentPairs.end(),
              [](const PairedOverlap & i, const PairedOverlap & j) {
      return i.combinedScore > j.combinedScore;
    });
    uint32_t actualSpeciesTaxID = taxDB.getTaxIDAtRank(actualTaxID, "species");
    uint32_t actualGenusTaxID = taxDB.getTaxIDAtRank(actualTaxID, "genus");
    auto it2 = taxIDAndResult.find(actualTaxID);
    if (it2 == taxIDAndResult.end()) {
      Result temp;
      temp.taxID = actualTaxID;
      temp.speciesTaxID = taxDB.getTaxIDAtRank(actualTaxID, "species");
      temp.genusTaxID = taxDB.getTaxIDAtRank(actualTaxID, "genus");
      it2 = taxIDAndResult.insert({
        temp.taxID, temp
      }).first;
    }
    bool species = false;
    bool genus = false;
    uint32_t topScore = 0;
    for (auto& overlap : readPair.alignmentPairs) {
      if (overlap.combinedScore >= topScore) topScore = overlap.combinedScore;
      if (taxDB.getTaxIDAtRank(
              genbankIndex.entries[overlap.entryPosInArray].taxonomyID,
              "species") == actualSpeciesTaxID) {
        if (overlap.combinedScore == topScore) containsCorrectSpecies = true;
        species = true;
        genus = true;
        break;
      } else if (taxDB.getTaxIDAtRank(
                     genbankIndex.entries[overlap.entryPosInArray].taxonomyID,
                     "genus") == actualGenusTaxID) {
        genus = true;
      }
    }
    if (species) {
      std::cout << readPair.r1PosInArray << "\t";
      if (!containsCorrectSpecies)
        std::cout << "Incorrect\n";
      else
        std::cout << "Correct \n";
      std::cout << "Correct species\t" << actualSpeciesTaxID << "\n";
      for (auto& overlap : readPair.alignmentPairs) {
        std::cout << genbankIndex.entries[overlap.entryPosInArray].taxonomyID
                  << "\t" << overlap.combinedScore << "\n";
      }
    }
    if (species)
      it2->second.numSpecies++;
    else
      numIncorrect++;
    if (genus) it2->second.numGenus++;
  }
  std::cout << "Incorrect\t" << numIncorrect << "\n";
  for (auto& result : taxIDAndResult) {
    std::cout << taxDB.getScientificName(result.second.taxID) << "\t"
              << result.second.numSpecies << "\t" << result.second.numGenus
              << "\n";
  }
}
void evaluateResults(std::string fileName, std::string locusToTaxIDFileName,
                     const TaxonomyDB taxDB, const std::string filename) {
  std::ofstream out(filename);
  class Result {
   public:
    uint32_t taxID = 0;
    uint32_t speciesTaxID = 0;
    uint32_t genusTaxID = 0;
    uint32_t numSubSpecies = 0;
    uint32_t numSpecies = 0;
    uint32_t numGenus = 0;
    //    uint32_t numIncorrect=0;
  };
  std::unordered_map<uint32_t, Result> taxIDAndResult;
  std::unordered_map<std::string, uint32_t> locusToTaxID;
  unsigned numIncorrect = 0;
  //  std::unordered_set<uint32_t> taxIDs;
  //  GenbankIndex index = getIndexFromBoostSerial("NCBIBacGenomes_Serialized");
  //  for(auto & entry : index.entries){
  //    taxIDs.insert(entry.taxonomyID);
  //  }
  std::string line;
  std::ifstream inFile2(locusToTaxIDFileName);
  if (!inFile2.good()) throw std::runtime_error("unable to open map file");
  while (safeGetline(inFile2, line)) {
    std::vector<std::string> tokens = tokenise(line, "\t");
    if (tokens.size() < 2) break;
    uint32_t taxID = std::stoi(tokens[1]);
      //    std::cout<<tokens[0]<<"\t"<<taxID<<"\n";
    locusToTaxID.insert({
      tokens[0], taxID
    });
  }

  std::ifstream inFile(fileName);
  if (!inFile.good()) throw std::runtime_error("unable to open results file");
  while (safeGetline(inFile, line)) {
    //    std::cout<<line<<"\n";

    std::vector<std::string> tokens = tokenise(line, "\t");
    if (tokens.size() < 2) continue;
    std::string seqID = tokens[0];
    uint32_t predictedTaxID = std::stoi(tokens[1]);
    //    std::cout<<"predicted
    // taxid\t"<<predictedTaxID<<"\t"<<taxDB.getScientificName(predictedTaxID)<<"\n";
    auto dotPos = seqID.find_first_of(". \t");
    std::string locus = seqID.substr(0, dotPos);
    uint32_t actualTaxID = 0;
    auto it = locusToTaxID.find(locus);
    if (it == locusToTaxID.end()) {
      std::cout << "Not found\t" << locus << "\n";
    } else {
      actualTaxID = it->second;
      //      std::cout<<"actual
      // taxid\t"<<actualTaxID<<"\t"<<taxDB.getScientificName(actualTaxID)<<"\n";
    }

    uint32_t actualSpeciesTaxID = taxDB.getTaxIDAtRank(actualTaxID, "species");
    uint32_t actualGenusTaxID = taxDB.getTaxIDAtRank(actualTaxID, "genus");
    uint32_t predictedSpeciesTaxID =
        taxDB.getTaxIDAtRank(predictedTaxID, "species");
    uint32_t predictedGenusTaxID =
        taxDB.getTaxIDAtRank(predictedTaxID, "genus");

    auto it2 = taxIDAndResult.find(actualTaxID);
    if (it2 == taxIDAndResult.end()) {
      //      std::cout<<"inserting\n";
      Result temp;
      temp.taxID = actualTaxID;
      temp.speciesTaxID = taxDB.getTaxIDAtRank(actualTaxID, "species");
      temp.genusTaxID = taxDB.getTaxIDAtRank(actualTaxID, "genus");
      it2 = taxIDAndResult.insert({
        temp.taxID, temp
      }).first;
    }
    bool isSubspecies = false;
    if (predictedSpeciesTaxID) {
      //        if(predictedSpeciesTaxID!=predictedTaxID){
      //          if(predictedTaxID==actualTaxID){
      //            it2->second.numSubSpecies++;
      //          }
      //          else
      //            continue;
      //        }
      if (predictedSpeciesTaxID == actualSpeciesTaxID)
        it2->second.numSpecies++;
      else
        numIncorrect++;
    }
    if (predictedGenusTaxID) {
      if (predictedGenusTaxID == actualGenusTaxID) it2->second.numGenus++;
      //        else
      //          continue;
    }
  }
  std::cout << "Incorrect\t" << numIncorrect << "\n";
  for (auto& result : taxIDAndResult) {
    //    std::cout<<std::boolalpha<<(taxIDs.find(result.second.taxID)!=taxIDs.end())<<"\t";
    out << taxDB.getScientificName(result.second.taxID) << "\t" <<
        //        result.second.taxID<<"\t"<<
        //        result.second.numSubSpecies<<"\t"<<
        result.second.numSpecies << "\t" << result.second.numGenus << "\n";
  }
}
}

#endif /* COMPARISON_H_ */
