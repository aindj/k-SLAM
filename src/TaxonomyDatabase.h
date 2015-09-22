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
#ifndef TAXONOMYDATABASE_H_
#define TAXONOMYDATABASE_H_
#include <unordered_map>
#include "sequenceTools.h"
#include <string>
namespace SLAM {
class TaxonomyEntry {
 public:
  uint32_t taxonomyID = 0;
  uint32_t parentTaxonomyID = 0;
  std::string scientificName;
  std::string rank;
  inline bool operator==(const TaxonomyEntry& other) const {
    return this->taxonomyID == other.taxonomyID &&
           this->parentTaxonomyID == other.parentTaxonomyID &&
           this->scientificName == other.scientificName;
  }
  TaxonomyEntry* parent = nullptr;
  std::vector<TaxonomyEntry*> children;
  unsigned numReadsAligned = 0;
  unsigned numReadsAlignedToChildren = 0;
  bool used = false;
  uint64_t genomeSize = 0;
  uint64_t genomeSizeOfChildren = 0;
  uint64_t numBelow = 0;
};

class TaxonomyDB {
 public:
  TaxonomyDB(const std::string inFileName);
  TaxonomyDB() {}
  ;
  std::unordered_map<uint32_t, TaxonomyEntry> taxIDsAndEntries;
  void parseNamesDump(const std::string namesDumpFileName);
  void parseNodesDump(const std::string nodesDumpFileName);
  uint32_t getTaxIDAtRank(const uint32_t taxID, const std::string& rank) const;
  std::string getScientificName(const uint32_t taxID) const;
  std::string getRank(const uint32_t taxID) const;
  uint32_t getLowestCommonAncestor(const std::vector<uint32_t>& taxIDs) const;
  uint32_t getParentTaxID(const uint32_t taxID) const;
  std::string getLineage(uint32_t taxonomyID) const;
  std::string getMetaPhlAnLineage(uint32_t taxonomyID) const;
  char* getIndexFileName(const uint32_t hostTaxID) const;
  void readTaxonomyIndex(const std::string inFileName);
  void writeTaxonomyIndex(const std::string outFileName,
                          const std::string namesDumpFileName,
                          const std::string nodesDumpFileName);
  bool isSubSpecies(uint32_t taxonomyID) const;
  int isBelowInTree(uint32_t upper, uint32_t lower) const;
  void createPointers();
  //  bool isArtificial(uint32_t taxID);
  //  bool isEndogenousRetrovirus(uint32_t taxID);
  //  bool isBetweenSpeciesAndGenus(uint32_t taxID);
  //  bool isPlasmid(uint32_t taxID);
  //  std::string getSpecies(uint32_t taxonomyID);
  //  uint32_t getSpeciesTaxID(uint32_t taxonomyID);
  //  bool isVirus(uint32_t taxonomyID);
  //  void testNearestNeighbour();
  //  void print();
};
void TaxonomyDB::createPointers() {
  for (auto& tax : taxIDsAndEntries) {
    auto parentIt = taxIDsAndEntries.find(tax.second.parentTaxonomyID);
    if (parentIt != taxIDsAndEntries.end()) {
      tax.second.parent = &(parentIt->second);
      parentIt->second.children.push_back(&tax.second);
    }
  }
}
TaxonomyDB::TaxonomyDB(const std::string inFileName) {
  log("Building taxonomy index");
  readTaxonomyIndex(inFileName);
  createPointers();
  log("Built a taxonomy tree with " + std::to_string(taxIDsAndEntries.size()) +
      " nodes");
}

void TaxonomyDB::parseNodesDump(const std::string nodesDumpFileName) {
  std::ifstream nodesDumpFile(nodesDumpFileName);
  if (!nodesDumpFile.is_open())
    throw std::runtime_error("unable to open nodes file");
  std::string line;
  while (nodesDumpFile.good()) {
    getline(nodesDumpFile, line);
    std::vector<std::string> tokens = tokenise(line, "\t|");
    if (tokens.size() > 2) {
      TaxonomyEntry newEntry;
      newEntry.taxonomyID = stoi(tokens[0]);
      newEntry.parentTaxonomyID = stoi(tokens[1]);
      newEntry.rank = tokens[2];
      auto entryIt = taxIDsAndEntries.insert({
        newEntry.taxonomyID, newEntry
      });
      if (!entryIt.second) {
        entryIt.first->second.taxonomyID = newEntry.taxonomyID;
        newEntry.parentTaxonomyID = stoi(tokens[1]);
      }
    }
  }
}
void TaxonomyDB::parseNamesDump(const std::string namesDumpFileName) {
  std::ifstream namesDumpFile(namesDumpFileName);
  if (!namesDumpFile.is_open())
    throw std::runtime_error("unable to open names file");
  std::string line;
  while (namesDumpFile.good()) {
    getline(namesDumpFile, line);
    std::vector<std::string> tokens = tokenise(line, "|");
    for (auto& token : tokens) {
      if (token.size() > 1) {
        if (token[0] == '\t') token.erase(0, 1);
        if (token[token.size() - 1] == '\t') token.erase(token.size() - 1, 1);
      }
    }
    if (tokens.size() > 3) {
      TaxonomyEntry newEntry;
      newEntry.taxonomyID = stoi(tokens[0]);
      //			for(auto & token : tokens)
      //				std::cout<<token<<"\n";
      if (tokens[3] == "scientific name") {
        //		std::cout<<"Found\n";
        newEntry.scientificName = tokens[1];
        //		std::cout<<newEntry.scientificName<<"\n";
      } else
        continue;
      auto entryIt = taxIDsAndEntries.insert({
        newEntry.taxonomyID, newEntry
      });
      if (!entryIt.second) {
        entryIt.first->second.scientificName = newEntry.scientificName;
      }
    }
  }
}

void TaxonomyDB::writeTaxonomyIndex(const std::string outFileName,
                                    const std::string namesDumpFileName,
                                    const std::string nodesDumpFileName) {
  parseNodesDump(nodesDumpFileName);
  parseNamesDump(namesDumpFileName);
  std::ofstream outFile(outFileName);
  for (auto& entry : taxIDsAndEntries) {
    outFile << entry.first << "\n" << entry.second.parentTaxonomyID << "\n"
            << entry.second.scientificName << "\n" << entry.second.rank << "\n";
  }
  outFile.close();
}

void TaxonomyDB::readTaxonomyIndex(const std::string inFileName) {
  std::ifstream inFile(inFileName);
  if (!inFile.is_open())
    throw std::runtime_error("unable to open taxonomy index file");
  for (std::string line; getline(inFile, line);) {
    TaxonomyEntry newEntry;
    newEntry.taxonomyID = stoi(line);
    getline(inFile, line);
    newEntry.parentTaxonomyID = stoi(line);
    getline(inFile, line);
    newEntry.scientificName = line;
    getline(inFile, line);
    newEntry.rank = line;
    taxIDsAndEntries.insert({
      newEntry.taxonomyID, newEntry
    });
  }
}

uint32_t TaxonomyDB::getLowestCommonAncestor(
    const std::vector<uint32_t>& taxIDs) const {
  if (taxIDs.size() == 0) {
    return 0;
  }
  std::vector<std::vector<uint32_t> > paths;
  for (auto& taxID : taxIDs) {
    bool good = true;
    std::vector<uint32_t> path;
    uint32_t tempTaxID = taxID;
    while (tempTaxID != 0) {
      path.push_back(tempTaxID);
      tempTaxID = getParentTaxID(tempTaxID);
    }
    if (good) paths.push_back(path);
  }
  if (paths.size() == 0) {
    return 0;
  }
  for (auto& path : paths)
    std::reverse(path.begin(), path.end());
  std::sort(paths.begin(), paths.end(),
            [](std::vector<uint32_t> i, std::vector<uint32_t> j) {
    return i.size() < j.size();
  });
  uint32_t consensus = 0;
  for (unsigned i = 0; i < paths[0].size(); i++) {
    uint32_t temp = 0;
    for (auto& path : paths) {
      if (temp == 0)
        temp = path[i];
      else if (temp != path[i]) {
        return consensus;
      }
    }
    consensus = temp;
  }
  return consensus;
}

uint32_t TaxonomyDB::getParentTaxID(const uint32_t taxID) const {
  auto entry = taxIDsAndEntries.find(taxID);
  if (entry != taxIDsAndEntries.end() && entry->second.parentTaxonomyID != 1)
    return entry->second.parentTaxonomyID;
  else
    return 0;
}

std::string TaxonomyDB::getScientificName(const uint32_t taxID) const {
  auto entry = taxIDsAndEntries.find(taxID);
  if (entry != taxIDsAndEntries.end()) {
    return entry->second.scientificName;
  } else
    return std::string();
}

std::string TaxonomyDB::getRank(const uint32_t taxID) const {
  auto entry = taxIDsAndEntries.find(taxID);
  if (entry != taxIDsAndEntries.end()) {
    return entry->second.rank;
  } else
    return std::string();
}

std::string TaxonomyDB::getLineage(uint32_t taxonomyID) const {
  std::string lineage;
  while (true) {
    // 131567 = Cellular organisms
    if (taxonomyID != 131567) {
      if (lineage.size()) lineage.insert(0, "; ");
      lineage.insert(0, getScientificName(taxonomyID));
      if (getRank(taxonomyID) == "species") lineage.clear();
    }
    taxonomyID = getParentTaxID(taxonomyID);
    if (taxonomyID == 0) {
      if (lineage.size()) lineage.append(".");
      break;
    }
  }
  return lineage;
}
std::string TaxonomyDB::getMetaPhlAnLineage(uint32_t taxonomyID) const {
  std::string rank = getRank(taxonomyID);
  if (rank == "superphylum") return std::string();
  std::string lineage;
  while (true) {
    // 131567 = Cellular organisms
    if (taxonomyID != 131567) {
      std::string rank = getRank(taxonomyID);
      if (rank == "species") {
        lineage.insert(0, "|s__");
        lineage.insert(4, getScientificName(taxonomyID));
      } else if (rank == "genus") {
        lineage.insert(0, "|g__");
        lineage.insert(4, getScientificName(taxonomyID));
      } else if (rank == "family") {
        lineage.insert(0, "|f__");
        lineage.insert(4, getScientificName(taxonomyID));
      } else if (rank == "order") {
        lineage.insert(0, "|o__");
        lineage.insert(4, getScientificName(taxonomyID));
      } else if (rank == "class") {
        lineage.insert(0, "|c__");
        lineage.insert(4, getScientificName(taxonomyID));
      } else if (rank == "phylum") {
        lineage.insert(0, "|p__");
        lineage.insert(4, getScientificName(taxonomyID));
      } else if (rank == "superkingdom") {
        lineage.insert(0, "k__");
        lineage.insert(3, getScientificName(taxonomyID));
      }
    }
    taxonomyID = getParentTaxID(taxonomyID);
    if (taxonomyID == 0) {
      break;
    }
  }
  std::replace(lineage.begin(), lineage.end(), ' ', '_');
  return lineage;
}

uint32_t TaxonomyDB::getTaxIDAtRank(const uint32_t taxID,
                                    const std::string& rank) const {
  auto entry = taxIDsAndEntries.find(taxID);
  while (entry != taxIDsAndEntries.end() &&
         entry->second.parentTaxonomyID != 1) {
    if (entry->second.rank == rank) {
      return entry->second.taxonomyID;
    } else
      entry = taxIDsAndEntries.find(entry->second.parentTaxonomyID);
  }
  return 0;
}
int TaxonomyDB::isBelowInTree(uint32_t upper, uint32_t lower) const {
  auto entry = taxIDsAndEntries.find(lower);
  unsigned level = 0;
  while (entry != taxIDsAndEntries.end() &&
         entry->second.parentTaxonomyID != 1) {
    if (entry->first == upper) {
      return level;
    } else {
      entry = taxIDsAndEntries.find(entry->second.parentTaxonomyID);
      level++;
    }
  }
  return -1;
}
bool TaxonomyDB::isSubSpecies(uint32_t taxonomyID) const {
  bool isSubSpecies = false;
  auto entry = taxIDsAndEntries.find(taxonomyID);
  int numLevels = 0;
  while (entry != taxIDsAndEntries.end() &&
         entry->second.parentTaxonomyID != 1) {
    if (entry->second.rank == "species") {
      if (numLevels > 0) {
        isSubSpecies = true;
      }
      break;
    } else
      entry = taxIDsAndEntries.find(entry->second.parentTaxonomyID);
    numLevels++;
  }
  return isSubSpecies;
}

//
//
// uint32_t TaxonomyDB::getParentTaxID(const uint32_t taxID) const {
//  auto entry = taxIDsAndEntries.find(taxID);
//  if (entry != taxIDsAndEntries.end() && entry->second.parentTaxonomyID != 1)
//    return entry->second.parentTaxonomyID;
//  else
//    return 0;
//}
//
//
// char* TaxonomyDB::getIndexFileName(const uint32_t hostTaxID) const {
//  uint32_t orderTaxID = getTaxIDAtRank(hostTaxID, "order");
//  //is primate?
//  if (orderTaxID == 9443) return "primateindex";
//  //is rodent?
//  if (orderTaxID == 9989) return "rodentindex.gz";
//  uint32_t classTaxID = getTaxIDAtRank(hostTaxID, "class");
//  //is mammal?
//  if (classTaxID == 40674) return "mammalindex.gz";
//  //is bacterium?
//  uint32_t superkingdomTaxID = getTaxIDAtRank(hostTaxID, "superkingdom");
//  if (superkingdomTaxID == 2) return "bacindex.gz";
//
//  //todo is plant?
//  //todo is vertebrate?
//  //todo is invertebrate?
//  return nullptr;
//}
//
//
//
// void TaxonomyDB::print() {
//  for (auto& taxonomyEntryIt : taxIDsAndEntries) {
//    auto& taxonomyEntry = taxonomyEntryIt.second;
//    if (taxonomyEntry.rank == "species") {
//      uint32_t genusTaxID = getTaxIDAtRank(taxonomyEntry.taxonomyID, "genus");
//      std::cout << getSpecies(taxonomyEntry.taxonomyID) << "\t"
//                << taxonomyEntry.scientificName << "\t"
//                << taxonomyEntry.taxonomyID << "\tGenus:\t"
//                << getScientificName(genusTaxID) << "\t" << genusTaxID <<
// "\n";
//    }
//    continue;
//    while (true) {
//      uint32_t parentTaxID = taxonomyEntry.parentTaxonomyID;
//      std::cout << taxonomyEntry.rank << "\t" << taxonomyEntry.scientificName
//                << "; ";  //"\t"<<parentTaxID<<"\t";
//      auto taxIt = taxIDsAndEntries.find(parentTaxID);
//      if (taxIt != taxIDsAndEntries.end() && parentTaxID != 1) {
//        taxonomyEntry = taxIt->second;
//      } else {
//        std::cout << "\n";
//        break;
//      }
//    }
//
//    //std::cout<<taxonomyEntryIt.second.taxonomyID<<"\t"<<taxonomyEntryIt.second.parentTaxonomyID<<"\t"<<taxonomyEntryIt.second.scientificName<<"\n";
//  }
//}
// bool TaxonomyDB::isBetweenSpeciesAndGenus(uint32_t taxID) {
//  bool aboveSpecies = true;
//  bool belowGenus = false;
//  auto entry = taxIDsAndEntries.find(taxID);
//  bool isFirst = true;
//  while (entry != taxIDsAndEntries.end() &&
//         entry->second.parentTaxonomyID != 1) {
//    if (entry->second.rank == "species") {
//      aboveSpecies = false;
//    } else if (entry->second.rank == "genus") {
//      if (!isFirst) belowGenus = true;
//      break;
//    }
//    entry = taxIDsAndEntries.find(entry->second.parentTaxonomyID);
//    isFirst = false;
//  }
//  return (aboveSpecies && belowGenus);
//}
// bool testEquality(std::vector<uint32_t>& inVector) {
//  if (inVector.size()) {
//    uint32_t first = inVector[0];
//    for (auto& element : inVector) {
//      if (element != first) return false;
//    }
//    return true;
//  }
//  return false;
//}
// void TaxonomyDB::testNearestNeighbour() {
//  std::vector<uint32_t> taxIDs;
//  for (auto& entry : taxIDsAndEntries) {
//    taxIDs.push_back(entry.first);
//  }
//  for (int i = 0; i < 10000; i++) {
//    std::vector<uint32_t> randomTaxIDs;
//    unsigned int numIDs = 4;  //2;//rand()%3+1;
//    for (unsigned int j = 0; j < numIDs; j++) {
//      uint32_t randomTaxID = taxIDs[rand() % (taxIDs.size())];
//      randomTaxIDs.push_back(randomTaxID);
//      //
// std::cout<<randomTaxID<<"\t"<<getScientificName(randomTaxID)<<"\n";
//    }
//    uint32_t nearestNeighbour = findNearestNeighbour(randomTaxIDs);
//    if ((nearestNeighbour != 0) && (nearestNeighbour != 2759) &&
//        (nearestNeighbour != 131567)) {
//      std::cout << "Random IDs\n";
//      for (auto& randomID : randomTaxIDs) {
//        std::cout << randomID << "\t" << getScientificName(randomID) << "\n";
//      }
//      std::cout << "Nearest neighbour\t" << nearestNeighbour << "\t"
//                << getScientificName(nearestNeighbour) << "\n";
//    }
//  }
//}
// bool TaxonomyDB::isVirus(uint32_t taxID) {
//  auto entry = taxIDsAndEntries.find(taxID);
//  while (entry != taxIDsAndEntries.end()) {
//    if (entry->second.rank == "superkingdom") {
//      return (entry->second.taxonomyID == 10239);
//    } else if (entry->second.parentTaxonomyID != 1)
//      entry = taxIDsAndEntries.find(entry->second.parentTaxonomyID);
//    else
//      return false;
//  }
//  return false;
//}
// bool TaxonomyDB::isPlasmid(uint32_t taxID) {
//  auto entry = taxIDsAndEntries.find(taxID);
//  while (entry != taxIDsAndEntries.end()) {
//    if (entry->second.taxonomyID == 36549) {
//      return true;
//    } else if (entry->second.parentTaxonomyID != 1)
//      entry = taxIDsAndEntries.find(entry->second.parentTaxonomyID);
//    else
//      return false;
//  }
//  return false;
//}
// bool TaxonomyDB::isArtificial(uint32_t taxID) {
//  auto entry = taxIDsAndEntries.find(taxID);
//  while (entry != taxIDsAndEntries.end()) {
//    if (entry->second.taxonomyID == 81077) {
//      return true;
//    } else if (entry->second.parentTaxonomyID != 1)
//      entry = taxIDsAndEntries.find(entry->second.parentTaxonomyID);
//    else
//      return false;
//  }
//  return false;
//}
// bool TaxonomyDB::isEndogenousRetrovirus(uint32_t taxID) {
//  auto entry = taxIDsAndEntries.find(taxID);
//  while (entry != taxIDsAndEntries.end()) {
//    if (entry->second.taxonomyID == 206037) {
//      return true;
//    } else if (entry->second.parentTaxonomyID != 1)
//      entry = taxIDsAndEntries.find(entry->second.parentTaxonomyID);
//    else
//      return false;
//  }
//  return false;
//}
// std::string TaxonomyDB::getSpecies(uint32_t taxID) {
//  uint32_t speciesID = 0;
//  auto entry = taxIDsAndEntries.find(taxID);
//  while (entry != taxIDsAndEntries.end() &&
//         entry->second.parentTaxonomyID != 1) {
//    if (entry->second.rank == "species") {
//      speciesID = entry->second.taxonomyID;
//      break;
//    } else
//      entry = taxIDsAndEntries.find(entry->second.parentTaxonomyID);
//  }
//  if (speciesID != 0) {
//    std::string speciesName = getScientificName(speciesID);
//    if (speciesName.size()) {
//      size_t endOfFirstWord = speciesName.find(' ');
//      if (endOfFirstWord != std::string::npos &&
//          endOfFirstWord < (speciesName.size() - 1)) {
//        return speciesName.substr(endOfFirstWord + 1, std::string::npos);
//      }
//    }
//  }
//  return std::string();
//}
//
// uint32_t TaxonomyDB::getSpeciesTaxID(uint32_t taxonomyID) {
//  uint32_t speciesID = 0;
//  auto entry = taxIDsAndEntries.find(taxonomyID);
//  while (entry != taxIDsAndEntries.end() &&
//         entry->second.parentTaxonomyID != 1) {
//    if (entry->second.rank == "species") {
//      speciesID = entry->second.taxonomyID;
//      break;
//    } else
//      entry = taxIDsAndEntries.find(entry->second.parentTaxonomyID);
//  }
//  return speciesID;
//}
}
#endif /* TAXONOMYDATABASE_H_ */
