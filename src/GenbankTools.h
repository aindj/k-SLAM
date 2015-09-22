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
#ifndef GENBANKTOOLS_H_
#define GENBANKTOOLS_H_
#include <string>
#include <mutex>
#include <parallel/algorithm>
#include "sequenceTools.h"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <unordered_map>
#include <inttypes.h>
#include "TaxonomyDatabase.h"
#include <boost/progress.hpp>
/*
 * Contains a variety of classes and functions
 * which deal with files in the genbank format
 */
/*
 * Note about Boost::serialization:
 * The database of genbank entries is stored on disk using
 * boost serialization for speed.
 */
namespace SLAM {
/*
 * Represents a coding sequence
 * Start/Stop: the zero indexed start/stop position
 * Complement: the complement field as listed in the genbank entry
 */
class CDS {
 public:
  CDS(const uint32_t start, const uint32_t stop, const bool complement)
      : start(start), stop(stop), complement(complement) {}
  ;
  CDS() {}
  ;
  uint32_t start = 0;
  uint32_t stop = 0;
  bool complement = false;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &start;
    ar &stop;
    ar &complement;
  }
};
/*
 * Represents a gene, with fields as listed in the Genbank File Format
 */
class Gene {
 public:
  Gene(const std::string &geneName, const std::string &locusTag,
       const std::string &proteinID, const std::string &product,
       const std::string &referenceSequence, const CDS &codingSequence)
      : geneName(geneName),
        locusTag(locusTag),
        proteinID(proteinID),
        product(product),
        referenceSequence(referenceSequence),
        codingSequence(codingSequence) {}
  ;
  Gene() {}
  ;
  // todo this may not be the correct equality thing
  inline bool operator==(const Gene &other) const {
    if (this->proteinID.size() == 0 && other.proteinID.size() == 0)
      return this->geneName == other.geneName;
    if (this->proteinID == other.proteinID)
      return this->product == other.product;
    else
      return this->proteinID == other.proteinID;
  }
  ;
  std::string geneName;
  std::string locusTag;
  std::string proteinID;
  std::string product;
  std::string referenceSequence;
  uint32_t geneID=0;
  CDS codingSequence;
  int count = 1;
  bool isTRNA = false;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &geneName;
    ar &locusTag;
    ar &proteinID;
    ar &product;
    ar &referenceSequence;
    ar &geneID;
    ar &codingSequence;
  }
};

// todo change this if gene equality changes
/*
 * Sorts based on gene name then product then protein id.
 */
struct geneSort {
  bool operator()(const Gene &i, const Gene &j) const {
    if (i.proteinID.size() == 0 && j.proteinID.size() == 0)
      return i.geneName < j.geneName;
    if (i.proteinID == j.proteinID)
      return i.product < j.product;
    else
      return i.proteinID < j.proteinID;
  }
};

/*
 * Represents an entry in the genbank database
 * speciesTaxID: species of the entry, found by traversing
 * 		taxonomy tree until species level reached
 * isPlasmid: Whether the entry is from a plasmid sequence
 * 		(this is inferred from DEFINITION line of entry,
 * 		 see parse function for more details)
// * is16S: Same as isPlasmid but for 16S rRNA
 */
class GenbankEntry {
 public:
  std::string bases;
  uint32_t taxonomyID = 0;
  uint32_t speciesTaxID = 0;
  uint32_t genbankID = 0;
  bool isPlasmid = false;
  bool is16S = false;
  std::string locusTag;
  std::string definition;
  std::string organismName;
  std::string taxonomy;
  std::string strain;
  std::vector<Gene> genes;
  std::vector<Gene> getGenesInRange(const uint32_t startPos,
                                    const uint32_t endPos) const;
  const Gene *getGene(const int32_t startPos, const int32_t endPos) const;
  void writeIndexEntry(std::ofstream &outFile);
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &bases;
    ar &taxonomyID;
    ar &genbankID;
    ar &isPlasmid;
    ar &is16S;
    ar &locusTag;
    ar &genes;
  }
};

/*
 * Returns gene which has the biggest overlap with the range
 * sequence[startPos]->sequence[endPos] inclusive
 */
const Gene *GenbankEntry::getGene(const int32_t startPos,
                                  const int32_t endPos) const {
  const Gene *bestMatch = nullptr;
  int32_t largestOverlap = 0;
  for (auto &gene : genes) {
    auto &cds = gene.codingSequence;
    int32_t numBasesOverlap = 0;
    numBasesOverlap =
        std::min<int>(endPos, cds.stop) - std::max<int>(startPos, cds.start);
    if (numBasesOverlap > largestOverlap) {
      bestMatch = &gene;
      largestOverlap = numBasesOverlap;
    }
  }
  return bestMatch;
}
/*
 * Represents a collection of genbank entries
 */
class GenbankIndex {
 public:
  std::vector<GenbankEntry> entries;
  template <typename KMerInt, unsigned K>
  void getKMers(std::vector<KMerAndData<KMerInt, K>> &kMers,
                unsigned gap) const;
  void buildIndex(std::vector<std::string> indexFileNames, uint32_t taxID);
//  void parseFASTA(const std::string fastaFileName);
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &entries;
  }
  void writeIndexToBoostSerial(const std::string serialFileName) {
    std::ofstream ofs(serialFileName);
    boost::archive::text_oarchive oa(ofs);
    oa << *this;
  }
  GenbankIndex getReducedDatabase();
};
/*
 * Adds kMers (separated by 'gap') from all genbank entries to a vector
 */
template <typename KMerInt, unsigned K>
void GenbankIndex::getKMers(std::vector<KMerAndData<KMerInt, K>> &kMers,
                            unsigned gap) const {
  log("Getting k-mers from index");
  size_t numKMers = kMers.size();
  getKMers_parallel(entries, kMers, true, gap);
  log("Obtained " + std::to_string(kMers.size() - numKMers) +
      " k-mers from index");
}
///*
// * Creates a GenbankIndex from a FASTA file and a file which maps locuses to
// * taxonomyIDs. The file should be in the format LOCUS tab TaxID
// */
inline void createIndexFromFASTA(std::vector<std::string> fastaFileNames,
                                const std::string outFileName) {
  GenbankIndex index;
  std::ofstream outFile(outFileName);
  boost::progress_display show_progress(fastaFileNames.size());
  for (auto &fileName : fastaFileNames) {
    std::string line;
    log("Parsing\t" + fileName);
    std::ifstream in;
    in.open(fileName);
    if (!in.good()) throw std::runtime_error("unable to open FASTA file");
    GenbankEntry newEntry;
    while (safeGetline(in, line)) {
      if (line.size() == 0) continue;
      if (line[0] == '>') {
        if (newEntry.bases.size()) {
          index.entries.push_back(std::move(newEntry));
        }
        newEntry.bases.clear();
        newEntry.locusTag.clear();
        auto spacePos=line.find(' ');
        if((spacePos != std::string::npos) && (spacePos != 0))
          newEntry.locusTag=line.substr(1,spacePos-1);
      } else {
        newEntry.bases.append(line);
      }
    }
    if (newEntry.bases.size()) {
      index.entries.push_back(std::move(newEntry));
    }
    ++show_progress;
  }
  for (auto &entry : index.entries) {
    inPlaceConvertToUpperCase(entry.bases);
  }
  index.writeIndexToBoostSerial(outFileName);
}
//void GenbankIndex::parseFASTA(const std::string fastaFileName,
//                              const std::string taxMapFileName) {
//  std::string line;
//  std::ifstream taxMap;
//  log("Reading taxonomy map file: " + std::string(taxMapFileName));
//  taxMap.open(taxMapFileName);
//  if (!taxMap.good()) throw std::runtime_error("unable to open taxonomy file");
//  std::unordered_map<std::string, uint32_t> locusToTaxIDMap;
//  while (safeGetline(taxMap, line)) {
//    auto tabPos = line.find('\t');
//    if (tabPos == std::string::npos) continue;
//    locusToTaxIDMap.insert({
//      line.substr(0, tabPos), std::stoi(line.substr(tabPos))
//    });
//  }
//
//  std::ifstream in;
//  log("Reading FASTA file: " + std::string(fastaFileName));
//  in.open(fastaFileName);
//  if (!in.good()) throw std::runtime_error("unable to open FASTA file");
//  GenbankEntry newEntry;
//  while (safeGetline(in, line)) {
//    if (line.size() == 0) continue;
//    if (line[0] == '>') {
//      if (newEntry.bases.size() && newEntry.genbankID) {
//        entries.push_back(std::move(newEntry));
//      }
//      newEntry.bases.clear();
//      newEntry.genbankID = 0;
//      newEntry.locusTag.clear();
//      std::vector<std::string> tokens = tokenise(line, "|");
//      if (tokens.size() < 5) continue;
//      try {
//        newEntry.genbankID = std::stoi(tokens[2]);
//      }
//      catch (...) {
//        continue;
//      }
//      newEntry.locusTag = tokens[4];
//    } else {
//      newEntry.bases.append(line);
//    }
//  }
//  if (newEntry.bases.size() && newEntry.genbankID) {
//    entries.push_back(std::move(newEntry));
//  }
//  for (auto &entry : entries) {
//    inPlaceConvertToUpperCase(entry.bases);
//  }
//}

void printEntry(const GenbankEntry &entry);
GenbankIndex getIndexFromBoostSerial(const std::string serialFileName);
void parseSection(const std::string &field, GenbankEntry &entry);
void createIndexFromGBFF(std::vector<std::string> genbankFileNames,
                         const std::string outFileName);

inline void printEntry(const GenbankEntry &entry) {
  std::cout << "TaxID\t" << entry.taxonomyID << "\n"
            << "GID\t" << entry.genbankID << "\n"
            << "Locus\t" << entry.locusTag << "\n"
            << "Length\t" << entry.bases.length() << "\n";
  for (auto &gene : entry.genes) {
    std::cout << "Gene\t" << gene.geneName << "\n";
    std::cout << "locus\t" << gene.locusTag << "\n";
    std::cout << "product\t" << gene.product << "\n";
    std::cout << "protein\t" << gene.proteinID << "\n";
    std::cout << "cdsb\t" << gene.codingSequence.start << "\n";
    std::cout << "cdse\t" << gene.codingSequence.stop << "\n";
  }
}
/*
 * Reads index from the boost serialization output. This is how the database
 * is initialized
 */
inline GenbankIndex getIndexFromBoostSerial(const std::string serialFileName) {
  log("Building index from serial file\t " + std::string(serialFileName));
  GenbankIndex index;
  std::ifstream ifs(serialFileName);
  if (!ifs.good()) throw std::runtime_error("unable to open index file");
  boost::archive::text_iarchive ia(ifs);
  ia >> index;
  return index;
}
/*
 * Reads a section of a genbank format entry
 */
inline void parseSection(const std::string &field, GenbankEntry &entry) {
  auto start = find_if(field.begin(), field.end(), [](const char & i) {
    return i != ' ';
  });
  if (start == field.end()) return;
  size_t startPos = 0;
  size_t endPos = 0;
  auto stop = find_if(start, field.end(), [](const char & i) {
    return i == ' ';
  });
  std::string tag(start, stop);
  start = find_if(stop, field.end(), [](const char & i) {
    return i != ' ';
  });
  if (tag == "VERSION") {
    stop = find_if(start, field.end(), [](const char & i) {
      return i == ' ';
    });
    entry.locusTag = std::string(start, stop);
    start = find_if(stop, field.end(), [](const char & i) {
      return isdigit(i);
    });
    stop = field.end();
    try {
      entry.genbankID = stoul(std::string(start, stop));
    }
    catch (...) {
    }
  } else if (tag == "DEFINITION") {
    entry.definition = std::string(start, field.end());
  } else if (tag == "source") {
    startPos = field.find("/db_xref=\"taxon:");
    endPos = field.find('\"', startPos);
    if (startPos != std::string::npos && endPos != std::string::npos) {
      startPos += 16;
      if (startPos < field.size()) try {
          entry.taxonomyID = stoul(field.substr(startPos, endPos - startPos));
        }
      catch (...) {
      }
    }
  } else if (tag == "CDS" || tag == "tRNA" || tag == "gene") {
    Gene gene;
      //Coding Sequence
    start = find_if(start, field.end(), [](const char & i) {
      return isdigit(i);
    });
    stop = find_if(start, field.end(), [](const char & i) {
      return !isdigit(i);
    });
    try {
      gene.codingSequence.start = stoul(std::string(start, stop));
    }
    catch (...) {
    }
    start = find_if(stop, field.end(), [](const char & i) {
      return isdigit(i);
    });
    stop = find_if(start, field.end(), [](const char & i) {
      return !isdigit(i);
    });
    try {
      gene.codingSequence.stop = stoul(std::string(start, stop));
    }
    catch (...) {
    }

    //Product
    startPos = field.find("/product=\"");
    if (startPos != std::string::npos) {
      startPos += 10;
      endPos = field.find("\"", startPos);
      if (endPos != std::string::npos) {
        if (startPos < field.size())
          gene.product = field.substr(startPos, endPos - startPos);
      }
    }

    //Protein ID
    startPos = field.rfind("/protein_id=\"");
    if (startPos != std::string::npos) {
      startPos += 13;
      endPos = field.find("\"", startPos);
      if (endPos != std::string::npos) {
        if (startPos < field.size())
          gene.proteinID = field.substr(startPos, endPos - startPos);
      }
    }

    //Locus tag
    startPos = field.find("/locus_tag=\"");
    if (startPos != std::string::npos) {
      startPos += 12;
      endPos = field.find("\"", startPos);
      if (endPos != std::string::npos) {
        if (startPos < field.size())
          gene.locusTag = field.substr(startPos, endPos - startPos);
      }
    }

    //GeneID
    startPos = field.find("GeneID:");
    if (startPos != std::string::npos) {
      startPos += 7;
      endPos = field.find("\"", startPos);
      if (endPos != std::string::npos) {
        if (startPos < field.size())
          gene.geneID = std::stoul(field.substr(startPos, endPos - startPos));
      }
    }

    //Gene name
    startPos = field.find("/gene=\"");
    if (startPos != std::string::npos) {
      startPos += 7;
      endPos = field.find('\"', startPos);
      if (endPos != std::string::npos) {
        if (startPos < field.size())
          gene.geneName = field.substr(startPos, endPos - startPos);
      }
    }
    gene.referenceSequence = entry.locusTag;
    entry.genes.push_back(gene);
  } else if (tag.size() && isdigit(tag[0])) {
    for (; start != field.end(); start++) {
      if (*start != ' ') entry.bases.push_back(toupper(*start));
    }
  }
}
/*
 * Builds a Genbank index from raw files in Genbank Flat File format and
 * produces a SLAM database using boost serialize
 */
inline void createIndexFromGBFF(std::vector<std::string> genbankFileNames,
                                const std::string outFileName) {
  TaxonomyDB taxDB("taxDB");
  GenbankIndex index;
  std::ofstream outFile(outFileName);
  boost::progress_display show_progress(genbankFileNames.size());
  for (auto &fileName : genbankFileNames) {
    log("Parsing\t" + fileName);
    std::ifstream genbankFile;
    genbankFile.open(fileName);
    if (!genbankFile.good())
      throw std::runtime_error("unable to open index file");
    std::string line, section;
    GenbankEntry entry;
    while (getline(genbankFile, line)) {
      if (line.size() == 0) continue;
      auto startPos = line.find_first_not_of(' ');
      if (startPos < 12) {
        parseSection(section, entry);
        section = line;
        if (line == "ORIGIN") {
          continue;
        } else if (line == "//") {
          std::sort(entry.genes.begin(), entry.genes.end(),
                    [](const Gene & i, const Gene & j) {
            if (i.codingSequence.start == j.codingSequence.start)
              return i.proteinID.size() > j.proteinID.size();
            return i.codingSequence.start < j.codingSequence.start;
          });
          auto it = std::unique(entry.genes.begin(), entry.genes.end(),
                                [](const Gene & i, const Gene & j) {
            return i.codingSequence.start == j.codingSequence.start;
          });
          entry.genes.resize(std::distance(entry.genes.begin(), it));
          index.entries.push_back(entry);
          entry = GenbankEntry();
        }
      } else if (startPos == std::string::npos)
        continue;
      else if (startPos > 0) {
        section.append(line.substr(startPos - 1, std::string::npos));
      }
    }
    ++show_progress;
  }
  index.writeIndexToBoostSerial(outFileName);
}
//std::vector<Gene> GenbankEntry::getGenesInRange(const uint32_t startPos,
//                                                const uint32_t endPos) const {
//  std::vector<Gene> intersectedGenes;
//  for (auto &gene : genes) {
//    auto &cds = gene.codingSequence;
//    if ((endPos >= cds.start) && (startPos <= cds.stop)) {
//      intersectedGenes.push_back(gene);
//    }
//  }
//  return intersectedGenes;
//}
///*
// * Produces a new GenbankIndex with one strain per species
// * This function chooses the strain with the longest sequence length
// * to pick the most complete assemblies
// */
//GenbankIndex GenbankIndex::getReducedDatabase() {
//  GenbankIndex newIndex;
//  struct allEntriesForTax {
//    uint32_t taxID = 0;
//    uint32_t speciesTaxID = 0;
//    uint32_t length = 0;
//    std::vector<GenbankEntry *> entries;
//  };
//  std::cout << "Num entries\t" << entries.size() << "\n";
//  std::unordered_map<uint32_t, allEntriesForTax> entriesByTaxID;
//  for (auto &entry : entries) {
//    auto it = entriesByTaxID.find(entry.taxonomyID);
//    if (it == entriesByTaxID.end()) {
//      allEntriesForTax temp;
//      temp.taxID = entry.taxonomyID;
//      temp.speciesTaxID = entry.speciesTaxID;
//      temp.entries.push_back(&entry);
//      temp.length = entry.bases.size();
//      entriesByTaxID.insert({temp.taxID, temp});
//    } else {
//      it->second.entries.push_back(&entry);
//      it->second.length += entry.bases.size();
//    }
//  }
//  std::vector<allEntriesForTax> entriesByTaxIDVec;
//  for (auto &el : entriesByTaxID) {
//    entriesByTaxIDVec.push_back(el.second);
//  }
//  std::sort(entriesByTaxIDVec.begin(), entriesByTaxIDVec.end(),
//            [](const allEntriesForTax &i, const allEntriesForTax &j) {
//    if (i.speciesTaxID == j.speciesTaxID)
//      return i.length > j.length;
//    else
//      return i.speciesTaxID < j.speciesTaxID;
//  });
//  auto end =
//      std::unique(entriesByTaxIDVec.begin(), entriesByTaxIDVec.end(),
//                  [](const allEntriesForTax &i, const allEntriesForTax &j) {
//        return i.speciesTaxID == j.speciesTaxID;
//      });
//  entriesByTaxIDVec.resize(std::distance(entriesByTaxIDVec.begin(), end));
//  for (auto &tax : entriesByTaxIDVec) {
//    for (auto &entry : tax.entries) {
//      newIndex.entries.push_back(*entry);
//    }
//  }
//  return newIndex;
//}
/////*
//// * Builds a "GenbankIndex" from a list of file names
//// * The input files used here are of the type produced by the parse function
//// * If the taxID variable is set then only entries with this taxonomy ID are
//// * stored
//// */
//void GenbankIndex::buildIndex(std::vector<std::string> indexFileNames,
//                              uint32_t taxID = 0) {
//  log("Building index");
//  size_t count = 0;
//  std::ifstream in;
//  for (auto indexFileName : indexFileNames) {
//    log("Reading index file: " + indexFileName);
//    in.open(indexFileName);
//    if (!in.good()) throw std::runtime_error("unable to open index file");
//    GenbankEntry tempEntry;
//    while (in.good()) {
//      std::string line;
//      getline(in, line);
//      if (line.size() == 0) continue;
//      tempEntry.genbankID = std::stoi(line);
//      getline(in, line);
//      tempEntry.taxonomyID = std::stoi(line);
//      if (taxID != 0) {
//        if (tempEntry.taxonomyID != taxID) {
//          while (line != "//") getline(in, line);
//          continue;
//        }
//      }
//      getline(in, line);
//      tempEntry.speciesTaxID = std::stoi(line);
//      getline(in, line);
//      tempEntry.locusTag = line;
//      getline(in, line);
//      tempEntry.organismName = line;
//      getline(in, line);
//      tempEntry.taxonomy = line;
//      getline(in, line);
//      tempEntry.strain = line;
//      getline(in, line);
//      tempEntry.isPlasmid = std::stoi(line);
//      getline(in, line);
//      tempEntry.genes.clear();
//      while (line == "G") {
//        Gene gene;
//        getline(in, line);
//        gene.codingSequence.complement = std::stoi(line);
//        getline(in, line);
//        gene.codingSequence.start = std::stoi(line);
//        getline(in, line);
//        gene.codingSequence.stop = std::stoi(line);
//        if (gene.codingSequence.start != 0 && gene.codingSequence.stop != 0) {
//          // convert to zero indexed
//          gene.codingSequence.start--;
//          gene.codingSequence.stop--;
//        }
//        getline(in, line);
//        gene.geneName = line;
//        getline(in, line);
//        gene.locusTag = line;
//        getline(in, line);
//        gene.proteinID = line;
//        getline(in, line);
//        gene.product = line;
//        getline(in, line);
//        gene.referenceSequence = tempEntry.locusTag;
//        getline(in, line);
//        tempEntry.genes.push_back(gene);
//      }
//      tempEntry.bases = std::move(line);
//      getline(in, line);
//      count++;
//      if (tempEntry.bases.size() != 0) {
//        entries.push_back(std::move(tempEntry));
//      }
//      // todo obv remove this
//      //            if (count > 0) break;
//    }
//    in.close();
//  }
//  //  log("Upper");
//  //  for(auto & entry : entries){
//  //  inPlaceConvertToUpperCase(entry.bases);
//  //  }
//  log("Constructed an index of " + std::to_string(entries.size()) +
//      " Genbank entries");
//}
/*
 * Converts a genbank entry in flat file format into a GenbankEntry in SLAM
 * format
 */
//GenbankEntry convertLinesToEntry(std::vector<std::string> &lines,
//                                 const TaxonomyDB &taxDB) {
//  GenbankEntry entry;
//  // std::vector<std::string> geneNames;
//  std::string gene, product, proteinID, locusTag;
//  CDS codingSequence;
//  // std::pair<unsigned int,unsigned int> CDSRange={0,0};
//  for (unsigned int i = 0; i < lines.size(); i++) {
//    std::vector<std::string> tokens = tokenise(lines[i], " =\":");
//    if (tokens.size() > 0) {
//      if (tokens[0] == "CDS") {
//        gene.clear();
//        product.clear();
//        proteinID.clear();
//        locusTag.clear();
//        std::vector<std::string> tempTokens = tokenise(lines[i], " ().<>");
//        if (tempTokens.size() > 1) {
//          if (tempTokens[1] == "complement" && tempTokens.size() > 3) {
//            try {
//              codingSequence.complement = true;
//              codingSequence.start = stoi(tempTokens[2]);
//              codingSequence.stop = stoi(tempTokens[3]);
//            }
//            catch (...) {
//              codingSequence = CDS(0, 0, 0);
//              continue;
//            }
//          } else if (tempTokens.size() > 2) {
//            try {
//              codingSequence.complement = true;
//              codingSequence.start = stoi(tempTokens[1]);
//              codingSequence.stop = stoi(tempTokens[2]);
//            }
//            catch (...) {
//              codingSequence = CDS(0, 0, 0);
//              continue;
//            }
//          }
//        }
//      } else if (tokens[0] == "/protein_id") {
//        if (tokens.size() > 1) {
//          proteinID = tokens[1];
//        }
//      }
//      else if (tokens[0] == "/locus_tag") {
//        if (tokens.size() > 1) {
//          locusTag = tokens[1];
//        }
//      } else if (tokens[0] == "/gene") {
//        if (tokens.size() > 1) {
//          gene = tokens[1];
//        }
//      } else if (tokens[0] == "DEFINITION") {
//        for (auto &token : tokens) {
//          if (token == "16S" || token == "23S") entry.is16S = true;
//          if (token == "plasmid" || token == "Plasmid") {
//            entry.isPlasmid = true;
//          }
//        }
//      } else if (tokens[0] == "/product") {
//        product.clear();
//        for (unsigned int j = 1; j < tokens.size(); j++) {
//          if (j != 1) product.append(" ");
//          product.append(tokens[j]);
//        }
//        for (unsigned int j = i + 1; j < lines.size(); j++) {
//          std::vector<std::string> tempTokens = tokenise(lines[j], " \"");
//          if ((lines[j].substr(0, 5) == "     ") && tempTokens.size() &&
//              (tempTokens[0][0] != '/')) {
//            for (unsigned int k = 0; k < tempTokens.size(); k++) {
//              product.append(" ");
//              product.append(tempTokens[k]);
//            }
//          } else
//            break;
//        }
//      } else if (tokens[0] == "ORGANISM") {
//        entry.organismName.clear();
//        for (unsigned int j = 1; j < tokens.size(); j++) {
//          if (j != 1) entry.organismName.append(" ");
//          entry.organismName.append(tokens[j]);
//        }
//        for (unsigned int j = i + 1; j < lines.size(); j++) {
//          std::vector<std::string> tempTokens = tokenise(lines[j], " \".");
//          if (tempTokens.size()) {
//            if (tempTokens[0] == "REFERENCE" || tempTokens[0] == "COMMENT")
//              break;
//            for (unsigned int k = 0; k < tempTokens.size(); k++) {
//              if (!(j == (i + 1) && k == 0)) entry.taxonomy.append(" ");
//              entry.taxonomy.append(tempTokens[k]);
//            }
//          }
//        }
//      } else if (tokens[0] == "/translation") {
//        entry.genes.push_back(Gene(gene, locusTag, proteinID, product,
//                                   entry.locusTag, codingSequence));
//      } else if (tokens[0] == "VERSION") {
//        if (tokens.size() > 3) {
//          entry.locusTag = tokens[1];
//          entry.genbankID = stoi(tokens[3]);
//        }
//      } else if (tokens[0] == "/db_xref") {
//        if (tokens.size() > 2) {
//          if (tokens[1] == "taxon" && entry.taxonomyID == 0) {
//            uint32_t taxID = stoi(tokens[2]);
//            entry.taxonomyID = taxID;
//            uint32_t speciesTaxID = taxDB.getTaxIDAtRank(taxID, "species");
//            if (speciesTaxID)
//              entry.speciesTaxID = speciesTaxID;
//            else
//              std::cout << "Error" << std::endl;
//          }
//        }
//      } else if (tokens[0] == "/strain") {
//        auto equalsPos = lines[i].find('=');
//        if (lines[i].size() > (equalsPos + 2)) {
//          entry.strain = lines[i].substr(equalsPos + 2, std::string::npos);
//          entry.strain.resize(entry.strain.size() - 1);
//        }
//      } else if (tokens[0] == "ORIGIN" && lines[i].size() >= 6 &&
//                 lines[i].substr(0, 6) == "ORIGIN") {
//        i++;
//        for (; i < lines.size(); i++) {
//          std::vector<std::string> tempTokens = tokenise(lines[i], " ");
//          for (unsigned int j = 1; j < tempTokens.size(); j++) {
//            entry.bases.append(tempTokens[j]);
//          }
//        }
//      }
//    }
//  }
//  return entry;
//}
//std::unordered_map<uint32_t, uint32_t> buildTaxIDToSpeciesTaxIDMap(
//    const std::string taxDBFileName) {
//  std::unordered_map<uint32_t, uint32_t> taxIDMap;
//  std::ifstream taxDBFile;
//  taxDBFile.open(taxDBFileName);
//  std::string line;
//  while (!getline(taxDBFile, line).eof()) {
//    std::vector<std::string> tokens = tokenise(line, " \t");
//    taxIDMap.insert({stoi(tokens[2]), stoi(tokens[1])});
//  }
//  taxDBFile.close();
//  return taxIDMap;
//}
//void createIndexFromGBFF(std::vector<std::string> genbankFileNames,
//                         const std::string outFileName) {
//  TaxonomyDB taxDB("taxDB");
//  GenbankIndex index;
//  std::ofstream outFile(outFileName);
//  //  std::unordered_map<uint32_t, uint32_t> taxIDMap =
//  //      buildTaxIDToSpeciesTaxIDMap("/storage/taxid");
////  boost::progress_display show_progress(genbankFileNames.size());
//  for (auto &fileName : genbankFileNames) {
//    //    std::cout << "Parsing\t" << fileName << std::endl;
//    log("Parsing\t" + fileName);
//    std::ifstream genbankFile;
//    genbankFile.open(fileName);
//    if (!genbankFile.good())
//      throw std::runtime_error("unable to open index file");
//    std::vector<std::string> lines;
//    while (genbankFile.good()) {
//      std::string temp;
//      getline(genbankFile, temp);
//      if (temp == "//") {
//        GenbankEntry tempEntry = convertLinesToEntry(lines, taxDB);
////
///std::cout<<tempEntry.taxonomyID<<"\t"<<tempEntry.speciesTaxID<<std::endl;
//        lines.clear();
//        //        if((tempEntry.genbankID!=0) && (!tempEntry.is16S))
//        //todo do again and remove this
////        if (!tempEntry.is16S) {
//          inPlaceConvertToUpperCase(tempEntry.bases);
////          printEntry(tempEntry);
////          index.entries.push_back(tempEntry);
//          //          tempEntry.writeIndexEntry(outFile);
////        }
//      } else {
//        lines.push_back(temp);
//      }
//    }
////    ++show_progress;
//  }
//  index.writeIndexToBoostSerial(outFileName);
//}
}

#endif /* GENBANKTOOLS_H_ */
