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
#ifndef METAGENOMICRESULTS_H_
#define METAGENOMICRESULTS_H_

#include "GenbankTools.h"
#include <algorithm>
#include "Overlap.h"
#include "PairedOverlap.h"
#include "GenbankTools.h"
#include "TaxonomyDatabase.h"
#include <inttypes.h>
namespace SLAM {
/*
 * Represents all of the reads and genes assigned to a taxonomyID
 */
class IdentifiedTaxonomy {
 public:
  uint32_t taxonomyID = 0;
  std::string strain;
  std::vector<std::string> reads;
  std::vector<Gene> genes;
//  std::vector<const GenbankEntry *> entries;
  void addGene(const Gene &gene);
  void addGenes(const std::vector<Gene> &genes);
  void addRead(const std::string &read);
};

void IdentifiedTaxonomy::addRead(const std::string &read) {
  reads.push_back(read);
}

void IdentifiedTaxonomy::addGene(const Gene &newGene) {
  if (std::find(genes.begin(), genes.end(), newGene) == genes.end())
    genes.push_back(newGene);
}
void IdentifiedTaxonomy::addGenes(const std::vector<Gene> &genes) {
  for (auto &gene : genes)
    addGene(gene);
}

void sortResults(std::vector<IdentifiedTaxonomy> &resultsVec);
void writeAbbreviatedResultsFile(std::vector<IdentifiedTaxonomy> &results,
                                 const std::string outFileName,
                                 const TaxonomyDB &taxDB);
void writeResults(std::vector<IdentifiedTaxonomy> &results, std::ostream &out,
                  const TaxonomyDB &taxDB, const unsigned numReads);
template <typename FASTQType, typename ITType>
std::vector<
    IdentifiedTaxonomy> convertAlignmentsToIdentifiedTaxonomies_parallel(
    const ITType begin, const ITType end, const std::vector<FASTQType> &reads,
    const GenbankIndex &index, const TaxonomyDB &taxDB);
std::vector<IdentifiedTaxonomy> combineTaxonomies(
    std::vector<IdentifiedTaxonomy> &identifiedTaxonomies);
template <typename ITType>
IdentifiedTaxonomy combineRangeOfIdentifiedTaxonomy(ITType begin, ITType end);
std::string getXML(const Gene &gene);
std::string getXML(const IdentifiedTaxonomy &entry,
                   const unsigned int totalNumReads, const TaxonomyDB &taxDB);
template <typename FASTQType>
IdentifiedTaxonomy getResultFromPairedOverlaps(
    ReadPairAndOverlaps &overlaps, const GenbankIndex &index,
    const std::vector<FASTQType> &reads, const TaxonomyDB &taxDB,
    std::unordered_map<uint32_t, uint32_t> &taxIDsAndFreq, bool final);
/*
 * Takes a set of paired alignments and infers taxonomy.
 * The taxonomy is then inferred by making a vector of every taxonomy
 * from the overlapping genbank entries and then finding the lowest common
 * ancestor in the taxonomy tree.
 * For each overlap, the best match gene is added to a vector. The vector is
 * then iterated over to find only unique genes.
 */
template <typename FASTQType>
inline IdentifiedTaxonomy getResultFromPairedOverlaps(
    ReadPairAndOverlaps &overlaps, const GenbankIndex &index,
    const std::vector<FASTQType> &reads, const TaxonomyDB &taxDB) {
  IdentifiedTaxonomy result;
  auto &pairedAlignments = overlaps.alignmentPairs;
  if (pairedAlignments.size() == 0) return result;
  std::vector<uint32_t> taxIDs;
  auto &read = reads[overlaps.r1PosInArray];
  for (auto &alignment : pairedAlignments) {
    auto &entry = index.entries[alignment.entryPosInArray];
    taxIDs.push_back(entry.taxonomyID);
    const Gene *gene = entry.getGene(alignment.refStart,alignment.refEnd);
    if (gene != nullptr) {
      result.genes.push_back(*gene);
    }
  }
  std::sort(result.genes.begin(), result.genes.end(), geneSort());
  auto it = std::unique(result.genes.begin(), result.genes.end());
  result.genes.resize(std::distance(result.genes.begin(), it));
  result.addRead(read.sequenceIdentifier);
  result.taxonomyID = taxDB.getLowestCommonAncestor(taxIDs);
  return result;
}

/*
 * Takes a range of "IdentifiedTaxonomy"s with the same tax ids and
 * combines them into one result
 */
template <typename ITType>
inline IdentifiedTaxonomy combineRangeOfIdentifiedTaxonomy(ITType begin,
                                                           ITType end) {
  IdentifiedTaxonomy taxonomy = *begin;
  for (auto tax = begin; tax != end; tax++) {
    if (tax == begin) continue;
    taxonomy.genes
        .insert(taxonomy.genes.end(), tax->genes.begin(), tax->genes.end());
    taxonomy.reads
        .insert(taxonomy.reads.end(), tax->reads.begin(), tax->reads.end());
  }
  std::sort(taxonomy.genes.begin(), taxonomy.genes.end(), geneSort());
  auto first = taxonomy.genes.begin();
  auto last = taxonomy.genes.end();
  if (first != last) {
    auto result = first;
    while (++first != last) {
      if (!(*result == *first)) {
        *(++result) = *first;
      } else
        result->count++;
    }
    ++result;
    taxonomy.genes.resize(std::distance(taxonomy.genes.begin(), result));
  }
  return taxonomy;
}

/*
 * Takes a set of identified taxonomies (one from every read) and combines
 * them such that there is just one IdentifiedTaxonomy per taxonomyID
 */
inline std::vector<IdentifiedTaxonomy> combineTaxonomies(
    std::vector<IdentifiedTaxonomy> &identifiedTaxonomies) {
  log("Combining taxonomies");
  __gnu_parallel::sort(
      identifiedTaxonomies.begin(), identifiedTaxonomies.end(),
      [](const IdentifiedTaxonomy & i, const IdentifiedTaxonomy & j) {
    return i.taxonomyID < j.taxonomyID;
  });
  std::vector<IdentifiedTaxonomy> combinedTaxonomies;
  if (identifiedTaxonomies.size() == 0) return combinedTaxonomies;
  uint32_t testTaxID = 0;
  auto start = identifiedTaxonomies.begin();
  for (auto tax = identifiedTaxonomies.begin();
       tax != identifiedTaxonomies.end(); tax++) {
    if (tax == identifiedTaxonomies.begin()) continue;
    if (tax->taxonomyID != testTaxID) {
      if (testTaxID != 0)
        combinedTaxonomies.push_back(
            combineRangeOfIdentifiedTaxonomy(start, tax));
      testTaxID = tax->taxonomyID;
      start = tax;
    }
  }
  if ((start != identifiedTaxonomies.end()) && (start->taxonomyID != 0))
    combinedTaxonomies.push_back(
        combineRangeOfIdentifiedTaxonomy(start, identifiedTaxonomies.end()));
  return combinedTaxonomies;
}

/*
 * Takes all overlaps for each read and infers taxonomy
 * See getResultFromPairedOverlaps for more info
 */
template <typename FASTQType, typename ITType>
inline std::vector<IdentifiedTaxonomy>
convertAlignmentsToIdentifiedTaxonomies_parallel(
    const ITType begin, const ITType end, const std::vector<FASTQType> &reads,
    const GenbankIndex &index, const TaxonomyDB &taxDB) {
  log("Converting alignments to metagenomic results");
  std::vector<IdentifiedTaxonomy> identifiedTaxonomies(
      std::distance(begin, end));
  typedef typename std::iterator_traits<ITType>::value_type value_type;

  __gnu_parallel::transform(begin, end, identifiedTaxonomies.begin(),
                            [&](value_type & i) {
    return getResultFromPairedOverlaps(i, index, reads, taxDB);
  });
  return identifiedTaxonomies;
}

/*
 * Turns all subspecies taxonomies into their parent species
 */
void convertToSpeciesLevel(std::vector<IdentifiedTaxonomy> &results,
                           const TaxonomyDB &taxDB) {
  for (auto &result : results) {
    if (taxDB.isSubSpecies(result.taxonomyID)) {
      result.taxonomyID = taxDB.getTaxIDAtRank(result.taxonomyID, "species");
    }
  }
}
/*
 * Produces an xml file containing the results
 */
inline void writeResults(std::vector<IdentifiedTaxonomy> &results,
                         std::ostream &out, const TaxonomyDB &taxDB,
                         const unsigned numReads) {
  log("Writing results file");
  sortResults(results);
  unsigned numAligned = 0;
  //  convertToSpeciesLevel(results, taxDB);
  for (auto &result : results) {
    out << getXML(result, numReads, taxDB);
    numAligned += result.reads.size();
  }

  //  std::cout << "Num aligned\t" << numAligned << "\n";
}
inline void writePerReadResults(std::vector<IdentifiedTaxonomy> &results,
                                const std::string outFileName) {
  log("Writing results file");
  std::ofstream out(outFileName);
  for (auto &result : results) {
    for(auto & read : result.reads){
      out<<read<<"\t"<<result.taxonomyID<<"\n";
    }
  }
}
inline void writeAbbreviatedResultsFile(
    std::vector<IdentifiedTaxonomy> &results, const std::string outFileName,
    const TaxonomyDB &taxDB, const unsigned numReads) {
  log("Writing results file");
  std::ofstream outFile(outFileName);
  sortResults(results);
  for (auto &result : results) {
    outFile << taxDB.getScientificName(result.taxonomyID) << "\t"
            << result.reads.size() * 100.0 / numReads << "\n";
  }

  //  std::cout << "Num aligned\t" << numAligned << "\n";
}

/*
 * Sorts results so that the output file is the same for repeated executions
 */
inline void sortResults(std::vector<IdentifiedTaxonomy> &resultsVec) {
  std::sort(resultsVec.begin(), resultsVec.end(),
            [](const IdentifiedTaxonomy & i, const IdentifiedTaxonomy & j) {
    if (i.reads.size() == j.reads.size())
      return i.taxonomyID < j.taxonomyID;
    else
      return i.reads.size() > j.reads.size();
  });
  for (auto &entry : resultsVec) {
    std::sort(entry.reads.begin(), entry.reads.end());
    std::sort(entry.genes.begin(), entry.genes.end(),
              [](const Gene & i, const Gene & j) {
      if (i.count == j.count) {
        if (i.codingSequence.start == j.codingSequence.start) {
          return i.locusTag < j.locusTag;
        }
        return i.codingSequence.start < j.codingSequence.start;
      } else
        return i.count > j.count;
    });
  }
}
inline std::string correctXML(const std::string & input){
  std::string output;
  for(auto c:input){
    switch(c){
      case '<':
        output.append("&lt;");
        break;
      case '>':
        output.append("&gt;");
        break;
      case '&':
        output.append("&amp;");
        break;
      case '\'':
        output.append("&apos;");
        break;
      case '\"':
        output.append("&quot;");
        break;
      default:
        output.push_back(c);
        break;
    }
  }
  return output;
}
inline std::string getXML(const Gene &gene) {
  std::string geneXML = "    <gene ";
  geneXML += "protein=\"";
  geneXML += correctXML(gene.proteinID);
  geneXML += "\" locus=\"";
  geneXML += correctXML(gene.locusTag);
  geneXML += "\" product=\"";
  geneXML += correctXML(gene.product);
  geneXML += "\" GeneID=\"";
  geneXML += std::to_string(gene.geneID);
  geneXML += "\" reference=\"";
  geneXML += correctXML(gene.referenceSequence);
  geneXML += "\" numReads=\"";
  geneXML += std::to_string(gene.count);
  geneXML += "\" cdsStart=\"";
  geneXML += std::to_string(gene.codingSequence.start);
  geneXML += "\" cdsEnd=\"";
  geneXML += std::to_string(gene.codingSequence.stop);
  geneXML += "\">";
  geneXML += correctXML(gene.geneName);
  geneXML += "</gene>";
  return geneXML;
}

inline std::string getXML(const IdentifiedTaxonomy &entry,
                          const unsigned int totalNumReads,
                          const TaxonomyDB &taxDB) {

  std::string outXML;
  outXML += "<taxon>\n";

  outXML += "  <abundance numReads=\"";
  outXML += std::to_string(entry.reads.size());
  outXML += "\">";
  outXML += std::to_string(entry.reads.size() * 100.0 / totalNumReads);
  outXML += "</abundance>\n";

  outXML += "  <taxonomyID>";
  outXML += std::to_string(entry.taxonomyID);
  outXML += "</taxonomyID>\n";
  outXML += "  <lineage>";
  outXML += correctXML(taxDB.getLineage(entry.taxonomyID));
  outXML += "</lineage>\n";

//  std::string rank = taxDB.getRank(entry.taxonomyID);
//  uint32_t genus = taxDB.getTaxIDAtRank(entry.taxonomyID, "genus");

  outXML += "  <name>";
  outXML += correctXML(taxDB.getScientificName(entry.taxonomyID));
  outXML += "</name>\n";

  outXML += "  <genes>\n";
  for (auto &gene : entry.genes) {
    outXML += getXML(gene);
    outXML += "\n";
  }
  outXML += "  </genes>\n";

  outXML += "  <reads>\n";
  for (auto &read : entry.reads) {
    outXML += "    <read>";
    outXML += correctXML(read);
    outXML += "</read>\n";
  }
  outXML += "  </reads>\n";
  outXML += "</taxon>\n";
  return outXML;
}
void fillInNumReadsAligned(
    const std::vector<IdentifiedTaxonomy> &identifiedTaxonomies,
    TaxonomyDB &taxDB, GenbankIndex &index) {
  for (auto &tax : identifiedTaxonomies) {
    auto entryIt = taxDB.taxIDsAndEntries.find(tax.taxonomyID);
    if (entryIt != taxDB.taxIDsAndEntries.end()) {
      entryIt->second.numReadsAligned = tax.reads.size();
    }
  }
  for (auto &entry : index.entries) {
    auto it = taxDB.taxIDsAndEntries.find(entry.taxonomyID);
    if (it != taxDB.taxIDsAndEntries.end()) {
      it->second.genomeSize += entry.bases.size();
    }
  }
  for (auto taxIt = taxDB.taxIDsAndEntries.begin();
       taxIt != taxDB.taxIDsAndEntries.end(); taxIt++) {
    TaxonomyEntry *tax = &taxIt->second;
    if (tax->used) continue;
    auto entry = tax;
    uint64_t runningTotal = 0;
    uint64_t genomeSizeOfChildren = 0;
    uint32_t numBelow = 0;
    while (entry != nullptr) {
      entry->numReadsAlignedToChildren += runningTotal;
      entry->genomeSizeOfChildren += genomeSizeOfChildren;
      entry->numBelow += numBelow;
      if (!entry->used) runningTotal += entry->numReadsAligned;
      if (!entry->used) genomeSizeOfChildren += entry->genomeSize;
      if (!entry->used && entry->genomeSize) numBelow++;
      entry->used = true;
      if (entry->taxonomyID == 1) {
        break;
      }
      entry = entry->parent;
    }
  }
  std::vector<std::pair<std::string, double>> taxAndAbundance;
  for (auto &tax : taxDB.taxIDsAndEntries) {
    if (tax.second.numReadsAligned + tax.second.numReadsAlignedToChildren !=
        0) {
      if (tax.second.taxonomyID == 1) {
        taxAndAbundance.push_back({
          "Root", (tax.second.numReadsAligned +
                   tax.second.numReadsAlignedToChildren)
        });
      }
      std::string lineage = taxDB.getMetaPhlAnLineage(tax.second.taxonomyID);
      if (lineage.size()) {
        double scaleFactor = 1;
        //        if(tax.second.genomeSize && tax.second.numBelow==0){
        //          scaleFactor=1.0/tax.second.genomeSize;
        //        }
        //        else if(tax.second.genomeSizeOfChildren &&
        // tax.second.numBelow
        // &&
        //            tax.second.genomeSize==0){
        //          scaleFactor=1.0*tax.second.numBelow/tax.second.genomeSizeOfChildren;
        //        }
        //        taxAndAbundance.push_back({
        //          lineage, (tax.second.numReadsAligned +
        //                    tax.second.numReadsAlignedToChildren) *
        // scaleFactor
        //        });
        if (tax.second.genomeSize != 0)
          taxAndAbundance.push_back({
            lineage, tax.second.numReadsAligned * 1.0 / tax.second.genomeSize
          });
      }
    }
  }
  std::sort(taxAndAbundance.begin(), taxAndAbundance.end(),
            [](const std::pair<std::string, double> & i,
               const std::pair<std::string, double> & j) {
    return i.second > j.second;
  });
  double max = 0;
  if (taxAndAbundance.size()) {
    max = taxAndAbundance[0].second;
    for (auto &el : taxAndAbundance) {
      if (el.first != "Root")
        std::cout << el.first << "\t" << el.second * 100.0 / max << "\n";
    }
  }
}
void writePerReadResults(
    const std::vector<IdentifiedTaxonomy> &identifiedTaxonomies,
    std::ostream &output) {
  log("Writing per read results");
  for (auto &tax : identifiedTaxonomies) {
    if (tax.reads.size())
      output << tax.reads[0] << "\t" << tax.taxonomyID << "\n";
  }
}
}
#endif /* METAGENOMICRESULTS_H_ */
