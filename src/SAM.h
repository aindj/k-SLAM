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
#ifndef SAM_H_
#define SAM_H_
#include <fstream>
#include <vector>
#include "GenbankTools.h"
#include "Overlap.h"
#include <inttypes.h>
namespace SLAM {
class SequenceDifference {
 public:
  std::string cigar;
  std::string MD;
  uint32_t NM = 0;
  double logProbability = 0;
};
std::vector<double> getLogMatchLookupTable() {
  std::vector<double> table;
  table.push_back(std::log10(1.0 - std::pow(10.0, 1.0 / -10.0)));
  for (int i = 1; i < 100; i++) {
    table.push_back(std::log10(1.0 - std::pow(10.0, i / -10.0)));
  }
  return table;
}
std::vector<double> getLogMismatchLookupTable() {
  std::vector<double> table;
  table.push_back(1 / -10.0);
  for (int i = 1; i < 100; i++) {
    table.push_back(i / -10.0);
  }
  return table;
}
std::string recreateRef(const std::string &query, const std::string &cigar,
                        const std::string &MD) {
  std::string ref;
  std::string strLength;
  int queryPos = 0;
  for (int i = 0; i < cigar.length(); i++) {
    if (isdigit(cigar[i]))
      strLength.push_back(cigar[i]);
    else {
      int length = stoi(strLength);
      char operation = cigar[i];
      for (int j = 0; j < length; j++) {
        switch (operation) {
          case 'M':
            ref.push_back(query[queryPos]);
            queryPos++;
            break;
          case 'D':
            ref.push_back('-');
            break;
          case 'I':
            queryPos++;
            break;
          case 'S':
            //            ref.push_back('-');
            queryPos++;
            break;
          default:
            break;
        }
      }
      strLength.clear();
    }
  }
  strLength.clear();
  queryPos = 0;
  for (int i = 0; i < MD.length(); i++) {
    if (isdigit(MD[i]))
      strLength.push_back(MD[i]);
    else if (MD[i] == '^') {

    } else {
      int length = 0;
      if (strLength.size()) length = stoi(strLength);
      strLength.clear();
      queryPos += length;
      ref[queryPos] = MD[i];
      queryPos++;
    }
  }
  return ref;
}
template <typename FASTQType>
SequenceDifference getCigarAndMD(const Overlap &overlap,
                                 const std::vector<FASTQType> &reads,
                                 const GenbankIndex &genbankIndex) {
  const static auto matchTable = getLogMatchLookupTable();
  const static auto misMatchTable = getLogMismatchLookupTable();
  SequenceDifference sequenceDiff;
  std::vector<std::string> MDcomponents;
  auto &ref = genbankIndex.entries[overlap.entryPosInArray].bases;
  auto query =
      overlap.revComp ? reverseComplement(reads[overlap.readPosInArray].bases)
                      : reads[overlap.readPosInArray].bases;
  auto quality = reads[overlap.readPosInArray].quality;
  if (overlap.revComp) {
    std::reverse(quality.begin(), quality.end());
  }
  if (!overlap.alignment.cigar) return sequenceDiff;
  int refPos = overlap.alignment.ref_begin;
  int queryPos = 0;
  if (overlap.alignment.query_begin > 0) {
    auto length = overlap.alignment.query_begin;
    sequenceDiff.cigar.append(std::to_string(length));
    sequenceDiff.cigar.push_back('S');
    queryPos += length;
  }
  for (auto cigarEl = overlap.alignment.cigar;
       cigarEl < overlap.alignment.cigar + overlap.alignment.cigarLen;
       cigarEl++) {
    int numMatch = 0;
    uint32_t length = *cigarEl >> 4;
    sequenceDiff.cigar.append(std::to_string(length));
    uint32_t operation = *cigarEl & 0b1111;
    std::string temp;
    switch (operation) {
      case 0:  // M
        sequenceDiff.cigar.push_back('M');
        for (int i = 0; i < length; i++) {
//          if (refPos >= ref.size())
//            std::cout << "refPos>=ref.size()\t" << refPos << "\t" << ref.size()
//                      << std::endl;
//          if (queryPos >= query.size())
//            std::cout << "queryPos>=query.size()\t" << queryPos << "\t"
//                      << query.size() << std::endl;
          if (ref[refPos] == query[queryPos]) {
            numMatch++;
            sequenceDiff.logProbability += matchTable[quality[queryPos] - 33];
          } else {
            sequenceDiff.NM++;
            if (numMatch) MDcomponents.push_back(std::to_string(numMatch));
            MDcomponents.push_back(std::string(1, ref[refPos]));
            sequenceDiff.logProbability +=
                misMatchTable[quality[queryPos] - 33];
            numMatch = 0;
          }
          refPos++;
          queryPos++;
        }
        if (numMatch) {
          MDcomponents.push_back(std::to_string(numMatch));
        }
        break;
      case 1:  // I
        sequenceDiff.cigar.push_back('I');
        sequenceDiff.NM += length;
        queryPos += length;
        break;
      case 2:  // D
        sequenceDiff.cigar.push_back('D');
        MDcomponents.push_back("^");
        //        std::cout << "skipping over in ref ";
        for (int i = 0; i < length; i++) {
//          if (refPos >= ref.size())
//            std::cout << "refPos>=ref.size()\t" << refPos << "\t" << ref.size()
//                      << std::endl;
          temp.push_back(ref[refPos]);
          sequenceDiff.NM++;
          refPos++;
        }
        MDcomponents.push_back(temp);
        break;
      default:
        break;
    }
  }
  int end = query.size() - overlap.alignment.query_end - 1;
  if (end > 0) {
    auto length = end;
    sequenceDiff.cigar.append(std::to_string(length));
    sequenceDiff.cigar.push_back('S');
    queryPos += length;
  }
  bool ambiguous = false;
  for (auto it = MDcomponents.begin(); it != MDcomponents.end();) {
    if (*it == "^") {
      sequenceDiff.MD.append(*it);
      it++;
      sequenceDiff.MD.append(*it);
      ambiguous = true;
      it++;
    } else if (isdigit((*it)[0])) {
      int runningTotal = 0;
      while (it != MDcomponents.end() && isdigit((*it)[0])) {
        runningTotal += std::stoi(*it);
        it++;
      }
      sequenceDiff.MD.append(std::to_string(runningTotal));
      ambiguous = false;
    } else {
      if (ambiguous) {
        sequenceDiff.MD.append("0");
        ambiguous = false;
      }
      sequenceDiff.MD.append(*it);
      it++;
    }
    if (it == MDcomponents.end()) break;
  }
//  std::cout << sequenceDiff.cigar << "\t" << sequenceDiff.MD << std::endl;
//  std::cout << "q\t" << query << std::endl << "r\t";
//  for (int i = overlap.alignment.ref_begin; i <= overlap.alignment.ref_end; i++)
//    std::cout << ref[i];
//  std::cout << std::endl;
//  //  std::cout<<"r2\t"<<ref.substr(overlap.alignment.ref_begin,overlap.alignment.ref_end+1-overlap.alignment.ref_begin)<<std::endl;
//  std::cout << "rec\t" << recreateRef(query, sequenceDiff.cigar,
//                                      sequenceDiff.MD) << std::endl;
//
//  std::string refRec = recreateRef(query, sequenceDiff.cigar, sequenceDiff.MD);
//  //  auto it = std::remove(refRec.begin(),refRec.end(),'-');
//  //  refRec.resize(std::distance(refRec.begin(),it));
//  std::cout << "rec\t" << refRec << std::endl;
//  std::cout << std::boolalpha
//            << (refRec == ref.substr(overlap.alignment.ref_begin,
//                                     overlap.alignment.ref_end + 1 -
//                                         overlap.alignment.ref_begin))
//            << std::endl;
  return sequenceDiff;
}
class SAMEntry {
 public:
  template <typename FASTQType>
  void init(const Overlap &overlap, const std::vector<FASTQType> &reads,
            const GenbankIndex &genbankIndex);
  std::string getEntry();
  uint16_t getFlag();
  std::string qname;
  std::string rname;
  uint32_t pos = 0;
  uint8_t mapq = 255;
  std::string cigar="*";
  std::string rnext = "=";
  uint32_t pnext = 0;
  int32_t tlen = 0;
  std::string seq;
  std::string qual;
  //flags
  bool multipleSegments = false;
  bool allSegmentsAligned = false;
  bool thisSegmentUnmapped = false;
  bool nextSegmentUnmapped = false;
  bool revComp = false;
  bool nextRevComp = false;
  bool first = false;
  bool secondary = true;
  //Optional
  std::string MD;   //Mismatching positions
  uint16_t AS;      //Alignment score
  uint32_t NM;      //Edit distance
  //SLAM specific
  uint16_t XS;      //Score assigned by SLAM
  uint32_t XO = 0;  //num hits
  uint32_t XT = 0;  //taxonomy id
  std::string XG;   //gene
  std::string XP;   //protein id
  std::string XR;   //product
  std::string XA;   //alternate hits
  double prob = 0;
};
std::string SAMEntry::getEntry() {
  std::string out;
  out +=
      qname + '\t' + std::to_string(getFlag()) + '\t' + rname + '\t' +
      std::to_string(pos) + '\t' + std::to_string(mapq) + '\t' + (reportCigar ? cigar : "*") + '\t' +
      rnext + '\t' + std::to_string(pnext) + '\t' + std::to_string(tlen) +
      '\t' +
      //seq
      '*'
      + '\t' +
      //qual;
      '*';
  if(thisSegmentUnmapped)
    return out;
  if(reportCigar){
    out.append("\tMD:Z:");
    out += MD;
  }
  out.append("\tAS:i:");
  out += std::to_string(AS) + '\t' + "XS:i:" + std::to_string(XS) + '\t' +
      "NM:i:" + std::to_string(NM) + '\t' + "X0:i:" +
      std::to_string(XO);
  if (XT!=0) out += "\tXT:i:" + std::to_string(XT);
  if (XG.size()) out += "\tXG:Z:" + XG;
  if (XP.size()) out += "\tXP:Z:" + XP;
  if (XR.size()) out += "\tXR:Z:" + std::string("\"") + XR + std::string("\"");
//  if (XA.size() && SAMXA) out += "\tXA:Z:" + XA;
  return out;
}
uint16_t SAMEntry::getFlag() {
  uint16_t flag = 0;
  if (multipleSegments) flag |= 0x1;
  if (allSegmentsAligned) flag |= 0x2;
  if (thisSegmentUnmapped) flag |= 0x4;
  if (nextSegmentUnmapped) flag |= 0x8;
  if (revComp) flag |= 0x10;
  if (nextRevComp) flag |= 0x20;
  if(pairedData){
    if (first)
      flag |= 0x40;
    else
      flag |= 0x80;
  }
  if (secondary) flag |= 0x100;
  return flag;
}
std::string decodeFlag(uint16_t flag){
  std::string flagString;
  if (flag & 0x1) flagString.append("multipleSegments_");
  if (flag & 0x2) flagString.append("allSegmentsAligned_");
  if (flag & 0x4) flagString.append("thisSegmentUnmapped_");
  if (flag & 0x8) flagString.append("nextSegmentUnmapped_");
  if (flag & 0x10) flagString.append("revComp_");
  if (flag & 0x20) flagString.append("nextRevComp_");
  if (flag & 0x40)
    flagString.append("r1_");
  if (flag & 0x80)
    flagString.append("r2_");
  if (flag & 0x100) flagString.append("secondary_");
  return flagString;
}
template <typename FASTQType>
void SAMEntry::init(const Overlap &overlap, const std::vector<FASTQType> &reads,
                    const GenbankIndex &genbankIndex) {
  auto &entry = genbankIndex.entries[overlap.entryPosInArray];
  auto sequenceDiff = getCigarAndMD(overlap, reads, genbankIndex);
  cigar = sequenceDiff.cigar;
  MD = sequenceDiff.MD;
  NM = sequenceDiff.NM;
  prob = std::pow(10, sequenceDiff.logProbability);
  rname = entry.locusTag;
  pos = overlap.alignment.ref_begin + 1;
  AS = overlap.alignment.sw_score;
}
template <typename FASTQType>
std::pair<SAMEntry, SAMEntry> getSAMFromPair(
    const PairedOverlap &alignmentPair, const std::vector<FASTQType> &reads,
    const GenbankIndex &index) {
  auto &entry = index.entries[alignmentPair.entryPosInArray];
  SAMEntry r1;
  SAMEntry r2;
  r1.first = true;
  r2.first = false;
  const Gene *gene =
      entry.getGene(alignmentPair.refStart, alignmentPair.refEnd);
  if (gene != nullptr) {
    r1.XG = gene->geneName;
    r2.XG = gene->geneName;
    r1.XP = gene->proteinID;
    r2.XP = gene->proteinID;
    r1.XR = gene->product;
    r2.XR = gene->product;
  }
  r1.XT = entry.taxonomyID;
  r2.XT = entry.taxonomyID;
  bool conventionalSequence = true;
  bool bothAligned = alignmentPair.hasR1 && alignmentPair.hasR2;
  if (pairedData) {
    r1.multipleSegments = true;
    r2.multipleSegments = true;
  }
  if (bothAligned) {
    r1.allSegmentsAligned = true;
    r2.allSegmentsAligned = true;
    conventionalSequence = alignmentPair.r1Overlap.alignment.ref_begin <
                           alignmentPair.r2Overlap.alignment.ref_begin;
    if (alignmentPair.r1Overlap.revComp) {
      r1.revComp = true;
      r2.nextRevComp = true;
    }
    if (alignmentPair.r2Overlap.revComp) {
      r2.revComp = true;
      r1.nextRevComp = true;
    }
  } else if (alignmentPair.hasR1) {
    r1.nextSegmentUnmapped = true;
    r2.thisSegmentUnmapped = true;
    if (alignmentPair.r1Overlap.revComp) r1.revComp = true;
  } else if (alignmentPair.hasR2) {
    r2.nextSegmentUnmapped = true;
    r1.thisSegmentUnmapped = true;
    if (alignmentPair.r2Overlap.revComp) r2.revComp = true;
  }
  if (alignmentPair.hasR1)
    r1.init(alignmentPair.r1Overlap, reads, index);
  if (alignmentPair.hasR2)
    r2.init(alignmentPair.r2Overlap, reads, index);
  r1.pnext = r2.pos;
  r2.pnext = r1.pos;
  if (!alignmentPair.hasR1){
    r1.rname=r2.rname;
    r1.pos=r2.pos;
    r2.pnext=r2.pos;
    r1.pnext=r2.pos;
  }
  if (!alignmentPair.hasR2){
    r2.rname=r1.rname;
    r2.pos=r1.pos;
    r1.pnext=r1.pos;
    r2.pnext=r1.pos;
  }
  if(!pairedData){
    r1.rnext="*";
    r1.pnext=0;
    r1.nextSegmentUnmapped=false;
  }
  int32_t tlen = alignmentPair.refEnd - alignmentPair.refStart + 1;
  if(!(alignmentPair.hasR1||alignmentPair.hasR2))
    tlen=0;
  if (!conventionalSequence) tlen *= -1;
  r1.tlen = tlen;
  r2.tlen = tlen * -1;
  r1.XS = alignmentPair.combinedScore;
  r2.XS = alignmentPair.combinedScore;
  return { r1, r2 };
}
template <typename FASTQType>
void writeSAMOutputPairs(std::ofstream &outFile, ReadPairAndOverlaps &read,
                         const std::vector<FASTQType> &reads,
                         const GenbankIndex &genbankIndex);

void appendXA(const SAMEntry &entry, std::string &XA) {
  XA +=entry.rname + ',' + (entry.revComp ? '-' : '+') + std::to_string(entry.pos) + ',' +
      (reportCigar ? entry.cigar : "*") + ',' + std::to_string(entry.NM) + ';';
}
template <typename FASTQType>
inline void writeSAMOutputPairs(std::ofstream &outFile,
                                ReadPairAndOverlaps &read,
                                const std::vector<FASTQType> &reads,
                                const GenbankIndex &genbankIndex) {
  std::sort(read.alignmentPairs.begin(), read.alignmentPairs.end(),
            [](const PairedOverlap & i, const PairedOverlap & j) {
      return i.combinedScore > j.combinedScore;
  });
  std::vector<std::pair<SAMEntry, SAMEntry> > SAMPairs;
  uint32_t r1NumHits = 0;
  uint32_t r2NumHits = 0;
  for (auto &alignmentPair : read.alignmentPairs) {
    if (alignmentPair.hasR1) r1NumHits++;
    if (alignmentPair.hasR2) r2NumHits++;
    SAMPairs.push_back(getSAMFromPair(alignmentPair, reads, genbankIndex));
    if(SAMPairs.size()>=numSAMAlignments)
      break;
  }
  std::string r1XA, r2XA;
  auto primary = SAMPairs.begin();
  double r1SumProb = 0;
  double r2SumProb = 0;
  auto &r1Read = reads[read.r1PosInArray];
  auto &r2Read = reads[read.r2PosInArray];
  for (auto SAMPair = SAMPairs.begin(); SAMPair != SAMPairs.end(); SAMPair++) {
    SAMPair->first.qname=r1Read.sequenceIdentifier;
    SAMPair->second.qname=r2Read.sequenceIdentifier;
    SAMPair->first.seq = SAMPair->first.revComp ? reverseComplement(r1Read.bases) : r1Read.bases;
    SAMPair->first.qual = r1Read.quality.size() ? r1Read.quality : "*";
    if(SAMPair->first.revComp)
      std::reverse(SAMPair->first.qual.begin(),SAMPair->first.qual.end());
    SAMPair->second.seq = SAMPair->second.revComp ? reverseComplement(r2Read.bases) : r2Read.bases;
    SAMPair->second.qual = r2Read.quality.size() ? r2Read.quality : "*";
    if(SAMPair->second.revComp)
      std::reverse(SAMPair->second.qual.begin(),SAMPair->second.qual.end());
    r1SumProb += SAMPair->first.prob;
    r2SumProb += SAMPair->second.prob;
    SAMPair->first.XO = r1NumHits;
    SAMPair->second.XO = r2NumHits;
    if (SAMPair == primary) continue;
    if (!SAMPair->first.thisSegmentUnmapped)
      appendXA(SAMPair->first, r1XA);
    if (!SAMPair->second.thisSegmentUnmapped)
      appendXA(SAMPair->second, r2XA);
  }
  primary->first.secondary = false;
  primary->second.secondary = false;
  primary->first.XA = r1XA;
  primary->second.XA = r2XA;
  for (auto SAMPair2 = SAMPairs.begin(); SAMPair2 != SAMPairs.end(); SAMPair2++) {
    double temp = 1.0 - SAMPair2->first.prob / r1SumProb;
    if (temp <= 0.00001) temp = 0.00001;
    double temp2 = 1.0 - SAMPair2->second.prob / r2SumProb;
    if (temp2 <= 0.00001) temp2 = 0.00001;
    SAMPair2->first.mapq = ceil(-10.0 * std::log10(temp));
    SAMPair2->second.mapq = ceil(-10.0 * std::log10(temp2));
//    SAMPair2->first.mapq = floor(-10.0 * std::log10(temp) + 0.5);
//    SAMPair2->second.mapq = floor(-10.0 * std::log10(temp2) + 0.5);
//    if (!SAMPair2->first.thisSegmentUnmapped)
//      outFile << decodeFlag(SAMPair2->first.getFlag())<<"\t"<<SAMPair2->first.getEntry() << "\n";
      outFile << SAMPair2->first.getEntry() << "\n";
//    if (!SAMPair2->second.thisSegmentUnmapped)
      if(pairedData)
        outFile << SAMPair2->second.getEntry() << "\n";
//        outFile << decodeFlag(SAMPair2->second.getFlag())<<"\t"<<SAMPair2->second.getEntry() << "\n";
      if(SAMXA)
        break;
  }
}
inline std::string getHeader(const GenbankIndex & index){
  std::string header;
  header += "@HD\tVN:1.0\tSO:unsorted\n";
  for(auto & entry : index.entries){
    header += "@SQ\tSN:";
    header += entry.locusTag;
    header += "\tLN:";
    header += std::to_string(entry.bases.size());
    if(entry.taxonomyID){
      header += "\tSP:";
      header += std::to_string(entry.taxonomyID);
    }
    header += "\n";
  }
  header += "@PG\tID:SLAM\tPN:SLAM\tVN:1.0\tCL:\"";
  header += commandLine;
  header += "\"\n";
  return header;
}
}

#endif /* SAM_H_ */
