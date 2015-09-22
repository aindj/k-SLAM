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
#ifndef METAGENOMICFASTQSEQUENCE_H_
#define METAGENOMICFASTQSEQUENCE_H_
#include <string>
namespace SLAM {

/*
 * Represents a FASTQ sequence but with all the associated metadata
 * needed for metagenomic analysis
 */
class MetagenomicFASTQSequence : public FASTQSequence {
 public:
  MetagenomicFASTQSequence(std::string seqID, std::string inBases,
                           std::string inQuality)
      : FASTQSequence(seqID, inBases, inQuality) {}
  ;
  MetagenomicFASTQSequence() {}
  ;
  bool isHost = false;
  bool isLowComplexity = false;
};
}
#endif /* METAGENOMICFASTQSEQUENCE_H_ */
