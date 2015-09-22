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
#ifndef GLOBALS_H_
#define GLOBALS_H_
#include <inttypes.h>

// Type used to store kMer
typedef uint64_t KMerInt;
// The default length of a kMer
const unsigned k = 32;
// Smith Waterman alignment parameters
uint32_t match = 0;
uint32_t misMatch = 0;
uint32_t gapOpen = 0;
uint32_t gapExtend = 0;
uint32_t scoreThreshold = 0;

uint32_t numSAMAlignments=10;

bool performPseudoAssembly = true;
bool reportCigar = false;
bool pairedData = true;
bool SAMXA=false;
bool justAlign=false;
std::string commandLine;
//uint32_t numGBKMers=100;
double scoreFractionThreshold=0.95;


#endif /* GLOBALS_H_ */
