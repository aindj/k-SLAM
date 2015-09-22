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
#ifndef LOOKUPTABLE_H_
#define LOOKUPTABLE_H_
#include <fstream>
#include <stdexcept>
#include <memory>
#include <cstring>

namespace SLAM {
/*
 * A class implementing a simple boolean lookup table
 * The table stores "size" on/off bits
 * The type "T" is the unsigned integer type of the
 * array index (recommended to use uint32 for efficiency)
 * The table's state can be read/written from/to a file
 * The datatype used to store the array can be modified
 * using the tableType typedef
 * The number of allocated "tableType"s is stored in
 * the sizeOfArray variable
 */
typedef uint32_t tableType;
constexpr unsigned numBits = sizeof(tableType) * 8;
constexpr tableType tableOne = 1;

template <typename T> class LookupTable {
  tableType *table;
  T max = 0;
  T sizeOfArray = 0;

 public:
  LookupTable(const T size);
  LookupTable(const T size, const char *inFileName);
  ;
  ~LookupTable() { delete[] table; }
  bool lookup(const T position) const;
  void set(const T position);
  void writeToFile(const char *fileName) const;
  bool operator==(const LookupTable<T> &rhs) const {
    if (rhs.sizeOfArray != this->sizeOfArray) return false;
    for (unsigned i = 0; i < sizeOfArray; i++) {
      if (this->table[i] != rhs.table[i]) return false;
    }
    return true;
  }
};

/*
 * Initializes a zeroed array of "size" bits
 */
template <typename T> LookupTable<T>::LookupTable(const T max) : max(max) {
  //  if (size == 0) throw std::runtime_error("Lookup Table Size 0");
  sizeOfArray = 1 + max / numBits;
  table = new tableType[sizeOfArray]();
}

/*
 * Sets bit at location "position" to 1
 */
template <typename T> void LookupTable<T>::set(const T position) {
  table[position / numBits] |= (tableOne << (position % numBits));
}

/*
 * Returns content of bit at "position"
 */
template <typename T>
inline bool LookupTable<T>::lookup(const T position) const {
  return table[position / numBits] & (tableOne << (position % numBits));
}

/*
 * Writes memory of "table" to file
 */
template <typename T>
void LookupTable<T>::writeToFile(const char *fileName) const {
  std::ofstream outFile(fileName, std::ios::out | std::ios::binary);
  outFile.write((char *)table, sizeof(tableType) * sizeOfArray);
  outFile.close();
}

/*
 * Initialises table from file
 * Note that saved table must have the same tableType, T and size
 */
template <typename T>
LookupTable<T>::LookupTable(const T max, const char *fileName)
    : max(max) {
  std::ifstream inFile(fileName,
                       std::ios::in | std::ios::binary | std::ios::ate);
  if (!inFile.good()) throw std::runtime_error("unable to open table file");
  std::ifstream::pos_type fileSize = inFile.tellg();
  sizeOfArray = fileSize / sizeof(tableType);
  table = new uint32_t[sizeOfArray];
  inFile.seekg(0, std::ios::beg);
  inFile.read((char *)table, sizeof(tableType) * sizeOfArray);
  inFile.close();
}
}

#endif /* LOOKUPTABLE_H_ */
