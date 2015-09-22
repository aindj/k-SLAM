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

/*
 * Contains tools for file read and string manipulation
 */
#ifndef SEQUENCETOOLS_H_
#define SEQUENCETOOLS_H_
#include <string>
#include <fstream>
#include <stdexcept>
#include <iostream>
#include <cstdio>
#include <algorithm>
#include <chrono>
#include <iomanip>
namespace SLAM {

template <typename T> double executionTime(const T startTime);
void log(const std::string &line);
inline void inPlaceConvertToUpperCase(std::string &input) {
  for (auto &c : input) {
    c = toupper(c);
  }
}
/*
 * Function used to read a line from a file which is in either Windows or Unix
 * format (newlines)
 */
inline std::istream &safeGetline(std::istream &is, std::string &t) {
  t.clear();

  // The characters in the stream are read one-by-one using a std::streambuf.
  // That is faster than reading them one-by-one using the std::istream.
  // Code that uses streambuf this way must be guarded by a sentry object.
  // The sentry object performs various tasks,
  // such as thread synchronization and updating the stream state.

  std::istream::sentry se(is, true);
  std::streambuf *sb = is.rdbuf();

  for (;;) {
    int c = sb->sbumpc();
    switch (c) {
      case '\n':
        return is;
      case '\r':
        if (sb->sgetc() == '\n') sb->sbumpc();
        return is;
      case EOF:
        // Also handle the case when the last line has no line ending
        if (t.empty()) is.setstate(std::ios::eofbit);
        return is;
      default:
        t += (char) c;
    }
  }
}
/*
* Returns the reverse complement of the input string
*/
inline std::string reverseComplement(const std::string &forward) {
  std::string rc = forward;
  std::reverse(rc.begin(), rc.end());
  for_each(rc.begin(), rc.end(), [](char & i) {
    switch (i) {
      case 'A':
        i = 'T';
        break;
      case 'C':
        i = 'G';
        break;
      case 'T':
        i = 'A';
        break;
      case 'G':
        i = 'C';
        break;
    }
  });
  return rc;
}
inline void inPlaceReverseComplement(std::string &bases) {
  std::reverse(bases.begin(), bases.end());
  for_each(bases.begin(), bases.end(), [](char & i) {
    switch (i) {
      case 'A':
        i = 'T';
        break;
      case 'C':
        i = 'G';
        break;
      case 'T':
        i = 'A';
        break;
      case 'G':
        i = 'C';
        break;
    }
  });
}
std::vector<std::string> tokenise(const std::string &line,
                                  const std::string &delimiters) {
  std::vector<std::string> tokens;
  // Skip delimiters at beginning.
  std::string::size_type lastPos = line.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  std::string::size_type pos = line.find_first_of(delimiters, lastPos);
  while (std::string::npos != pos || std::string::npos != lastPos) {
    // Found a token, add it to the vector.
    tokens.push_back(line.substr(lastPos, pos - lastPos));
    // Skip delimiters.  Note the "not_of"
    lastPos = line.find_first_not_of(delimiters, pos);
    // Find next "non-delimiter"
    pos = line.find_first_of(delimiters, lastPos);
  }
  return tokens;
}
class Date {
 public:
  Date(const unsigned day, const unsigned month, const unsigned year)
      : day(day), month(month), year(year) {}
  ;
  Date() {}
  ;
  unsigned day = 0;
  unsigned month = 0;
  unsigned year = 0;
  bool operator<(const Date &other) const {
    if (this->year == other.year) {
      if (this->month == other.month)
        return this->day < other.day;
      else
        return this->month < other.month;
    } else
      return this->year < other.year;
  }
};
class Log {
  std::ofstream logFile;
  decltype(std::chrono::system_clock::now()) startTime;

 public:
  Log(const std::string logFileName);
  void write(const std::string &entry);
  void reset() { startTime = std::chrono::system_clock::now(); }
};
Log::Log(const std::string logFileName) {
  logFile.open(logFileName);
  if (!logFile.good()) throw std::runtime_error("unable to open log file");
  startTime = std::chrono::system_clock::now();
  logFile.setf(std::ios::fixed, std::ios::floatfield);
  logFile.precision(2);
}

void Log::write(const std::string &entry) {
  logFile << "[t = " << executionTime(startTime) << "s]\t" << entry
          << std::endl;
}
void log(const std::string &line) {
  static Log L("log.txt");
  if (line == "reset") L.reset();
  L.write(line);
}
/*
 * Returns number of seconds elapsed since the start of execution
 */
template <typename T> inline double executionTime(const T startTime) {
  const auto elapsedTime = std::chrono::system_clock::now() - startTime;
  return std::chrono::duration_cast<std::chrono::milliseconds>(elapsedTime)
             .count() * 1.0 / 1000;
}
}

/*
 * Converts a number into binary
 */
template <typename T> inline std::string getBinary(T input) {
  std::string out;
  unsigned length = sizeof(T) * 8;
  out.resize(length);
  for (int i = length - 1; i >= 0; i--) {
    out[i] = (input & 0b1) ? '1' : '0';
    input >>= 1;
  }
  return out;
}
// uint32_t singleGapAlign(const std::string &ref, const std::string &query,
//                        int32_t relPos) {
//  struct Seed {
//    uint32_t startPos = 0;
//    uint32_t endPos = 0;
//    int32_t score = 0;
//    int32_t offset = 0;
//  };
//  int32_t maxScore = 0;
//  int32_t runningScore = 0;
//  std::vector<Seed> seeds;
//  int32_t score = 0;
//  Seed seed;
//  for (int32_t i = 0; i < query.size(); i++) {
//    int32_t refIndex = relPos + i;
//    if (refIndex < 0) continue;
//    if (refIndex >= ref.size()) break;
//    if (ref[refIndex] == query[i])
//      runningScore += 2;
//    else
//      runningScore -= 3;
//    if (runningScore > maxScore) {
//      maxScore = runningScore;
//      seed.endPos = i;
//    }
//  }
//  for (int32_t i = seed.endPos - 1; i >= 0; i--) {
//    int32_t refIndex = relPos + i;
//    if (refIndex < 0) break;
//    if (refIndex >= ref.size()) continue;
//    if (ref[refIndex] == query[i])
//      runningScore += 2;
//    else
//      runningScore -= 3;
//    if (runningScore > maxScore) {
//      maxScore = runningScore;
//      seed.startPos = i;
//    }
//  }
//  return score > 0 ? score : 0;
//}
// uint32_t ungappedAlign(const std::string &ref, const std::string &query,
//                       int32_t relPos) {
//  int32_t score = 0;
//  for (int32_t i = 0; i < query.size(); i++) {
//    int32_t refIndex = relPos + i;
//    if (refIndex < 0) continue;
//    if (refIndex >= ref.size()) break;
//    if (ref[refIndex] == query[i])
//      score += 2;
//    else
//      score -= 3;
//  }
//  //  std::cout<<"Score\t"<<score<<"\n";
//  return score > 0 ? score : 0;
//}
// uint32_t ungappedAlignRevComp(const std::string &ref, const std::string
// &query,
//                              int32_t relPos) {
//  int32_t score = 0;
//  for (int32_t i = 0; i < query.size(); i++) {
//    int32_t refIndex = relPos + i;
//    if (refIndex < 0) continue;
//    if (refIndex >= ref.size()) break;
//    switch (ref[refIndex]) {
//      case 'A':
//        if (query[query.size() - 1 - i] == 'T')
//          score += 2;
//        else
//          score -= 3;
//        break;
//      case 'C':
//        if (query[query.size() - 1 - i] == 'G')
//          score += 2;
//        else
//          score -= 3;
//        break;
//      case 'T':
//        if (query[query.size() - 1 - i] == 'A')
//          score += 2;
//        else
//          score -= 3;
//        break;
//      case 'G':
//        if (query[query.size() - 1 - i] == 'C')
//          score += 2;
//        else
//          score -= 3;
//        break;
//      default:
//        score -= 3;
//        break;
//    }
//  }
//  //  std::cout<<"Score\t"<<score<<"\n";
//  return score > 0 ? score : 0;
//}

#endif /* SEQUENCETOOLS_H_ */
