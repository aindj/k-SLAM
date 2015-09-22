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
#ifndef PARALLELTOOLS_H_
#define PARALLELTOOLS_H_

#include <omp.h>
#include <vector>
namespace SLAM {
void forceParallel();
template <typename RAIter, typename Functor>
std::vector<RAIter> getStartPositions(const RAIter begin, const RAIter end,
                                      Functor breakAt);
template <typename RAIter, typename Functor, typename Functor2, typename T>
void parallelize(const RAIter begin, const RAIter end, T &returnContainer,
                 Functor functor, Functor2 breakAt);
template <typename RAIter, typename Functor, typename Functor2>
void parallelForEachWithSplit(const RAIter begin, const RAIter end,
                              Functor functor, Functor2 breakAt);
/*
 * Tells the gnu parallel libstdc++ to always use parallel mode, this overrides
 * the default serial mode of the for_each for ranges <1000 elements
 */
inline void forceParallel() {
  __gnu_parallel::_Settings s;
  s.algorithm_strategy = __gnu_parallel::force_parallel;
  __gnu_parallel::_Settings::set(s);
}

/*
 * Tries to break up a range for paralellisation
 * Only breaks a range between i and i+1 if breakAt(i,i+1)=true
 */
template <typename RAIter, typename Functor>
std::vector<RAIter> getStartPositions(const RAIter begin, const RAIter end,
                                      Functor breakAt) {
  int numThreads = omp_get_max_threads();
  auto numElements = std::distance(begin, end);
  size_t numElementsPerThread = numElements / numThreads;
  std::vector<RAIter> startPositions { begin }
  ;
  for (int i = 1; i < numThreads; i++) {
    if (numElementsPerThread == 0) {
      startPositions.push_back(begin);
      continue;
    }
    auto tempBegin = begin + numElementsPerThread * i - 1;
    while (!breakAt(tempBegin, (tempBegin + 1))) {
      tempBegin++;
    }
    startPositions.push_back(++tempBegin);
  }
  return startPositions;
}
template <typename RAIter>
void mergeVectors(std::vector<std::vector<RAIter> > &inVectors,
                  std::vector<RAIter> &returnVec) {
//  log("merge");
  size_t totalSize = 0;
  for (auto &vec : inVectors)
    totalSize += vec.size();
  returnVec.reserve(totalSize);
  for (auto &vec : inVectors)
    std::move(vec.begin(), vec.end(), std::back_inserter(returnVec));
//  log(std::to_string(totalSize));
}

/*
 * Function template for parallelisation.
 * This function is intended to parallelise functions which operate on a range
 * and return a container of results describing that range.
 * Parallelise takes as its input a range and a functor to apply on that range,
 * returning a container of results. The range is split such that the function
 * breakAt evaluates to true at the break-points between subranges
 * Example using a range of integers and a function which returns the
 * consecutive repeated elements:
 * std::vector<int> vec{1,1,2,3,4,4,5,5,5,6};
 * auto getConsecutiveRepeats [](std::vector<int>::iterator begin,
 * 								    std::vector<int>::iterator
 * end)
 * 								    {
 * 								    	std::vector<int>
 * consecutiveRepeats;
 * 								    	for(auto
 * it=begin;it!=(end-1);it++)
 * 								    		if(*it==*(it+1))
 * 								    			consecutiveRepeats.push_back(*it);
 * 								    }
 * auto notEqual [](std::vector<int>::iterator i,
 * 					std::vector<int>::iterator j){
 * 					return (*i)!=(*j);
 * 					}
 * std::vector<int> returnVec =
 * parallelize(vec.begin(),vec.end(),getConsecutiveRepeats,notEqual);
 * In this function, to get the correct repeated elements, the ranges must not
 * be split between equal consecutive elements
 */
template <typename RAIter, typename Functor, typename Functor2, typename T>
void parallelize(const RAIter begin, const RAIter end, T &returnContainer,
                 Functor functor, Functor2 breakAt) {
  int numThreads = omp_get_max_threads();
  std::vector<T> returnVec;
  returnVec.resize(numThreads);
  std::vector<std::thread *> threads;
  std::vector<RAIter> startPositions = getStartPositions(begin, end, breakAt);
  for (unsigned i = 0; i < startPositions.size(); i++) {
    auto tempBegin = startPositions[i];
    auto tempEnd =
        i == (startPositions.size() - 1) ? end : startPositions[i + 1];
    threads.push_back(
        new std::thread([i, tempBegin, tempEnd, &returnVec, &functor]() {
      returnVec[i] = functor(tempBegin, tempEnd);
    }));
  }

  for (auto &thread : threads)
    thread->join();
  mergeVectors(returnVec, returnContainer);
//  log("done");
}

/*
 * Applies "functor" to all elements in a range, only splitting the range
 * where breakAt=true (see "getStartPositions"
 */
template <typename RAIter, typename Functor, typename Functor2>
void parallelForEachWithSplit(const RAIter begin, const RAIter end,
                              Functor functor, Functor2 breakAt) {
  std::vector<std::thread *> threads;
  std::vector<RAIter> startPositions = getStartPositions(begin, end, breakAt);
  for (unsigned i = 0; i < startPositions.size(); i++) {
    auto tempBegin = startPositions[i];
    auto tempEnd =
        i == (startPositions.size() - 1) ? end : startPositions[i + 1];
    threads.push_back(new std::thread([tempBegin, tempEnd, &functor, i]() {
      functor(tempBegin, tempEnd);
    }));
  }
  for (auto &thread : threads)
    thread->join();
}
}
#endif /* PARALLELTOOLS_H_ */
