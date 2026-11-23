/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_tools_MergeVectorTools_h
#define __PLUMED_tools_MergeVectorTools_h

#include "small_vector/small_vector.h"
#include <vector>
#include <algorithm>
#include <type_traits>

namespace PLMD {

namespace mergeVectorTools {

/// Merge sorted vectors.
/// Takes a vector of pointers to containers and merge them.
/// Containers should be already sorted.
/// The content is appended to the result vector.
/// Optionally, uses a priority_queue implementation.
template<class C>
static void mergeSortedVectors(const C* const* vecs, std::size_t size, std::vector<typename C::value_type> & result) {

  /// local class storing the range of remaining objects to be pushed
  class Entry {
    typename C::const_iterator fwdIt,endIt;

  public:
    explicit Entry(C const& v) : fwdIt(v.begin()), endIt(v.end()) {}
    /// check if this vector still contains something to be pushed
    bool empty() const {
      return fwdIt == endIt;
    }
    // TODO: revert "typename C::value_type" to "auto": nvc++ and icpc seems to do not deduce automatically the return type
    const typename C::value_type & top() const {
      return *fwdIt;
    }
    /// to allow using a priority_queu, which selects the highest element.
    /// we here (counterintuitively) define < as >
    bool operator< (Entry const& rhs) const {
      return top() > rhs.top();
    }
    void next() {
      ++fwdIt;
    };
  };

  // preprocessing
  {
    std::size_t maxsize=0;
    for(unsigned i=0; i<size; i++) {
      // find the largest
      maxsize=std::max(maxsize,vecs[i]->size());
    }
    // this is just to decrease the number of reallocations on push_back
    result.reserve(maxsize);
    // if vectors are empty we are done
    if(maxsize==0) {
      return;
    }
  }

  // A few vectors (e.g. the same large group requested by a few actions) are merged two at a time,
  // which is much faster than the heap below. The cost grows with the number of vectors, so many
  // vectors (e.g. many actions with a few atoms each) use the heap.
  constexpr std::size_t maxPairwise=8;
  if(size<=maxPairwise) {
    using T=typename C::value_type;
    // acc holds the elements merged so far, sorted and without repetitions
    std::vector<T> acc,tmp;
    for(std::size_t i=0; i<size; i++) {
      tmp.clear();
      tmp.reserve(acc.size()+vecs[i]->size());
      auto a=acc.cbegin();
      const auto aend=acc.cend();
      auto b=vecs[i]->begin();
      const auto bend=vecs[i]->end();
      auto push=[&tmp](const T& val) {
        if(tmp.empty() || val > tmp.back()) {
          tmp.push_back(val);
        }
      };
      while(a!=aend && b!=bend) {
        if(*b < *a) {
          push(*b);
          ++b;
        } else {
          push(*a);
          ++a;
        }
      }
      for(; a!=aend; ++a) {
        push(*a);
      }
      for(; b!=bend; ++b) {
        push(*b);
      }
      acc.swap(tmp);
    }
    result.insert(result.end(),acc.begin(),acc.end());
    return;
  }

  // start
  // heap containing the ranges to be pushed
  // avoid allocations when it's small
  gch::small_vector<Entry,32> heap;

  {
    for(unsigned i=0; i<size; i++) {
      if(!vecs[i]->empty()) {
        // add entry at the end of the array
        heap.emplace_back(Entry(*vecs[i]));
      }
    }
    // make a sorted heap
    std::make_heap(std::begin(heap),std::end(heap));
  }

  // first iteration, to avoid an extra "if" in main iteration
  // explanations are below
  {
    std::pop_heap(std::begin(heap), std::end(heap));
    auto & tmp=heap.back();
    const auto val=tmp.top();
    // here, we append inconditionally
    result.push_back(val);
    tmp.next();
    if(tmp.empty()) {
      heap.pop_back();
    } else {
      std::push_heap(std::begin(heap), std::end(heap));
    }
  }

  while(!heap.empty()) {
    // move entry with the smallest element to the back of the array
    std::pop_heap(std::begin(heap), std::end(heap));
    // entry
    auto & tmp=heap.back();
    // element
    const auto val=tmp.top();
    // if the element is larger than the current largest element,
    // push it to result
    if(val > result.back()) {
      result.push_back(val);
    }
    // move forward the used entry
    tmp.next();
    // if this entry is exhausted, remove it from the array
    if(tmp.empty()) {
      heap.pop_back();
    }
    // otherwise, sort again the array
    else {
      std::push_heap(std::begin(heap), std::end(heap));
    }
  }
}

template<typename T, typename = void>
struct has_size_and_data : std::false_type {};

template<typename T>
struct has_size_and_data<T, std::void_t<decltype(std::declval<T>().size()), decltype(std::declval<T>().data())>> : std::true_type {};

template<class C, class D>
auto mergeSortedVectors(C& vecs, std::vector<D> & result) -> typename std::enable_if<has_size_and_data<C>::value, void>::type {
  mergeSortedVectors(vecs.data(),vecs.size(),result);
}

}
}

#endif

