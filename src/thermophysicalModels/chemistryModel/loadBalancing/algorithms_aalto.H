#pragma once

#include <iterator>





namespace Foam{


///
///@brief Converts an iterator to a reverse iterator 
///
///@tparam InputIt input iterator type
///@param i the iterator to convert
///@return std::reverse_iterator<Iter> a reverse iterator 
///        with ++ and -- operators shifting to left and right, respectively.
///
template<class InputIt>
std::reverse_iterator<InputIt> make_reverse(InputIt i){
    return std::reverse_iterator<InputIt>(i);
}

///
///@brief Computest the number of sequential elements for which the predicate is true.
///
///@tparam InputIt Random access iterator
///@tparam UnaryPredicate Unary predicate type
///@param first iterator to the first element to test
///@param last  iterator to the last element to test
///@param predicate the predicate to test
///@return std::iterator_traits<InputIt>::difference_type the number of sequential elements
///
template<class InputIt, class UnaryPredicate>
typename std::iterator_traits<InputIt>::difference_type count_while(InputIt first, InputIt last, UnaryPredicate predicate){

    typename std::iterator_traits<InputIt>::difference_type ret = 0;
    for (; first != last; ++first){
      if(!predicate(*first)){
          break;
      }
      ret++; 
  }
  return ret;

}

}