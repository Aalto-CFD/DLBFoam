#include <iterator>

namespace Foam{


template<class InputIt, class UnaryPredicate>
typename std::iterator_traits<InputIt>::difference_type count_while(InputIt first, InputIt last, UnaryPredicate predicate){

    typename std::iterator_traits<InputIt>::difference_type ret = 0;
    for (auto it = first; it != last; ++it){
      if(!predicate(*it)){
          break;
      }
      ret++; 
  }
  return ret;

}

}