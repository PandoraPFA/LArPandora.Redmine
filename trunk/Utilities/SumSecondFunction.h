/// SumSecondFunction
/// 28-Jul-2009 William Seligman <seligman@nevis.columbia.edu>

/// When using STL, maps, and the std::accumulate function, there's one
/// function I keep coding over and over again: accumulate the second
/// member of a pair.
///
/// To save myself time (and others who know and use STL), here's a
/// complete STL-compatible implementation of that function.  To use
/// it, assume you have an object:
///
///    std::map<K,V> myMap;
///
/// To sum all the V's in the map:
///
/// V sum = std::accumulate( myMap.begin(),myMap.end(),V(),SumSecondFunction<K,V>() );
///
/// (Yes, I know there are other, better ways to do this, if one has
/// access to BOOST.  Unfortunately, we're not supposed to use BOOST in
/// LArSoft, since as of Jul-2009 it's not universally installed on
/// FNAL machines.)
///
#ifndef Utilities_SumSecondFunction_h
#define Utilities_SumSecondFunction_h

#include <functional>

namespace util {

  template < typename _Key, typename _Value, typename _BinaryOperation = std::plus<_Value> >
  class SumSecondFunction
    : public std::binary_function< _Value,
				   std::pair<_Key, _Value>,
				   _Value >
  {
  public:
    const _Value operator() ( const _Value& value,
			      const std::pair<_Key, _Value>& entry ) const
    {
      return _BinaryOperation()( value, entry.second );
    }
  };

} // namespace util

#endif // Utilities_SumSecondFunction_h


