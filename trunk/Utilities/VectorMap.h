/// $Id: SConscript,v 1.3 2010/05/18 21:47:28 kutschke Exp $
/// VectorMap
/// 17-Apr-2008 William Seligman <seligman@nevis.columbia.edu>
///
/// This class is an implementation of a concept discussion in
/// "Effective STL" by Scott Meyers:
///
/// STL maps are useful because their contents are always sorted, so
/// they're effective for fast searches.  However, in almost every
/// other respect vectors are superior: They take up less space, and
/// they use random-access iterators.
///
/// This class implements "sorted vector maps," that is, an STL-style
/// map implemented as a sorted STL vector of pairs.  I've done my best
/// to implement all aspects of the std::map interface in this class,
/// with some additions; if you've defined the following:
///
/// VectorMap<key_type, data_type> svm;
///
/// - svm(i) will return the "i-th" value in the map; that is, "i" is a
///   numeric index instead of a key.  (Note the use of parenthesis
///   instead of square brackets.)  This is a boon to physicists, most
///   of whom couldn't tell an iterator from a hole in the wall.
///
/// - svm.Key(i) will return the "i-th" key in the map.
///
/// - svm.Data(i) will return the same result as svm(i).
///
/// - svm[key_type] will now return the corresponding data_type in both
///   const and non-const contexts.  However, if you ask for svm[key]
///   and the key isn't in the map, and you're in a const context, the
///   routine will throw an out-of-range exception.
///
/// INCREDIBLY IMPORTANT NOTE: The "key type" of a VectorMap
/// cannot be a "const" type (unlike maps); it won't even compile.
/// When you do an insert, the underlying vector has to move things
/// around within its list, and it uses the assignment operator=() (or
/// "vector{i+1)=vector(i)" if you like).  You can't do that if either
/// the key or the data is const.
///
/// As with a map, there's no way to insert items at a specific
/// location in a VectorMap.  The insertion methods (including
/// operator[]) all operate on a sorted sequence according to the key.
/// Because of this, insertions take a long time.
///
/// However, for our processing, this doesn't matter much; for almost
/// all our maps, we typically have:
///
/// - Initialization, where the time doesn't matter (e.g., tracks in a
///   Monte Carlo).
///
/// - Access, where efficient or "simple" access to the map's contents
///   are important.  In general, we access a map many, many more times
///   that we create one.
///
/// - After we create/initialize a map, we never change its contents.
///
/// For this usage, a sorted vector is generally superior to a map.
///
/// This class just implements the equivalent of an STL map, not a
/// multimap, set, nor a multiset.  If there's a need, I may create
/// additional classes.
///
///
/// Is there any map feature that's not implemented?  Yes:
///
/// - equal_range in a const context (which causes some weird ROOT
///   dictionary problem); this isn't likely to be used for a map
///   anyway (multimaps or multisets would be a different story).
///
/// Advanced implementation note: Depending on the application, it
/// might be possible to speed up this class by using "lazy
/// evaluation"; that is, we wouldn't actually sort the vector until
/// the user actually tries to access its contents.  I'm not going to
/// do this, because:
///
/// A) I don't think my programming skills are up to the task.
///
/// B) In the primary application for which I plan to use this class
///    (Monte-Carlo particle tracks), we're performing at least one
///    search after every insert; lazy evaluation wouldn't be much of a
///    speed improvement.
///
/// Feb-2011 WGS: VectorMap mostly looks like a map, but there are some
/// memory-management issues that relate to it being a vector. Include
/// the vector-based routines reserve() and capacity().

#ifndef Utilities_VectorMap_H
#define Utilities_VectorMap_H

#include <cstddef>
#include <vector>
#include <map>
#include <functional>
#include <algorithm>

namespace util {

  // Start with the same template terms as used for a map (copied from
  // GNU C++'s stl_map.h).  In general, if you see variables that begin
  // with underscores ("_"), it means I copied the code directly from
  // GNU C++, either from stl_map.h or stl_vector.h.
  template < typename _Key, typename _Tp, typename _Compare = std::less<_Key> > 
  class VectorMap {
  public:
    // Define lots of handy types.

    typedef _Key                                       key_type;
    typedef _Tp                                        mapped_type;
    typedef std::pair<_Key, _Tp>                       value_type;
    typedef _Compare                                   key_compare;
    typedef std::allocator< std::pair<_Key, _Tp> >     allocator_type;

    typedef std::vector<value_type>                    vector_type;

    typedef typename vector_type::pointer              pointer;
    typedef typename vector_type::const_pointer        const_pointer;
    typedef typename vector_type::reference            reference;
    typedef typename vector_type::const_reference      const_reference;
    typedef typename vector_type::iterator             iterator;
    typedef typename vector_type::const_iterator       const_iterator;
    typedef std::reverse_iterator<const_iterator>      const_reverse_iterator;
    typedef std::reverse_iterator<iterator>            reverse_iterator;
    typedef size_t                                     size_type;
    typedef ptrdiff_t                                  difference_type;

  public:

    // The "value_compare" class is an internal utility class for
    // VectorMap.  It extends the "_Compare" template entry (the
    // name of a class that compares the keys) to a binary operator that
    // compares two entries in the "sortedVectorMap" private member.

    // Note that because of the context-based multiple definitions of
    // operator(), this class cannot inherit from the STL
    // binary_function template.  This means that it's not "adaptable";
    // e.g., you can't use it as an argument to bind2nd.

    // If you're getting confused by all this, just think of it as a
    // fancy way of defining "less than" and ignore it.
    class value_compare
    {
      friend class VectorMap<_Key, _Tp, _Compare>;
    protected:

      value_compare(_Compare __c = _Compare())
	: comp(__c) { }
      _Compare GetCompare() const {  return comp; }
      _Compare comp;
    
    public:

      // From an idea suggested by item 23 in "Effective STL" by Scott
      // Meyers:
      bool operator()(const value_type& __x, 
		      const value_type& __y) const
      { return keyLess(__x.first, __y.first); }
      bool operator()(const value_type& __x, 
		      const key_type& __y) const
      { return keyLess(__x.first, __y); }
      bool operator()(const key_type& __x, 
		      const value_type& __y) const
      { return keyLess(__x, __y.first); }
    private:
      bool keyLess(const key_type& __x, 
		   const key_type& __y) const
      { return comp(__x, __y); }
    };

  private:
    // The vector that contains the sorted pair<Key,Value> entries.
    vector_type sortedVectorMap; // The sorted <key,data> pairs.

    // The actual key-comparison object.
    value_compare valueCompare; //! Don't write this to the ROOT file.
  
  public:
    // After copying a lot of stuff from GNU C++, here's where I start
    // implementing some methods on my own.  I'm trying to be complete,
    // and implementing all the methods that STL map supports, even if I
    // don't think I'll ever use them.

    allocator_type get_allocator() const
    {
      return sortedVectorMap.get_allocator();
    }
 
    iterator begin()
    {
      return sortedVectorMap.begin();
    }

    const_iterator begin() const
    {
      return sortedVectorMap.begin();
    }

    iterator end()
    {
      return sortedVectorMap.end();
    }

    const_iterator end() const
    {
      return sortedVectorMap.end();
    }

    reverse_iterator rbegin()
    {
      return sortedVectorMap.rbegin();
    }

    const_reverse_iterator rbegin() const
    {
      return sortedVectorMap.rbegin();
    }

    reverse_iterator rend()
    {
      return sortedVectorMap.rend();
    }

    const_reverse_iterator rend() const
    {
      return sortedVectorMap.rend();
    }

    bool empty() const
    {
      return sortedVectorMap.empty();
    }

    size_t size() const
    {
      return sortedVectorMap.size();
    }

    size_t max_size() const
    {
      return sortedVectorMap.max_size();
    }

    mapped_type& operator[](const key_type& key)
    {
      // Do a binary search for the key.
      iterator i = lower_bound(key);

      // If the key is not found, or i->first is less than key, then
      // the key is not found.
      if (i == end() || key_comp()(key, (*i).first))
	// Insert this key into the map, with a default value.  Thanks
	// to the lower_bound call above, i already points to correct
	// place to insert the value in the sorted vector.
	i = sortedVectorMap.insert(i, value_type(key, mapped_type() ));

      return (*i).second;
    }

    // Not part of STL, even in a GNU extension, but something I always
    // wanted: operator[] in a const context.  Derived from the at()
    // method below.
    const mapped_type& operator[](const key_type& __k) const
    {
      const_iterator __i = lower_bound(__k);
      if (__i == end() || key_comp()(__k, (*__i).first))
	std::__throw_out_of_range(__N("Utilities/VectorMap::operator[]"));
      return (*__i).second;
    }

    // at(), equivalent to operator[] except that it can throw an
    // exception.
    mapped_type& at(const key_type& __k)
    {
      iterator __i = lower_bound(__k);
      if (__i == end() || key_comp()(__k, (*__i).first))
	std::__throw_out_of_range(__N("Utilities/VectorMap::at"));
      return (*__i).second;
    }

    const mapped_type& at(const key_type& __k) const
    {
      const_iterator __i = lower_bound(__k);
      if (__i == end() || key_comp()(__k, (*__i).first))
	std::__throw_out_of_range(__N("Utilities/VectorMap::at"));
      return (*__i).second;
    }

    // One of the key members of this adapted class.  Attempt to insert
    // the entry (a pair<key,value>) into the map.  Since this is a
    // unique insert, the map will only be changed if the key is not
    // already present.  In the returned pair, the iterator points to
    // the map entry that contains the key; the bool indicates whether
    // the insertion succeeded or failed.  Note the combination of a
    // binary search (lower_bound) with an insert-in-the-middle
    // ("sortedVectorMap.insert()"); that what's makes this method so
    // slow.
    std::pair<iterator,bool> insert(const value_type& entry)
    {
      const key_type& key = entry.first;
      iterator i = lower_bound(key);
      if (i == end() || key_comp()(key, (*i).first))
	{
	  // The entry was not found.  In that case, lower_bound has
	  // already found the correct point at which we want to insert
	  // the entry to maintain the sort.
	  i = sortedVectorMap.insert(i, entry);
	  return std::make_pair( i, true);
	}

      // If we get here, the entry was found.  Return failure.
      return std::make_pair( i, false );
    }

    // In this form of insert(), the iterator in the argument is
    // supposed to give us a hint about where to insert the item.
    // Actually, this hint is useless (since lower_bound doesn't take a
    // hint <heh, heh>).  So just do a regular insert instead.
    iterator insert(iterator, const value_type& entry)
    { 
      std::pair<iterator,bool> result = insert(entry);
      return result.first;
    }

    // Mass insertion.
    template <typename _InputIterator>
    void insert(_InputIterator __first, _InputIterator __last)
    { 
      for ( ; __first != __last; ++__first)
	insert(*__first);
    }

    void erase(iterator __position)
    { 
      sortedVectorMap.erase(__position);
    }

    // Erases all the entries with the key, and returns the number of
    // erasures.
    size_type erase(const key_type& key)
    { 
      iterator i = find(key);
      if ( i == end() ) return 0;
      erase(i);
      return 1;
    }

    // Erase a range.
    void erase(iterator __first, iterator __last)
    { 
      sortedVectorMap.erase(__first, __last); 
    }

    // Swap two maps. For VectorMap, this is pretty simple: use
    // the standard vector mechanism for swapping the vector portion of
    // the maps, then swap the value_compare objects (if any) by the
    // usual "temp" method.
    void swap(VectorMap& other)
    { 
      sortedVectorMap.swap(other.sortedVectorMap);
      value_compare temp(valueCompare);
      valueCompare = other.valueCompare;
      other.valueCompare = temp;
    }

    void clear()
    { 
      sortedVectorMap.clear(); 
    }

    // Returns the key-comparison object used for this VectorMap.
    key_compare key_comp() const
    { 
      return valueCompare.GetCompare(); 
    }

    // Returns the value-comparison object, which just compares the
    // entry.first of the two pairs.
    value_compare value_comp() const
    { 
      return valueCompare; 
    }

    iterator find(const key_type& key)
    { 
      iterator i = lower_bound(key); 
      if (i == end() || key_comp()(key, (*i).first))
	return end();

      return i;
    }

    const_iterator find(const key_type& key) const
    { 
      const_iterator i = lower_bound(key); 
      if (i == end() || key_comp()(key, (*i).first))
	return end();

      return i;
    }

    // Count the number of entries with a given key (which will be 0 or
    // 1 for a map).
    size_type count(const key_type& __x) const
    { 
      return find(__x) == end() ? 0 : 1; 
    }

    iterator lower_bound(const key_type& __x)
    { 
      return std::lower_bound(begin(),end(),__x,valueCompare); 
    }

    const_iterator lower_bound(const key_type& __x) const
    { 
      return std::lower_bound(begin(),end(),__x,valueCompare); 
    }

    iterator upper_bound(const key_type& __x)
    { 
      return std::upper_bound(begin(),end(),__x,valueCompare); 
    }

    const_iterator upper_bound(const key_type& __x) const
    { 
      return std::upper_bound(begin(),end(),__x,valueCompare); 
    }

    std::pair<iterator, iterator> equal_range(const key_type& key)
    { 
      return std::equal_range(begin(),end(),key,valueCompare); 
    }

    // The following code does not compile when processed by ROOT's
    // dictionary.  For now, this is not a big deal; no one is likely to
    // use equal_range on a map anyway.  If we ever have to implement a
    // multimap or multiset using the same scheme, this could be a
    // problem.
    /*
      std::pair<const_iterator, const_iterator> equal_range(const key_type& key) const
      { 
      return std::equal_range(begin(),end(),key,valueCompare); 
      }
    */

    // My own little extras, as described near the top of this header
    // file's comments:

    const mapped_type& operator()(const size_type& index) const
    {
      return sortedVectorMap[index].second;
    }

    const mapped_type& Data(const size_type& index) const
    {
      return sortedVectorMap[index].second;
    }

    const key_type& Key(const size_type& index) const
    {
      return sortedVectorMap[index].first;
    }


    // Vector-based memory management.
    void reserve( size_type i )
    {
      sortedVectorMap.reserve(i);
    }
    size_type capacity()
    {
      return sortedVectorMap.capacity();
    }

  
    // Friend definitions for comparison operators.
    template <typename _K1, typename _T1, typename _C1>
    friend bool operator== (const VectorMap<_K1, _T1, _C1>&,
			    const VectorMap<_K1, _T1, _C1>&);

    template <typename _K1, typename _T1, typename _C1>
    friend bool operator< (const VectorMap<_K1, _T1, _C1>&,
			   const VectorMap<_K1, _T1, _C1>&);

  public:

  };

} // namespace util

namespace util {

  // Comparison operators.
  template <typename _Key, typename _Tp, typename _Compare>
  inline bool operator==(const VectorMap<_Key, _Tp, _Compare>& __x,
			 const VectorMap<_Key, _Tp, _Compare>& __y)
  { 
    return __x.sortedVectorMap == __y.sortedVectorMap; 
  }

  template <typename _Key, typename _Tp, typename _Compare>
  inline bool operator<(const VectorMap<_Key, _Tp, _Compare>& __x,
			const VectorMap<_Key, _Tp, _Compare>& __y)
  { 
    return std::lexicographical_compare(__x.sortedVectorMap.begin(),
					__x.sortedVectorMap.end(),
					__y.sortedVectorMap.begin(),
					__y.sortedVectorMap.end(),
					__x.valueCompare);
  }

  /// Based on operator==
  template <typename _Key, typename _Tp, typename _Compare>
  inline bool operator!=(const VectorMap<_Key, _Tp, _Compare>& __x,
			 const VectorMap<_Key, _Tp, _Compare>& __y)
  { 
    return !(__x == __y); 
  }

  /// Based on operator<
  template <typename _Key, typename _Tp, typename _Compare>
  inline bool operator>(const VectorMap<_Key, _Tp, _Compare>& __x,
			const VectorMap<_Key, _Tp, _Compare>& __y)
  { 
    return __y < __x; 
  }

  /// Based on operator<
  template <typename _Key, typename _Tp, typename _Compare>
  inline bool operator<=(const VectorMap<_Key, _Tp, _Compare>& __x,
			 const VectorMap<_Key, _Tp, _Compare>& __y)
  { 
    return !(__y < __x); 
  }

  /// Based on operator<
  template <typename _Key, typename _Tp, typename _Compare>
  inline bool operator>=(const VectorMap<_Key, _Tp, _Compare>& __x,
			 const VectorMap<_Key, _Tp, _Compare>& __y)
  { 
    return !(__x < __y); 
  }

  /// See VectorMap::swap().
  template <typename _Key, typename _Tp, typename _Compare>
  inline void swap(VectorMap<_Key, _Tp, _Compare>& __x,
		   VectorMap<_Key, _Tp, _Compare>& __y)
  { 
    __x.swap(__y); 
  }

} // namespace util

#endif // Utilities_VectorMap_H
