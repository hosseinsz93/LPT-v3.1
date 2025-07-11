// (C) Copyright David Abrahams 2002.
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// Boost versions of
//
//    std::iterator_traits<>::iterator_category
//    std::iterator_traits<>::difference_type
//    std::distance()
//
// ...for all compilers and iterators
//
// Additionally, if X is a pointer
//    std::iterator_traits<X>::pointer

// Otherwise, if partial specialization is supported or X is not a pointer
//    std::iterator_traits<X>::value_type
//    std::iterator_traits<X>::pointer
//    std::iterator_traits<X>::reference
//
// See http://www.boost.org for most recent version including documentation.

// Revision History
// 04 Mar 2001 - More attempted fixes for Intel C++ (David Abrahams)
// 03 Mar 2001 - Put all implementation into namespace
//               boost::detail::iterator_traits_. Some progress made on fixes
//               for Intel compiler. (David Abrahams)
// 02 Mar 2001 - Changed BOOST_MSVC to BOOST_MSVC_STD_ITERATOR in a few
//               places. (Jeremy Siek)
// 19 Feb 2001 - Improved workarounds for stock MSVC6; use yes_type and
//               no_type from type_traits.hpp; stopped trying to remove_cv
//               before detecting is_pointer, in honor of the new type_traits
//               semantics. (David Abrahams)
// 13 Feb 2001 - Make it work with nearly all standard-conforming iterators
//               under raw VC6. The one category remaining which will fail is
//               that of iterators derived from std::iterator but not
//               boost::iterator and which redefine difference_type.
// 11 Feb 2001 - Clean away code which can never be used (David Abrahams)
// 09 Feb 2001 - Always have a definition for each traits member, even if it
//               can't be properly deduced. These will be incomplete types in
//               some cases (undefined<void>), but it helps suppress MSVC errors
//               elsewhere (David Abrahams)
// 07 Feb 2001 - Support for more of the traits members where possible, making
//               this useful as a replacement for std::iterator_traits<T> when
//               used as a default template parameter.
// 06 Feb 2001 - Removed useless #includes of standard library headers
//               (David Abrahams)

#ifndef ITERATOR_DWA122600_HPP_
#    define ITERATOR_DWA122600_HPP_

#if 0
#    include <boost/config.hpp>
#    include <iterator>
#endif


// This is the case where everything conforms except BOOST_NO_STD_ITERATOR_TRAITS

namespace boost { namespace detail {

// Rogue Wave Standard Library fools itself into thinking partial
// specialization is missing on some platforms (e.g. Sun), so fails to
// supply iterator_traits!
template <class Iterator>
struct iterator_traits
{
    typedef typename Iterator::value_type value_type;
    typedef typename Iterator::reference reference;
    typedef typename Iterator::pointer pointer;
    typedef typename Iterator::difference_type difference_type;
    typedef typename Iterator::iterator_category iterator_category;
};

template <class T>
struct iterator_traits<T*>
{
    typedef T value_type;
    typedef T& reference;
    typedef T* pointer;
    typedef std::ptrdiff_t difference_type;
    typedef std::random_access_iterator_tag iterator_category;
};

template <class T>
struct iterator_traits<T const*>
{
    typedef T value_type;
    typedef T const& reference;
    typedef T const* pointer;
    typedef std::ptrdiff_t difference_type;
    typedef std::random_access_iterator_tag iterator_category;
};

}} // namespace boost::detail


namespace boost { namespace detail {

namespace iterator_traits_
{
  template <class Iterator, class Difference>
  struct distance_select
  {
      static Difference execute(Iterator i1, const Iterator i2, ...)
      {
          Difference result = 0;
          while (i1 != i2)
          {
              ++i1;
              ++result;
          }
          return result;
      }

      static Difference execute(Iterator i1, const Iterator i2, std::random_access_iterator_tag*)
      {
          return i2 - i1;
      }
  };
} // namespace boost::detail::iterator_traits_

template <class Iterator>
inline typename iterator_traits<Iterator>::difference_type
distance(Iterator first, Iterator last)
{
    typedef typename iterator_traits<Iterator>::difference_type diff_t;
    typedef typename ::boost::detail::iterator_traits<Iterator>::iterator_category iterator_category;
    
    return iterator_traits_::distance_select<Iterator,diff_t>::execute(
        first, last, (iterator_category*)0);
}

}}


#endif // ITERATOR_DWA122600_HPP_
