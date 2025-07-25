// --------------------------------------------------
//
// (C) Copyright Chuck Allison and Jeremy Siek 2001 - 2002.
// (C) Copyright Gennaro Prota                 2003 - 2004.
//
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)
//
// -----------------------------------------------------------

//  See http://www.boost.org/libs/dynamic_bitset for documentation.


#ifndef BOOST_DETAIL_DYNAMIC_BITSET_HPP
#define BOOST_DETAIL_DYNAMIC_BITSET_HPP

#include <cstddef> // for std::size_t
#if 0
#include "boost/config.hpp"
#include "boost/detail/workaround.hpp"
//#include "boost/static_assert.hpp" // gps
#endif


namespace boost {

  namespace detail {

    // Gives (read-)access to the object representation
    // of an object of type T (3.9p4). CANNOT be used
    // on a base sub-object
    //
    template <typename T>
    inline const unsigned char * object_representation (T* p)
    {
        return static_cast<const unsigned char *>(static_cast<const void *>(p));
    }


    // ------- count function implementation --------------

    namespace dynamic_bitset_count_impl {

    typedef unsigned char byte_type;

    enum mode { access_by_bytes, access_by_blocks };

    template <mode> struct mode_to_type {};

    // the table: wrapped in a class template, so
    // that it is only instantiated if/when needed
    //
    template <bool dummy_name = true>
    struct count_table { static const byte_type table[]; };

    template <>
    struct count_table<false> { /* no table */ };


     const unsigned int table_width = 8;
     template <bool b>
     const byte_type count_table<b>::table[] =
     {
       // Automatically generated by GPTableGen.exe v.1.0
       //
     0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
     1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
     1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
     2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
     1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
     2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
     2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
     3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8
     };


     // overload for access by bytes
     //

     template <typename Iterator>
     inline std::size_t do_count(Iterator first, std::size_t length,
                                 int /*dummy param*/,
                                 mode_to_type<access_by_bytes>* )
     {
         std::size_t num = 0;
         if (length)
         {
             const byte_type * p = object_representation(&*first);
             length *= sizeof(*first);

              do {
                 num += count_table<>::table[*p];
                 ++p;
                 --length;

             } while (length);
         }

         return num;
     }


     // overload for access by blocks
     //
     template <typename Iterator, typename ValueType>
     inline std::size_t do_count(Iterator first, std::size_t length, ValueType,
                                 mode_to_type<access_by_blocks>*)
     {
         std::size_t num = 0;
         while (length){

             ValueType value = *first;
             while (value) {
                 num += count_table<>::table[value & ((1u<<table_width) - 1)];
                 value >>= table_width;
             }

             ++first;
             --length;
         }

         return num;
     }


    } // dynamic_bitset_count_impl
    // -------------------------------------------------------


    // Some library implementations simply return a dummy
    // value such as
    //
    //   size_type(-1) / sizeof(T)
    //
    // from vector<>::max_size. This tries to get out more
    // meaningful info.
    //
    template <typename T>
    typename T::size_type vector_max_size_workaround(const T & v) {

      typedef typename T::allocator_type allocator_type;

      const typename allocator_type::size_type alloc_max =
                                                  v.get_allocator().max_size();
      const typename T::size_type container_max = v.max_size();

      return alloc_max < container_max?
                    alloc_max :
                    container_max;
    }

    // for static_asserts
    template <typename T>
    struct dynamic_bitset_allowed_block_type {
        enum { value = T(-1) > 0 }; // ensure T has no sign
    };

    template <>
    struct dynamic_bitset_allowed_block_type<bool> {
        enum { value = false };
    };


  } // namespace detail

} // namespace boost

#endif // include guard

