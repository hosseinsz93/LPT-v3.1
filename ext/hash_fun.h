// 'struct hash' from SGI -*- C++ -*-

// Copyright (C) 2001, 2002, 2003, 2004 Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2, or (at your option)
// any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along
// with this library; see the file COPYING.  If not, write to the Free
// Software Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307,
// USA.

// As a special exception, you may use this file as part of a free software
// library without restriction.  Specifically, if other files instantiate
// templates or use macros or inline functions from this file, or you compile
// this file and link it with other files to produce an executable, this
// file does not by itself cause the resulting executable to be covered by
// the GNU General Public License.  This exception does not however
// invalidate any other reasons why the executable file might be covered by
// the GNU General Public License.

/*
 * Copyright (c) 1996-1998
 * Silicon Graphics Computer Systems, Inc.
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.  Silicon Graphics makes no
 * representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 *
 *
 * Copyright (c) 1994
 * Hewlett-Packard Company
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.  Hewlett-Packard Company makes no
 * representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 *
 */

/** @file ext/hash_fun.h
 *  This file is a GNU extension to the Standard C++ Library (possibly
 *  containing extensions from the HP/SGI STL subset).
 */

#ifndef _MKHASH_HASH_FUN_H
#define _MKHASH_HASH_FUN_H 1

#include <cstddef>
#include <string>

//namespace __gnu_cxx
namespace MKHash
{
  template<class _Key>
    struct hash { };

  inline std::size_t
  __stl_hash_string(const char* __s)
  {
    unsigned long __h = 0;
    for ( ; *__s; ++__s)
      __h = 5 * __h + *__s;
    return std::size_t(__h);
  }

  template<>
    struct hash<char*>
    {
      std::size_t
      operator()(const char* __s) const
      { return __stl_hash_string(__s); }
    };

  template<>
    struct hash<const char*>
    {
      std::size_t
      operator()(const char* __s) const
      { return __stl_hash_string(__s); }
    };

  template<>
    struct hash<std::string>
    {
      std::size_t
      operator()(std::string const & __s) const
      { return __stl_hash_string(__s.c_str()); }
    };

  template<>
    struct hash<const std::string>
    {
      std::size_t
      operator()(std::string const & __s) const
      { return __stl_hash_string(__s.c_str()); }
    };

  template<>
    struct hash<char>
    { 
      std::size_t
      operator()(char __x) const
      { return __x; }
    };

  template<>
    struct hash<unsigned char>
    { 
      std::size_t
      operator()(unsigned char __x) const
      { return __x; }
    };

  template<>
    struct hash<signed char>
    {
      std::size_t
      operator()(unsigned char __x) const
      { return __x; }
    };

  template<>
    struct hash<short>
    {
      std::size_t
      operator()(short __x) const
      { return __x; }
    };

  template<>
    struct hash<unsigned short>
    {
      std::size_t
      operator()(unsigned short __x) const
      { return __x; }
    };

  template<>
    struct hash<int>
    { 
      std::size_t 
      operator()(int __x) const 
      { return __x; }
    };

  template<>
    struct hash<unsigned int>
    { 
      std::size_t
      operator()(unsigned int __x) const
      { return __x; }
    };

  template<>
    struct hash<long>
    {
      std::size_t
      operator()(long __x) const
      { return __x; }
    };

  template<>
    struct hash<unsigned long>
    {
      std::size_t
      operator()(unsigned long __x) const
      { return __x; }
    };

} // namespace MKHash (old __gnu_cxx)

#endif
