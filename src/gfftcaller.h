/***************************************************************************
 *   Copyright (C) 2009 by Volodymyr Myrnyy                                *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 ***************************************************************************/

#ifndef __caller_h
#define __caller_h

#include "loki/Typelist.h"

/** \file
    \brief Caller classes
*/


namespace GFFT {


/** \class {GFFT::Caller}
\brief Calls member function apply() in every type of the typelist

This metaprogram goes through all types in TList and call
function apply(T*). All the classes in TList must 
include function apply(T*) with one parameter
as a pointer of type T*.
*/

template<class TList>
struct Caller;

template<class Head, class Tail>
struct Caller<Loki::Typelist<Head,Tail> > {
   template<typename T>
   void apply(T* data) {
      obj_.apply(data);
      next_.apply(data);
   }
private:
   Head obj_;
   Caller<Tail> next_;
};

template<>
struct Caller<Loki::NullType> {
   template<typename T>
   void apply(T*) { }
};


/** \class {GFFT::StaticCaller}
\brief Calls static member function apply() in every type of the typelist

This metaprogram goes through all types in TList and call
static function apply(T*). All the classes in TList must 
include static function apply(T*) with one parameter
as a pointer of type T*.
*/
template<class TList>
struct StaticCaller;

template<class Head, class Tail>
struct StaticCaller<Loki::Typelist<Head,Tail> > {
   template<typename T>
   static void apply(T* data) {
      Head::apply(data);
      StaticCaller<Tail>::apply(data);
   }
};

template<>
struct StaticCaller<Loki::NullType> {
   template<typename T>
   static void apply(T*) { }
};

} // namespace GFFT

#endif /*__caller_h*/
