/***************************************************************************
 *   Copyright (C) 2007-2009 by Volodymyr Myrnyy                           *
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

#ifndef __finit_h
#define __finit_h

/** \file
    \brief Initialization of object factory
*/

#include "loki/Typelist.h"

template<class TList>
struct FactoryInit;

/// Object factory initialization class
/** Register all classes from TypeList in
    object factory Fact
*/

template<class H, class T>
struct FactoryInit<Loki::Typelist<H,T> > {
   template<class Fact>
   static void apply(Fact& f) {
      f.Register(H::ID,H::Create);
      FactoryInit<T>::apply(f);
   }
};

template<>
struct FactoryInit<Loki::NullType> {
   template<class Fact>
   static void apply(Fact&) { }
};


#endif /*__finit_h*/
