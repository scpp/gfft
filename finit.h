/***************************************************************************
 *   Copyright (C) 2007 by Volodymyr Myrnyy                                *
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

#include "loki/Typelist.h"

////// template class FactoryInit
// object factory initialization class
// register all the classes from TList in
// object factory Fact

template<class TList>
struct FactoryInit;

template<class H, class T>
struct FactoryInit<Loki::Typelist<H,T> > {
   template<class Fact>
   static void apply(Fact& f) {
      f.Register(H::id,H::Create);
      FactoryInit<T>::apply(f);
   }
};

template<>
struct FactoryInit<Loki::NullType> {
   template<class Fact>
   static void apply(Fact&) { }
};

#endif /*__finit_h*/
