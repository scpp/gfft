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

/** \file
    \brief Definition of holders for static integer types
*/

#ifndef __sint_h
#define __sint_h

/// Integer number metacontainer.
/**     Integer N is wrapped into container class to handle integers and other
        compile-time number-classes using the same operation classes specializing
        them for particular number-container.
 \param N an integer number
*/
template<int N>
struct SInt {
   enum { Value = N };
};

#define STATIC_INTEGER_CLASS(Type, Name) \
template<Type N>                         \
struct s_##Name {                        \
   typedef Type value_type;              \
   static const Type value = N;          \
};

STATIC_INTEGER_CLASS(int, int)
STATIC_INTEGER_CLASS(unsigned int, uint)
STATIC_INTEGER_CLASS(long, long)
STATIC_INTEGER_CLASS(unsigned long, ulong)

#undef STATIC_INTEGER_CLASS

#endif
