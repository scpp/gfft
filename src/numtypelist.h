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

#ifndef __numtypelist_h
#define __numtypelist_h

#include "loki/Typelist.h"

#include "sint.h"

using namespace Loki;

/// \class Max
/// \brief Calculates the maximum value of a numlist
/// \param NList a numtypelist
/// \return an integer that is the maximum value of NList \n
/// Max<NList>::value
////////////////////////////////////////////////////////////////////////////////

        template <class NList> struct Max;
        template <class Num>
        struct Max< Typelist<Num,NullType> >
        {
            typedef Num Result;
        };

        template <class Num, class Tail>
        struct Max< Typelist<Num, Tail> >
        {
        private:
            enum { temp = Max<Tail>::Value };
        public:
            typedef SInt<(temp > Num::Value) ? temp : Num::Value> Result;
        };

/// \class Min
/// \brief Calculates the minimum value of a numlist
/// \param NList a numlist
/// \return an integer that is the minimum value of NList \n
/// Min<NList>::value
////////////////////////////////////////////////////////////////////////////////

        template <class NList> struct Min;
        template <class Num>
        struct Min< Typelist<Num,NullType> >
        {
            typedef Num Result;
        };

        template <class Num, class Tail>
        struct Min< Typelist<Num, Tail> >
        {
        private:
            enum { temp = Min<Tail>::Value };
        public:
            typedef SInt<(temp < Num::Value) ? temp : Num::Value> Result;
        };

/// \class AddConst
/// \brief Adds a constant integer to all elements of a numlist
/// \param NList a numlist
/// \param Num a constant integer
/// \return a numlist which elements increased on Num comparing to NList \n
/// AddConst<NList,Num>::Result
////////////////////////////////////////////////////////////////////////////////

        template <class NList, class Num> struct AddConst;
        template <class Num>
        struct AddConst<NullType,Num>
        {
            typedef NullType Result;
        };

        template <class Head, class Tail, class Num>
        struct AddConst<Typelist<Head,Tail>,Num>
        {
            typedef Typelist<SInt<Head::Value + Num::Value>,
                    typename AddConst<Tail,Num>::Result> Result;
        };

/// \class AddAt
/// \brief Adds an integer to the element with given index of a numlist
/// \param NList a numlist
/// \param index a non-negative integer
/// \param Num an integer
/// \return a numlist in which the i-th element increased on Num comparing to NList
/// An out-of-bounds index results in a compile-time error \n
/// AddAt<NList,index,Num>::Result
////////////////////////////////////////////////////////////////////////////////

        template <class NList, unsigned int index, class Num> struct AddAt;

        template <class Head, class Tail, class Num>
        struct AddAt<Typelist<Head, Tail>, 0, Num>
        {
            typedef Typelist<SInt<Head::Value + Num::Value>, Tail> Result;
        };

        template <class Head, class Tail, unsigned int i, class Num>
        struct AddAt<Typelist<Head, Tail>, i, Num>
        {
            typedef Typelist<Head,
                    typename AddAt<Tail, i-1, Num>::Result>
                Result;
        };

/// \class Add
/// \brief Adds two numlists even then if their lengths are different
/// \param NList1 a numlist
/// \param NList2 a numlist
/// \return a numlist which elements are sums of correspondent in NList1 and NList2
/// the tail of the longest numlist is appended \n
/// Add<NList1,NList2>::Result
////////////////////////////////////////////////////////////////////////////////

        template <class NList1, class NList2> struct Add;

        template <>
        struct Add<NullType,NullType>
        {
            typedef NullType Result;
        };

        template <class Head, class Tail>
        struct Add<NullType,Typelist<Head,Tail> >
        {
            typedef Typelist<Head,Tail> Result;
        };

        template <class Head, class Tail>
        struct Add<Typelist<Head,Tail>,NullType>
        {
            typedef Typelist<Head,Tail> Result;
        };

        template <class Head1, class Tail1, class Head2, class Tail2>
        struct Add<Typelist<Head1,Tail1>,Typelist<Head2,Tail2> >
        {
            typedef Typelist<SInt<Head1::Value + Head2::Value>,
                    typename Add<Tail1,Tail2>::Result> Result;
        };


/// \class SubConst
/// \brief Subtracts a constant integer to all elements of a numlist
/// \param NList a numlist
/// \param Num a constant integer
/// \return a numlist which elements decreased on Num comparing to NList \n
/// SubConst<NList,Num>::Result
////////////////////////////////////////////////////////////////////////////////

        template <class NList, class Num> struct SubConst;

        template <class Num>
        struct SubConst<NullType,Num>
        {
            typedef NullType Result;
        };

        template <class Head, class Tail, class Num>
        struct SubConst<Typelist<Head,Tail>,Num>
        {
            typedef Typelist<SInt<Head::Value - Num::Value>,
                    typename SubConst<Tail,Num>::Result> Result;
        };

/// \class SubAt
/// \brief Subtracts an integer to the element with given index of a numlist
/// \param NList a numlist
/// \param index a non-negative integer
/// \param Num an integer
/// \return a numlist in which the i-th element is decreased on Num comparing to NList
/// An out-of-bounds index results in a compile-time error \n
/// SubAt<NList,index,Num>::Result
////////////////////////////////////////////////////////////////////////////////

        template <class NList, unsigned int index, class Num> struct SubAt;

        template <class Head, class Tail, class Num>
        struct SubAt<Typelist<Head, Tail>, 0, Num>
        {
            typedef Typelist<SInt<Head::Value - Num::Value>, Tail> Result;
        };

        template <class Head, class Tail, unsigned int i, class Num>
        struct SubAt<Typelist<Head, Tail>, i, Num>
        {
            typedef Typelist<Head,
                    typename SubAt<Tail, i-1, Num>::Result> Result;
        };

/// \class Sub
/// \brief Subtracts two numlists even then if their lengths are different
/// \param NList1 a numlist
/// \param NList2 a numlist
/// \return a numlist which elements are differences of correspondent in NList1 and NList2
/// the tail of the longest numlist is appended \n
/// Sub<NList1,NList2>::Result
////////////////////////////////////////////////////////////////////////////////

        template <class NList1, class NList2> struct Sub;

        template <>
        struct Sub<NullType,NullType>
        {
            typedef NullType Result;
        };

        template <class Head, class Tail>
        struct Sub<NullType,Typelist<Head,Tail> >
        {
            typedef Typelist<Head,Tail> Result;
        };

        template <class Head, class Tail>
        struct Sub<Typelist<Head,Tail>,NullType>
        {
            typedef Typelist<Head,Tail> Result;
        };

        template <class Head1, class Tail1, class Head2, class Tail2>
        struct Sub<Typelist<Head1,Tail1>,Typelist<Head2,Tail2> >
        {
            typedef Typelist<SInt<Head1::Value-Head2::Value>,
                    typename Sub<Tail1,Tail2>::Result> Result;
        };

/// \class MultConst
/// \brief Multiplies all elements of a numlist by a constant integer
/// \param NList a numlist
/// \param Num a constant integer
/// \return a numlist which elements multiplied by Num comparing to NList \n
/// MultConst<NList,Num>::Result
////////////////////////////////////////////////////////////////////////////////

        template <class NList, class Num> struct MultConst;

        template <class Num>
        struct MultConst<NullType,Num>
        {
            typedef NullType Result;
        };

        template <class Head, class Tail, class Num>
        struct MultConst<Typelist<Head,Tail>,Num>
        {
            typedef Typelist<SInt<Head::Value * Num::Value>,
                    typename MultConst<Tail,Num>::Result> Result;
        };

/// \class MultAt
/// \brief Multiplies an element with given index of a numlist by an integer
/// \param NList a numlist
/// \param index a non-negative integer
/// \param Num an integer
/// \return a numlist in which the i-th element multiplied by Num comparing to NList
/// An out-of-bounds index results in a compile-time error \n
/// MultAt<NList,index,Num>::Result
////////////////////////////////////////////////////////////////////////////////

        template <class NList, unsigned int index, class Num> struct MultAt;

        template <class Head, class Tail, class Num>
        struct MultAt<Typelist<Head, Tail>, 0, Num>
        {
            typedef Typelist<SInt<Head::Value * Num::Value>, Tail> Result;
        };

        template <class Head, class Tail, unsigned int i, class Num>
        struct MultAt<Typelist<Head, Tail>, i, Num>
        {
            typedef Typelist<Head,
                    typename MultAt<Tail, i-1, Num>::Result>
                Result;
        };

/// \class Mult
/// \brief Multiplies element by element of two numlists even with different lengths
/// \param NList1 a numlist
/// \param NList2 a numlist
/// \return a numlist which elements are products of correspondent in NList1 and NList2
/// The tail of the longest numlist is truncated \n
/// Mult<NList1,NList2>::Result
////////////////////////////////////////////////////////////////////////////////

        template <class NList1, class NList2> struct Mult;
        template <>
        struct Mult<NullType,NullType>
        {
            typedef NullType Result;
        };

        template <class Head, class Tail>
        struct Mult<NullType,Typelist<Head,Tail> >
        {
            typedef NullType Result;
        };

        template <class Head, class Tail>
        struct Mult<Typelist<Head,Tail>,NullType>
        {
            typedef Typelist<Head,Tail> Result;
        };

        template <class Head1, class Tail1, class Head2, class Tail2>
        struct Mult<Typelist<Head1,Tail1>,Typelist<Head2,Tail2> >
        {
            typedef Typelist<SInt<Head1::Value*Head2::Value>,
                    typename Add<Tail1,Tail2>::Result> Result;
        };

/// \class Sum
/// \brief Computes the sum of elements in a numlist
/// \param NList a numlist
/// \return an integer that is the sum of NList \n
/// Sum<NList>::value
////////////////////////////////////////////////////////////////////////////////

        template <class NList> struct Sum;
        template <> struct Sum<NullType>
        {
            typedef SInt<0> Result;
        };

        template <class Num, class Tail>
        struct Sum< Typelist<Num, Tail> >
        {
            typedef SInt<Num::Value + Sum<Tail>::Result::Value> Result;
        };

/// \class Compare
/// \brief Compares elements in the numlist NList1 with NList2
/// \param NList1 a numlist
/// \param NList2 a numlist
/// \return Positive value, if the numlist NList1 greater than NList2,
///         negative value otherwise. Returns zero, if NList1 and NList2 are equal.
///         (last element is the most significant) \n
/// Compare<NList1,NList2>::value
////////////////////////////////////////////////////////////////////////////////

        template <class NList1, class NList2> struct Compare;

        template <class H1, class T1, class H2, class T2>
        struct Compare<Typelist<H1,T1>,Typelist<H2,T2> >
        {
            enum { v = Compare<T1,T2>::Result::Value };
            typedef SInt<(v==0) ? (H1::Value-H2::Value) : v> Result;
        };

        template <class H, class T>
        struct Compare<Typelist<H,T>,NullType>
        {
            typedef SInt<1> Result;
        };

        template <class H, class T>
        struct Compare<NullType,Typelist<H,T> >
        {
            typedef SInt<-1> Result;
        };

        template <class H1, class H2>
        struct Compare<Typelist<H1,NullType>,Typelist<H2,NullType> >
        {
            typedef SInt<(H1::Value-H2::Value)> Result;
        };

/// \class Sort
/// \brief Sorts a numlist applying consequent the operation Min
/// \param NList a numlist
/// \return a numlist which is sorted NList \n
/// Sort<NList>::Result
////////////////////////////////////////////////////////////////////////////////

        template <class NList> struct Sort;
        template <> struct Sort<NullType>
        {
            typedef NullType Result;
        };

        template <class Head, class Tail>
        struct Sort<Typelist<Head,Tail> >
        {
        private:
            typedef typename Min<Typelist<Head,Tail> >::Result _Min;
            typedef typename TL::Replace<Tail,_Min,Head>::Result temp;

        public:
            typedef Typelist<_Min,typename Sort<temp>::Result> Result;
        };



#endif
