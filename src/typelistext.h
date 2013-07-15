/***************************************************************************
                          typelist.h  -  description
                             -------------------
    copyright            : (C) 2004 by Volodymyr Myrnyy (Vladimir Mirnyi)
 ***************************************************************************
 Permission to use, copy, modify, distribute and sell this software for any
 purpose is hereby granted without fee, provided that the above copyright
 notice appear in all copies and that both that copyright notice and this
 permission notice appear in supporting documentation.
 ***************************************************************************/
/*! \file
    \brief Some extentions to the template metaalgorithms on typelists from
           Typelist.h of the Loki library by Andrei Alexandrescu (see his book
	   "Modern C++ design" and <a href="http://www.moderncppdesign.com" target="_blank">www.moderncppdesign.com</a>).

	   - Another empty type UnitType is added
	   - Append algorithm is extended to understand UnitType
	   - Replace algorithm can now replace a type with a typelist (not only another type)
	   - Algorithm Repeat builds a typelist of type T repeated N times
	   - very important class Print forces compiler to print a typelist in the compiler
	     messages for debugging purposes
*/

#ifndef __typelist_h
#define __typelist_h

#include "loki/Typelist.h"

namespace Loki {

/// \brief Another empty class, which can also end a typelist
class UnitType {};

namespace TL {

/// \brief A complement of a class template Loki::TL::Append
/// to handle typelists, which end with UnitType
template <> struct Append<UnitType, UnitType>
{
   typedef UnitType Result;
};

template <class T> struct Append<UnitType, T>
{
   typedef Typelist<T,UnitType> Result;
};

template <class Head, class Tail>
struct Append<UnitType, Typelist<Head, Tail> >
{
   typedef Typelist<Head, Tail> Result;
};


/// \brief A complement of a class template Loki::TL::Replace, which substitutes
/// Typelist<Head1,Tail1> instead of the type T in Typelist<T, Tail>.
template <class T, class Tail, class Head1, class Tail1>
struct Replace<Typelist<T, Tail>, T, Typelist<Head1,Tail1> >
{
   typedef Typelist<Head1,typename Append<Tail1,Tail>::Result> Result;
};


/// \brief Builds a typelist, which includes type T N times
template<class T, unsigned int N>
struct Repeat
{
   typedef Typelist<T,typename Repeat<T,N-1>::Result> Result;
};

template<class T>
struct Repeat<T,0>
{
   typedef NullType Result;
};


/// \class Print
/// \brief Prints a typelist as compiler messages at compile-time
///
/// This metaalgorithm is developed for the debugging purposes only.
/// If you type Print<SomeTypelist>, you can see the whole hierarchy of
/// SomeTypelist in the compiler messages.
/// Specialization of the end of a Typelist generates the error:
/// no type named 'Result' in struct ...
/// The compilation succeeds only if SomeTypelist is empty.
///////////////////////////////////////////////////////////
template<class TList> struct Print;

template<> struct Print<NullType> { };
template<> struct Print<UnitType> { };

template<class Head, class Tail>
struct Print<Typelist<Head,Tail> > {
   typedef typename Print<Tail>::Result Result;
};

/// \class Next
/// \brief Makes N steps into Typelist recursive structure and returns 
/// the next Typelist or Nulltype, if the end has been reached
///////////////////////////////////////////////////////////////
template<class TList, unsigned int N, unsigned int I=N>
struct Next;

template<class H, class T, unsigned int N, unsigned int I>
struct Next<Typelist<H,T>,N,I>
{
   typedef typename Next<T,N,I-1>::Result Result;
};

template<class H, class T, unsigned int N>
struct Next<Typelist<H,T>,N,0>
{
   typedef Typelist<H,T> Result;
};

template<unsigned int N, unsigned int I>
struct Next<NullType,N,I>
{
   typedef NullType Result;
};

template<unsigned int N>
struct Next<NullType,N,0>
{
   typedef NullType Result;
};

/// \class BitsetSelect
/// \brief Selects a subset from TList corresponding to the Bitset
///////////////////////////////////////////////////////////////
template<int Bitset, class TList>
struct BitsetSelect;

template<int Bitset, class H, class T>
struct BitsetSelect<Bitset,Typelist<H,T> > {
private:
   typedef typename BitsetSelect<Bitset/2,T>::Result next;
public:
   typedef typename Select<((Bitset%2)==0),
            next,Typelist<H,next> >::Result Result;
};

template<int Bitset>
struct BitsetSelect<Bitset,NullType> {
   typedef NullType Result;
};

template<class H, class T>
struct BitsetSelect<0,Typelist<H,T> > {
   typedef NullType Result;
};



template<class TList, unsigned int Index> 
struct EraseAt;

template<class Head, class Tail, unsigned int Index> 
struct EraseAt<Typelist<Head,Tail>,Index>
{
  typedef Typelist<Head, typename EraseAt<Tail,Index-1>::Result> Result;
};

template<class Head, class Tail> 
struct EraseAt<Typelist<Head,Tail>,0>
{
  typedef Tail Result;
};

template<unsigned int Index> 
struct EraseAt<NullType,Index>
{
  typedef NullType Result;
};


template<class TList, unsigned int N>
struct ShiftLeft : public Next<TList, N> {};

template<class TList, unsigned int N, class DefaultType, unsigned int I=N>
struct ShiftRight {
  typedef Typelist<DefaultType,typename ShiftRight<TList,N,DefaultType,I-1>::Result> Result;
};
  
template<class TList, unsigned int N, class DefaultType>
struct ShiftRight<TList,N,DefaultType,0> {
  typedef TList Result;
};


} //namespace TL
} //namespace Loki

#endif /* __typelist_h */
