/***************************************************************************
 *   Copyright (C) 2008 by Volodymyr Myrnyy                                *
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


#ifndef __gfftconf_h
#define __gfftconf_h

#include "gfft.h"

#include "Numlist.h"

#include <complex>

namespace GFFT {

template<class TList> struct Print;

template<> struct Print<Loki::NullType> { };

template<class Head, class Tail>
struct Print<Loki::Typelist<Head,Tail> > {
   typedef typename Print<Tail>::Result Result;
};

struct DOUBLE {
   enum { ID = 0 };
   typedef double ValueType;
};

struct FLOAT {
   enum { ID = 1 };
   typedef float ValueType;
};

struct COMPLEX_DOUBLE {
   enum { ID = 2 };
   typedef std::complex<double> ValueType;
};

struct COMPLEX_FLOAT {
   enum { ID = 3 };
   typedef std::complex<float> ValueType;
};

template<unsigned int N>
struct SInt {
   enum { Value = N, ID = N-1 };
};

typedef TYPELIST_2(DOUBLE,FLOAT) ValueTypeList;

typedef TYPELIST_2(COMPLEX,REAL) TransformTypeList;

typedef TYPELIST_2(FORWARD,BACKWARD) DirectionList;

typedef TYPELIST_2(Serial,OpenMP<2>) ParallelizationList;

typedef TYPELIST_2(INTIME,INFREQ) DecimationList;

typedef TYPELIST_4(DFT,IDFT,RDFT,IRDFT) NewTransformTypeList;


template<unsigned Begin, unsigned End>
struct GenNumList {
   typedef Loki::Typelist<SInt<Begin>,
      typename GenNumList<Begin+1,End>::Result> Result;
};

template<unsigned End>
struct GenNumList<End,End> {
   typedef Loki::NullType Result;
};

// template<class TList, int ID>
// struct DefineGFFT {
//    typedef typename TList::Tail::Tail::Tail::Head VType;
//    typedef GFFT<TList::Tail::Tail::Tail::Tail::Head::Value,VType,
//                 typename TList::Tail::Tail::Head,
//                 typename TList::Tail::Head,
//                 typename TList::Head,
//                 AbstractFFT<VType>,ID> Result;
// };

/// Takes types from TList to define GFFT class
template<class TList, int ID>
struct DefineGFFT {
   typedef typename TList::Tail::Head VType;
   typedef GFFT<TList::Head::Value,VType,
                typename TList::Tail::Tail::Head,
                typename TList::Tail::Tail::Tail::Head,
                typename TList::Tail::Tail::Tail::Tail::Head,
                typename TList::Tail::Tail::Tail::Tail::Tail::Head,
                AbstractFFT<typename VType::ValueType>,ID> Result;
};

/// Takes types from TList to define Transform class
template<class TList, int ID>
struct DefineTransform {
   typedef typename TList::Tail::Head VType;
   typedef Transform<typename TList::Head,VType,
                typename TList::Tail::Tail::Head,
                typename TList::Tail::Tail::Tail::Head,
                typename TList::Tail::Tail::Tail::Tail::Head,
                AbstractFFT<typename VType::ValueType>,ID> Result;
};

/// Generates all different combinations of parameters given 
/// in the two-dimensional compile-time array TList taking one
/// parameter from every TList's entry. The entry may be either 
/// a type or a Typelist.
template<class TList, class TLenList, 
         template<class,int> class DefTrans,
         class WorkingList=Loki::NullType, int ID=0>
struct ListGenerator;

// H is a simple type
template<class H, class Tail, int N, class NTail, 
         template<class,int> class DefTrans, class WorkingList, int ID>
struct ListGenerator<Loki::Typelist<H,Tail>,NL::Numlist<N,NTail>,DefTrans,WorkingList,ID> {
   typedef Loki::Typelist<H,WorkingList> WList;
   typedef typename ListGenerator<Tail,NTail,DefTrans,WList,(ID*N)+H::ID>::Result Result;
};

// Typelist is in the head
template<class H, class T, class Tail, int N, class NTail, 
         template<class,int> class DefTrans, class WorkingList, int ID>
struct ListGenerator<Loki::Typelist<Loki::Typelist<H,T>,Tail>,NL::Numlist<N,NTail>,DefTrans,WorkingList,ID> {
   typedef Loki::Typelist<H,WorkingList> WList;
   typedef typename ListGenerator<Loki::Typelist<T,Tail>,NL::Numlist<N,NTail>,DefTrans,WorkingList,ID>::Result L1;
   typedef typename ListGenerator<Tail,NTail,DefTrans,WList,(ID*N)+H::ID>::Result L2;
   typedef typename Loki::TL::Append<L1,L2>::Result Result;
};

template<class H, class Tail, int N, class NTail, 
         template<class,int> class DefTrans, class WorkingList, int ID>
struct ListGenerator<Loki::Typelist<Loki::Typelist<H,Loki::NullType>,Tail>,
                     NL::Numlist<N,NTail>,DefTrans,WorkingList,ID> {
   typedef Loki::Typelist<H,WorkingList> WList;
   typedef typename ListGenerator<Tail,NTail,DefTrans,WList,(ID*N)+H::ID>::Result Result;
};

template<template<class,int> class DefTrans, class WorkingList, int ID>
struct ListGenerator<Loki::NullType,NL::NullType,DefTrans,WorkingList,ID> {
   typedef Loki::Typelist<typename DefTrans<WorkingList,ID>::Result,Loki::NullType> Result;
};



template<class NList>
struct TranslateID;

template<int N, class T>
struct TranslateID<NL::Numlist<N,T> > {
   static unsigned int apply(const unsigned int* n) {
      return TranslateID<T>::apply(n+1)*N + *n;
   }
};

template<int N>
struct TranslateID<NL::Numlist<N,NL::NullType> > {
   static unsigned int apply(const unsigned int* n) {
      return *n;
   }
};



template<unsigned Begin, unsigned End,
class T          = ValueTypeList,
class TransType  = TransformTypeList,     // COMPLEX, REAL
class Direction  = DirectionList,
class Parall     = ParallelizationList,
class Decimation = INFREQ>        // INTIME, INFREQ
class GenList {
   typedef typename GenNumList<Begin,End>::Result NList;
   enum { L1 = Loki::TL::Length<NList>::value };
   enum { L2 = Loki::TL::Length<ValueTypeList>::value };
   enum { L3 = Loki::TL::Length<TransformTypeList>::value };
   enum { L4 = Loki::TL::Length<DirectionList>::value };
   enum { L5 = Loki::TL::Length<ParallelizationList>::value };
   enum { L6 = Loki::TL::Length<DecimationList>::value };
   typedef NUMLIST_6(L1,L2,L3,L4,L5,L6) LenList;
//   typedef TYPELIST_5(NList,T,TransType,Decimation,Direction) List;

//   typedef typename Loki::TL::Reverse<List>::Result RevList;  // ? fails
   typedef typename NL::Reverse<LenList>::Result RevLenList;

//    typedef NUMLIST_5(L5,L4,L3,L2,L1) RevLenList;
   typedef TYPELIST_6(Decimation,Parall,Direction,TransType,T,NList) RevList;

   typedef TranslateID<LenList> Trans;
public:
   typedef typename ListGenerator<RevList,RevLenList,DefineGFFT>::Result Result;

   static unsigned int trans_id(const unsigned int* n) {
      return Trans::apply(n);
   }
};

template<unsigned Begin, unsigned End,
class T         /* = ValueTypeList*/,        // has to be set explicitely because of the AbstractFFT<T>
class TransType  = NewTransformTypeList,     // DFT, IDFT, RDFT, IRDFT
class Parall     = ParallelizationList,
class Decimation = INFREQ>        // INTIME, INFREQ
class Generate {
   typedef typename GenNumList<Begin,End>::Result NList;
   enum { L1 = Loki::TL::Length<NList>::value };
   enum { L2 = Loki::TL::Length<ValueTypeList>::value };
   enum { L3 = Loki::TL::Length<NewTransformTypeList>::value };
   enum { L5 = Loki::TL::Length<ParallelizationList>::value };
   enum { L6 = Loki::TL::Length<DecimationList>::value };
   typedef NUMLIST_5(L1,L2,L3,L5,L6) LenList;
//   typedef TYPELIST_5(NList,T,TransType,Decimation,Direction) List;

//   typedef typename Loki::TL::Reverse<List>::Result RevList;  // ? fails
   typedef typename NL::Reverse<LenList>::Result RevLenList;

//    typedef NUMLIST_5(L5,L4,L3,L2,L1) RevLenList;
   typedef TYPELIST_5(Decimation,Parall,TransType,T,NList) RevList;

   typedef TranslateID<LenList> Trans;

   typedef AbstractFFT<typename T::ValueType> AbstrFFT;
   Loki::Factory<AbstrFFT, unsigned int> gfft;
  //  FactoryInit<List::Result>::apply(gfft);

public:
   typedef typename ListGenerator<RevList,RevLenList,DefineTransform>::Result Result;

   Generate() {
      FactoryInit<Result>::apply(gfft);
   }

   static unsigned int trans_id(const unsigned int* n) {
      unsigned int nn[5];
      for (int i=0; i<5; ++i) nn[i] = n[i];
      nn[0]--;
      return Trans::apply(nn);
   }

};

}  //namespace

#endif
