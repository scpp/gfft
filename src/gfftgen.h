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


#ifndef __gfftgen_h
#define __gfftgen_h

/** \file
    \brief Configuration and generation classes
*/

#include "sint.h"


/// Main namespace
namespace GFFT {

typedef unsigned int uint;

/** \class {GFFT::Transform}
\brief Generic Fast Fourier transform in-place class
\tparam Power2 defines transform length, which is 2^Power2
\tparam VType type of data element
\tparam Type type of transform: DFT, IDFT, RDFT, IRDFT
\tparam Dim dimension of transform, defined as SIntID<N>, N=1,2,...
\tparam Parall parallelization
\tparam Decimation in-time or in-frequency: INTIME, INFREQ
\tparam FactoryPolicy policy used to create an object factory. Don't define it explicitely, if unsure

Use this class only, if you need transform of a single fixed type and length.
Otherwise, rely to template class GenerateTransform
*/
template<class Power2,
class VType,
class Type,                    // DFT, IDFT, RDFT, IRDFT
class Dim,
class Parall,
class Decimation,              // INTIME, INFREQ
class FactoryPolicy=Empty,
uint IDN = Power2::ID>
class Transform:public FactoryPolicy {
   enum { P = Type::template Length<Power2>::Value };
   enum { N = 1<<P };
   typedef typename VType::ValueType T;
   typedef typename Parall::template Swap<P,T>::Result Swap;
   typedef typename Type::template Direction<N,T> Dir;
   typedef Separate<N,T,Dir::Sign> Sep;
   typedef typename Decimation::template List<N,T,Swap,Dir,Parall::NParProc>::Result TList;
   typedef typename Type::template Algorithm<TList,Sep>::Result Alg;

   Caller<Loki::Typelist<Parall,Alg> > run;
public:
   typedef VType ValueType;
   typedef Type TransformType;
   typedef Parall ParallType;
   typedef Decimation DecimationType;

   enum { ID = IDN, Len = N, PLen = P };

   static FactoryPolicy* Create() {
      return new Transform<Power2,VType,Type,Dim,Parall,Decimation,FactoryPolicy>();
   }

   void fft(T* data) {
      run.apply(data);
   }
};

/*!
\defgroup gr_groups GFFT parameter groups
\brief Classes that represent groups of parameters to define transform
*/

/// \brief Lists all acceptable and default value types
/// \ingroup gr_groups
struct ValueTypeGroup
{
  typedef TYPELIST_4(DOUBLE,FLOAT,COMPLEX_DOUBLE,COMPLEX_FLOAT) FullList;
  static const uint Length = 4;
  typedef DOUBLE Default;
  static const uint default_id = DOUBLE::ID;
};

/// \brief Lists all acceptable and default types of Fast Fourier transform
/// \ingroup gr_groups
struct TransformTypeGroup
{
  typedef TYPELIST_4(DFT,IDFT,RDFT,IRDFT) FullList;
  static const uint Length = 4;
  typedef TYPELIST_2(DFT,IDFT) Default;
  static const uint default_id = DFT::ID;
};

/// \brief Lists all acceptable and default parallelization methods
/// \ingroup gr_groups
struct ParallelizationGroup
{
  typedef TYPELIST_2(Serial,OpenMP<2>) FullList;
  static const uint Length = 2;
  typedef Serial Default;
  static const uint default_id = Serial::ID;
};

/// \brief Lists all acceptable and default decimation versions
/// \ingroup gr_groups
struct DecimationGroup
{
  typedef TYPELIST_2(INTIME,INFREQ) FullList;
  static const uint Length = 2;
  typedef INFREQ Default;
  static const uint default_id = INFREQ::ID;
};

/// Static unsigned integer class holder with additional definition of ID
template<unsigned int N>
struct SIntID : public s_uint<N> {
   static const uint ID = N-1;
};

/// Generates Typelist with types SIntID<N>, N = Begin,...,End
template<unsigned int Begin, unsigned int End>
struct GenNumList {
   typedef Loki::Typelist<SIntID<Begin>,
      typename GenNumList<Begin+1,End>::Result> Result;
};

template<unsigned End>
struct GenNumList<End,End> {
   typedef Loki::Typelist<SIntID<End>,Loki::NullType> Result;
};

/// Takes types from TList as parameters to define Transform class
/**
\tparam TList is Typelist containing template parameters to define Transform class
\tparam ID is a unique id passed to Transform class as well

This template class extracts all the parameters from TList directly without 
iterations and substitutes them into Transform class. Obviously,
all the types of substituted template parameters must be correct.
*/
template<class TList, uint ID>
struct DefineTransform {
   typedef typename TList::Tail::Head VType;
   typedef Transform<typename TList::Head, VType,
                typename TList::Tail::Tail::Head,
                typename TList::Tail::Tail::Tail::Head,
                typename TList::Tail::Tail::Tail::Tail::Head,
                typename TList::Tail::Tail::Tail::Tail::Tail::Head,
                AbstractFFT<typename VType::ValueType>,ID> Result;
};

/** \class {GFFT::ListGenerator}
    \brief Generates all different combinations of given parameters.
\tparam TList is one- or two-dimensional (an entry can be a Typelist too) TypeList.
\tparam TLenList is Typelist of \a s_uint<N>, where N is maximum possible 
        lengths of every Typelist in TList. This array of numbers is used 
        to comute unique ID for each generated unique set of parameters.
\tparam DefTrans is transform definition class. 
        DefineTransform class is substituted here, but also definitions of 
        other template classes with suited parameters are potentially possible.
\tparam WorkingList is the working Typelist, which accumulates a unique set of
        parameters. It is being passed as parameter to the class DefineTransform,
        when the set is complete. This parameter must not be set explicitely.
\tparam ID unique id for each set of parameters. This static constant is generated 
        at compile-time and must not be set explicitely.

The parameters are given in the two-dimensional compile-time array.
This metaprogram takes one parameter from every TList's entry
and generates Typelist of unique sets of parameters to define Transform.
The entry may be either a type or a Typelist.
*/
template<class TList, class TLenList, 
         template<class,uint> class DefTrans,
         class WorkingList=Loki::NullType, uint ID=0>
struct ListGenerator;

// H is a simple type
template<class H, class Tail, uint N, class NTail, 
         template<class,uint> class DefTrans, class WorkingList, uint ID>
struct ListGenerator<Loki::Typelist<H,Tail>,Loki::Typelist<s_uint<N>,NTail>,DefTrans,WorkingList,ID> {
   typedef Loki::Typelist<H,WorkingList> WList;
   typedef typename ListGenerator<Tail,NTail,DefTrans,WList,(ID*N)+H::ID>::Result Result;
};

// Typelist is in the head
template<class H, class T, class Tail, uint N, class NTail, 
         template<class,uint> class DefTrans, class WorkingList, uint ID>
struct ListGenerator<Loki::Typelist<Loki::Typelist<H,T>,Tail>,Loki::Typelist<s_uint<N>,NTail>,DefTrans,WorkingList,ID> {
   typedef Loki::Typelist<H,WorkingList> WList;
   typedef typename ListGenerator<Loki::Typelist<T,Tail>,Loki::Typelist<s_uint<N>,NTail>,DefTrans,WorkingList,ID>::Result L1;
   typedef typename ListGenerator<Tail,NTail,DefTrans,WList,(ID*N)+H::ID>::Result L2;
   typedef typename Loki::TL::Append<L1,L2>::Result Result;
};

template<class H, class Tail, uint N, class NTail, 
         template<class,uint> class DefTrans, class WorkingList, uint ID>
struct ListGenerator<Loki::Typelist<Loki::Typelist<H,Loki::NullType>,Tail>,
                     Loki::Typelist<s_uint<N>,NTail>,DefTrans,WorkingList,ID> {
   typedef Loki::Typelist<H,WorkingList> WList;
   typedef typename ListGenerator<Tail,NTail,DefTrans,WList,(ID*N)+H::ID>::Result Result;
};

template<template<class,uint> class DefTrans, class WorkingList, uint ID>
struct ListGenerator<Loki::NullType,Loki::NullType,DefTrans,WorkingList,ID> {
   typedef Loki::Typelist<typename DefTrans<WorkingList,ID>::Result,Loki::NullType> Result;
};



template<class NList>
struct TranslateID;

template<uint N, class T>
struct TranslateID<Loki::Typelist<s_uint<N>,T> > {
   static unsigned int apply(const unsigned int* n) {
      return TranslateID<T>::apply(n+1)*N + *n;
   }
};

template<uint N>
struct TranslateID<Loki::Typelist<s_uint<N>,Loki::NullType> > {
   static unsigned int apply(const unsigned int* n) {
      return *n;
   }
};

template <typename IdentifierType, class AbstractProduct>
struct TransformFactoryError
{
   struct Exception : public std::exception 
   {
      const char* what() const throw() { 
        return "Requested transform is not compiled! Check your instantiation of GenerateTransform class!"; 
      }
   };
   static AbstractProduct* OnUnknownType(IdentifierType) {
      throw Exception();
   }
};


/// MAIN CLASS TO USE! Generates a set of transform classes
/**
\tparam Begin defines minimum transform length as a power of two (2^Begin)
\tparam End defines maximum transform length as a power of two (2^End)
\tparam T type of data element
\tparam TransType type of transform: DFT, IDFT, RDFT, IRDFT
\tparam Dim dimension of transform, defined as SIntID<N>, N=1,2,...
        By now, only one-dimensional transforms are implemented.
\tparam Parall parallelization method
\tparam Decimation in-time or in-frequency: INTIME, INFREQ
\tparam FactoryPolicy policy used to create an object factory. Don't define it explicitely, if unsure

This generator class makes possible to generate a set of necessary transforms.
The first three template parameters: minimum and maximum power of two and value type
must be defined. Further parameters have default values and may be omitted.
Default values for template parameters are taken from the corresponding group-classes
(see \ref gr_groups ), from returned type Default. To specify another set 
of transforms than the default one, you have to specify every template parameter either as
a Typelist or a single parameter class (see \ref gr_params ).

In the following example we define forward complex-valued set of transforms
of the length from 2^10 to 2^15 using double precision:
\code
typedef GenerateTransform<10, 15, GFFT::DOUBLE, GFFT::DFT> TransformSet;
\endcode

The next example defines set of forward and backward real-valued transforms of single recision,
where code for single and parallelized for two threads is generated:
\code
typedef TYPELIST_2(GFFT::RDFT, GFFT::IRDFT) RealTransforms;
typedef TYPELIST_2(GFFT::Serial, GFFT::OpenMP<2>) ParallSet;
typedef GenerateTransform<10, 15, GFFT::FLOAT, RealTransforms, SIntID<1>, ParallSet> TransformSet;
\endcode

You may use return type FullList of the group-classes to specify all types of a parameter. 
The next example returns forward and backward, complex- and real-valued transforms using double precision:
\code
typedef GenerateTransform<10, 15, GFFT::DOUBLE, GFFT::TransformTypeGroup::FullList> TransformSet;
\endcode

\sa Transform, ListGenerator
*/
template<unsigned Begin, unsigned End,
class T         /* = ValueTypeList*/,        // has to be set explicitely because of the AbstractFFT<T>
class TransType  = TransformTypeGroup::Default,     // DFT, IDFT, RDFT, IRDFT
class Dim        = SIntID<1>,
class Parall     = ParallelizationGroup::Default,
class Decimation = DecimationGroup::Default>        // INTIME, INFREQ
class GenerateTransform {
   typedef typename GenNumList<Begin,End>::Result NList;
   enum { L1 = Loki::TL::Length<NList>::value };
   enum { L2 = Loki::TL::Length<ValueTypeGroup::FullList>::value };
   enum { L3 = Loki::TL::Length<TransformTypeGroup::FullList>::value };
   enum { L4 = 1 };
   enum { L5 = Loki::TL::Length<ParallelizationGroup::FullList>::value };
   enum { L6 = Loki::TL::Length<DecimationGroup::FullList>::value };
   typedef TYPELIST_6(s_uint<L1>,s_uint<L2>,s_uint<L3>,s_uint<L4>,s_uint<L5>,s_uint<L6>) LenList;

   typedef typename Loki::TL::Reverse<LenList>::Result RevLenList;

   typedef TYPELIST_6(Decimation,Parall,Dim,TransType,T,NList) RevList;

   typedef TranslateID<LenList> Translate;

public:
   typedef typename ListGenerator<RevList,RevLenList,DefineTransform>::Result Result;
   typedef AbstractFFT<typename T::ValueType> ObjectType;

   Loki::Factory<ObjectType,uint,ObjectType*(*)(),TransformFactoryError> factory;

   GenerateTransform() {
      FactoryInit<Result>::apply(factory);
   }

   ObjectType* CreateTransformObject(uint p, uint vtype_id, 
                                   uint trans_id = TransformTypeGroup::default_id, 
                                   uint dim = 1, 
                                   uint parall_id = ParallelizationGroup::default_id, 
                                   uint decim_id = DecimationGroup::default_id) 
   {
      uint n[] = {p-1, vtype_id, trans_id, dim-1, parall_id, decim_id};
      uint obj_id = Translate::apply(n);
//std::cout<<obj_id<<std::endl;
      return factory.CreateObject(obj_id);
   }

//    static unsigned int trans_id(const unsigned int* n) {
//       unsigned int nn[6];
//       for (int i=0; i<6; ++i) nn[i] = n[i];
//       nn[0]--;
//       return Translate::apply(nn);
//    }

};

}  //namespace

#endif
