/***************************************************************************
 *   Copyright (C) 2012 by Volodymyr Myrnyy                                *
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
#include "typelistgen.h"
#include "gfftparamgroups.h"


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
class FactoryPolicy = Empty,
uint IDN = Power2::ID>
class Transform2:public FactoryPolicy {
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
      return new Transform2<Power2,VType,Type,Dim,Parall,Decimation,FactoryPolicy>();
   }

   // in-place transform
   void fft(T* data) { run.apply(data); }

   // out-of-place transform
   void fft(const T* src, T* dst) { run.apply(src, dst); }
};


template<class N,
class VType,
class Type,                    // DFT, IDFT, RDFT, IRDFT
class Dim,
class Parall,
class Decimation,              // INTIME, INFREQ
class FactoryPolicy = Empty,
uint IDN = N::ID>
class Transform:public FactoryPolicy {
   //enum { P = Type::template Length<Power2>::Value };
   //enum { N = 1<<P };
   typedef typename VType::ValueType T;
   typedef typename Parall::template Swap<N::value,T>::Result Swap;
   typedef typename Type::template Direction<N::Value,T> Dir;
   typedef Separate<N::Value,T,Dir::Sign> Sep;
   typedef Caller<Loki::NullType> EmptySwap;
   typedef typename Decimation::template List<N::Value,T,Swap,Dir,Parall::NParProc>::Result TList;
   typedef typename Type::template Algorithm<TList,Sep>::Result Alg;

   Caller<Loki::Typelist<Parall,Alg> > run;
public:
   typedef VType ValueType;
   typedef Type TransformType;
   typedef Parall ParallType;
   typedef Decimation DecimationType;

   enum { ID = IDN, Len = N::Value };

   static FactoryPolicy* Create() {
      return new Transform<N,VType,Type,Dim,Parall,Decimation,FactoryPolicy>();
   }

   // in-place transform
   void fft(T* data) { run.apply(data); }

   // out-of-place transform
   void fft(const T* src, T* dst) { run.apply(src, dst); }
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
class GeneratePower2Transform {
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

   GeneratePower2Transform() {
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



template<class NList,
class T         /* = ValueTypeList*/,        // has to be set explicitely because of the AbstractFFT<T>
class TransType  = TransformTypeGroup::Default,     // DFT, IDFT, RDFT, IRDFT
class Dim        = SIntID<1>,
class Parall     = ParallelizationGroup::Default,
class Decimation = DecimationGroup::Default>        // INTIME, INFREQ
class GenerateTransform {
   //typedef typename GenNumList<Begin,End>::Result NList;
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

   ObjectType* CreateTransformObject(uint n, uint vtype_id, 
                                   uint trans_id = TransformTypeGroup::default_id, 
                                   uint dim = 1, 
                                   uint parall_id = ParallelizationGroup::default_id, 
                                   uint decim_id = DecimationGroup::default_id) 
   {
      uint narr[] = {n-1, vtype_id, trans_id, dim-1, parall_id, decim_id};
      uint obj_id = Translate::apply(narr);
//std::cout<<obj_id<<std::endl;
      return factory.CreateObject(obj_id);
   }

};

  
}  //namespace

#endif
