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

/// Static unsigned integer class holder with additional definition of ID
template<int_t N>
struct SIntID : public SInt<N> {
   static const uint ID = N-1;
};

template<int_t M, int_t P>
struct PowerHolder {
   static const int_t N = IPow<M,P>::value;
   static const int_t ID = N-1;
   static const int_t value = N;
};

template<int_t P>
struct PowerHolder<2,P> {
   static const int_t N = 1<<P;
   static const int_t ID = N-1;
   static const int_t value = N;
};

template<int_t P>
struct Power2holder : public PowerHolder<2,P> {};

template<int_t P>
struct Power3holder : public PowerHolder<3,P> {};


/// Generates Typelist with types Holder<N>, N = Begin,...,End
template<int_t Begin, int_t End, 
template<int_t> class Holder = SIntID>
struct GenNumList {
   typedef Loki::Typelist<Holder<Begin>,
      typename GenNumList<Begin+1,End,Holder>::Result> Result;
};

template<int_t End, 
template<int_t> class Holder>
struct GenNumList<End,End,Holder> {
   typedef Loki::Typelist<Holder<End>,Loki::NullType> Result;
};


template<class W1, int_t N, int Accuracy, class W, int_t Count, int_t I = 2>
struct __RootListLoop {
  typedef typename Simplify<SFraction<SInt<2*I>,SInt<N> > >::Result SF;
//typedef typename NL::Print<SF>::Result TTT;
  typedef typename GetNextRoot<SF::Numer::value,SF::Denom::value,W1,W,Accuracy>::Result WW;
  typedef Compute<typename WW::Re,Accuracy> CRe;
  typedef Compute<typename WW::Im,Accuracy> CIm;
  typedef typename __RootListLoop<W1,N,Accuracy,WW,Count,I+1>::Result Next;
  typedef Loki::Typelist<Pair<CRe,CIm>,Next> Result;
};

template<class W1, int_t N, int Accuracy, class W, int_t Count>
struct __RootListLoop<W1,N,Accuracy,W,Count,Count> {
  typedef Loki::NullType Result;
};

template<class RList>
struct GenerateSymmetricPart;

template<class H, class T>
struct GenerateSymmetricPart<Loki::Typelist<H,T> > {
  typedef typename H::first T1;
  typedef typename Negate<typename H::second::BigInt>::Result T2;
  typedef typename GenerateSymmetricPart<T>::Result Next;
  typedef Loki::Typelist<Pair<T1,T2>,Next> Result;
};

template<>
struct GenerateSymmetricPart<Loki::NullType> {
  typedef Loki::NullType Result;
};


template<int_t N, int S, int Accuracy>
class GetFirstRoot {
  //typedef typename SinPiDecimal<1,N,Accuracy>::Result Sin1;
  typedef typename SinPiDecimal<2,N,Accuracy>::Result Sin2;
  
  typedef typename Loki::Select<(S<0),Sin2,
          typename Negate<Sin2>::Result>::Result WI;
  //typedef typename FractionToDecimal<WI,Accuracy>::Result WIDec;
  
//   typedef typename Sub<SInt<1>,typename Mult<SInt<2>,
//           typename Mult<Sin1,Sin1>::Result>::Result>::Result WR;
  typedef typename CosPiDecimal<2,N,Accuracy>::Result WR;
  //typedef typename FractionToDecimal<WR,Accuracy>::Result WRDec;
  
public:
  typedef MComplex<WR,WI> Result;
};



template<int_t N, int S, int Accuracy>
class GenerateRootList 
{
  typedef typename GetFirstRoot<N,-S,Accuracy>::Result W1;
  
public:
  typedef Compute<typename W1::Re,Accuracy> CRe;
  typedef Compute<typename W1::Im,Accuracy> CIm;
  typedef Loki::Typelist<Pair<CRe,CIm>,typename __RootListLoop<W1,N,Accuracy,W1,(N%2==0) ? N/2 : N/2+1>::Result> FirstHalf;
  typedef typename Loki::Select<(N%2==0),
          typename Loki::TL::Append<FirstHalf,typename GenerateRootList<2,S,Accuracy>::Result>::Result,FirstHalf>::Result FirstHalfMod;
//  typedef typename Loki::TL::Reverse<typename GenerateSymmetricPart<FirstHalf>::Result>::Result SecondHalf;

//public:
//  typedef typename Loki::TL::Append<FirstHalfMod,SecondHalf>::Result Result;
  typedef FirstHalf Result;
};

template<int S, int Accuracy>
class GenerateRootList<2,S,Accuracy> {
  typedef Compute<SInt<-1>,Accuracy> CRe;
  typedef Compute<SInt<0>, Accuracy> CIm;
public:
  typedef Loki::Typelist<Pair<CRe,CIm>,Loki::NullType> Result;
};

template<int S, int Accuracy>
class GenerateRootList<1,S,Accuracy> {
public:
  typedef Loki::NullType Result;
};

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
Otherwise, rely on template class GenerateTransform
*/
template<class N,
class VType,
class Type,                    // DFT, IDFT, RDFT, IRDFT
class Dim,
class Parall,
class Decimation,              // INTIME, INFREQ
class FactoryPolicy = Empty,
uint IDN = N::ID>
class Transform : public FactoryPolicy 
{
   
   typedef typename VType::ValueType T;
   typedef typename Parall::template Swap<N::value,T>::Result Swap;
   typedef typename Type::template Direction<N::value,T> Dir;
   typedef Separate<N::value,T,Dir::Sign> Sep;
   typedef Caller<Loki::NullType> EmptySwap;
   typedef typename Factorization<N, SInt>::Result NFact;

   //typedef typename GenerateRootList<N::value,Dir::Sign,2>::Result RootList;
   static const int Accuracy = 2;
   typedef typename GetFirstRoot<N::value,Dir::Sign,Accuracy>::Result W1;
   
   typedef typename Decimation::template List<N::value,NFact,T,Swap,Dir,Parall::NParProc,W1>::Result TList;
   typedef typename Type::template Algorithm<TList,Sep>::Result Alg;
   
   Caller<Loki::Typelist<Parall,Alg> > run;
   
   T* buf;
   
public:
   typedef VType ValueType;
   typedef Type TransformType;
   typedef Parall ParallType;
   typedef Decimation DecimationType;

   enum { ID = IDN, Len = N::value };

   static FactoryPolicy* Create() {
      return new Transform<N,VType,Type,Dim,Parall,Decimation,FactoryPolicy>();
   }

   Transform() 
   {
     if (Decimation::ID == 1)
       buf = new T[2*N::value];
   }
 
   ~Transform() 
   {
     if (Decimation::ID == 1)
       delete [] buf;
   }

   // in-place transform
   void fft(T* data) { run.apply(data); }

   // out-of-place transform
   void fft(const T* src, T* dst) { run.apply(src, dst, buf); }
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
   static unsigned int apply(const int_t* n) {
      return TranslateID<T>::apply(n+1)*N + *n;
   }
};

template<uint N>
struct TranslateID<Loki::Typelist<s_uint<N>,Loki::NullType> > {
   static unsigned int apply(const int_t* n) {
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
\tparam NList Typelist containing transform lengths
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
template<class NList,
class T         /* = ValueTypeList*/,        // has to be set explicitely because of AbstractFFT<T>
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

   Loki::Factory<ObjectType,int_t,ObjectType*(*)(),TransformFactoryError> factory;

   GenerateTransform() {
      FactoryInit<Result>::apply(factory);
   }

   ObjectType* CreateTransformObject(int_t n, int_t vtype_id, 
                                   int_t trans_id = TransformTypeGroup::default_id, 
                                   int_t dim = 1, 
                                   int_t parall_id = ParallelizationGroup::default_id, 
                                   int_t decim_id = DecimationGroup::default_id) 
   {
      int_t narr[] = {n-1, vtype_id, trans_id, dim-1, parall_id, decim_id};
      int_t obj_id = Translate::apply(narr);
//std::cout<<obj_id<<std::endl;
      return factory.CreateObject(obj_id);
   }

};

  
}  //namespace

#endif
