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
class Place,              // IN_PLACE, OUT_OF_PLACE
class FactoryPolicy = Empty,
uint IDN = N::ID>
class Transform;


template<class N,
class VType,
class Type,                    // DFT, IDFT, RDFT, IRDFT
class Dim,
class Parall,
class FactoryPolicy,
uint IDN>
class Transform<N,VType,Type,Dim,Parall,IN_PLACE,FactoryPolicy,IDN> : public FactoryPolicy 
{
   typedef typename VType::ValueType T;
   typedef typename Parall::template Swap<N::value,T>::Result Swap;
   typedef typename Type::template Direction<N::value,T> Dir;
   typedef Separate<N::value,T,Dir::Sign> Sep;
   typedef Caller<Loki::NullType> EmptySwap;
   typedef typename Factorization<N, SInt>::Result NFact;
   //typedef typename IN_PLACE::template Check<NFact>::Result Checked;

   //typedef typename GenerateRootList<N::value,Dir::Sign,2>::Result RootList;
   static const int Accuracy = 2;
   typedef typename GetFirstRoot<N::value,Dir::Sign,Accuracy>::Result W1;
   
   typedef typename IN_PLACE::template List<N::value,NFact,T,Swap,Dir,Parall::NParProc,W1>::Result TList;
   typedef typename Type::template Algorithm<TList,Sep>::Result Alg;
   
   Caller<Loki::Typelist<Parall,Alg> > run;
   
public:
   typedef VType ValueType;
   typedef Type TransformType;
   typedef Parall ParallType;
   typedef IN_PLACE PlaceType;

   enum { ID = IDN, Len = N::value };

   static FactoryPolicy* Create() {
      return new Transform<N,VType,Type,Dim,Parall,IN_PLACE,FactoryPolicy,ID>();
   }

   Transform() { }
    ~Transform() { }

//   in-place transform
   void fft(T* data) { run.apply(data); }
};


template<class N,
class VType,
class Type,                    // DFT, IDFT, RDFT, IRDFT
class Dim,
class Parall,
class FactoryPolicy,
uint IDN>
class Transform<N,VType,Type,Dim,Parall,OUT_OF_PLACE,FactoryPolicy,IDN> : public FactoryPolicy 
{
   typedef typename VType::ValueType T;
   //typedef typename Parall::template Swap<N::value,T>::Result Swap;
   typedef typename Type::template Direction<N::value,T> Dir;
   typedef Separate<N::value,T,Dir::Sign> Sep;
   typedef Caller<Loki::NullType> EmptySwap;
   typedef typename Factorization<N, SInt>::Result NFact;

   //typedef typename GenerateRootList<N::value,Dir::Sign,2>::Result RootList;
   static const int Accuracy = 2;
   typedef typename GetFirstRoot<N::value,Dir::Sign,Accuracy>::Result W1;
   
   typedef typename OUT_OF_PLACE::template List<N::value,NFact,T,EmptySwap,Dir,Parall::NParProc,W1>::Result TList;
   typedef typename Type::template Algorithm<TList,Sep>::Result Alg;
   
   Caller<Loki::Typelist<Parall,Alg> > run;
   
   T* buf;
   
public:
   typedef VType ValueType;
   typedef Type TransformType;
   typedef Parall ParallType;
   typedef OUT_OF_PLACE PlaceType;

   enum { ID = IDN, Len = N::value };

   static FactoryPolicy* Create() {
      return new Transform<N,VType,Type,Dim,Parall,OUT_OF_PLACE,FactoryPolicy>();
   }

   Transform() 
   {
       buf = new T[2*N::value];
   }
 
   ~Transform() 
   {
       delete [] buf;
   }

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
   typedef typename TList::Tail::Tail::Head TransformType;
   typedef typename TList::Tail::Tail::Tail::Tail::Tail::Head Place;
   typedef typename Place::template Interface<typename VType::ValueType>::Result Abstract;
   
   typedef Transform<typename TList::Head, VType, TransformType,
                typename TList::Tail::Tail::Tail::Head,
                typename TList::Tail::Tail::Tail::Tail::Head,
                Place,Abstract,ID> Result;
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
class Place      = PlaceGroup::Default>             // IN_PLACE, OUT_OF_PLACE
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

   typedef TYPELIST_6(Place,Parall,Dim,TransType,T,NList) RevList;

   typedef TranslateID<LenList> Translate;

public:
   typedef typename ListGenerator<RevList,RevLenList,DefineTransform>::Result Result;
   typedef typename Place::template Interface<typename T::ValueType>::Result ObjectType;

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
