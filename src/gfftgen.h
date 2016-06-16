/***************************************************************************
 *   Copyright (C) 2012-2015 by Vladimir Mirnyy                            *
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

#include "Singleton.h"

/// Main namespace
namespace GFFT {

  
/** \class {GFFT::Transform}
\brief Generic Fast Fourier transform class
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
id_t IDN = N::ID>
class Transform 
{
   typedef typename VType::ValueType T;
   
   typedef typename Parall::template Factor<N>::Result NFactor;
   //static const long_t NR = PrecomputeRoots;
   //typedef ExtractFactor<NR, NFactor> EF;
   //static const long_t NN = N::value/NFactor::Head::first::value;
   //typedef Loki::Typelist<Pair<typename NFactor::Head::first, long_<NFactor::Head::second::value-1> >, typename NFactor::Tail> NFactNext;
   
   //typedef Loki::SingletonHolder<RootsHolder<NR,typename EF::Result,VType,Type::Sign> > Twiddles;

   typedef typename Type::template Algorithm<N::value,NFactor,VType,Parall,Place>::Result Alg;
   
   typedef typename Place::template Interface<typename VType::ValueType>::Result ReturnType;
   typedef typename Place::template Function<Caller<Loki::Typelist<Parall,Alg> >, T> ExecType;
   
public:
   typedef VType ValueType;
   typedef Type TransformType;
   typedef Parall ParallType;
   typedef Place PlaceType;

   typedef ExecType Instance;

   enum { ID = IDN };
   static const long_t Len = N::value;

   static ReturnType* Create() {
     return new ExecType();
   }

   Transform() { }
    ~Transform() { }
};


/// Takes types from TList as parameters to define Transform class
/**
\tparam TList is Typelist containing template parameters to define Transform class
\tparam ID is a unique id passed to Transform class as well

This template class extracts all the parameters from TList directly without 
iterations and substitutes them into Transform class. Obviously,
all the types of substituted template parameters must be correct.
*/
template<class TList, id_t ID>
struct DefineTransform {
   typedef typename TList::Tail::Head VType;
   typedef typename TList::Tail::Tail::Head TransformType;
   typedef typename TList::Tail::Tail::Tail::Tail::Tail::Head Place;
//    typedef typename Place::template Interface<typename VType::ValueType>::Result Abstract;
   
   typedef Transform<typename TList::Head, VType, TransformType,
                typename TList::Tail::Tail::Tail::Head,
                typename TList::Tail::Tail::Tail::Tail::Head,
                Place,ID> Result;
};



template<class NList>
struct TranslateID;

template<ulong_t N, class T>
struct TranslateID<Loki::Typelist<ulong_<N>,T> > {
   static ulong_t apply(const ulong_t* n) {
      return TranslateID<T>::apply(n+1)*N + *n;
   }
};

template<ulong_t N>
struct TranslateID<Loki::Typelist<ulong_<N>,Loki::NullType> > {
   static ulong_t apply(const ulong_t* n) {
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
class Dim        = ulong_<1>,
class Parall     = ParallelizationGroup::Default,
class Place      = PlaceGroup::Default>             // IN_PLACE, OUT_OF_PLACE
class GenerateTransform {
   //typedef typename GenNumList<Begin,End>::Result NList;
   static const ulong_t L1 = Loki::TL::Length<NList>::value;
   static const ulong_t L2 = Loki::TL::Length<ValueTypeGroup::FullList>::value;
   static const ulong_t L3 = Loki::TL::Length<TransformTypeGroup::FullList>::value;
   static const ulong_t L4 = 1;
   static const ulong_t L5 = Loki::TL::Length<ParallelizationGroup::FullList>::value;
   static const ulong_t L6 = Loki::TL::Length<PlaceGroup::FullList>::value;
   typedef TYPELIST_6(ulong_<L1>,ulong_<L2>,ulong_<L3>,ulong_<L4>,ulong_<L5>,ulong_<L6>) LenList;

   typedef typename Loki::TL::Reverse<LenList>::Result RevLenList;

   typedef TYPELIST_6(Place,Parall,Dim,TransType,T,NList) RevList;

   typedef TranslateID<LenList> Translate;

public:
   typedef typename ListGenerator<RevList,RevLenList,DefineTransform>::Result Result;
   typedef typename Place::template Interface<typename T::ValueType>::Result ObjectType;
   typedef Place PlaceType;

   Loki::Factory<ObjectType,ulong_t,ObjectType*(*)(),TransformFactoryError> factory;

   GenerateTransform() {
      FactoryInit<Result>::apply(factory);
   }

   ObjectType* CreateTransformObject(ulong_t n, ulong_t vtype_id,
                                     ulong_t trans_id = TransformTypeGroup::Default::ID,
                                     ulong_t dim = 1,
                                     ulong_t parall_id = ParallelizationGroup::Default::ID,
                                     ulong_t place_id = PlaceGroup::Default::ID)
   {
      ulong_t narr[] = {n-1, vtype_id, trans_id, dim-1, parall_id, place_id};
      ulong_t obj_id = Translate::apply(narr);
      return factory.CreateObject(obj_id);
   }

};

  
}  //namespace

#endif
