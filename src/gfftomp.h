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

#ifndef __gfftomp_h
#define __gfftomp_h

/** \file
    \brief Parallelization template classes using %OpenMP standard
*/

#include "gfftalg.h"

#include <omp.h>

namespace GFFT {

/** \var static const uint SwitchToOMP
This static constant defines FFT length for that 
OpenMP parallelization is switched on. The overhead is
too large for transforms with smaller sizes.
*/
static const unsigned int SwitchToOMP = (1<<12);

/** \class {GFFT::InTimeOMP}
\brief %OpenMP parallelized Danielson-Lanczos section of the decimation-in-time FFT version.
\tparam NThreads is number of threads
\tparam N current transform length
\tparam T value type of the data array
\tparam S sign of the transform: 1 - forward, -1 - backward
\tparam C condition to ensure that (N>NThreads) and (N>=SwitchToOMP), otherwise 
        parallelization is meaningless and sequential implementation InTime
        is inherited.

Comparing to sequential implementation in template class InTime, this class
runs apply() function of both instances of the half length (N/2) in the separated
threads and so on until NThreads has become equal 1. Then the sequential version
in template class InTime is inherited.
\sa InFreqOMP, InTime, InFreq
*/
template<unsigned int NThreads, unsigned long N, typename T, int S, bool C=((N>NThreads) && (N>=SwitchToOMP))>
class InTimeOMP;

template<unsigned int NThreads, unsigned long N, typename T, int S>
class InTimeOMP<NThreads,N,T,S,true> {
   typedef typename TempTypeTrait<T>::Result LocalVType;
   InTimeOMP<NThreads/2,N/2,T,S> next;
   static const long IN = N;
public:
   void apply(T* data) {

      LocalVType wpr,wpi,t,wtemp,tempr,tempi,wr,wi;

      #pragma omp parallel shared(data,wpr,wpi,t) private(wtemp,tempr,tempi,wr,wi)
      {
        #pragma omp sections
        {
          #pragma omp section
          next.apply(data);

          #pragma omp section
          next.apply(data+N);
        }

      t = Sin<N,1,LocalVType>::value();
      wpr = -2.0*t*t;
      wpi = -S*Sin<N,2,LocalVType>::value();
      wr = 1.0;
      wi = 0.0;
      long i,chunk = N/2;

      #pragma omp for schedule(static,chunk)
      for (i=0; i<IN; i+=2) {
        tempr = data[i+N]*wr - data[i+N+1]*wi;
        tempi = data[i+N]*wi + data[i+N+1]*wr;
        data[i+N] = data[i]-tempr;
        data[i+N+1] = data[i+1]-tempi;
        data[i] += tempr;
        data[i+1] += tempi;

        wtemp = wr;
        wr += wr*wpr - wi*wpi;
        wi += wi*wpr + wtemp*wpi;
      }
      } // parallel*/
   }
};

template<unsigned long N, typename T, int S>
class InTimeOMP<1,N,T,S,true> : public InTime<N,T,S> { };

template<unsigned int NThreads, unsigned long N, typename T, int S>
class InTimeOMP<NThreads,N,T,S,false> : public InTime<N,T,S> { };


/** \class {GFFT::InFreqOMP}
\brief %OpenMP parallelized Danielson-Lanczos section of the decimation-in-time FFT version.
\tparam NThreads is number of threads
\tparam N current transform length
\tparam T value type of the data array
\tparam S sign of the transform: 1 - forward, -1 - backward
\tparam C condition to ensure that (N>NThreads) and (N>=SwitchToOMP), otherwise 
        parallelization is meaningless and sequential implementation InFreq
        is inherited.

Comparing to sequential implementation in template class InFreq, this class
runs apply() function of both instances of the half length (N/2) in the separated
threads and so on until NThreads has become equal 1. Then the sequential version
in template class InTime is inherited.
\sa InFreqOMP, InTime, InFreq
*/
template<unsigned int NThreads, unsigned long N, typename T, int S, 
bool C=((N>NThreads) && (N>=SwitchToOMP))>
class InFreqOMP;

template<unsigned int NThreads, unsigned long N, typename T, int S>
class InFreqOMP<NThreads,N,T,S,true> {
   typedef typename TempTypeTrait<T>::Result LocalVType;
   InFreqOMP<NThreads/2,N/2,T,S> next;
   static const long IN = N;
public:
   void apply(T* data) {

      LocalVType wtemp,tempr,tempi,wr,wi,wpr,wpi;

      #pragma omp parallel shared(data,wpr,wpi) private(wtemp,tempr,tempi,wr,wi)
      {

      wtemp = Sin<N,1,LocalVType>::value();
      wpr = -2.0*wtemp*wtemp;
      wpi = -S*Sin<N,2,LocalVType>::value();
      wr = 1.0;
      wi = 0.0;
      long i,chunk = N/2;

      #pragma omp for schedule(static,chunk)
      for (i=0; i<IN; i+=2) {
        tempr = data[i] - data[i+N];
        tempi = data[i+1] - data[i+N+1];
        data[i] += data[i+N];
        data[i+1] += data[i+N+1];
        data[i+N] = tempr*wr - tempi*wi;
        data[i+N+1] = tempi*wr + tempr*wi;

        wtemp = wr;
        wr += wr*wpr - wi*wpi;
        wi += wi*wpr + wtemp*wpi;
      }

      #pragma omp sections
      {
        #pragma omp section
        next.apply(data);
        #pragma omp section
        next.apply(data+N);
      }
      } // parallel
   }
};

template<unsigned long N, typename T, int S>
class InFreqOMP<1,N,T,S,true> : public InFreq<N,T,S> { };

template<unsigned int NThreads, unsigned long N, typename T, int S>
class InFreqOMP<NThreads,N,T,S,false> : public InFreq<N,T,S> { };




// doesn't work
template<unsigned int NThreads, unsigned long N, typename T>
class GFFTswapOMP {
public:
   void apply(T* data) {
     unsigned long i,m,j=0;
     unsigned long n = N;
     #pragma omp parallel firstprivate(data) private(i,j,m)
     {
       #pragma omp for ordered
         for (i=0; i<2*n-1; i+=2) {
           if (j>i) {
             std::swap(data[j], data[i]);
             std::swap(data[j+1], data[i+1]);
           }
           m = n;
           while (m>=2 && j>=m) {
             j -= m;
             m >>= 1;
           }
           j += m;
         }
     }
   }
};

/** \class {GFFT::GFFTswap2OMP}
\brief Binary reordering paralelized by %OpenMP
\tparam NThreads is number of threads
\tparam P current transform length as power of two
\tparam T value type of the data array
\tparam I is internal counter
\tparam C condition to ensure that (N>NThreads) and (N>=SwitchToOMP), otherwise 
        parallelization is meaningless and the sequential implementation GFFTswap2
        is inherited.
*/
template<unsigned int NThreads, unsigned int P, typename T,
unsigned int I=0, bool C=(((1<<P)>NThreads) && ((1<<P)>=SwitchToOMP))>
class GFFTswap2OMP;

template<unsigned int NThreads, unsigned int P, typename T, unsigned int I>
class GFFTswap2OMP<NThreads,P,T,I,true> {
   static const unsigned long BN = 1<<(I+1);
   static const unsigned long BR = 1<<(P-I);
   GFFTswap2OMP<NThreads/2,P,T,I+1> next;
public:
   void apply(T* data, const unsigned long n=0, const unsigned long r=0) {
     #pragma omp parallel shared(data)
     {
       #pragma omp sections
       {
         #pragma omp section
         next.apply(data,n,r);

         #pragma omp section
         next.apply(data,n|BN,r|BR);
       }
     }
   }
};

template<unsigned int NThreads, unsigned int P, typename T>
class GFFTswap2OMP<NThreads,P,T,P,true> {
public:
   void apply(T* data, const unsigned long n, const unsigned long r) {
      if (n>r) {
        swap(data[n],data[r]);
        swap(data[n+1],data[r+1]);
      }
   }
};

template<unsigned int P, typename T, unsigned int I>
class GFFTswap2OMP<1,P,T,I,true> : public GFFTswap2<P,T,I> { };

template<unsigned int P, typename T>
class GFFTswap2OMP<1,P,T,P,true> : public GFFTswap2<P,T,P> { };

template<unsigned int NThreads, unsigned int P, typename T, unsigned int I>
class GFFTswap2OMP<NThreads,P,T,I,false> : public GFFTswap2<P,T,I> { };



} //namespace

#endif
