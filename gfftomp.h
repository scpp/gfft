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

#ifndef __gfftomp_h
#define __gfftomp_h

#include "gfftalg.h"


namespace DFT {

/// OpenMP parallelized Danielson-Lanczos section of the decimation-in-time FFT version.
template<unsigned NThreads, unsigned N, typename T, int S, bool C=((N>NThreads) && (N>4))>
class InTimeOMP;

template<unsigned NThreads, unsigned N, typename T, int S>
class InTimeOMP<NThreads,N,T,S,true> {
   typedef typename TempTypeTrait<T>::Result LocalVType;
   InTimeOMP<NThreads/2,N/2,T,S> next;
   static const int IN = N;
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
      int i,chunk = N/2;

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

template<unsigned N, typename T, int S>
class InTimeOMP<1,N,T,S,true> : public InTime<N,T,S> { };

template<unsigned NThreads, unsigned N, typename T, int S>
class InTimeOMP<NThreads,N,T,S,false> : public InTime<N,T,S> { };


/// Danielson-Lanczos section of the decimation-in-frequency FFT version
template<unsigned NThreads, unsigned N, typename T, int S, bool C=((N>NThreads) && (N>4))>
class InFreqOMP;

template<unsigned NThreads, unsigned N, typename T, int S>
class InFreqOMP<NThreads,N,T,S,true> {
   typedef typename TempTypeTrait<T>::Result LocalVType;
   InFreqOMP<NThreads/2,N/2,T,S> next;
public:
   void apply(T* data) {

      LocalVType wtemp,tempr,tempi,wr,wi,wpr,wpi;

      #pragma omp parallel shared(data,wpr,wpi,t) private(wtemp,tempr,tempi,wr,wi)
      {
      wtemp = Sin<N,1,LocalVType>::value();
      wpr = -2.0*wtemp*wtemp;
      wpi = -S*Sin<N,2,LocalVType>::value();
      wr = 1.0;
      wi = 0.0;
      int i,chunk = N/2;

      #pragma omp for schedule(static,chunk)
      for (i=0; i<N; i+=2) {
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

template<unsigned N, typename T, int S>
class InTimeOMP<1,N,T,S,true> : public InTime<N,T,S> { };

template<unsigned NThreads, unsigned N, typename T, int S>
class InTimeOMP<NThreads,N,T,S,false> : public InTime<N,T,S> { };


// doesn't work
template<unsigned NThreads, unsigned N, typename T>
class GFFTswapOMP {
public:
   void apply(T* data) {
     unsigned int i,m,j=0;
     unsigned int n = N;
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

//-------------------------------------------------

template<unsigned NThreads, unsigned P, typename T,
unsigned I=0, bool C=(NThreads<(1<<P))>
class GFFTswap2OMP;

template<unsigned NThreads, unsigned P, typename T,
unsigned I>
class GFFTswap2OMP<NThreads,P,T,I,true> {
   enum { BN = 1<<(I+1), BR = 1<<(P-I) };
   GFFTswap2OMP<NThreads/2,P,T,I+1> next;
public:
   void apply(T* data, const unsigned n=0, const unsigned r=0) {
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

template<unsigned NThreads, unsigned P, typename T>
class GFFTswap2OMP<NThreads,P,T,P,true> {
public:
   void apply(T* data, const unsigned n, const unsigned r) {
      if (n>r) {
        swap(data[n],data[r]);
        swap(data[n+1],data[r+1]);
      }
   }
};

template<unsigned P, typename T, unsigned I>
class GFFTswap2OMP<1,P,T,I,true> : public GFFTswap2<P,T,I> { };

template<unsigned P, typename T>
class GFFTswap2OMP<1,P,T,P,true> : public GFFTswap2<P,T,P> { };

template<unsigned NThreads, unsigned P, typename T, unsigned I>
class GFFTswap2OMP<NThreads,P,T,I,false> : public GFFTswap2<P,T,I> { };



} //namespace

#endif
