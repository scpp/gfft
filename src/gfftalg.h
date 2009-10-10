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

#ifndef __gfftalg_h
#define __gfftalg_h

/** \file
    \brief General algorithms and short-radix FFT specifications
*/

#include "metafunc.h"


namespace GFFT {


template<typename T>
inline void _spec2(T* data) {
      T tr = data[2];
      T ti = data[3];
      data[2] = data[0]-tr;
      data[3] = data[1]-ti;
      data[0] += tr;
      data[1] += ti;
}


template<typename T>
struct TempTypeTrait;

template<>
struct TempTypeTrait<float> {
   typedef double Result;
};

template<>
struct TempTypeTrait<double> {
   typedef long double Result;
};

template<typename T,
template<typename> class Complex>
struct TempTypeTrait<Complex<T> > {
   typedef typename TempTypeTrait<T>::Result Result;
};

template<typename T, typename A,
template<typename,typename> class Complex>
struct TempTypeTrait<Complex<T,A> > {
   typedef typename TempTypeTrait<T>::Result Result;
};

// template<typename T, typename A,
// template<typename,typename> class Complex>
// struct TempTypeTrait<Complex<T,A> > {
//    typedef T Result;
// };


/// Danielson-Lanczos section of the decimation-in-time FFT version
template<unsigned N, typename T, int S>
class InTime {
   typedef typename TempTypeTrait<T>::Result LocalVType;
   InTime<N/2,T,S> next;
public:
   void apply(T* data) {

      LocalVType wtemp,tempr,tempi,wr,wi,wpr,wpi,t;

      next.apply(data);
      next.apply(data+N);

//    Change dynamic calculation to the static one
//      wtemp = sin(S*M_PI/N);
      t = Sin<N,1,LocalVType>::value();
      wpr = -2.0*t*t;
//      wpi = -sin(2*M_PI/N);
      wpi = -S*Sin<N,2,LocalVType>::value();
      wr = 1.0;
      wi = 0.0;
      for (unsigned int i=0; i<N; i+=2) {
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
   }
};

/*
/// Specialization for N=4, decimation-in-time
template<typename T, int S>
class InTime<4,T,S> {
public:
   void apply(T* data) {
      T tr = data[2];
      T ti = data[3];
      data[2] = data[0]-tr;
      data[3] = data[1]-ti;
      data[0] += tr;
      data[1] += ti;
      tr = data[6];
      ti = data[7];
      data[6] = S*(data[5]-ti);
      data[7] = S*(tr-data[4]);
      data[4] += tr;
      data[5] += ti;

      tr = data[4];
      ti = data[5];
      data[4] = data[0]-tr;
      data[5] = data[1]-ti;
      data[0] += tr;
      data[1] += ti;
      tr = data[6];
      ti = data[7];
      data[6] = data[2]-tr;
      data[7] = data[3]-ti;
      data[2] += tr;
      data[3] += ti;
   }
};
*/

/// Specialization for N=2, decimation-in-time
template<typename T, int S>
class InTime<2,T,S> {
public:
   void apply(T* data) { _spec2(data); }
};

/// Specialization for N=1, decimation-in-time
template<typename T, int S>
class InTime<1,T,S> {
public:
   void apply(T* data) { }
};

/// Danielson-Lanczos section of the decimation-in-frequency FFT version
template<unsigned N, typename T, int S>
class InFreq {
   typedef typename TempTypeTrait<T>::Result LocalVType;
   InFreq<N/2,T,S> next;
public:
   void apply(T* data) {

      LocalVType wtemp,tempr,tempi,wr,wi,wpr,wpi;
//    Change dynamic calculation to the static one
//      wtemp = sin(M_PI/N);
      wtemp = Sin<N,1,LocalVType>::value();
      wpr = -2.0*wtemp*wtemp;
//      wpi = -sin(2*M_PI/N);
      wpi = -S*Sin<N,2,LocalVType>::value();
      wr = 1.0;
      wi = 0.0;
      for (unsigned i=0; i<N; i+=2) {
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

      next.apply(data);
      next.apply(data+N);
   }
};

/// Specialization for N=4, decimation-in-frequency
// template<typename T, int S>
// class InFreq<4,T,S> {
// public:
//    void apply(T* data) {
//       T tr = data[4];
//       T ti = data[5];
//       data[4] = data[0]-tr;
//       data[5] = data[1]-ti;
//       data[0] += tr;
//       data[1] += ti;
//       tr = data[6];
//       ti = data[7];
//       data[6] = S*(data[3]-ti);
//       data[7] = S*(tr-data[2]);
//       data[2] += tr;
//       data[3] += ti;
// 
//       tr = data[2];
//       ti = data[3];
//       data[2] = data[0]-tr;
//       data[3] = data[1]-ti;
//       data[0] += tr;
//       data[1] += ti;
//       tr = data[6];
//       ti = data[7];
//       data[6] = data[4]-tr;
//       data[7] = data[5]-ti;
//       data[4] += tr;
//       data[5] += ti;
//    }
// };

/// Specialization for N=2, decimation-in-frequency
template<typename T, int S>
class InFreq<2,T,S> {
public:
   void apply(T* data) { _spec2(data); }
};

/// Specialization for N=1, decimation-in-frequency
template<typename T, int S>
class InFreq<1,T,S> {
public:
   void apply(T* data) { }
};


/// Binary reordering of array elements
template<unsigned N, typename T>
class GFFTswap {
public:
   void apply(T* data) {
     unsigned int m,j=0;
     for (unsigned int i=0; i<2*N-1; i+=2) {
        if (j>i) {
            std::swap(data[j], data[i]);
            std::swap(data[j+1], data[i+1]);
        }
        m = N;
        while (m>=2 && j>=m) {
            j -= m;
            m >>= 1;
        }
        j += m;
     }
   }
};


//-------------------------------------------------

template<unsigned P, typename T=double,
unsigned I=0>
class GFFTswap2 {
   enum { BN = 1<<(I+1), BR = 1<<(P-I) };
   GFFTswap2<P,T,I+1> next;
public:
   void apply(T* data, unsigned n=0, unsigned r=0) {
      next.apply(data,n,r);
      next.apply(data,n|BN,r|BR);
   }
};

template<unsigned P, typename T>
class GFFTswap2<P,T,P> {
public:
   void apply(T* data, unsigned n=0, unsigned r=0) {
      if (n>r) {
        swap(data[n],data[r]);
        swap(data[n+1],data[r+1]);
      }
   }
};


/// Reordering of data for real-valued transforms
template<unsigned N, typename T, int S>
class Separate {
   typedef typename TempTypeTrait<T>::Result LocalVType;
   static const int M = (S==1) ? 2 : 1;
public:
   void apply(T* data) {
      unsigned int i,i1,i2,i3,i4;
      LocalVType wtemp,wr,wi,wpr,wpi;
      LocalVType h1r,h1i,h2r,h2i,h3r,h3i;
      wtemp = Sin<2*N,1,LocalVType>::value();
      wpr = -2.*wtemp*wtemp;
      wpi = -S*Sin<N,1,LocalVType>::value();
      wr = 1.+wpr;
      wi = wpi;
      for (i=1; i<N/2; ++i) {
        i1 = i+i;
        i2 = i1+1;
        i3 = 2*N-i1;
        i4 = i3+1;
        h1r = 0.5*(data[i1]+data[i3]);
        h1i = 0.5*(data[i2]-data[i4]);
        h2r = S*0.5*(data[i2]+data[i4]);
        h2i =-S*0.5*(data[i1]-data[i3]);
        h3r = wr*h2r - wi*h2i;
        h3i = wr*h2i + wi*h2r;
        data[i1] = h1r + h3r;
        data[i2] = h1i + h3i;
        data[i3] = h1r - h3r;
        data[i4] =-h1i + h3i;

        wtemp = wr;
        wr += wr*wpr - wi*wpi;
        wi += wi*wpr + wtemp*wpi;
      }
      h1r = data[0];
      data[0] = M*0.5*(h1r + data[1]);
      data[1] = M*0.5*(h1r - data[1]);

      if (N>1) data[N+1] = -data[N+1];
   }
};


}  //namespace DFT

#endif /*__gfftalg_h*/
