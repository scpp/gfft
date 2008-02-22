/***************************************************************************
 *   Copyright (C) 2007 by Volodymyr Myrnyy                                *
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
/*

 Generic simple and efficient Fast Fourier Transform (FFT) implementation
 using template metaprogramming

 ***************************************************************************/

#ifndef __gfft_h
#define __gfft_h

#include "loki/Typelist.h"
#include "loki/Factory.h"
#include "loki/Singleton.h"
//#include "loki/Functor.h"

#include "finit.h"
#include "metafunc.h"
#include "dl.h"

namespace DFT {


////// template class AbstractFFT
// abstract base class to build an object factory

template<typename T>
class AbstractFFT {
public:
   virtual void fft(T*) = 0;
};

////// template class Empty
// abstract empty base class to pass instead of AbstractFFT
// (if object factory is not needed)
// and to avoid a virtual function call penalty

class Empty { };


//static Loki::Factory<AbstractFFT<double>,unsigned int> gfft;

//typedef Loki::Factory<AbstractFFT<double>,unsigned int,Loki::Functor<AbstractFFT<double>*> > GFFT_FactoryD;
//typedef Loki::Factory<AbstractFFT<double>,unsigned int> GFFT_FactoryD;

//typedef Loki::Factory<AbstractFFT<float>,unsigned int> GFFT_FactoryF;


////// template class GFFT
// generic fast Fourier transform
// with decimation-in-time

template<unsigned P, typename T=double,
class FactoryPolicy=Empty>
class GFFTt:public FactoryPolicy {
   enum { N = 1<<P };
   GFFTswap<N,T> scramble;
   DLTime<N,T> recursion;
public:
   enum { id = P };
   static FactoryPolicy* Create() {
      return new GFFTt<P,T,FactoryPolicy>();
   }
   void fft(T* data) {
      scramble.apply(data);
      recursion.apply(data);
   }
};

////// template class GFFT
// generic fast Fourier transform
// with decimation-in-frequency

template<unsigned P, typename T=double,
class FactoryPolicy=Empty>
class GFFTf:public FactoryPolicy {
   enum { N = 1<<P };
   GFFTswap<N,T> scramble;
   DLFreq<N,T> recursion;
public:
   enum { id = P };
   static FactoryPolicy* Create() {
      return new GFFTf<P,T,FactoryPolicy>();
   }
   void fft(T* data) {
      recursion.apply(data);
      scramble.apply(data);
   }
};

////// template class IGFFT
// inverse generic fast Fourier transform
// with decimation-in-time

template<unsigned P, typename T=double,
class FactoryPolicy=Empty>
class IGFFTt:public FactoryPolicy {
   enum { N = 1<<P };
   GFFTswap<N,T> scramble;
   DLTime<N,T,true> recursion;
public:
   enum { id = P };
   static FactoryPolicy* Create() {
      return new IGFFTt<P,T,FactoryPolicy>();
   }
   void fft(T* data) {
      scramble.apply(data);
      recursion.apply(data);
      for (T* i=data; i<data+2*N; ++i) *i/=N;
   }
};

////// template class IGFFT
// inverse generic fast Fourier transform
// with decimation-in-frequency

template<unsigned P, typename T=double,
class FactoryPolicy=Empty>
class IGFFTf:public FactoryPolicy {
   enum { N = 1<<P };
   GFFTswap<N,T> scramble;
   DLFreq<N,T,true> recursion;
public:
   enum { id = P };
   static FactoryPolicy* Create() {
      return new IGFFTf<P,T,FactoryPolicy>();
   }
   void fft(T* data) {
      recursion.apply(data);
      scramble.apply(data);
      for (T* i=data; i<data+2*N; ++i) *i/=N;
   }
};

////// Separate
template<unsigned N, typename T=double, bool isInverse=false>
class Separate {
   enum { S = isInverse ? -1 : 1 };
   enum { M = isInverse ? 1 : 2 };
public:
   void apply(T* data) {
      unsigned int i,i1,i2,i3,i4;
      T wtemp,tempr,tempi,wr,wi,wpr,wpi;
      T h1r,h1i,h2r,h2i;
      wtemp = Sin<2*N,1,T>::value();
      wpr = -2.0*wtemp*wtemp;
      wpi = -S*Sin<N,1,T>::value();
      wr = 1.0+wpr;
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
        data[i1] = h1r + wr*h2r - wi*h2i;
        data[i2] = h1i + wr*h2i + wi*h2r;
        data[i3] = h1r - wr*h2r + wi*h2i;
        data[i4] =-h1i + wr*h2i + wi*h2r;

        wtemp = wr;
        wr += wr*wpr - wi*wpi;
        wi += wi*wpr + wtemp*wpi;
      }
      h1r = data[0];
      data[0] = M*0.5*(h1r + data[1]);
      data[1] = M*0.5*(h1r - data[1]);
   }
};

////// template class RVGFFTf
// Real-valued generic fast Fourier transform in-place
// with decimation-in-frequency

template<unsigned P, typename T=double,
class FactoryPolicy=Empty>
class RVGFFTf:public FactoryPolicy {
   enum { N = 1<<P };
/*   GFFTswap<N,T> scramble;
   DLFreq<N,T> recursion;*/
   GFFTf<P,T,Empty> gfft;
   Separate<N,T> sep;
public:
   enum { id = P };
   static FactoryPolicy* Create() {
      return new RVGFFTf<P,T,FactoryPolicy>();
   }
   void fft(T* data) {
      gfft.fft(data);
      sep.apply(data);
   }
};

////// template class IRVGFFTf
// inverse Real-valued generic fast Fourier transform in-place
// with decimation-in-frequency

template<unsigned P, typename T=double,
class FactoryPolicy=Empty>
class IRVGFFTf:public FactoryPolicy {
   enum { N = 1<<P };
/*   GFFTswap<N,T> scramble;
   DLFreq<N,T> recursion;*/
   IGFFTf<P,T,Empty> igfft;
   Separate<N,T,true> sep;
public:
   enum { id = P };
   static FactoryPolicy* Create() {
      return new IRVGFFTf<P,T,FactoryPolicy>();
   }
   void fft(T* data) {
      sep.apply(data);
      igfft.fft(data);
   }
};


////// template class GFFTList
// create a Typelist of given FFT classes to pass it
// to the initialization class FactoryInit

template<
template<unsigned,class,class> class FFT,
unsigned Begin, unsigned End,
typename T=double,
class FactoryPolicy=AbstractFFT<T> >
struct GFFTList {
   typedef Loki::Typelist<FFT<Begin,T,FactoryPolicy>,
           typename GFFTList<FFT,Begin+1,End,T,
                    FactoryPolicy>::Result> Result;
};

template<
template<unsigned,class,class> class FFT,
unsigned End, typename T, class FactoryPolicy>
struct GFFTList<FFT,End,End,T,FactoryPolicy> {
   typedef Loki::NullType Result;
};


// GFFT Singleton holds a GFFT object factory

template<
template<unsigned,class,class> class FFT,
unsigned Min, unsigned Max, typename T>
class GFFT_Singleton
: public Loki::SingletonHolder<Loki::Factory<AbstractFFT<T>,unsigned int> > {
   typedef Loki::Factory<AbstractFFT<T>,unsigned int> GFFT_Factory;
   typedef Loki::SingletonHolder<GFFT_Factory> Base;
   typedef FactoryInit<typename GFFTList<FFT,Min,Max,T>::Result> FInit;

   //Protection
   //GFFT_Singleton();
public:
   static GFFT_Factory& Instance() {
      FInit::apply(Base::Instance());
      return Base::Instance();
   }
};


}  // namespace DFT

#endif /*__gfft_h*/
