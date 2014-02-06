/***************************************************************************
 *   Copyright (C) 2006-2014 by Vladimir Mirnyy                            *
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

#ifndef __gffswap_h
#define __gffswap_h

/** \file
    \brief Reordering of data for FFT
*/


namespace GFFT {

using namespace MF;



/// Binary reordering of array elements
/*!
\tparam M^P length of the data
\tparam T value type

This is C-like implementation. It has been written very 
similar to the one presented in the book 
"Numerical recipes in C++".
\sa GFFTswap2
*/
template<uint_t M, uint_t P, typename T>
class SwapNR {
  static const uint_t N = IPow<M,P>::value;
public:
   void apply(T* data) {
     int_t m, j = 0;
     for (int_t i=0; i<2*N-1; i+=2) {
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


/// Binary reordering of array elements
/*!
\tparam N length of the data
\tparam T value type

This is second version of binary reordering.
It is based on template class recursion,
where member function apply() is called twice recursively
for the left and the right half-tree
building the parameters n and r, which are the
indexes of the exchanged data values.
\sa SwapNR, GFFTswap2OMP
*/
template<uint_t M, uint_t P, typename T, uint_t I=0>
class GFFTswap2 {
   static const int_t BN = IPow<M,I>::value;
   static const int_t BR = IPow<M,P-I-1>::value;
   GFFTswap2<M,P,T,I+1> next;
public:
   void apply(T* data, int_t n=0, int_t r=0) {
     const int_t qn = n/BN;
     const int_t rn = n%BN;
     const int_t qr = r/BR;
     const int_t rr = r%BR;
     for (uint_t i = 0; i < M; ++i) { 
       next.apply(data,(qn+i)*BN+rn,(qr+i)*BR+rr);
     }
   }
};

template<uint_t M, uint_t P, typename T>
class GFFTswap2<M,P,T,P> { 
public:
   void apply(T* data, int_t n=0, int_t r=0) {
      if (n>r) {
	const int_t n2 = 2*n;
	const int_t r2 = 2*r;
        std::swap(data[n2],data[r2]);
        std::swap(data[n2+1],data[r2+1]);
      }
   }
};

template<uint_t P, typename T, uint_t I>
class GFFTswap2<2,P,T,I> {
   static const int_t BN = 1<<(I+1);
   static const int_t BR = 1<<(P-I);
   GFFTswap2<2,P,T,I+1> next;
public:
   void apply(T* data, int_t n=0, int_t r=0) {
      next.apply(data,n,r);
      next.apply(data,n|BN,r|BR);
   }
};

template<uint_t P, typename T>
class GFFTswap2<2,P,T,P> {
public:
   void apply(T* data, int_t n=0, int_t r=0) {
      if (n>r) {
        std::swap(data[n],data[r]);
        std::swap(data[n+1],data[r+1]);
      }
   }
};

/// Binary reordering of array elements
/*!
\tparam M^P length of the data
\tparam T value type
*/
template<uint M, uint P, typename T,
template<typename> class Complex>
class SwapNR<M,P,Complex<T> > 
{
  static const uint N = IPow<M,P>::value;
public:
   void apply(Complex<T>* data) {
     int_t m,j=0;
     for (int_t i=0; i<N; ++i) {
        if (j>i) {
            std::swap(data[j], data[i]);
        }
        m = N/2;
        while (m>=1 && j>=m) {
            j -= m;
            m >>= 1;
        }
        j += m;
     }
   }
};

template<uint_t M, uint_t P, typename T, uint_t I,
template<typename> class Complex>
class GFFTswap2<M,P,Complex<T>,I> {
   static const int_t BN = IPow<M,I>::value;
   static const int_t BR = IPow<M,P-I-1>::value;
   GFFTswap2<M,P,Complex<T>,I+1> next;
public:
   void apply(Complex<T>* data, int_t n=0, int_t r=0) {
     const int_t qn = n/BN;
     const int_t rn = n%BN;
     const int_t qr = r/BR;
     const int_t rr = r%BR;
     for (uint_t i = 0; i < M; ++i) { 
       next.apply(data,(qn+i)*BN+rn,(qr+i)*BR+rr);
     }
   }
};

template<uint_t M, uint_t P, typename T,
template<typename> class Complex>
class GFFTswap2<M,P,Complex<T>,P> { 
public:
   void apply(Complex<T>* data, int_t n=0, int_t r=0) {
      if (n>r) 
        std::swap(data[n],data[r]);
   }
};

template<uint_t P, typename T, uint_t I,
template<typename> class Complex>
class GFFTswap2<2,P,Complex<T>,I> {
   static const int_t BN = 1<<I;
   static const int_t BR = 1<<(P-I-1);
   GFFTswap2<2,P,Complex<T>,I+1> next;
public:
   void apply(Complex<T>* data, int_t n=0, int_t r=0) {
      next.apply(data,n,r);
      next.apply(data,n|BN,r|BR);
   }
};

template<uint_t P, typename T,
template<typename> class Complex>
class GFFTswap2<2,P,Complex<T>,P> {
public:
   void apply(Complex<T>* data, int_t n=0, int_t r=0) {
      if (n>r)
        swap(data[n],data[r]);
   }
};


/// Reordering of data for real-valued transforms
/*!
\tparam N length of the data
\tparam T value type
\tparam S sign of the transform: 1 - forward, -1 - backward
*/
template<int_t N, typename T, int S>
class Separate {
   typedef typename TempTypeTrait<T>::Result LocalVType;
   static const int M = (S==1) ? 2 : 1;
public:
   void apply(T* data) {
      int_t i,i1,i2,i3,i4;
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
   
   void apply(const T*, T*) 
   { 
      std::cout << "Not implemented!" << std::endl;
      exit(1);
   }
};


/// Reordering of data for real-valued transforms
/*!
\tparam N length of the data
\tparam T value type
\tparam S sign of the transform: 1 - forward, -1 - backward
*/
template<int_t N, typename T, int S,
template<typename> class Complex>
class Separate<N,Complex<T>,S> {
   typedef typename TempTypeTrait<T>::Result LocalVType;
   typedef Complex<LocalVType> LocalComplex;
   static const int M = (S==1) ? 2 : 1;
public:
   void apply(Complex<T>* data) {
      int_t i,i1;
      LocalComplex h1,h2,h3;
      LocalVType wtemp = Sin<2*N,1,LocalVType>::value();
      LocalComplex wp(-2.*wtemp*wtemp,-S*Sin<N,1,LocalVType>::value());
      LocalComplex w(1.+wp.real(),wp.imag());

      for (i=1; i<N/2; ++i) {
        i1 = N-i;
        h1 = Complex<LocalVType>(static_cast<LocalVType>(0.5*(data[i].real()+data[i1].real())),
                                 static_cast<LocalVType>(0.5*(data[i].imag()-data[i1].imag())));
        h2 = Complex<LocalVType>(static_cast<LocalVType>( S*0.5*(data[i].imag()+data[i1].imag())),
                                 static_cast<LocalVType>(-S*0.5*(data[i].real()-data[i1].real())));
        h3 = w*h2;
        data[i] = h1 + h3;
        data[i1]= h1 - h3;
        data[i1] = Complex<T>(data[i1].real(), -data[i1].imag());

        w += w*wp;
      }
      wtemp = data[0].real();
      data[0] = Complex<T>(M*0.5*(wtemp + data[0].imag()), M*0.5*(wtemp - data[0].imag()));

      data[N/2] = Complex<T>(data[N/2].real(), -data[N/2].imag());
   }

   void apply(const Complex<T>*, Complex<T>*) 
   { 
      std::cout << "Not implemented!" << std::endl;
      exit(1);
   }
};

// Policy for a definition of forward FFT
template<int_t N, typename T>
struct Forward {
   enum { Sign = 1 };
   void apply(T*) { }
   void apply(const T*, T*) { }
   void apply(const T*, T*, T*) { }
};

template<int_t N, typename T,
template<typename> class Complex>
struct Forward<N,Complex<T> > {
   enum { Sign = 1 };
   void apply(Complex<T>*) { }
   void apply(const Complex<T>*, Complex<T>*) { }
   void apply(const Complex<T>*, Complex<T>*, Complex<T>*) { }
};

// Policy for a definition of backward FFT
template<int_t N, typename T>
struct Backward {
   enum { Sign = -1 };
   void apply(T* data) {
      for (T* i=data; i<data+2*N; ++i) *i/=N;
   }
   void apply(const T*, T* dst) { apply(dst); }
   void apply(const T*, T* dst, T*) { apply(dst); }
};

template<int_t N, typename T,
template<typename> class Complex>
struct Backward<N,Complex<T> > {
   enum { Sign = -1 };
   void apply(Complex<T>* data) {
      for (int_t i=0; i<N; ++i) {
        data[i]/=N;
      }
   }
   void apply(const Complex<T>*, Complex<T>* dst) { apply(dst); }
   void apply(const Complex<T>*, Complex<T>* dst, Complex<T>*) { apply(dst); }
};


}  //namespace DFT

#endif /*__gfftalg_h*/
