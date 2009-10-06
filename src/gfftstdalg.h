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

#ifndef __gfftstdalg_h
#define __gfftstdalg_h

/** \file
    \brief General algorithms and short-radix FFT specializations
           using standard classes complex<> and vector<>
*/


namespace GFFT {


template<typename T,
template<typename> class Complex>
inline void _spec2(Complex<T>* data) {
      Complex<T> t(data[1]);
      data[1] = data[0]-t;
      data[0] += t;
}

/// Danielson-Lanczos section of the decimation-in-time FFT version
/// if data is given as array of complex<>
template<unsigned N, typename T, int S,
template<typename> class Complex>
class InTime<N,Complex<T>,S> {
   typedef typename TempTypeTrait<T>::Result LocalVType;
   static const unsigned N2 = N/2;
   InTime<N2,Complex<T>,S> next;
public:
   void apply(Complex<T>* data) {
      next.apply(data);
      next.apply(data+N2);

      LocalVType wtemp = Sin<N,1,LocalVType>::value();
      Complex<T> temp;
      Complex<LocalVType> w(1.,0.);
      Complex<LocalVType> wp(-2.0*wtemp*wtemp,-S*Sin<N,2,LocalVType>::value());

      for (unsigned i=0; i<N2; ++i) {
        // rewritten componentwise because of the different types of the components
        temp.real() = data[i+N2].real()*w.real() - data[i+N2].imag()*w.imag();
        temp.imag() = data[i+N2].real()*w.imag() + data[i+N2].imag()*w.real();
        //temp = data[i+N2]*Complex<T>(w);
        data[i+N2] = data[i]-temp;
        data[i] += temp;

        w += w*wp;
      }
   }
};

/// Specialization for N=4, decimation-in-time
template<typename T, int S,
template<typename> class Complex>
class InTime<4,Complex<T>,S> {
public:
   void apply(Complex<T>* data) {
      Complex<T> t(data[1]);
      data[1] = data[0]-t;
      data[0] += t;
      t = data[3];
      data[3] = Complex<T>(S*(data[2].imag()-t.imag()),S*(t.real()-data[2].real()));
      data[2] += t;

      t = data[2];
      data[2] = data[0]-t;
      data[0] += t;
      t = data[3];
      data[3] = data[1]-t;
      data[1] += t;
   }
};

/// Specialization for N=2, decimation-in-time
template<typename T, int S,
template<typename> class Complex>
class InTime<2,Complex<T>,S> {
public:
   void apply(Complex<T>* data) { _spec2(data); }
};

/// Specialization for N=1, decimation-in-time
template<typename T, int S,
template<typename> class Complex>
class InTime<1,Complex<T>,S> {
public:
   void apply(Complex<T>* data) { }
};

/// Danielson-Lanczos section of the decimation-in-frequency FFT version
template<unsigned N, typename T, int S,
template<typename> class Complex>
class InFreq<N,Complex<T>,S> {
   typedef typename TempTypeTrait<T>::Result LocalVType;
   static const unsigned N2 = N/2;
   InFreq<N2,Complex<T>,S> next;
public:
   void apply(Complex<T>* data) {

      LocalVType wtemp = Sin<N,1,LocalVType>::value();
      Complex<T> temp;
      Complex<LocalVType> w(1.,0.);
      Complex<LocalVType> wp(-2.0*wtemp*wtemp,-S*Sin<N,2,LocalVType>::value());

      for (unsigned i=0; i<N2; ++i) {
        temp = data[i]-data[i+N2];
        data[i] += data[i+N2];
        // rewritten componentwise because of the different types of the components
        data[i+N2].real() = temp.real()*w.real() - temp.imag()*w.imag();
        data[i+N2].imag() = temp.imag()*w.real() + temp.real()*w.imag();
        //data[i+N2] = temp*w;

        w += w*wp;
      }

      next.apply(data);
      next.apply(data+N2);
   }
};

/// Specialization for N=4, decimation-in-frequency
template<typename T, int S,
template<typename> class Complex>
class InFreq<2,Complex<T>,S> {
public:
   void apply(Complex<T>* data) { _spec2(data); }
};

/// Specialization for N=1, decimation-in-frequency
template<typename T, int S,
template<typename> class Complex>
class InFreq<1,Complex<T>,S> {
public:
   void apply(Complex<T>* data) { }
};

/// Binary reordering of array elements
template<unsigned N, typename T,
template<typename> class Complex>
class GFFTswap<N,Complex<T> > {
public:
   void apply(Complex<T>* data) {
     unsigned int m,j=0;
     for (unsigned int i=0; i<N; ++i) {
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

/// Reordering of data for real-valued transforms
template<unsigned N, typename T, int S,
template<typename> class Complex>
class Separate<N,Complex<T>,S> {
   typedef typename TempTypeTrait<T>::Result LocalVType;
   typedef Complex<LocalVType> LocalComplex;
   static const int M = (S==1) ? 2 : 1;
public:
   void apply(Complex<T>* data) {
      unsigned int i,i1;
      LocalComplex h1,h2,h3;
      LocalVType wtemp = Sin<2*N,1,LocalVType>::value();
      LocalComplex wp(-2.*wtemp*wtemp,-S*Sin<N,1,LocalVType>::value());
      LocalComplex w(1.+wp.real(),wp.imag());

      for (i=1; i<N/2; ++i) {
        i1 = N-i;
        h1.real() = 0.5*(data[i].real()+data[i1].real());
        h1.imag() = 0.5*(data[i].imag()-data[i1].imag());
        h2.real() = S*0.5*(data[i].imag()+data[i1].imag());
        h2.imag() =-S*0.5*(data[i].real()-data[i1].real());
        h3 = w*h2;
        data[i] = h1 + h3;
        data[i1]= h1 - h3;
        data[i1].imag() = -data[i1].imag();

        w += w*wp;
      }
      wtemp = data[0].real();
      data[0].real() = M*0.5*(wtemp + data[0].imag());
      data[0].imag() = M*0.5*(wtemp - data[0].imag());

      data[N/2].imag() = -data[N/2].imag();
   }
};

} //namespace

#endif
