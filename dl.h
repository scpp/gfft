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

#ifndef __dl_h
#define __dl_h

namespace DFT {

///// template class DLTime
// Danielson-Lanczos section of the decimation-in-time
// FFT version

template<unsigned N, typename T=double, bool isInverse=false>
class DLTime {
   DLTime<N/2,T,isInverse> next;
   enum { S = isInverse ? -1 : 1 };
public:
   void apply(T* data) {
      next.apply(data);
      next.apply(data+N);

      T wtemp,tempr,tempi,wr,wi,wpr,wpi;
//    Change dynamic calculation to the static one
//      wtemp = sin(S*M_PI/N);
      wtemp = Sin<N,1,T>::value();
      wpr = -2.0*wtemp*wtemp;
//      wpi = -sin(2*M_PI/N);
      wpi = -S*Sin<N,2,T>::value();
      wr = 1.0;
      wi = 0.0;
      for (unsigned i=0; i<N; i+=2) {
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

template<typename T, bool isInverse>
class DLTime<4,T,isInverse> {
   enum { S = isInverse ? -1 : 1 };
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

template<typename T, bool isInverse>
class DLTime<2,T,isInverse> {
public:
   void apply(T* data) {
      T tr = data[2];
      T ti = data[3];
      data[2] = data[0]-tr;
      data[3] = data[1]-ti;
      data[0] += tr;
      data[1] += ti;
   }
};

template<typename T, bool isInverse>
class DLTime<1,T,isInverse> {
public:
   void apply(T* data) { }
};

///// template class DLFreq
// Danielson-Lanczos section of the decimation-in-frequency
// FFT version

template<unsigned N, typename T=double, bool isInverse=false>
class DLFreq {
   DLFreq<N/2,T,isInverse> next;
   enum { S = isInverse ? -1 : 1 };
public:
   void apply(T* data) {

      T wtemp,tempr,tempi,wr,wi,wpr,wpi;
//    Change dynamic calculation to the static one
//      wtemp = sin(M_PI/N);
      wtemp = Sin<N,1,T>::value();
      wpr = -2.0*wtemp*wtemp;
//      wpi = -sin(2*M_PI/N);
      wpi = -S*Sin<N,2,T>::value();
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

template<typename T, bool isInverse>
class DLFreq<4,T,isInverse> {
   enum { S = isInverse ? -1 : 1 };
public:
   void apply(T* data) {
      T tr = data[4];
      T ti = data[5];
      data[4] = data[0]-tr;
      data[5] = data[1]-ti;
      data[0] += tr;
      data[1] += ti;
      tr = data[6];
      ti = data[7];
      data[6] = S*(data[3]-ti);
      data[7] = S*(tr-data[2]);
      data[2] += tr;
      data[3] += ti;

      tr = data[2];
      ti = data[3];
      data[2] = data[0]-tr;
      data[3] = data[1]-ti;
      data[0] += tr;
      data[1] += ti;
      tr = data[6];
      ti = data[7];
      data[6] = data[4]-tr;
      data[7] = data[5]-ti;
      data[4] += tr;
      data[5] += ti;
   }
};

template<typename T, bool isInverse>
class DLFreq<2,T,isInverse> {
public:
   void apply(T* data) {
      T tr = data[2];
      T ti = data[3];
      data[2] = data[0]-tr;
      data[3] = data[1]-ti;
      data[0] += tr;
      data[1] += ti;
   }
};

template<typename T, bool isInverse>
class DLFreq<1,T,isInverse> {
public:
   void apply(T* data) { }
};

///////////////////////////

template<unsigned N, typename T=double>
class GFFTswap {
public:
   void apply(T* data) {
     int i,m,j=1;
     for (i=1; i<2*N; i+=2) {
        if (j>i) {
            std::swap(data[j-1], data[i-1]);
            std::swap(data[j], data[i]);
        }
        m = N;
        while (m>=2 && j>m) {
            j -= m;
            m >>= 1;
        }
        j += m;
     }
   }
};

}  //namespace DFT

#endif /*__dl_h*/
