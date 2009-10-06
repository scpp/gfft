// Some programs to test grom Numerical Recipes in C++


// Replaces data[1..2*nn] by its discrete Fourier transform.  data[] is a
// complex array
template<typename T>
void four1(T* data, unsigned long nn, const int isign)
{
    unsigned long n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    T tempr, tempi;

    n = nn << 1;
    j=1;

    for (i=1; i < n; i+= 2) {                                       //   1
        if (j > i)                                                  //   2
        {                                                           //   3
            std::swap(data[j-1], data[i-1]);                                 //   4
            std::swap(data[j], data[i]);                             //   5
        }                                                           //   6
        m = n >> 1;                                                 //   7
        while (m >= 2 && j > m) {                                   //   8
            j -= m;                                                 //   9
            m >>= 1;                                                //  10
        }                                                           //  11
        j += m;                                                     //  12
    };                                                              //  13
                                                                    //  14
    mmax=2;                                                         //  15
                                                                    //  16
    while (n > mmax) {                                              //  17
        istep = mmax << 1;                                          //  18
        theta = isign*(-6.28318530717959/mmax);                      //  19
        wtemp = sin(0.5*theta);                                     //  20
        wpr = -2.0*wtemp*wtemp;                                     //  21
        wpi = sin(theta);                                           //  22
//std::cout<<wpr<<" "<<wpi<<std::endl;
        wr = 1.0;                                                   //  23
        wi = 0.0;                                                   //  24
                                                                    //  25
       for (m=1; m < mmax; m += 2) {                               //  26
                                                                    //  27
            for (i=m; i <= n; i += istep) {                         //  28
                                                                    //  29
                j=i+mmax;                                           //  30
                tempr = wr*data[j-1] - wi*data[j];                  //  31
                tempi = wr * data[j] + wi*data[j-1];                //  32
                                                                    //  33
                                                                    //  34
                data[j-1] = data[i-1] - tempr;                          //  35
                data[j] = data[i] - tempi;                      //  36
                data[i-1] += tempr;                                   //  37
                data[i] += tempi;                                 //  38
            }                                                       //  39
              wtemp=wr;                                             //  40
              wr = wr*wpr-wi*wpi+wr;                                //  41
              wi=wi*wpr+wtemp*wpi+wi;                               //  42
        }                                                           //  43
        mmax=istep;                                                 //  44
    }                                                               //  45
}

template<typename T>
void realft(T* data, unsigned long nn, const int isign)
{
      unsigned int n,i,i1,i2,i3,i4;
      T c2,theta,wtemp,tempr,tempi,wr,wi,wpr,wpi,c1=0.5;
      T h1r,h1i,h2r,h2i;

      n = nn << 1;
      theta = -6.28318530717959/n;
      if (isign==1) {
         c2 = -0.5;
         four1(data,nn,1);
      } else {
         c2 = 0.5;
         theta = -theta;
      }
      wtemp = sin(0.5*theta);
      wpr = -2.0*wtemp*wtemp;
      wpi = sin(theta);
      wr = 1.0+wpr;
      wi = wpi;
      for (i=1; i<(n>>2); ++i) {
        i1 = i+i;
        i2 = i1+1;
        i3 = n-i1;
        i4 = i3+1;
        h1r = c1*(data[i1]+data[i3]);
        h1i = c1*(data[i2]-data[i4]);
        h2r =-c2*(data[i2]+data[i4]);
        h2i = c2*(data[i1]-data[i3]);
        data[i1] = h1r + wr*h2r - wi*h2i;
        data[i2] = h1i + wr*h2i + wi*h2r;
        data[i3] = h1r - wr*h2r + wi*h2i;
        data[i4] =-h1i + wr*h2i + wi*h2r;

        wtemp = wr;
        wr += wr*wpr - wi*wpi;
        wi += wi*wpr + wtemp*wpi;
      }
      if (isign==1) {
         h1r = data[0];
         data[0] = h1r + data[1];
         data[1] = h1r - data[1];
      } else {
         h1r = data[0];
         data[0] = c1*(h1r + data[1]);
         data[1] = c1*(h1r - data[1]);
         four1(data,nn,-1);
      }

}

