#ifndef TIMER_H
#define TIMER_H

//#include "standard.h"
//#include <unistd.h>

#include <cstdlib>
#include <cstdio>
#include <iostream>
//#include <csignal>
//#include <climits>
#include <ctime>
//#include <sys/time.h>
//#include <cerrno>
#include <cmath>

#define TIMING_INTERVAL 1

using namespace std;

typedef unsigned long ulong;

void ConvertTime(const double msec, int &hour, short &min, float &sec)
{
      ulong dt=(ulong)msec;
      sec = (dt%60) + (msec-dt);
      dt/=60;
      min=dt%60;
      dt/=60;
      hour=dt;
}

class progClock {
  double all_time;
  clock_t cur_time;
  time_t  abs_time;
public:
  progClock() { start(); }

  void start() {
      cur_time=clock();
      abs_time=time(0);
      all_time=0;
  };

  ~progClock() {};

  double GetTime(char *ch)
    {
      double d;
      clock_t tmp=clock();
      //if (sign(tmp)==sign(cur_time))
	       d = (double)(fabs((double)tmp)-fabs((double)cur_time))/CLOCKS_PER_SEC;
      /*else
	    if (tmp<0)
	       d = (double)(2*double(LONG_MAX)+tmp-cur_time)/CLOCKS_PER_SEC;
	    else
	       d = (double)(tmp+fabs(cur_time))/CLOCKS_PER_SEC;*/
      cur_time=tmp;
      all_time+=d;
      cout<<ch<<all_time<<" sec  ( ";
      int hr;
      short min;
      float sec;
      ConvertTime(all_time, hr, min, sec);
      cout<<hr<<"stu "<<min<<"min "<<sec<<"sec )"<<endl;
      return all_time;
    };

  double GetTime() {
      double d;
      clock_t tmp=clock();
      //if (sign(tmp)==sign(cur_time))
	       d = (abs(tmp)-abs(cur_time))/(double)CLOCKS_PER_SEC;
      /*else {
	       if (tmp<0)
	          d = (double)(2*double(LONG_MAX)+tmp-cur_time)/CLOCKS_PER_SEC;
	       else
    	      d = (double)(tmp+fabs(cur_time))/CLOCKS_PER_SEC;
      }*/
      cur_time=tmp;
      all_time+=fabs(d);
      return all_time;
  }

  ulong GetAbsTime(char *ch) {
      time_t t=time(0)-abs_time;
      //double d=t*1000;
      int hr;
      short min;
      float sec;
      ConvertTime(all_time, hr, min, sec);
      cout<<ch<<t<<" sec  ( ";
      cout<<hr<<"stu "<<min<<"min "<<sec<<"sec )"<<endl;
      return t;
  }

};

#endif
