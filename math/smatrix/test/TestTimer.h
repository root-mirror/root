#ifndef TESTTIMER_H
#define TESTTIMER_H

// simple class to measure time

#include "TStopwatch.h"
#include <string>

namespace ROOT {

   namespace Math{

      namespace test {

#ifdef REPORT_TIME
         void reportTime( std::string s, double time);
#endif

         void printTime(TStopwatch & time, std::string s) {
            int pr = std::cout.precision(8);
            std::cout << s << "\t" << " time = " << time.RealTime() << "\t(sec)\t"
            //    << time.CpuTime()
            << std::endl;
            std::cout.precision(pr);
         }



         class Timer {

         public:

            Timer(const std::string & s = "") : fName(s), fTime(0)
            {
               fWatch.Start();
            }
            Timer(double & t, const std::string & s = "") : fName(s), fTime(&t)
            {
               fWatch.Start();
            }

            ~Timer() {
               fWatch.Stop();
               printTime(fWatch,fName);
#ifdef REPORT_TIME
               // report time
               reportTime(fName, fWatch.RealTime() );
#endif
               if (fTime) *fTime += fWatch.RealTime();
            }


         private:

            std::string fName;
            double * fTime;
            TStopwatch fWatch;

         };
      }

   }
}

#endif
