// Timer and memory functions.

#define INCLUDE_MEMORY_CODE 1
//#define OUTPUT_MEMORY_PLOTPOINT

#ifndef INCLUDE_MEMORY_CODE
#include "MemoryMonitor.h"
#endif

#include <iostream>
#include <limits>
#include <malloc.h>
#include <math.h>
#include <string>
#include <sys/resource.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>

#include "timer.h"

#if defined INCLUDE_MEMORY_CODE || OUTPUT_MEMORY_PLOTPOINT
#include <sstream>
#include <fstream>
#include <execinfo.h>
#include <string.h>
#include <sys/stat.h>


/* static int timerPointCount = 1; */
/* static bool outputGnuplotTimerPoints = false; */
#endif

static double elapsedWallTime(struct timeval t1)
{
  struct timeval t2;
  gettimeofday(&t2, NULL);
    
  struct timeval dt;
  dt.tv_usec = t2.tv_usec - t1.tv_usec;
  if (dt.tv_usec < 0)
    {
      dt.tv_usec += 1000000;
      dt.tv_sec = t2.tv_sec - t1.tv_sec - 1;
    }
  else
    dt.tv_sec = t2.tv_sec - t1.tv_sec;
  return dt.tv_sec + dt.tv_usec * 1.0e-6;
}


static double elapsedCPUTime(clock_t c1)
{
  clock_t c = clock();
  // if elapsed time is too large it can overflow so we compute delta in
  // double precision and if negative add the max amount of clock ticks.
  double delta = double(c) - double(c1);
  if (delta < 0.0)
    delta += double(std::numeric_limits<clock_t>::max());
  delta /= double(CLOCKS_PER_SEC);
  return delta;
}

static double cumulativeCPUTime()
{
  clock_t c = clock();
  // if time is too large it can overflow so we use
  // double precision and if negative add the max amount of clock ticks.
  double delta = double(c);
  if (delta < 0.0)
    delta += double(std::numeric_limits<clock_t>::max());
  delta /= double(CLOCKS_PER_SEC);
  return delta;
}

namespace Util
{
  struct Timer::TimerData
  { 
    timeval t1;
    clock_t c1;
  };
  
  
  Timer::Timer(ClockType clockmode, bool printmem)
    : _clock(clockmode), _printmem(printmem)
  {
#if ! INCLUDE_MEMORY_CODE
    _printmem = printmem;
#endif
    _data = new Timer::TimerData;
    this->reset();
  }

  Timer::~Timer()
  {
    delete _data;
  }

  void Timer::reset()
  {
    gettimeofday(&_data->t1, NULL);
    _data->c1 = clock();
  }

  double Timer::getElapsedWallTime() const
  {
    return elapsedWallTime(_data->t1);
  }

  std::string Timer::getElapsedCPUTimeStr() const
  {
    double time = getElapsedWallTime();

    // Calculate seconds
    double sec = fmod(time, 60.0);

    // Calculate minutes
    int min = round(time - sec) / 60;

    // If min is zero, return seconds
    if (min == 0)
    {
      char buff[128];
      sprintf(buff, "%.2f", sec);
      return std::string(buff);
    }

    // Calculate hours
    int hour = (min - (min % 60)) / 60;
    min %= 60;

    // If hour is zero, return min and sec.
    if (hour == 0)
    {
      char buff[128];
      sprintf(buff, "%d:%05.2f", min, sec);
      return std::string(buff);
    }

    // Calculate days
    int day = (hour - (hour % 24)) / 24;
    hour %= 24;

    // If day is zero, return hour, min and sec.
    if (day == 0)
    {
      char buff[128];
      sprintf(buff, "%d:%02d:%05.2f", hour, min, sec);
      return std::string(buff);
    }

    // return day, hour, min and sec.
    char buff[128];
    sprintf(buff, "%d:%02d:%02d:%05.2f", day, hour, min, sec);
    return std::string(buff);
  }

  double Timer::getElapsedCPUTime() const
  {
    return elapsedCPUTime(_data->c1);
  }
  
  double Timer::getCumulativeCPUTime() const
  {
    return cumulativeCPUTime();
  }
  
#if ! defined INCLUDE_MEMORY_CODE
  double Timer::getMemoryUsage() const
  {
    MemoryMonitor memory;
    return memory.getVirtualMemoryUsage(MemoryMonitor::current);
  }


  void Timer::outputMemoryPlotpoint(int reportingLevel, char const* label) const
  // reportingLevel: number of levels up to report from the backtrace
  //    i.e. 0 : this routine
  //         1 : the routine calling this one
  //         2 : the routine two levels up
  //         etc.
  //
  {
#ifdef OUTPUT_MEMORY_PLOTPOINT
    char const* memMonitorPath = getenv("MEMORYMONITOR_PATH");
    if (memMonitorPath && *memMonitorPath) {
      MemoryMonitor memory;
      char const* hostName = getenv("HOSTNAME");
      int imemory = int(memory.getVirtualMemoryUsage(MemoryMonitor::current));
      int icputime = (int)(this->getCumulativeCPUTime());
      pid_t pid = getpid();
      // get the process name
      char cmd[100];
      sprintf(cmd, "ps -o comm= -p %d",pid);
      FILE *fp;      
      fp = popen(cmd, "r");  // this creates other star-ccm+ processes
      char process[1000];
      process[0] = '\0';
      if (fp) {
        if (fgets(process, sizeof(process), fp) != NULL) {
          if (process[strlen(process) - 1] == '\n')
            process[strlen(process) - 1] = '\0';
        }
        pclose(fp);
      }

      // see if the gnuplot file has been written
      if (!outputGnuplotTimerPoints) {      
        std::ostringstream gnuname;
        gnuname.str("");
        if (hostName && *hostName)
          gnuname << memMonitorPath << "/" << hostName << "/" << process
                  << "_" << pid << ".gnuplot";
        else
          gnuname << memMonitorPath << "/" << process
                  << "_" << pid << ".gnuplot";
        // see if the gnuplot file exists yet
        struct stat stFileInfo;
        int intStat = stat(gnuname.str().c_str(),&stFileInfo);
        if (intStat == 0) {
          std::ofstream gnuplotFile(gnuname.str().c_str(), std::ios_base::app);
          if (gnuplotFile.is_open()) {
            gnuplotFile.close();
            outputGnuplotTimerPoints = true;
          }
        }
      }

      // output the point to the gnuplot file
      if (outputGnuplotTimerPoints) {
        std::ostringstream gnuname;
        gnuname.str("");
        if (hostName && *hostName) 
          gnuname << memMonitorPath <<  "/" << hostName << "/" << process
                  << "_" << pid << ".gnuplot";
        else
          gnuname << memMonitorPath << "/" << process
                  << "_" << pid << ".gnuplot";
        std::ofstream gnuplotFile(gnuname.str().c_str(), std::ios_base::app);
        if (gnuplotFile.is_open()) {
          gnuplotFile << "set label \"" << timerPointCount << "\" at " <<  icputime << ","
                      << imemory << " point ps 1" <<std::endl;
          gnuplotFile.close();
        }
      }

      // output the label to the plotpointlabels file
      std::ostringstream labelname;
      labelname.str("");
      if (hostName && *hostName) 
        labelname << memMonitorPath  <<  "/" << hostName << "/" << process
                  << "_" << pid << ".plotpointlabels";
      else
        labelname << memMonitorPath << "/" << process
                  << "_" << pid << ".plotpointlabels";
      std::ofstream backtraceFile(labelname.str().c_str(), std::ios_base::app);
      if (backtraceFile.is_open()) {
        // use the provided label
        if (label) {
          backtraceFile << timerPointCount << " (" <<  icputime << ","
                        << imemory << "):  " << label << std::endl;
        }
        // otherwise use backtrace to find the name of the routine which
        // called the timer
        else {      
          int kSize=reportingLevel+1;
          if (kSize > 11) kSize = 11;
          void *traceSymbols[11];
          int num = backtrace(traceSymbols,kSize);
          char **list = backtrace_symbols(traceSymbols,num);
          if (num == kSize) {
            // strip the path from the sting
            char const* s = list[reportingLevel];
            for (char const* p = s; *p; ++p)
	      if (*p == '/')
	        s = p+1;
            backtraceFile << timerPointCount << " (" <<  icputime << ","
                          << imemory << "):  " << s << std::endl;
          }
          if (list) free(list);
        }
        backtraceFile.close();
      }
      timerPointCount++;
    }
#endif
  }
#endif

  std::ostream &operator<<(std::ostream &s, Timer const &t)
  {
#if ! defined INCLUDE_MEMORY_CODE
    double mb = 0;
    if (t.getPrintMemory()) 
      {
        MemoryMonitor memory;
        mb = memory.getVirtualMemoryUsage(MemoryMonitor::current);
      }
#endif

    // Set to fixed format with 2 digits after decimal point
    std::streamsize oldprecision = s.precision();
    s.setf(std::ios_base::fixed,std::ios_base::floatfield);
    s.precision(2);

    if (t.getClockType() & Timer::CPU)
      s << "CPU Time: " << t.getElapsedCPUTime();

    if (t.getClockType() == Timer::Both)
      s << ", ";

    if (t.getClockType() & Timer::Wall)
      s << "Wall Time: " << t.getElapsedWallTime();

#if ! defined INCLUDE_MEMORY_CODE
    if (mb > 0.0)
      s << ", Memory: " << mb << " MB";
#endif

    // Set back to default (general) format
    s.setf(std::ios_base::fmtflags(0),std::ios_base::floatfield);
    s.precision(oldprecision);

#if ! defined INCLUDE_MEMORY_CODE && defined OUTPUT_MEMORY_PLOTPOINT
    t.outputMemoryPlotpoint(2);
#endif
    
    return s;
  }         
}
