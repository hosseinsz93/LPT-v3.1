#ifndef _UTIL_TIMER_H
#define _UTIL_TIMER_H

#include <ostream>

namespace Util
{
  /*! The Timer class is used to keep track of CPU and/or wall time for 
    output at runtime.  */
  class Timer
  {
    public:
    enum ClockType {CPU = 0x01, Wall = 0x02, Both = 0x03};
    /*! Construct the timer.  The current time will be encoded at the time
      of construction.  All subsequent calls to getElapsed... and prints
      through the << operator will report times since the construction of
      the Timer, unless the reset() method is called. 

      The \a clock argument specifies whether CPU time, Wall time, or both
      should be reported. */
    // I am changing the default to be Both: Aly -- 11/02/2004
    Timer(ClockType clock = Both, bool printmem = false);
    
    /*! Destructor. */
    ~Timer();
        
    /*! Set the clock type to report CPU time, Wall time, or both. This
      only affects how the timer is printed to ostreams.  */
    void setClockType(ClockType clock) {
      _clock = clock;
    }
    /*! Return the current clock type. */
    ClockType getClockType() const {return _clock;}

    /*! Enable/disable printing of current memory usage after time. */      
    void setPrintMemory(bool val) {
      _printmem = val;
    }
    /*! Return status of memory printing flag. */
    bool getPrintMemory() const {return _printmem;}

    /*! Reset the clock.  This will record the current time as of the reset()
      call.  All subsequent calls to getElapsed... and prints using << will
      report elapsed time since the most recent reset() call. */
    void reset();

    /*! Get the elapsed wall time since construction of the timer (or most
      recent reset() call).  This can be called even if the clock is in CPU
      mode. */
    double getElapsedWallTime() const;
    std::string getElapsedCPUTimeStr() const;
    /*! Get the elapsed CPU time since construction of the timer (or most
      recent reset() call).  This can be called even if the clock is in wall
      mode. */
    double getElapsedCPUTime() const;
    double getCumulativeCPUTime() const;
    double getMemoryUsage() const;
    void outputMemoryPlotpoint(int reportingLevel, char const* label=NULL) const;
  private:
    ClockType _clock;
    bool _printmem;
    
    struct TimerData;    
    TimerData *_data;
  };

  std::ostream & operator<<(std::ostream &s, Timer const &t);
}

#endif


// Automatic setting of emacs local variables.
// Local Variables:
// mode: C++
// tab-width: 8
// End:
