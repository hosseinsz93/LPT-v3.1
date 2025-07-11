#include "GeomCheckInterrupts.h"
#if DEBUG_INTERRUPTS
#include "GeomCpuTimer.h"
#include "GeomPrintStackTrace.h"
#endif

namespace Geom
{

  CheckInterrupts* CheckInterrupts::s_checker = 0;
  long CheckInterrupts::s_count = 0;

#if DEBUG_INTERRUPTS
  bool CheckInterrupts::s_checkedTime = false;
  long CheckInterrupts::s_totalCount = 0;       // total count
  long CheckInterrupts::s_totalChecks = 0;      // total interrupt checks

  // max update interval in nano-seconds, default is 1 secs
  static CpuTime_t maximumUpdateInterval = 1000000000;
  static CpuTime_t lastUpdateCheck = 0;
  static CpuTime_t prevUpdateCheck = 0;
  static CpuTime_t maxInterval = 0;             // max interval between checks

  inline static CpuTime_t updateTime()
  {
    prevUpdateCheck = lastUpdateCheck;
    lastUpdateCheck = CpuTimeDetails().currentTime();
    if (prevUpdateCheck == 0)
      prevUpdateCheck = lastUpdateCheck;
    return lastUpdateCheck - prevUpdateCheck;
  }

  void CheckInterrupts::reportInterval(double delta, ConstWhereRef where)
  {
    std::cout << "**** INTERRUPT INTERVAL EXCEEDED by " << delta << " seconds";
    if (where.file || where.func) {
      std::cout << ", at ";
      if (where.file)
        std::cout << where.file;
      if (where.file && where.func)
        std::cout << ':';
      if (where.func)
        std::cout << where.func;
      if (where.line)
        std::cout << ':' << where.line;
    }
    std::cout << std::endl;
    if (!where.file && !where.func)
      MK::printStackTrace(2,10);
  }

  void CheckInterrupts::checkTime(long count, ConstWhereRef where)
  {
    if (!s_checkedTime) {
      s_checkedTime = true;
      s_totalCount += count;
      ++s_totalChecks;
      CpuTime_t delta = updateTime();
      if (maxInterval < delta)
        maxInterval = delta;
      if (maximumUpdateInterval < delta)
        reportInterval(CpuTimeDetails().toSeconds(delta), where);
    }
  }
#endif

  long CheckInterrupts::getTotalCount()
  {
#if DEBUG_INTERRUPTS
    return s_totalCount;
#else
    return 0;
#endif
  }

  long CheckInterrupts::getTotalChecks()
  {
#if DEBUG_INTERRUPTS
    return s_totalChecks;
#else
    return 0;
#endif
  }

  void CheckInterrupts::setMaxInterval(double secs)
  {
#if DEBUG_INTERRUPTS
    maximumUpdateInterval = CpuTimeDetails().getResolution() * secs;
#endif
  }

  double CheckInterrupts::getMaxInterval()
  {
#if DEBUG_INTERRUPTS
    return CpuTimeDetails().toSeconds(maxInterval);
#else
    return 0;
#endif
  }

}
