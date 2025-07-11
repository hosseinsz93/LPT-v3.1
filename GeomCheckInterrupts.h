#ifndef _GEOM_INTERRUPTS
#define _GEOM_INTERRUPTS

// set DEBUG_INTERRUPTS to 1 to report on excessive time intervals
// between interrupt checks
#define DEBUG_INTERRUPTS 0

namespace Geom
{

  // Check for interrupt and abort without need for a Task object.
  // Reduces the number of calls to updateCurrentTask,
  // useful in loops where a check doesn't need to be done every iteration,
  // but the loop size can be large so some checks are desirable.
  // interval is how often an actual check is done.  The proper count is
  // maintained between instantiations and with nesting or recursion.
  class CheckInterrupts {
    struct CheckTime;

  public:
    // standard check interval frequencies
    enum Frequency {
      everytime	= 1,
      every10	= 10,
      veryhigh	= 100,
      high	= 1000,
      med	= 10000,
      low	= 100000,
      verylow	= 1000000
    };

    // for reporting location of function call when using DEBUG_INTERRUPTS
    struct Where;
    typedef Where const& ConstWhereRef;

    CheckInterrupts(long interval=everytime) :
      _count(restoreCount()),
      _interval(interval) {
      pushChecker(this);
    }
    ~CheckInterrupts() {
      popChecker();
      storeCount(_count);
    }

    //! perform interrupt check
    bool forceCheck() {
      return forceCheck(Where());
    }

    bool forceCheck(ConstWhereRef where) {
      CheckTime check(_count, where);
      _count = 0;
      return true;
    }
    //! perform interrupt check if interval exceeded, true if checked
    //! A Task object must be active for this call to ever abort.
    bool check() {
      return check(Where());
    }

    bool check(ConstWhereRef where) {
      return ((_interval <= ++_count) ? forceCheck(where) : false);
    }

    //! force interrupt check using current (top) checker in stack
    static bool forceCheckCurrent(ConstWhereRef where=Where()) {
      return (s_checker ? s_checker->forceCheck(where) : false);
    }
    //! perform interrupt check using current (top) checker in stack
    static bool checkCurrent(ConstWhereRef where=Where()) {
      return (s_checker ? s_checker->check(where) : false);
    }

    //! change check interval
    void setInterval(long interval)	{ _interval = interval; }
    //! get check interval
    long getInterval() const		{ return _interval; }
    //! get current count
    long getCount() const		{ return _count; }

    // these only work if DEBUG_INTERRUPTS is non-zero
    static long getTotalCount();
    static long getTotalChecks();
    static void setMaxInterval(double secs);
    static double getMaxInterval();

    struct Where {
#if DEBUG_INTERRUPTS
      const char* const file;
      const char* const func;
      int line;
      Where() : file(0), func(0), line(0) {}
      Where(const char* const fi, int li)
        : file(fi), func(0), line(li) {}
      Where(const char* const fi, const char* const fu, int li)
        : file(fi), func(fu), line(li) {}
#endif
    };

  private:
#if DEBUG_INTERRUPTS
    struct CheckTime {
      CheckTime(long count, ConstWhereRef where) { checkTime(count, where); }
      ~CheckTime() { s_checkedTime = false; }
    };
    static bool s_checkedTime;		// has time been checked?
    static long s_totalCount;           // total count
    static long s_totalChecks;          // total interrupt checks
    static void checkTime(long count=1, ConstWhereRef where=Where());
    static void reportInterval(double delta, ConstWhereRef where);
#else
    struct CheckTime {
      CheckTime(long, ConstWhereRef) {}
    };
    static void checkTime(long count=0, ConstWhereRef=Where()) {}
#endif
    friend class Task;

    CheckInterrupts* _checker;		// previous on checker stack
    long _count;			// local count
    long _interval;			// check interval
    static CheckInterrupts* s_checker;	// top of checker stack
    static long s_count;		// saved count

    long restoreCount() const {
      return (s_checker ? s_checker->_count : s_count);
    }
    void storeCount(long count) {
      if (s_checker)
	s_checker->_count = count;
      else
	s_count = count;
    }
    void pushChecker(CheckInterrupts* checker) {
      _checker  = s_checker;
      s_checker = checker;
    }
    void popChecker() {
      s_checker = _checker;
    }
  };

}

#if DEBUG_INTERRUPTS
// To use the Where object, when calling a function that takes this type,
// put WHEREAMI in the arg list, like so: interrupt.check(WHEREAMI);
// When a report is given it will then give the location of the call.
#include "mk/common/FunctionName.h"
#  define WHEREAMI CheckInterrupts::Where(__FILE__, FUNCTION_NAME, __LINE__)
#else
#  define WHEREAMI
#endif


#endif
