#ifndef _GEOM_COMPILE_TIME_ASSERT
#define _GEOM_COMPILE_TIME_ASSERT

// compile-time assertion macro

template <bool p>
struct _CompileTimeAssert
{};

template <>
struct _CompileTimeAssert<true>
{
  static inline void failedAssertion() {}
};

#define CompileTimeAssert(test) _CompileTimeAssert<(test)>::failedAssertion()

#endif
