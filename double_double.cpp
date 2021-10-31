#include "double_double.hpp"
#include <cmath>

namespace MandelMath {

// https://www.davidhbailey.com/dhbsoftware/
// https://www.davidhbailey.com/dhbsoftware/qd-2.3.23.tar.gz
// https://www.davidhbailey.com/dhbsoftware/dqfun-v01.tar.gz

#define X86

void fpu_fix_start(unsigned int *old_cw) {
#ifdef X86
#if defined(_WIN32) && !defined(__GNUC__)
#ifdef __BORLANDC__
  /* Win 32 Borland C */
  unsigned short cw = _control87(0, 0);
  _control87(0x0200, 0x0300);
  if (old_cw) {
    *old_cw = cw;
  }
#else
  /* Win 32 MSVC */
  unsigned int cw = _control87(0, 0);
  _control87(0x00010000, 0x00030000);
  if (old_cw) {
    *old_cw = cw;
  }
#endif
#else
  /* Linux */
  volatile unsigned short cw, new_cw;
  //_FPU_GETCW(cw);
  asm volatile ("fnstcw %0":"=m" (cw));


  new_cw = (cw & ~0x0300) | 0x0200;
  //_FPU_SETCW(new_cw);
  asm volatile ("fldcw %0": :"m" (new_cw));

  if (old_cw) {
    *old_cw = cw;
  }
#endif
#endif
}

void fpu_fix_end(unsigned int *old_cw) {
#ifdef X86
#if defined(_WIN32) && !defined(__GNUC__)

#ifdef __BORLANDC__
  /* Win 32 Borland C */
  if (old_cw) {
    unsigned short cw = (unsigned short) *old_cw;
    _control87(cw, 0xFFFF);
  }
#else
  /* Win 32 MSVC */
  if (old_cw) {
    _control87(*old_cw, 0xFFFFFFFF);
  }
#endif

#else
  /* Linux */
  if (old_cw) {
    int cw;
    cw = *old_cw;
    //_FPU_SETCW(cw);
    asm volatile ("fldcw %0": :"m" (cw));
  }
#endif
#endif
}

#define _QD_SPLITTER 134217729.0               // = 2^27 + 1
#define _QD_SPLIT_THRESH 6.69692879491417e+299 // = 2^996

/*********** Basic Functions ************/
/* Computes fl(a+b) and err(a+b).  Assumes |a| >= |b|. */
inline void dd_real::quick_two_sum(double a, double b)
{
  hi = a + b;
  lo = b + (a - hi);
}

/* Computes fl(a-b) and err(a-b).  Assumes |a| >= |b| */
inline double quick_two_diff(double a, double b, double &err) {
  double s = a - b;
  err = (a - s) - b;
  return s;
}

/* Computes fl(a+b) and err(a+b).  */
inline void dd_real::two_sum(double a, double b)
{
  hi = a + b;
  double bb = hi - a;
  lo = (a - (hi - bb)) + (b - bb);
}

/* Computes fl(a-b) and err(a-b).  */
inline double two_diff(double a, double b, double &err) {
  double s = a - b;
  double bb = s - a;
  err = (a - (s - bb)) - (b + bb);
  return s;
}

/* Computes high word and lo word of a */
inline void dd_real::split(double a)
{
  double temp;
  if (a > _QD_SPLIT_THRESH || a < -_QD_SPLIT_THRESH) {
    a *= 3.7252902984619140625e-09;  // 2^-28
    temp = _QD_SPLITTER * a;
    hi = temp - (temp - a);
    lo = a - hi;
    hi *= 268435456.0;          // 2^28
    lo *= 268435456.0;          // 2^28
  } else {
    temp = _QD_SPLITTER * a;
    hi = temp - (temp - a);
    lo = a - hi;
  }
}

/* Computes fl(a*b) and err(a*b). */
inline void dd_real::two_prod(double a, double b)
{
  dd_real aa;
  dd_real bb;
  hi = a * b;
  aa.split(a);
  bb.split(b);
  lo = ((aa.hi * bb.hi - hi) + aa.hi * bb.lo + aa.lo * bb.hi) + aa.lo * bb.lo;
}

/* Computes fl(a*a) and err(a*a).  Faster than the above method. */
inline void dd_real::two_sqr(double a)
{
  dd_real aa;
  hi = a * a;
  aa.split(a);
  lo = ((aa.hi * aa.hi - hi) + 2.0 * aa.hi * aa.lo) + aa.lo * aa.lo;
}

//hack if you don't want to mess with FPU control word... but is a nightmare
#define MAKE_DOUBLE( tmp, val ) ( tmp=(val), *(volatile double *)&tmp )

void dd_real::lshift(int exp)
{
  //lolwut "On many implementations, std::ldexp is less efficient than multiplication or division
  //        by a power of two using arithmetic operators"
  //should I just adjust exponent in *(int *)&hi or what
  hi=ldexp(hi, exp);
  lo=ldexp(lo, exp);
}

void dd_real::add_double(double h2)
{
  dd_real se;
  se.two_sum(hi, h2);
  se.lo += lo;
  quick_two_sum(se.hi, se.lo);
}

void dd_real::add(double h2, double l2)
{
//dqfun.f90#dqadd is the same as dd_inline.h#dd_real::sloppy_add
  dd_real se;

  se.two_sum(hi, h2);
  se.lo += (lo + l2);
  quick_two_sum(se.hi, se.lo);

//1998 inline.h#operator +(const doubledouble& x is the same as dd_inline.h#dd_real::ieee_add
/*
  dd_real s;
  dd_real t;

  s.two_sum(hi, h2);
  t.two_sum(lo, l2);
  s.lo += t.hi;
  s.quick_two_sum(s.hi, s.lo);
  s.lo += t.lo;
  quick_two_sum(s.hi, s.lo);
*/
}

void dd_real::mul(double h2, double l2)
{
  dd_real p;

  p.two_prod(hi, h2);
  p.lo += (hi * l2 + lo * h2);
  quick_two_sum(p.hi, p.lo);
/*
  double hx, tx, hy, ty, C, c;
  double th, tl;
  C = MAKE_DOUBLE(th, double_double_split*hi);
  hx = MAKE_DOUBLE(tl, C-hi);
  c = MAKE_DOUBLE(th, double_double_split*h2);
  hx = MAKE_DOUBLE(tl, C-hx);
  tx = MAKE_DOUBLE(th, hi-hx);
  hy = MAKE_DOUBLE(tl, c-h2);
  C = MAKE_DOUBLE(th, hi*h2);
  hy = MAKE_DOUBLE(tl, c-hy);
  ty = MAKE_DOUBLE(th, h2-hy);
  c = MAKE_DOUBLE(tl, ((((hx*hy-C)+hx*ty)+tx*hy)+tx*ty)+(hi*l2+lo*h2));

  hi = MAKE_DOUBLE(th, C+c);
  hx = MAKE_DOUBLE(th, C-hi);
  lo = MAKE_DOUBLE(tl, c+hx);
*/
}

void dd_real::sqr()
{
  dd_real p;
  dd_real s;
  p.two_sqr(hi);
  p.lo += 2.0 * hi * lo;
  p.lo += lo * lo;
  quick_two_sum(p.hi, p.lo);
}

/*
inline doubledouble recip(const doubledouble& y) {
  x86_FIX
  double  hc, tc, hy, ty, C, c, U, u;
  C = 1.0/y.h();
  c = doubledouble::Split*C;
  hc =c-C;
  u = doubledouble::Split*y.h();
  hc = c-hc; tc = C-hc; hy = u-y.h(); U = C*y.h(); hy = u-hy; ty = y.h()-hy;
  u = (((hc*hy-U)+hc*ty)+tc*hy)+tc*ty;
  c = ((((1.0-U)-u))-C*y.l())/y.h();
  doubledouble z; z.hi = C+c; z.lo = double(C-z.hi)+c;
  END_x86_FIX
  return z;
}

inline doubledouble operator /(const doubledouble& x,const doubledouble& y ) {
  x86_FIX
  double hc, tc, hy, ty, C, c, U, u;
  C = x.hi/y.hi; c = doubledouble::Split*C; hc =c-C;  u = doubledouble::Split*y.hi; hc = c-hc;
  tc = C-hc; hy = u-y.hi; U = C * y.hi; hy = u-hy; ty = y.hi-hy;
  u = (((hc*hy-U)+hc*ty)+tc*hy)+tc*ty;
  c = ((((x.hi-U)-u)+x.lo)-C*y.lo)/y.hi;
  doubledouble z; z.hi = C+c; z.lo = double(C-z.hi)+c;
  END_x86_FIX
  return z;
}
*/

#undef MAKE_DOUBLE

} // namespace MandelMath
