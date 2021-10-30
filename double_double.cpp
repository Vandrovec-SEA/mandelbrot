#include "double_double.hpp"
#include <cmath>

namespace MandelMath {

constexpr double double_double_split=134217729.0L; // 2^27+1, for IEEE double

#define MAKE_DOUBLE( tmp, val ) ( tmp=(val), *(volatile double *)&tmp )

void double_double::add_double(double h2)
{
  double H, h, S, s, e;
  double th;
  S = MAKE_DOUBLE(th, hi+h2);
  e = MAKE_DOUBLE(th, S-h2);
  s = MAKE_DOUBLE(th, S-e);
  s = MAKE_DOUBLE(th, (hi-e)+(h2-s));
  H = MAKE_DOUBLE(th, S+(s+lo));
  h = MAKE_DOUBLE(th, (s+lo)+(S-H));

  hi = MAKE_DOUBLE(th, H+h);
  e = MAKE_DOUBLE(th, H-hi);
  lo = MAKE_DOUBLE(th, h+e);
}

void double_double::add(double h2, double l2)
{
  double H, h, T, t, S, s, e, f;
  double th, tl;
  S = MAKE_DOUBLE(th, hi+h2);
  T = MAKE_DOUBLE(tl, lo+l2);
  e = MAKE_DOUBLE(th, S-hi);
  f = MAKE_DOUBLE(tl, T-lo);
  s = MAKE_DOUBLE(th, S-e);
  t = MAKE_DOUBLE(tl, T-f);
  s = MAKE_DOUBLE(th, (h2-e)+(hi-s));
  t = MAKE_DOUBLE(tl, (l2-f)+(lo-t));
  e = MAKE_DOUBLE(th, s+T);
  H = MAKE_DOUBLE(tl, S+e);
  h = MAKE_DOUBLE(th, e+(S-H));
  f = MAKE_DOUBLE(tl, t+h);

  hi = MAKE_DOUBLE(th, H+f);
  e  = MAKE_DOUBLE(th, H-hi);
  lo = MAKE_DOUBLE(tl, f+e);
}

void double_double::mul(double h2, double l2)
{
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
