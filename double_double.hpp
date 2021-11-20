#ifndef MANDELMATH_DOUBLE_DOUBLE_HPP
#define MANDELMATH_DOUBLE_DOUBLE_HPP

namespace MandelMath {

void fpu_fix_start(unsigned int *old_cw);
void fpu_fix_end(unsigned int *old_cw);

class dd_real {
protected:
  inline void quick_two_sum(double a, double b);
  inline void two_sum(double a, double b);
  inline void split(double a);
  inline void two_prod(double a, double b);
  inline void two_sqr(double a);
public:
  dd_real(): hi(0), lo(0) { }
  double hi;
  double lo;
  void assign(const dd_real &src) { hi=src.hi; lo=src.lo; }
  void chs() { hi=-hi; lo=-lo; }
  void lshift(int exp); //*=2^exp
  void add_double(double h2);
  void add(double h2, double l2);
  void mul(double h2, double l2);
  void sqr();
  void recip();
  void sqrt();
  bool isequal(const dd_real *other);
  bool is0();
  bool isle(const dd_real *other);
  bool isle0();
  bool isl1();
};

} // namespace MandelMath

#endif // MANDELMATH_NUMBER_DOUBLE_HPP
