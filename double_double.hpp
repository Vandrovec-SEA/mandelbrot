#ifndef MANDELMATH_DOUBLE_DOUBLE_HPP
#define MANDELMATH_DOUBLE_DOUBLE_HPP

namespace MandelMath {

class double_double {
public:
  double_double(): hi(0), lo(0) { }
  double hi;
  double lo;
  void chs() { hi=-hi; lo=-lo; }
  void add_double(double h2);
  void add(double h2, double l2);
  void mul(double h2, double l2);
};

} // namespace MandelMath

#endif // MANDELMATH_NUMBER_DOUBLE_HPP
