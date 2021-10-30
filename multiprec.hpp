#ifndef MANDELMATH_MULTIPREC_HPP
#define MANDELMATH_MULTIPREC_HPP


namespace MandelMath {

class multiprec
{
public:
  multiprec();
  void set(double val);
  void chs();
  void frac_pos();
  void add_double(double x);
  void mul_double(double x);
  void add(multiprec *other);
  void mul(multiprec *other);
  int round();
};


} // namespace MandelMath

#endif // MANDELMATH_MULTIPREC_HPP
