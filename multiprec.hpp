#ifndef MANDELMATH_MULTIPREC_HPP
#define MANDELMATH_MULTIPREC_HPP


namespace MandelMath {

class multiprec
{
public:
  multiprec();
  void set(double val);
  void chs();
  void lshift(int shoft);
  void frac_pos();
  void add_double(double x);
  void mul_double(double x);
  void add(multiprec *other);
  void mul(multiprec *other);
  void sqr();
  void recip();
  void sqrt();
  int round();
  bool isequal(const multiprec *other);
  bool is0();
  bool isle(const multiprec *other);
  bool isle0();
  bool isl1();
  double toDouble();
};


} // namespace MandelMath

#endif // MANDELMATH_MULTIPREC_HPP
