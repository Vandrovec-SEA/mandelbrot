#include "multiprec.hpp"
#include "MandelMath.hpp" //dbgPoint

namespace MandelMath {

multiprec::multiprec()
{
  dbgPoint();
}

void multiprec::set(double val)
{
  (void)val;
}

void multiprec::chs()
{

}

void multiprec::lshift(int shoft)
{

}

void multiprec::frac_pos()
{

}

void multiprec::add_double(double x)
{
  (void)x;
}

void multiprec::mul_double(double x)
{
  (void)x;
}

void multiprec::add(multiprec *other)
{
  (void)other;
}

void multiprec::mul(multiprec *other)
{
  (void)other;
}

void multiprec::sqr()
{
  mul(this);
}

int multiprec::round()
{
  return 0;
}

double multiprec::toDouble()
{
  return 0;
}

} // namespace MandelMath
