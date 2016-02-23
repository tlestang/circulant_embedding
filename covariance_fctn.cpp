#include <cmath>

double f(int x, int y)
{
  double d, a;
  double sigma = 32.0;
  d = x*x + y*y;
  //a = (1.0/(sqrt(2.0*M_PI)*sigma))*exp(-d/(2*sigma*sigma));
  a = exp(-d/(2*sigma*sigma));
  return a;
}
