#ifndef L2L3res
#define L2L3res

#include <cmath>
#include <iostream>

double getJEC_L2L3res(double jetPt)
{
  double correction = 1;
  //2014 Winter L2L3 Residual Correction
    const double x = jetPt;
    correction *= 1./(0.9773-0.0442*100./3.*(TMath::Max(0.0,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.0,1.03091-0.051154*pow(208.0,-0.154227)))-0.090*((-2.36997+0.413917*log(x))/x-(-2.36997+0.413917*log(208.0))/208));
    return jetPt*correction; 
}
#endif
