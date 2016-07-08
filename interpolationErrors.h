#include "TMath.h"

//x is number of gluon jets at 2 tev
//y is " at 5
//z is " at 7
//p2 is pp FF at 2 tev 
//p7 is pp FF at 7 tev
//a 'd' prefix means error
double getInterpolationError(double x, double dx, double y, double dy, double z, double dz, double p2, double dp2, double p7, double dp7)
{
  double errorDueToP2 = TMath::Power((y-z)*dp2/(x-z),2);
  double errorDueToP7 = TMath::Power((x-y)*dp7/(x-z),2);
  double errorDueToX  = TMath::Power((p2-p7)*(z-y)*dx/((x-z)*(x-z)),2);
  double errorDueToY  = TMath::Power((p2-p7)*dy/(x-z),2);
  double errorDueToZ  = TMath::Power((p2-p7)*(x-y)*dz/((x-z)*(x-z)),2);
  double total = TMath::Power(errorDueToP2 + errorDueToP7 + errorDueToX + errorDueToY + errorDueToZ,0.5);

  return total;
}
