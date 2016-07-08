
double get7tevPt(double perp, double pseudorapidity, int nVtx)
{
  double perp_scale = 1;

  switch(nVtx){
          case 1:
          #include "jec_pp7tev_pu1.h"
          break;
          case 2:
          #include "jec_pp7tev_pu2.h"
          break;
          case 3:
          #include "jec_pp7tev_pu3.h"
          break;
          case 4:
          #include "jec_pp7tev_pu4.h"
          break;
          case 5:
          #include "jec_pp7tev_pu5.h"
          break;
          case 6:
          #include "jec_pp7tev_pu6.h"
          break;
          case 7:
          #include "jec_pp7tev_pu7.h"
          break;
          case 8:
          #include "jec_pp7tev_pu8.h"
          break;
          case 9:
          #include "jec_pp7tev_pu9.h"
          break;
          case 10:
          #include "jec_pp7tev_pu10.h"
          break;
          case 11:
          #include "jec_pp7tev_pu11.h"
          break;
          case 12:
          #include "jec_pp7tev_pu12.h"
          break;
        }

 return perp*perp_scale; 
}
