#include <cmath>
#include <iostream>

double getJEC_2nd(double perp, double pseudorapidity, const char * mode)
{
  double perp_scale = 1;
  
  if(strcmp(mode,"pp2")==0) 
  {  
    #include "jec_2nd_pp_276.h"
  }
  else if(strcmp(mode,"pp7")==0)
  { 
    #include "jec_2nd_pp_7.h"
  }
  /*else if(strcmp(mode,"pPb5")==0 || strcmp(mode,"pp5")==0)
  { 
    #include "jec_2nd_ppb.h"
  }
  else if(strcmp(mode,"Pbp5")==0)
  { 
    #include "jec_2nd_pbp.h"
  }*/

  if(strcmp(mode,"pp7")==0)
  {
    return perp*perp_scale*1.01; 
  }
  else
  {
    return perp*perp_scale;
  }
}
