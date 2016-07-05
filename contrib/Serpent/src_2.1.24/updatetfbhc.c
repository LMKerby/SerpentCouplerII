/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : iteratetfb.c                                   */
/*                                                                           */
/* Created:       2012/08/20 (VVa)                                           */
/* Last modified: 2015/05/06 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Updates the temperature dependent heat conductivities for    */
/*              the temperature feedback calculations                        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "UpdateTFBhc:"

void UpdateTFBhc(long tfb)
{
  long reg, mattype;
  int i;
  double r0, r1, c0, c1, c2, N, h, simpsonsum, T;
  double Tin, Tout, bu, gad;

  /* Check mode */

  if ((long)RDB[tfb + TFB_MODE] < 1)
    return;

  /* Loop over regions to calculate thermal expansion */
  reg = (long)RDB[tfb + TFB_PTR_REG_LIST];
  CheckPointer(FUNCTION_NAME, "(reg)", DATA_ARRAY, reg);

  while(reg > VALID_PTR)
  {

    /* Get material type */
    mattype = RDB[reg + TFB_REG_MAT_TYPE];

    /* Read basic data concerning this region */
    r0=RDB[reg + TFB_REG_R0];
    r1=RDB[reg + TFB_REG_R1];
    c0=RDB[reg + TFB_REG_ITER_C0];
    c1=RDB[reg + TFB_REG_ITER_C1];
    c2=RDB[reg + TFB_REG_ITER_C2];

    /* Calculate the previous boundary temperatures */
    if(r0==0) Tin=c0*r0*r0+c2;
    else Tin=c0*r0*r0+c1*log(r0)+c2;

    Tout=c0*r1*r1+c1*log(r1)+c2;


    /* Update heat conductivities */

    if(mattype == TFB_MAT_TYPE_FUEL)
    {
      /* At the moment, only trivial burnups */
      bu=0.0;

      /* Trivial gadolinium */
      gad=0.0;

      /* FRAP model */
      /* Integration by the extended Simpson's rule */
      if(abs(Tin-Tout)< 0.1)
        WDB[reg + TFB_REG_HC] = (1/(0.0452+1.1599*gad+2.46E-4
        *Tin+0.00187*bu+(1-0.9*exp(-0.04*bu))*
        0.038*pow(bu,0.28)*1/(1+396*exp(-6380/Tin)))+
        3.5E9/pow(Tin,2)*exp(-16361/Tin))*0.01;
      else{
        /* Integration by the extended Simpson's rule */
        simpsonsum=0;
        T=Tout;
        N=2000;
        h=(Tin-Tout)/N;

        for(i=0; i<=N; i++){

          /* Calculate the next term in the Simpson's sum*/
          if(i==0 || i==N){
            simpsonsum+=1*(1/(0.0452+1.1599*gad+2.46E-4*T+0.00187*bu+(1-0.9*exp(-0.04*bu))*
          0.038*pow(bu,0.28)*1/(1+396*exp(-6380/T)))+
          3.5E9/pow(T,2)*exp(-16361/T))*0.01;
          }
          else if((i%2)==1){
            simpsonsum+=4*(1/(0.0452+1.1599*gad+2.46E-4*T+0.00187*bu+(1-0.9*exp(-0.04*bu))*
          0.038*pow(bu,0.28)*1/(1+396*exp(-6380/T)))+
          3.5E9/pow(T,2)*exp(-16361/T))*0.01;
          }
          else {
            simpsonsum+=2*(1/(0.0452+1.1599*gad+2.46E-4*T+0.00187*bu+(1-0.9*exp(-0.04*bu))*
          0.038*pow(bu,0.28)*1/(1+396*exp(-6380/T)))+
          3.5E9/pow(T,2)*exp(-16361/T))*0.01;
          }
          T+=h;
        }

      simpsonsum=simpsonsum*h/3;
      WDB[reg + TFB_REG_HC] = simpsonsum/(Tin-Tout);
      }
    }

    if(mattype == TFB_MAT_TYPE_ZIRCALOY)
    {
      /* MATPRO model*/
      if (abs(Tin-Tout)>0.1)
        WDB[reg + TFB_REG_HC] = (7.67E-9/4.0*(pow(Tin,4)-pow(Tout,4))-
            1.45E-5/3.0*(pow(Tin,3)-pow(Tout,3))+
            2.09E-2/2.0*(pow(Tin,2)-pow(Tout,2))+
            7.51*(Tin-Tout))/(Tin-Tout)*0.01;
      else
        WDB[reg + TFB_REG_HC] = (7.51+0.0209*Tin-0.0000145*pow(Tin,2)+.00000000767*pow(Tin,3))*0.01;

    }

    if(mattype == TFB_MAT_TYPE_HELIUM)
    {
      /* FRAP model*/
      if(abs(Tin-Tout)>0.1)
        WDB[reg + TFB_REG_HC] = (1/(0.7146+1)*2.531E-3*pow(Tin,0.7146+1)-1/(0.7146+1)*2.531E-3*pow(Tout,0.7146+1))/(Tin-Tout)*0.01;
      else
        WDB[reg + TFB_REG_HC] = 2.531E-3*pow(Tin,0.7146)*0.01;
    }

    reg = NextItem(reg);
  }

  return;
}
