/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : updatetfbgeom.c                                */
/*                                                                           */
/* Created:       2012/08/20 (VVa)                                           */
/* Last modified: 2012/08/21 (JLe)                                           */
/* Version:       2.1.8                                                      */
/*                                                                           */
/* Description: Updates the temperature dependent geometry dimensions for    */
/*              the temperature feedback calculations                        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "UpdateTFBgeom:"

void UpdateTFBgeom(long tfb)
{
  long reg, mat, mattype;
  int i;
  double r0, r1, c0, c1, c2, N, h, simpsonsum, r, gap0, gap1;
  double Nfuel, reloc, A, B, C, D, P, T, cum;

  /* Check mode */

  if ((long)RDB[tfb + TFB_MODE] != 1)
    return;

  gap0=0.0;
  gap1=0.0;
  cum=0.0;

  /* Loop over regions to calculate thermal expansion */
  reg = (long)RDB[tfb + TFB_PTR_REG_LIST];
  CheckPointer(FUNCTION_NAME, "(reg)", DATA_ARRAY, reg);

  /* This is used to calculate the number of fuel regions    */
  /* so that the relocation can be divided evenly among them */
  Nfuel=0;
  while(reg > VALID_PTR)
  {

    /* Get pointer to the material block */
    mat = (long)RDB[reg + TFB_REG_PTR_MAT];
    CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

    /* Get material type */
    mattype = RDB[reg + TFB_REG_MAT_TYPE];

    /* Read basic data concerning this region */
    r0=RDB[reg + TFB_REG_R0];
    r1=RDB[reg + TFB_REG_R1];
    c0=RDB[reg + TFB_REG_ITER_C0];
    c1=RDB[reg + TFB_REG_ITER_C1];
    c2=RDB[reg + TFB_REG_ITER_C2];

    /* Calculate the outer radius */

    /* Gas regions are a special case, since their dimensions */
    /* change according to the adjacent regions               */
    if(RDB[mat+MATERIAL_MDENS] < 0.01)
    {

       /* We'll calculate the change of the outer radius     */
       /* According to the material type of the next region  */
       if(NextItem(reg) > VALID_PTR)
       {
         mattype = RDB[NextItem(reg) + TFB_REG_MAT_TYPE];
         r0=r1; /* Expansion should not be integrated over gas */
       }

       /* Store the gas gap radii for relocation calculation */
       if(r0 != 0.0)
       {
        gap0=RDB[reg + TFB_REG_R0];
        gap1=RDB[reg + TFB_REG_R1];           
       }

       /* Reset the cumulative expansion, since it is not  */
       /* transmitted through compressable phases */
       cum=0.0;
    }

    /* Calculate the outer radius for fuel rings */
    if(mattype == TFB_MAT_TYPE_FUEL )
    { 

      /* FRAP-model */
      A=9.8E-6;
      B=2.61E-3;
      C=3.16E-1;
      D=1.32;

      /* Update the outer radius for gas rings         */
      /* based on the expansion of the inner radius of */
      /* the adjacent fuel ring */
      if(r1==r0)
      {
        T = c0*r1*r1+c1*log(r1)+c2;
        cum = cum + A*T-B+C*exp(-D/(1.380649E-4*T));
        WDB[reg + TFB_REG_ITER_R1] = RDB[reg + TFB_REG_R1] + 
           cum;
      }
      else /* Fuel rings*/
      {
        /* Integration by the extended Simpson's rule */
        N=2000;
        h=(r1-r0)/N;
        simpsonsum=0;
        r=r0;



        for(i=0; i<=N; i++){

          /* Calculate the temperature at this spot */
          if(r==0) T=c0*r*r+c2;
          else T = c0*r*r+c1*log(r)+c2;

          /* Calculate the next term in the Simpson's sum*/
          if(i==0 || i==N){
            simpsonsum+=1*(A*T-B+C*exp(-D/(1.380649E-4*T)));
          }
          else if((i%2)==1){
            simpsonsum+=4*(A*T-B+C*exp(-D/(1.380649E-4*T)));
          }
          else {
            simpsonsum+=2*(A*T-B+C*exp(-D/(1.380649E-4*T)));
          }

          r+=h;
        }

        simpsonsum=simpsonsum*h/3;

        /* Calculate the cumulative transition*/
        cum = cum + simpsonsum;
        WDB[reg + TFB_REG_ITER_R1] = RDB[reg + TFB_REG_R1] + cum;

        Nfuel = Nfuel + 1;
      }
    }

    /* Calculate the outer radius for zircaloy rings */
    if(mattype == TFB_MAT_TYPE_ZIRCALOY)
    { 
      /* MATPRO-09 model */
      A=6.721E-6;
      B=-2.0163E-3;

      /* This handles automatically gas rings and cumulative expansion */
      WDB[reg + TFB_REG_ITER_R1] = RDB[reg + TFB_REG_R1] +
         A*c1*(r1*(log(r1)-1)+r0)+A*c2*r1+B*r1;
    }

    /* Update the inner radius */

    if(r0 != 0.0)
      WDB[reg + TFB_REG_ITER_R0] = RDB[PrevItem(reg) + TFB_REG_ITER_R1];
    else
      WDB[reg + TFB_REG_ITER_R0] = r0;



    reg = NextItem(reg);
  }


  /* Calculate total pellet relocation */

  P = RDB[tfb + TFB_ITER_POW];

  if(P<200)
    reloc = (gap1-gap0)*0.3;
  else if(P<400)
    reloc = (gap1-gap0)*(0.28+(P-200)/200*0.05);
  else
    reloc = (gap1-gap0)*0.32;

  /* Initialize the cumulative relocation variable */
  cum = 0.0;

  /* Loop over the regions to add the effect of relocation */

  reg = (long)RDB[tfb + TFB_PTR_REG_LIST];
  CheckPointer(FUNCTION_NAME, "(reg)", DATA_ARRAY, reg);

  while(reg > VALID_PTR)
  {

    /* Get material type */
    mattype = RDB[reg + TFB_REG_MAT_TYPE];

    /* Update the outer radius of fuel rings */
    if(mattype == TFB_MAT_TYPE_FUEL)
    { 
      /* Update the cumulative relocation */
      cum=cum + reloc/Nfuel;

      /* Store the new outer radius*/
      WDB[reg + TFB_REG_ITER_R1] = RDB[reg + TFB_REG_ITER_R1] + cum;
    }


    /* Update the inner radius */
    if(RDB[reg + TFB_REG_R0] != 0.0)
    {
      WDB[reg + TFB_REG_ITER_R0] = RDB[PrevItem(reg) + TFB_REG_ITER_R1];
    }

    /* Do not allow closing of the gas regions */
    if(RDB[reg + TFB_REG_ITER_R0] >= RDB[reg + TFB_REG_ITER_R1])
    {
      WDB[reg + TFB_REG_ITER_R0] = RDB[reg + TFB_REG_ITER_R1] - 2E-4;
      WDB[PrevItem(reg) + TFB_REG_ITER_R1] = RDB[reg + TFB_REG_ITER_R0];
    }

    reg = NextItem(reg);
  }

  return;
}
