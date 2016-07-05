/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : otfsabscattering.c                             */
/*                                                                           */
/* Created:       2015/03/23 (TVi)                                           */
/* Last modified: 2015/05/25 (TVi)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Handles OTF S(a,b) scattering                                */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "OTFSabScattering:"

/*****************************************************************************/

void OTFSabScattering(long rea, double *E, double *u, double *v, double *w, 
		   long id)
{
  double mu, mut[2], E0, r, d1, d2, a, f, rnd[2], d1e[2], d2e[2], r_old;
  long law, ptr, l0, l1, l2, erg;
  long nc2, mode, ne, ne2, i, j, k, l, type, nuc, n;
  long sab, iso, ncol, mt, rea0;

  /* Avoid compiler warning */

  k = -1;
  l = -1;
  nc2 = -1;
  r = 0.0;
  r_old = -1.0;

  /***************************************************************************/
  
  /***** Remember some values before the collision ***************************/
  
  /* Initial energy */
  
  E0 = *E;
  
  /* Check initial energy */
  
  CheckValue(FUNCTION_NAME, "E0", "", *E, ZERO, INFTY);
  
  /***************************************************************************/

  /***** Get initial values and check ****************************************/

  /* Get nuclide and SAB pointers */

  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
  nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE]; 
  sab = (long)RDB[nuc + NUCLIDE_PTR_SAB]; 
  CheckPointer(FUNCTION_NAME, "(sab)", DATA_ARRAY, sab); 

  /* mt from 2004 / 2002 to 1004 / 1002 */ 

  mt = (long)RDB[rea + REACTION_MT]-1000;

  /******* Get first of the 2 S(a,b) nuclides between which to interpolate ***/
  /* (stored previously in OTFSabXS() ) */

  /* Collision number */

  ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
  ncol = (long)GetPrivateData(ptr, id);
  
  f = TestValuePair(sab + SAB_PTR_PREV_FRAC, ncol, id);

  if( ( sab = (long)TestValuePair(sab + SAB_PTR_PREV_SAB1, ncol, id)) 
      < VALID_PTR)
    Die(FUNCTION_NAME, "Previous sab data not available for nuclide %s", 
	GetText(nuc + NUCLIDE_PTR_NAME));
  
  /* Pre-sample random numbers for later use (same random numbers are used for 
     both temperatures */

  for(i=0; i<2; i++)
    rnd[i]=RandF(id);
  
  /*****************************************************************************/

  /******** Calculate mu and E for both temperatures by interpolating **********/
  
  for(j=0; j<2; j++) {

    iso = (long)RDB[sab + SAB_PTR_ISO];
    CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);
    
    /******* Find reactions ****************************************************/
    
    rea0 = (long)RDB[iso + NUCLIDE_PTR_REA];
    while(rea0 > VALID_PTR){
      
      if((long)RDB[rea0 + REACTION_MT]==mt)
	break;
      
      rea0 = NextItem(rea0);
    }
    
    /* Check that valid reactions were found */
    
    if(rea0 < VALID_PTR)
      Die(FUNCTION_NAME, "S(a,b) reaction %ld not found", mt);
 
    /* Get pointers to energy distributions */
    
    erg = (long)RDB[rea0 + REACTION_PTR_ERG];
    CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);
    
    /* Get energy distribution type */
    
    law = (long)RDB[erg + ERG_LAW];
    
    /* Get pointers to data */
    
    l0 = (long)RDB[erg + ERG_PTR_DATA];
    CheckPointer(FUNCTION_NAME, "(l0)", DATA_ARRAY, l0);
        
    /* Get scattering mode */
    
    mode = (long)RDB[l0++];
    
    /* Get pointers to incident energy grid */
    
    ptr = (long)RDB[l0++];
    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
        
    /***************************************************************************/
    
    /***** Sample energy and scattering angle **********************************/
    
    if (law == 1002)
      {
	
	/***********************************************************************/
	
	/***** S(a,b) elastic scattering ***************************************/
	
	/* Get number of cosines */
	
	nc2 = (long)RDB[l0++];
	
	/* Check mode */
	
	if (mode == 4)
	  {
	    /* Exact treatment, get index to energy grid */
	    
	    if ((i = GridSearch(ptr, E0)) < 0)
	      {
		/* Energy below or above grid, get grid size */
		
		ne = (long)RDB[ptr + ENERGY_GRID_NE];
		
		/* Avoid compiler warning */
		
		i = -1;
		
		/* Check incident energy */
		
		if (E0 < RDB[ptr + ENERGY_GRID_EMIN])
		  i = 0;
		else if (E0 > RDB[ptr + ENERGY_GRID_EMAX])
		  i = ne - 2;
		else
		  Die(FUNCTION_NAME, "law 1002: i = %E", i);
	      }
	    
	    /* Pointers to data */
	    
	    l1 = l0;
	    l2 = l0 + i + 1;
	    
	    /* Sample */
	    	    
	    d1 = rnd[0]*RDB[l2];
	    
	    /* Search */
	    
	    while (l2 != l1 + 1)
	      {	    
		n = (long)((l1 + l2)/2.0);
		
		if (d1 < RDB[n])
		  l2 = n;
		else
		  l1 = n;		
	      }
	    
	    /* Pointer to incident energy grid data */
	    
	    ptr = (long)RDB[ptr + ENERGY_GRID_PTR_DATA];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	    
	    /* Pointer to value */
	    
	    ptr = ptr + l2 - l0 - 1;
	    
	    /* Calculate mu */
	    
	    mut[j] = 1.0 - 2.0*RDB[ptr]/E0;
	    	    
	  }
      else if (mode == 3)
	{

	  /* Discrete cosines */
	  
	  if ((r = GridFactor(ptr, E0, id)) < 0.0)
	    {
	      /* Get number of energies */
	      
	      ne = (long)RDB[ptr + ENERGY_GRID_NE];			      

	      /* Avoid compiler warning */

	      i = -1;
	      r = 0.0;

	      /* Check energy */
	      
	      if (E0 < RDB[ptr + ENERGY_GRID_EMIN])
		{
		  i = 0;
		  r = 0.0;
		}
	      else if (E0 > RDB[ptr + ENERGY_GRID_EMAX])
		{
		  i = ne - 2;
		  r = 1.0;
		}
	      else
		Die(FUNCTION_NAME, "wtf?");
	    }
	  else
	    {
	      /* Separate integer and decimal parts */
	      
	      i = (long)r;
	      r = r - (double)i;
	    }
	  
	  /* Get pointers */
	  
	  l1 = l0 + i*nc2;
	  l2 = l0 + (i + 1)*nc2;
	  
	  /* Sample bin */

	  l = (long)(((double)nc2)*rnd[0]);
	  
	  /* Get cosine values */
	  
	  d1 = RDB[l1 + l];
	  d2 = RDB[l2 + l];
	  
	  /* Check values */
	  
	  CheckValue(FUNCTION_NAME, "d1 (law 1002)", "", d1, -1.0, 1.0);
	  CheckValue(FUNCTION_NAME, "d2 (law 1002)", "", d2, -1.0, 1.0);
	  
	  /* Interpolate mu value */
	  
	  mut[j] = d1 + r*(d2 - d1);
	}
      else
	Die(FUNCTION_NAME, "Invalid mode");

      /***********************************************************************/
      }
  else if (law == 1004)
    {
      
      /***** S(a,b) inelastic scattering *************************************/
      
      /* Get interpolation factor */
      /* HUOM! Tata voisi varmasti optimoida, koska energiagridit 
	 on lähes aina samat. Pitää vaan muistaa päivittää l0 -pointteria +=3 */

      if ((r = GridFactor(ptr, E0, id)) < 0.0)
	{
	  /* Get number of energies */
	  
	  ne = (long)RDB[ptr + ENERGY_GRID_NE];
	  
	  /* Avoid complier warning */
	  
	  i = -1;
	  
	  /* Check energy */
	  
	  if (E0 < RDB[ptr + ENERGY_GRID_EMIN])
	    {
	      i = 0;
	      r = 0.0;
	    }
	  else if (E0 > RDB[ptr + ENERGY_GRID_EMAX])
	    {
	      i = ne - 2;
	      r = 1.0;
	    }
	    else
	      Die(FUNCTION_NAME, "wtf?");
	  }
	else
	  {
	    /* Separate integer and decimal parts */
	    
	    i = (long)r;
	    r = r - (double)i;
	  }
	
      /* Check interpolation factor */
      
      CheckValue(FUNCTION_NAME, "r (law 1004)", "", r, 0.0, 1.0);
      
      /* Get number of secondary energies and cosines */
      
      ne2 = (long)RDB[l0++];
      nc2 = (long)RDB[l0++];
            
      /* Get secondary energy type and inelastic scattering mode */
      
      type = (long)RDB[l0++];
      
      /* Get index to energy bin (from sabcol-subroutine in MCNP source) */
      
      if (type == 2)
	Die(FUNCTION_NAME, 
	    "Continuous-energy S(a,b) treatment not yet implemented");
      else if (type == 0)
	k = (long)(rnd[0]*ne2);
      else if ((a = rnd[0]*(ne2 - 3.0)) > 1.0)
	k = (long)a + 1;
      else if (a > 0.6)
	k = ne2 - 2;
      else if (a > 0.5)
	k = ne2 - 1;
      else if (a > 0.1)
	k = 1;
      else
	k = 0;
	
      /* Check value */
	
      CheckValue(FUNCTION_NAME, "k", " (law 1004)", k, 0, ne2 - 1);
		
      /* Pointers to distributions */
	
      l1 = (long)RDB[l0 + i] + k*(nc2 + 1);
      l2 = (long)RDB[l0 + i + 1] + k*(nc2 + 1);
	
      /* Get energy values */
	
      d1e[j] = RDB[l1];
      d2e[j] = RDB[l2];
            
      /* Check cosine distribution mode */
      
      if (mode == 3)
	{
	  /* Sample from discrete cosines */
	  
	  l = (long)(((double)nc2)*rnd[1]);
	  
	  /* Get cosine values */
	  
	  d1 = RDB[l1 + 1 + l];
	  d2 = RDB[l2 + 1 + l];
	  
	  /* Check values */
	  
	  CheckValue(FUNCTION_NAME, "d1 (law 1004)", "", d1, -1.0, 1.0);
	  CheckValue(FUNCTION_NAME, "d2 (law 1004)", "", d2, -1.0, 1.0);
	  
	  /* Interpolate mu value */

	  mut[j] = d1 + r*(d2 - d1);
	  
	}
      else
	Die(FUNCTION_NAME, "Invalid inelastic scattering mode");
      
      /***********************************************************************/
    }
  else
    Die(FUNCTION_NAME, "Invalid S(a,b) mode");

    if(r_old > 0.0 && r_old != r) 
      Die(FUNCTION_NAME, "Energy grids differ in OTF S(a,b) interpolation (not yet implemented)"); 
      
    r_old = r; 
    sab = NextItem(sab);
  }

  /********* Interpolate mu and E in temperature *****************************/

  /* mu kay helposti */

  mu = (1-f)*mut[0] + f*mut[1];

  /* Monimutkaisehko lauseke energialle. Tässä oletetaan että r on sama 
     molemmille lämpötiloille (tämä tarkistetaan yllä). Käytännössä siis
     energiagridien on oltava samat. (onkohan ongelma ollenkaan Serpentissä?) */

  if(mt == 1004)
    *E = -(r-1)*d1e[0]*d1e[1]/(d1e[0]*f+d1e[1]*(1-f)) +
      r*d2e[0]*d2e[1]/(d2e[0]*f+d2e[1]*(1-f));

  /***************************************************************************/
  
  /* Check energy and cosines */

  CheckValue(FUNCTION_NAME, "E1", "", *E, ZERO, INFTY);
  CheckValue(FUNCTION_NAME, "r", "", *u**u+*v**v+*w**w - 1.0, -1E-5, 1E-5);

  /* Rotate direction cosines around a random azimuthal angle */
      
  AziRot(mu, u, v, w, id);      
}

/*****************************************************************************/

