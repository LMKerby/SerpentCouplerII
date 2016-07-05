/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processinterface.c                             */
/*                                                                           */
/* Created:       2012/02/14 (JLe)                                           */
/* Last modified: 2015/10/01 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Processes multi-physics interfaces                           */
/*                                                                           */
/* Comments: - Stats are allocated in allocinterfacestat.c                   */
/*                                                                           */
/*           - Polttoaineinterfacen aksiaalijako lisätty 3.4.2013            */ 
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessInterface:"

/*****************************************************************************/

void ProcessInterface(long update)
{
  long loc0, type, mat, mat0;

  /* Check pointer */

  if ((loc0 = (long)RDB[DATA_PTR_IFC0]) < VALID_PTR)
    return;

  fprintf(out, "Processing multi-physics interfaces...\n");

  /* Loop over interfaces */

  while (loc0 > VALID_PTR)
    {

      /* Get interface type */

      type = (long)RDB[loc0 + IFC_TYPE];

      /***** Processing for each interface type ************/

      switch(type)
	{
	case IFC_TYPE_PT_AVG:
	  
	  /***** Average of points *************************/

	  ProcessIFCPtAvg(loc0, update);

	  break;
	  /*************************************************/

	case IFC_TYPE_REG_MESH:
	  
	  /***** Regular mesh based distribution ***********/

	  ProcessIFCRegMesh(loc0, update);

	  break;
	  /*************************************************/

	case IFC_TYPE_FUNC:
	  
	  /***** Function based interface ******************/

	  ProcessIFCFunc(loc0, update);
	  
	  break;
	  /*************************************************/

	case IFC_TYPE_TET_MESH:
	  
	  /***** Unstructured tetrahedral mesh *************/
	  /* Can be OpenFOAM based */

	  ProcessIFCTetMesh(loc0, update);

	  break;
	  /*************************************************/

	case IFC_TYPE_FUEP:
	case IFC_TYPE_FPIP:
	  
	  /***** Fuel behavior interface *******************/

	  ProcessIFCFB(loc0, update);

	  break;
	  /*************************************************/

	default:
	  Die(FUNCTION_NAME, 
	      "Unknown interface type %ld in interface file: %s\n", 
	      type, GetText(loc0 + IFC_PTR_INPUT_FNAME));
      
	}

      /* Next interface */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Common stuff (JLe 1.10.2015 / 2.1.25) *******************************/

  /* Copy interface data to divided materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check pointer to parent */
      
      if ((mat0 = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) > VALID_PTR)
	{
	  /* Copy pointer */

	  WDB[mat + MATERIAL_PTR_IFC] = RDB[mat0 + MATERIAL_PTR_IFC];

	  /* Copy mode and limits */

	  WDB[mat + MATERIAL_TMS_TMIN] = RDB[mat0 + MATERIAL_TMS_TMIN];
	  WDB[mat + MATERIAL_TMS_TMAX] = RDB[mat0 + MATERIAL_TMS_TMAX];
	  WDB[mat + MATERIAL_TMS_MODE] = RDB[mat0 + MATERIAL_TMS_MODE];
	}

      /* Next material */

      mat = NextItem(mat);
    }  

  /***************************************************************************/

  fprintf(out, "OK.\n\n");

  /***************************************************************************/
}

/*****************************************************************************/
