#include "HeatToMoose.h"
#include "MooseMesh.h"
#include "UserObject.h"
/* LMK #define OLD_WAY_HEATTOMOOSE */
template<>

InputParameters validParams<HeatToMoose>()
{
  InputParameters params = validParams<GeneralUserObject>();

  return params;
}

// Constructor
HeatToMoose::HeatToMoose(const InputParameters & parameters) :
  GeneralUserObject(parameters)
{
}

Real HeatToMoose::heatValue(const int element_idx) const
{
  /* Return the heat production in this element */

  return _value_array[element_idx];
}

void HeatToMoose::initialize()
{
  int numE = _subproblem.mesh().nElem();

   std::cout << "Initializing heattomoose\n";
  _value_array = new Real[numE];

  for(int i = 0 ; i < numE ; i++)
  {
    _value_array[i] = 0.0;
  }
   std::cout << "Done\n";

}

void HeatToMoose::execute()
{
  int num_e = _subproblem.mesh().nElem();
  int id1, idx, eId;
  int outIdx;
  int tot_num_e;
  Real vol, value, relerr;
  std::string str;
  std::map <int, int> elemIdMap;
  std::ifstream ifcfile;
  const Elem *element;

  std::cout << "Executing heattomoose\n";

  /* Old way */

#ifdef OLD_WAY_HEATTOMOOSE

  /* Open interface file for reading */

  ifcfile.open("mooseifc.out0", std::ios::in);

  if (ifcfile.is_open())
  {

    /* Loop over the interface file and get values */

    while ( ifcfile >> id1 >> idx >> value >> relerr )
    {
      /* Store heat production to valuearray */

      _value_array[idx-1] = value;

      //	  std::cout << "Putting value " << value << " to " << idx-1 << std::endl;
    }

  }

#else

  /* New way (2016) */

  /* We'll first need to create a map from Serpent output indices */
  /* to element indices (we did this the other way around in      */
  /* ElementTransfer) */

  /* Get total number of elements in the mesh */

  tot_num_e = _subproblem.mesh().nElem();

  /* Loop over all elements to do the mapping */

  outIdx = 0;

  for(int i = 0 ; i < tot_num_e ; i++)
    {
      /* Get next element */

      element = _subproblem.mesh().elemPtr(i);

      /* Only map active elements! */

      if(!element->active())
	continue;

      /* Map element id to the outIdx */

      elemIdMap[outIdx] = element->id();

      /* Increment output index */

      outIdx++;
    }

  /* Let's move on to reading the results */

  ifcfile.open("mooseifc_new.m", std::ios::in);

  if (ifcfile.is_open())
    {

      /* Skip first line */

      getline(ifcfile,str);

      /* Loop over the interface file and get values */

      while (ifcfile >> id1 >> idx >> vol >> value >> relerr )
	{
	  /* Get corresponding elementId */

	  eId = elemIdMap[idx-1];

	  /* Store heat production to valuearray */

	  _value_array[eId] = value;

	  //	  std::cout << "Putting value " << value << " to " << eId << std::endl;
	}

    }

  ifcfile.close();

#endif

}

void HeatToMoose::finalize()
{

}
