#include "RunSerpent.h"
#include "header.h"
#include "locations.h"
#include "ElementTransfer.h"

template<>

InputParameters validParams<RunSerpent>()
{
  InputParameters params = validParams<GeneralUserObject>();
  params.addRequiredParam<UserObjectName>("transfer_user_object", "The name of the user object providing data transfer from MOOSE to Serpent.");

  return params;
}

// Constructor
RunSerpent::RunSerpent(const InputParameters & parameters) :
  GeneralUserObject(parameters),
  _element_transfer(getUserObject<ElementTransfer>("transfer_user_object"))
{
  _initialized = 0;
  _run = 0;
}

void RunSerpent::initialize()
{
  char **argumentti;
  int numarg = 4;

  if(_initialized == 0)
    {
      _element_transfer.printNodes();

      std::cout << "Allocating argumentti\n";

      argumentti = new char*[5];
      argumentti[0] = new char[80];
      argumentti[1] = new char[80];
      argumentti[2] = new char[80];
      argumentti[3] = new char[80];
      argumentti[4] = new char[80];
      //      argumentti[3] = new char[80];

      sprintf(argumentti[0],"input");
      sprintf(argumentti[2],"-omp");
      sprintf(argumentti[3],"3");
      //sprintf(argumentti[4],"-ext");
      sprintf(argumentti[1],"input");
      //      sprintf(argumentti[3],"-plot");
      std::cout << "Calling Cmain" << std::endl;
      Cmain(numarg,argumentti);
      std::cout << "Done\n";
      _initialized = 1;
    }
  else
    {
      _initialized = 1;
      std::cout << "Already initialized!" << std::endl;
    }
}

void RunSerpent::execute()
{

  std::cout << "At RunSerpent::execute()\n";

  /* Check mode and start calculation */

  if ((long)RDB[DATA_PTR_RIA0] > VALID_PTR)
    {
      /* RIA simulation */

      RIACycle();
    }
  else if ((long)RDB[DATA_BURNUP_CALCULATION_MODE] == YES)
    {
      /* Internal burnup mode, go to burnup cycle */

      BurnupCycle();
      /*
	UCBBurnupCycle();
      */
    }
  else
    {
      /* Single transport cycle */
      _element_transfer.printNodes();

//      UpdateInterface();

      PrepareTransportCycle();
      TransportCycle();
    }



}

void RunSerpent::finalize()
{

  if(_run==1)
  {
    /* Free memory */

    FreeMem();

    /* Finalize MPI */

    FinalizeMPI();
  }
}
