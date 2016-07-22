#ifndef RUNSERPENT_H
#define RUNSERPENT_H

// This is a general user object, so that we can do anything

#include "GeneralUserObject.h"
#include "ElementTransfer.h"

class RunSerpent;

template<>
InputParameters validParams<RunSerpent>();

// This class just runs serpent

class RunSerpent : public GeneralUserObject
{
public:
  // Constructor
  RunSerpent(const InputParameters & parameters);

  /**
   * This is called before execute so you can reset any internal data.
   */
  virtual void initialize();

  /**
   * Called on every "object" (like every element or node).
   * In this case, it is called at every quadrature point on every element.
   */
  virtual void execute();

  /**
   * Called _once_ after execute has been called for all "objects".
   */
  virtual void finalize();

/* LMK, added declarations for Cmain and UpdateInterface, 7/2016 */
//  int Cmain(int, char**);
//  void UpdateInterface();

private:

  int _initialized, _run;
  const ElementTransfer & _element_transfer;
};

#endif
