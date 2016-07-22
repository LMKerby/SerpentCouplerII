#ifndef ELEMENTTRANSFER_H
#define ELEMENTTRANSFER_H

// This is an elemental user object
#include "ElementIntegralVariablePostprocessor.h"
#include <iostream>
#include <fstream>

// Libmesh includes
#include "libmesh/mesh_tools.h"
//#include "libmesh/elem_type.h"

class ElementTransfer;

template<>
InputParameters validParams<ElementTransfer>();

// This class computes averages on each element and transfers them to Serpent

class ElementTransfer : public ElementIntegralVariablePostprocessor
{
public:
  // Constructor
  ElementTransfer(const InputParameters & parameters);

  void printNodes() const;

  /**
   * This is called before execute so you can reset any internal data.
   */
  virtual void initialize();

  /**
   * Called on every "object" (like every element or node).
   * In this case, it is called at every quadrature point on every element.
   */
  virtual void execute();

  /* Called to join another user object like this */
  virtual void threadJoin(const UserObject & y);

  /**
   * Called _once_ after execute has been called for all "objects".
   */
  virtual void finalize();

protected:

  Real *_integral_value_array;
  
  Real *_volume_value_array;

  Real *_average_value_array;
};

#endif
