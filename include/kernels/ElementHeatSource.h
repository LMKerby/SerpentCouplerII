#ifndef ELEMENTHEATSOURCE_H
#define ELEMENTHEATSOURCE_H

#include "Kernel.h"
#include "HeatToMoose.h"

//Forward Declarations
class ElementHeatSource;

template<>
InputParameters validParams<ElementHeatSource>();

class ElementHeatSource : public Kernel
{
public:

  ElementHeatSource(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  const HeatToMoose & _heat_to_moose;

};

#endif
