#include "ElementHeatSource.h"

template<>
InputParameters validParams<ElementHeatSource>()
{
  InputParameters params = validParams<Kernel>();
  params.set<Real>("value")=0.0;
  params.addRequiredParam<UserObjectName>("to_moose_object", "An user object that brings element-wise data to MOOSE");
  return params;
}

ElementHeatSource::ElementHeatSource(const InputParameters & parameters) :
    Kernel(parameters),
    _heat_to_moose(getUserObject<HeatToMoose>("to_moose_object"))
{
}

Real ElementHeatSource::computeQpResidual()
{
  return -_test[_i][_qp]*_heat_to_moose.heatValue(_current_elem->id())/(_current_elem->volume());
}
