#ifndef HEATTOMOOSE_H
#define HEATTOMOOSE_H

// This is a general user object, so that we can do anything

#include "GeneralUserObject.h"

class HeatToMoose;

template<>
InputParameters validParams<HeatToMoose>();

// This class just runs serpent

class HeatToMoose : public GeneralUserObject
{
public:
  // Constructor
  HeatToMoose(const InputParameters & parameters);

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

  Real heatValue(const int) const;

private:

  Real *_value_array;
};

#endif
