#ifndef SERPENTCOUPLERIIAPP_H
#define SERPENTCOUPLERIIAPP_H

#include "MooseApp.h"

class SerpentcoupleriiApp;

template<>
InputParameters validParams<SerpentcoupleriiApp>();

class SerpentcoupleriiApp : public MooseApp
{
public:
  SerpentcoupleriiApp(InputParameters parameters);
  virtual ~SerpentcoupleriiApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
};

#endif /* SERPENTCOUPLERIIAPP_H */
