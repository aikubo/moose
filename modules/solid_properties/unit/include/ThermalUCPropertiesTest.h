#pragma once

#include "MooseObjectUnitTest.h"
#include "ThermalUCProperties.h"

class ThermalUCPropertiesTest : public MooseObjectUnitTest
{
public:
  ThermalUCPropertiesTest() : MooseObjectUnitTest("SolidPropertiesApp") { buildObjects(); }

protected:
  void buildObjects()
  {
    InputParameters uo_pars1 = _factory.getValidParams("ThermalUCProperties");
    _fe_problem->addUserObject("ThermalUCProperties", "sp1", uo_pars1);
    _sp1 = &_fe_problem->getUserObject<ThermalUCProperties>("sp1");

    InputParameters uo_pars2 = _factory.getValidParams("ThermalUCProperties");
    uo_pars2.set<Real>("density") = 13000.0;
    _fe_problem->addUserObject("ThermalUCProperties", "sp2", uo_pars2);
    _sp2 = &_fe_problem->getUserObject<ThermalUCProperties>("sp2");
  }

  const ThermalUCProperties * _sp1;

  // model using a non-default density
  const ThermalUCProperties * _sp2;
};
