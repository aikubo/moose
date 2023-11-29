//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Moose.h"
#include "MooseInit.h"

namespace moose
{

/**
 * Create a MooseApp from command-line arguments.
 */
std::shared_ptr<MooseApp>
createMooseApp(const std::string & default_app_name, int argc, char * argv[]);

/**
 * Initialize, create and run a MooseApp
 */
template <typename DefaultAppType>
void
main(int argc, char * argv[])
{

  MooseInit init(argc, argv);

  DefaultAppType::registerApps();

  const auto default_app_name = MooseUtils::prettyCppType<DefaultAppType>();
  auto app = createMooseApp(default_app_name, argc, argv);

  app->run();
}

}
