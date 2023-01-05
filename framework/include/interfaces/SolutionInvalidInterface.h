//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
#include "Moose.h"
#include "SolutionInvalidity.h"
#include "FEProblemBase.h"

// Forward declarations
class MooseObject;

#define flagInvalidSolution(message)                                                               \
  do                                                                                               \
  {                                                                                                \
    static const auto __invalid_id = this->registerInvalidSolutionInternal(message);               \
    this->flagInvalidSolutionInternal(__invalid_id);                                               \
  } while (0)

/**
 * An interface to communicate the solutioninvalidity for one app to the materials
 */
class SolutionInvalidInterface
{
public:
  /**
   * A class to help communicating the solutioninvalidity for one app to the materials
   */
  SolutionInvalidInterface(MooseObject * const moose_object);

  void flagInvalidSolutionInternal(InvalidSolutionID _invalid_solution_id);

  InvalidSolutionID registerInvalidSolutionInternal(const std::string & message) const;

private:
  /// The MooseObject that owns this interface
  MooseObject & _si_moose_object;

  /// A reference to FEProblem base
  FEProblemBase & _si_problem;
};
