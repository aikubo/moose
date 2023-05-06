//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ScalarKernel.h"

#include <set>

// Forward Declarations
class GeneralizedPlaneStrainUOInterface;

class GeneralizedPlaneStrain : public ScalarKernel
{
public:
  static InputParameters validParams();

  GeneralizedPlaneStrain(const InputParameters & parameters);

  virtual void reinit(){};
  virtual void computeResidual();
  virtual void computeJacobian();

  const GeneralizedPlaneStrainUOInterface & _gps;
  const unsigned int _scalar_var_id;
  /// A pointer to the reference residual problem. This will be a nullptr if the problem type is not
  /// \p ReferenceResidualProblem
  const ReferenceResidualProblem * const _reference_residual_problem;
};
