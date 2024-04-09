//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "NavierStokesPhysicsBase.h"
#include "WCNSFVCoupledAdvectionPhysicsHelper.h"

/**
 * Creates all the objects needed to solve the Navier Stokes scalar advection equations
 */
class WCNSFVScalarAdvectionPhysics final : public NavierStokesPhysicsBase,
                                           public WCNSFVCoupledAdvectionPhysicsHelper
{
public:
  static InputParameters validParams();

  WCNSFVScalarAdvectionPhysics(const InputParameters & parameters);

  /// Get the names of the advected scalar quantity variables
  const std::vector<NonlinearVariableName> & getAdvectedScalarNames() const
  {
    return _passive_scalar_names;
  }

  /// Whether the physics is actually creating the scalar advection equations
  bool hasScalarEquations() const { return _has_scalar_equation; }

protected:
private:
  void addNonlinearVariables() override;
  void addInitialConditions() override;
  void addFVKernels() override;
  void addFVBCs() override;

  unsigned short getNumberAlgebraicGhostingLayersNeeded() const override;

  // /**
  //  * Functions adding kernels for the incompressible / weakly compressible scalar advection
  //  * equation
  //  * If the material properties are not constant, some of these can be used for
  //  * weakly-compressible simulations as well.
  //  */
  void addScalarTimeKernels();
  void addScalarDiffusionKernels();
  void addScalarAdvectionKernels();
  /// Equivalent of NSFVAction addScalarCoupledSourceKernels
  void addScalarSourceKernels();

  /// Functions adding boundary conditions for the incompressible simulation.
  /// These are used for weakly-compressible simulations as well.
  void addScalarInletBC();
  void addScalarWallBC();

  /// Names of the passive scalar variables
  std::vector<NonlinearVariableName> _passive_scalar_names;
  /// A boolean to help compatibility with the old Modules/NavierStokesFV syntax
  const bool _has_scalar_equation;

  /// Functors for the passive scalar sources. Indexing is scalar variable index
  std::vector<MooseFunctorName> _passive_scalar_sources;
  /// Functors for the passive scalar (coupled) sources. Inner indexing is scalar variable index
  std::vector<std::vector<MooseFunctorName>> _passive_scalar_coupled_sources;
  /// Coefficients multiplying for the passive scalar sources. Inner indexing is scalar variable index
  std::vector<std::vector<Real>> _passive_scalar_sources_coef;

  /// Passive scalar inlet boundary types
  MultiMooseEnum _passive_scalar_inlet_types;
  /// Functors describing the inlet boundary values. See passive_scalar_inlet_types for what the functors actually represent
  std::vector<std::vector<MooseFunctorName>> _passive_scalar_inlet_functors;
};
