mu=.01
rho=1

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = .1
    ymin = 0
    ymax = .1
    nx = 20
    ny = 20
  []
[]

[Modules]
  [NavierStokesFV]
    simulation_type = 'steady-state'
    compressibility = 'incompressible'
    porous_medium_treatment = false
    add_energy_equation = false

    density = 'rho'
    dynamic_viscosity = 'mu'
    gravity = '0 0 0'

    inlet_boundaries = 'top'
    momentum_inlet_types = 'fixed-velocity'
    momentum_inlet_function = 'velocity_x_inlet velocity_y_inlet'

    wall_boundaries = 'left right bottom'
    momentum_wall_types = 'noslip noslip noslip'
    pin_pressure = true

    pinned_pressure_type = average
    pinned_pressure_value = 0
  []
[]

[AuxVariables]
  [U]
    order = CONSTANT
    family = MONOMIAL
    fv = true
  []
[]

[AuxKernels]
  [mag]
    type = VectorMagnitudeAux
    variable = U
    x = vel_x
    y = vel_y
  []
[]

[Materials]
  [const]
    type = ADGenericFunctorMaterial
    prop_names = 'rho mu'
    prop_values = '${rho} ${mu}'
  []
[]

[Functions]
  [velocity_x_inlet]
    type = ParsedFunction
    value = '1'
  []
  [velocity_y_inlet]
    type = ParsedFunction
    value = '0'
  []
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_pc_type -sub_pc_factor_shift_type'
  petsc_options_value = 'asm      100                lu           NONZERO'
  nl_rel_tol = 1e-12
[]

[Outputs]
  exodus = true
  file_base = 'lid-driven_out'
[]
