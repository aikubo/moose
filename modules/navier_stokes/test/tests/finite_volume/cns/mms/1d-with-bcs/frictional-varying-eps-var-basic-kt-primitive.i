[GlobalParams]
  fp = fp
  limiter = 'central_difference'
  two_term_boundary_expansion = true
[]

[Mesh]
  [cartesian]
    type = GeneratedMeshGenerator
    dim = 1
    xmin = .1
    xmax = .6
    nx = 2
  []
[]

[Modules]
  [FluidProperties]
    [fp]
      type = IdealGasFluidProperties
    []
  []
[]

[Problem]
  kernel_coverage_check = false
  fv_bcs_integrity_check = false
[]

[Variables]
  [pressure]
    type = MooseVariableFVReal
  []
  [sup_vel_x]
    type = MooseVariableFVReal
  []
  [T_fluid]
    type = MooseVariableFVReal
  []
[]

[AuxVariables]
  [eps]
    type = MooseVariableFVReal
  []
[]

[ICs]
  [pressure]
    type = FunctionIC
    variable = pressure
    function = 'exact_p'
  []
  [sup_vel_x]
    type = FunctionIC
    variable = sup_vel_x
    function = 'exact_sup_vel_x'
  []
  [T_fluid]
    type = FunctionIC
    variable = T_fluid
    function = 'exact_T'
  []
  [eps]
    type = FunctionIC
    variable = eps
    function = 'eps'
  []
[]

[FVKernels]
  [mass_advection]
    type = PCNSFVKT
    variable = pressure
    eqn = "mass"
  []
  [mass_fn]
    type = FVBodyForce
    variable = pressure
    function = 'forcing_rho'
  []

  [momentum_x_advection]
    type = PCNSFVKT
    variable = sup_vel_x
    momentum_component = x
    eqn = "momentum"
  []
  [eps_grad]
    type = PNSFVPGradEpsilon
    variable = sup_vel_x
    momentum_component = 'x'
    epsilon_var = 'eps'
  []
  [drag]
    type = PNSFVMomentumFriction
    variable = sup_vel_x
    momentum_component = 'x'
    linear_coef_name = 1
    momentum_name = momentum
  []
  [momentum_fn]
    type = FVBodyForce
    variable = sup_vel_x
    function = 'forcing_rho_ud'
  []

  [fluid_energy_advection]
    type = PCNSFVKT
    variable = T_fluid
    eqn = "energy"
  []
  [energy_fn]
    type = FVBodyForce
    variable = T_fluid
    function = 'forcing_rho_et'
  []
[]

[FVBCs]
  [mass_left]
    variable = pressure
    type = PCNSFVKTBC
    boundary = left
    T_fluid_function = 'exact_T'
    superficial_velocity_function = 'exact_superficial_velocity'
    eqn = 'mass'
  []
  [momentum_left]
    variable = sup_vel_x
    type = PCNSFVKTBC
    boundary = left
    T_fluid_function = 'exact_T'
    superficial_velocity_function = 'exact_superficial_velocity'
    eqn = 'momentum'
    momentum_component = 'x'
  []
  [energy_left]
    variable = T_fluid
    type = PCNSFVKTBC
    boundary = left
    T_fluid_function = 'exact_T'
    superficial_velocity_function = 'exact_superficial_velocity'
    eqn = 'energy'
  []
  [mass_right]
    variable = pressure
    type = PCNSFVKTBC
    boundary = right
    eqn = 'mass'
    pressure_function = 'exact_p'
  []
  [momentum_right]
    variable = sup_vel_x
    type = PCNSFVKTBC
    boundary = right
    eqn = 'momentum'
    momentum_component = 'x'
    pressure_function = 'exact_p'
  []
  [energy_right]
    variable = T_fluid
    type = PCNSFVKTBC
    boundary = right
    eqn = 'energy'
    pressure_function = 'exact_p'
  []

  # help gradient reconstruction
  [pressure_right]
    type = FVFunctionDirichletBC
    variable = pressure
    function = exact_p
    boundary = 'right'
  []
  [sup_vel_x_left]
    type = FVFunctionDirichletBC
    variable = sup_vel_x
    function = exact_sup_vel_x
    boundary = 'left'
  []
  [T_fluid_left]
    type = FVFunctionDirichletBC
    variable = T_fluid
    function = exact_T
    boundary = 'left'
  []
  [eps]
    type = FVFunctionDirichletBC
    variable = eps
    function = eps
    boundary = 'left right'
  []
[]

[Materials]
  [var_mat]
    type = PorousPrimitiveVarMaterial
    pressure = pressure
    superficial_vel_x = sup_vel_x
    T_fluid = T_fluid
    porosity = porosity
  []
  # [porosity]
  #   type = PorosityVarMaterial
  #   porosity_var = 'eps'
  # []
  [porosity]
    type = PorosityVarMaterial
    porosity_function = 'eps'
  []
[]

[Functions]
[exact_rho]
  type = ParsedFunction
  value = '3.48788261470924*cos(x)'
[]
[forcing_rho]
  type = ParsedFunction
  value = '-3.83667087618017*sin(1.1*x)*cos(1.3*x) - 4.53424739912202*sin(1.3*x)*cos(1.1*x)'
[]
[exact_rho_ud]
  type = ParsedFunction
  value = '3.48788261470924*cos(1.1*x)*cos(1.3*x)'
[]
[forcing_rho_ud]
  type = ParsedFunction
  value = '(-(10.6975765229419*cos(1.5*x)/cos(x) - 0.697576522941849*cos(1.1*x)^2/cos(x)^2)*sin(x) + (10.6975765229419*sin(x)*cos(1.5*x)/cos(x)^2 - 1.3951530458837*sin(x)*cos(1.1*x)^2/cos(x)^3 + 1.53466835047207*sin(1.1*x)*cos(1.1*x)/cos(x)^2 - 16.0463647844128*sin(1.5*x)/cos(x))*cos(x))*cos(1.3*x) + 3.48788261470924*sin(x)*cos(1.1*x)^2*cos(1.3*x)/cos(x)^2 - 7.67334175236034*sin(1.1*x)*cos(1.1*x)*cos(1.3*x)/cos(x) - 4.53424739912202*sin(1.3*x)*cos(1.1*x)^2/cos(x) + 3.48788261470924*cos(1.1*x)'
[]
[exact_rho_et]
  type = ParsedFunction
  value = '26.7439413073546*cos(1.5*x)'
[]
[forcing_rho_et]
  type = ParsedFunction
  value = '1.0*(3.48788261470924*(3.06706896551724*cos(1.5*x)/cos(x) - 0.2*cos(1.1*x)^2/cos(x)^2)*cos(x) + 26.7439413073546*cos(1.5*x))*sin(x)*cos(1.1*x)*cos(1.3*x)/cos(x)^2 - 1.1*(3.48788261470924*(3.06706896551724*cos(1.5*x)/cos(x) - 0.2*cos(1.1*x)^2/cos(x)^2)*cos(x) + 26.7439413073546*cos(1.5*x))*sin(1.1*x)*cos(1.3*x)/cos(x) - 1.3*(3.48788261470924*(3.06706896551724*cos(1.5*x)/cos(x) - 0.2*cos(1.1*x)^2/cos(x)^2)*cos(x) + 26.7439413073546*cos(1.5*x))*sin(1.3*x)*cos(1.1*x)/cos(x) + 1.0*(-(10.6975765229419*cos(1.5*x)/cos(x) - 0.697576522941849*cos(1.1*x)^2/cos(x)^2)*sin(x) + (10.6975765229419*sin(x)*cos(1.5*x)/cos(x)^2 - 1.3951530458837*sin(x)*cos(1.1*x)^2/cos(x)^3 + 1.53466835047207*sin(1.1*x)*cos(1.1*x)/cos(x)^2 - 16.0463647844128*sin(1.5*x)/cos(x))*cos(x) - 40.1159119610319*sin(1.5*x))*cos(1.1*x)*cos(1.3*x)/cos(x)'
[]
[exact_T]
  type = ParsedFunction
  value = '0.0106975765229418*cos(1.5*x)/cos(x) - 0.000697576522941848*cos(1.1*x)^2/cos(x)^2'
[]
[exact_eps_p]
  type = ParsedFunction
  value = '3.48788261470924*(3.06706896551724*cos(1.5*x)/cos(x) - 0.2*cos(1.1*x)^2/cos(x)^2)*cos(x)*cos(1.3*x)'
[]
[exact_p]
  type = ParsedFunction
  value = '3.48788261470924*(3.06706896551724*cos(1.5*x)/cos(x) - 0.2*cos(1.1*x)^2/cos(x)^2)*cos(x)'
[]
[exact_sup_vel_x]
  type = ParsedFunction
  value = '1.0*cos(1.1*x)*cos(1.3*x)/cos(x)'
[]
[eps]
  type = ParsedFunction
  value = 'cos(1.3*x)'
[]
[exact_superficial_velocity]
  type = ParsedVectorFunction
  value_x = '1.0*cos(1.1*x)*cos(1.3*x)/cos(x)'
[]
[]

[Executioner]
  solve_type = NEWTON
  type = Transient
  num_steps = 1
  dtmin = 1
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_max_its = 50
  line_search = bt
  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-12
[]

[Outputs]
  exodus = true
  csv = true
[]

[Debug]
  show_var_residual_norms = true
[]

[Postprocessors]
  [h]
    type = AverageElementSize
    outputs = 'console csv'
    execute_on = 'timestep_end'
  []
  [L2pressure]
    type = ElementL2Error
    variable = pressure
    function = exact_p
    outputs = 'console csv'
    execute_on = 'timestep_end'
  []
  [L2sup_vel_x]
    variable = sup_vel_x
    function = exact_sup_vel_x
    type = ElementL2Error
    outputs = 'console csv'
    execute_on = 'timestep_end'
  []
  [L2T_fluid]
    variable = T_fluid
    function = exact_T
    type = ElementL2Error
    outputs = 'console csv'
    execute_on = 'timestep_end'
  []
[]
