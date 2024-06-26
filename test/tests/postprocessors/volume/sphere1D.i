# The volume of each block should be 3

[Mesh]
  file = sphere1D.e
  coord_type = RSPHERICAL
[]

[Functions]
  [fred]
    type = ParsedFunction
    expression = '200'
  []
[]

[AuxVariables]
  [constantVar]
    order = FIRST
    family = LAGRANGE
  []
[]

[Variables]
  [u]
    order = FIRST
    family = LAGRANGE
    initial_condition = 100
  []
[]

[AuxKernels]
  [fred]
    type = ConstantAux
    variable = constantVar
    block = 1
    value = 1
  []
[]

[ICs]
  [ic1]
    type = ConstantIC
    variable = constantVar
    value = 1
    block = 1
  []
[]

[Kernels]
  [heat_r]
    type = Diffusion
    variable = u
  []
[]

[BCs]
  [temps]
    type = FunctionDirichletBC
    variable = u
    boundary = 1
    function = fred
  []
[]

[Materials]
[]

[Executioner]
  type = Steady
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -snes_ls -ksp_gmres_restart'
  petsc_options_value = 'lu       basic                 101'
  line_search = 'none'
  nl_abs_tol = 1e-11
  nl_rel_tol = 1e-10
  l_max_its = 20
[]

[Postprocessors]
  [should_be_one]
    type = ElementAverageValue
    block = 1
    variable = constantVar
    execute_on = 'initial timestep_end'
  []
  [volume1]
    type = VolumePostprocessor
    block = 1
    execute_on = 'initial timestep_end'
  []
  [volume2]
    type = VolumePostprocessor
    block = 2
    execute_on = 'initial timestep_end'
  []
  [volume3]
    type = VolumePostprocessor
    block = 3
    execute_on = 'initial timestep_end'
  []
[]

[Outputs]
  exodus = true
[]
