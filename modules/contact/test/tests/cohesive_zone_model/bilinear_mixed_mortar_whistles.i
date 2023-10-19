[Mesh]
  [base]
    type = GeneratedMeshGenerator
    dim = 2
    xmax = 1.1
    ymax = 1
    xmin = -0.1
    nx = 1
    ny = 1
  []
  [rename_base]
    type = RenameBoundaryGenerator
    input = base
    old_boundary = 'top bottom left right'
    new_boundary = 'top_base bottom_base left_base right_base'
  []
  [base_id]
    type = SubdomainIDGenerator
    input = rename_base
    subdomain_id = 1
  []

  [top]
    type = GeneratedMeshGenerator
    dim = 2
    xmax = 1
    ymin = 1
    ymax = 2
    nx = 1
    ny = 1
  []
  [rename_top]
    type = RenameBoundaryGenerator
    input = top
    old_boundary = 'top bottom left right'
    new_boundary = '100 101 102 103'
  []
  [top_id]
    type = SubdomainIDGenerator
    input = rename_top
    subdomain_id = 2
  []
  [combined]
    type = MeshCollectionGenerator
    inputs = 'base_id top_id'
  []
  [top_node]
    type = ExtraNodesetGenerator
    coord = '0 2 0'
    input = combined
    new_boundary = top_node
  []
  [bottom_node]
    type = ExtraNodesetGenerator
    coord = '-0.1 0 0'
    input = top_node
    new_boundary = bottom_node
  []

  [secondary]
    type = LowerDBlockFromSidesetGenerator
    new_block_id = 10001
    new_block_name = 'secondary_lower'
    sidesets = 'top_base'
    input = bottom_node
  []
  [primary]
    type = LowerDBlockFromSidesetGenerator
    new_block_id = 10000
    sidesets = '101'
    new_block_name = 'primary_lower'
    input = secondary
  []

  #  [nodeset_bottom]
  # type = ExtraNodesetGenerator
  #
  #    []
  #  [sideset_bottom]
  # type = SideSetsFromNodeSetsGenerator
  #
  #    []
  # [secondary]
  #   type = LowerDBlockFromSidesetGenerator
  #   new_block_id = 10001
  #   new_block_name = 'secondary_lower'
  #   sidesets = '10'
  #   input = input_file
  # []
  # [primary]
  #   type = LowerDBlockFromSidesetGenerator
  #   new_block_id = 10000
  #   sidesets = '20'
  #   new_block_name = 'primary_lower'
  #   input = secondary
  # []
  patch_update_strategy = auto
  patch_size = 20
  allow_renumbering = false
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Modules]
  [TensorMechanics]
    [Master]
      generate_output = 'stress_yy'
      [all]
        strain = FINITE
        add_variables = true
        use_automatic_differentiation = true
        decomposition_method = TaylorExpansion
        generate_output = 'vonmises_stress'
      []
    []
  []
[]

[BCs]
  [fix_x]
    type = DirichletBC
    preset = true
    value = 0.0
    boundary = bottom_node
    variable = disp_x
  []

  [fix_top]
    type = DirichletBC
    preset = true
    boundary = 100
    variable = disp_x
    value = 0
  []

  [top]
    type = FunctionDirichletBC
    boundary = 100
    variable = disp_y
    function = 'if(t<=0.3,t,if(t<=0.6,0.3-(t-0.3),0.6-t))'
    preset = true
  []

  [bottom]
    type = DirichletBC
    boundary = bottom_base
    variable = disp_y
    value = 0
    preset = true
  []
[]

[Materials]
  [stress]
    type = ADComputeFiniteStrainElasticStress
  []
  [elasticity_tensor]
    type = ADComputeElasticityTensor
    fill_method = symmetric9
    C_ijkl = '1.684e5 0.176e5 0.176e5 1.684e5 0.176e5 1.684e5 0.754e5 0.754e5 0.754e5'
  []
  # [czm]
  #   type = BiLinearMixedModeTraction
  #   boundary = 'interface'
  #   penalty_stiffness = 1e6
  #   GI_c = 1e3
  #   GII_c = 1e2
  #   normal_strength = 1e4
  #   shear_strength = 1e3
  #   displacements = 'disp_x disp_y'
  #   eta = 2.2
  #   viscosity = 1e-3
  # []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient

  solve_type = 'NEWTON'
  line_search = none

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  automatic_scaling = true

  l_max_its = 2
  l_tol = 1e-14
  nl_max_its = 150
  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-10
  start_time = 0.0
  dt = 0.1
  end_time = 1.0
  dtmin = 0.1
[]

[Outputs]
  exodus = true
[]


[UserObjects]
  [friction_uo]
    type = PenaltySimpleCohesiveZoneModel
    primary_boundary = 101
    secondary_boundary = 'top_base'
    primary_subdomain = 10000
    secondary_subdomain = 10001
    disp_x = disp_x
    disp_y = disp_y
    friction_coefficient = 0.1 # with 2.0 works
    secondary_variable = disp_x
    penalty = 1e7
    czm_normal_stiffness = 1e4
    penalty_friction = 1e4
    czm_tangential_stiffness = 1e5
    use_physical_gap = true
    # unused
    czm_normal_strength = 1e3
    czm_tangential_strength = 1e3

  []
[]

[Constraints]
  [x]
    type = NormalMortarMechanicalContact
    primary_boundary = 101
    secondary_boundary = 'top_base'
    primary_subdomain = 10000
    secondary_subdomain = 10001
    secondary_variable = disp_x
    component = x
    use_displaced_mesh = true
    compute_lm_residuals = false
    weighted_gap_uo = friction_uo
  []
  [y]
    type = NormalMortarMechanicalContact
    primary_boundary = 101
    secondary_boundary = 'top_base'
    primary_subdomain = 10000
    secondary_subdomain = 10001
    secondary_variable = disp_y
    component = y
    use_displaced_mesh = true
    compute_lm_residuals = false
    weighted_gap_uo = friction_uo
  []
  [t_x]
    type = TangentialMortarMechanicalContact
    primary_boundary = 101
    secondary_boundary = 'top_base'
    primary_subdomain = 10000
    secondary_subdomain = 10001
    secondary_variable = disp_x
    component = x
    use_displaced_mesh = true
    compute_lm_residuals = false
    weighted_velocities_uo = friction_uo
  []
  [t_y]
    type = TangentialMortarMechanicalContact
    primary_boundary = 101
    secondary_boundary = 'top_base'
    primary_subdomain = 10000
    secondary_subdomain = 10001
    secondary_variable = disp_y
    component = y
    use_displaced_mesh = true
    compute_lm_residuals = false
    weighted_velocities_uo = friction_uo
  []
[]
