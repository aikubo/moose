[Mesh]
  [rmp]
    type = ReactorMeshParams
    dim = 3
    geom = "Square"
    assembly_pitch = 2.84126
    radial_boundary_id = 200
    axial_regions = '1.0'
    axial_mesh_intervals = '1'
    top_boundary_id = 201
    bottom_boundary_id = 202
  []

  [pin1]
    type = PinMeshGenerator
    reactor_params = rmp
    pin_type = 1
    pitch = 1.42063
    region_ids='1 2 5'
    quad_center_elements = true
    num_sectors = 2
    ring_radii = 0.2
    duct_halfpitch = 0.68
    mesh_intervals = '1 1 1'
  []

  [pin2]
    type = PinMeshGenerator
    reactor_params = rmp
    pin_type = 2
    pitch = 1.42063
    region_ids='2'
    quad_center_elements = true
    num_sectors = 2
    mesh_intervals = '2'
  []

  [pin3]
    type = PinMeshGenerator
    reactor_params = rmp
    pin_type = 3
    pitch = 1.42063
    region_ids='3 4'
    quad_center_elements = true
    num_sectors = 2
    ring_radii = 0.3818
    mesh_intervals = '1 1'
  []

  [amg1]
    type = AssemblyMeshGenerator
    assembly_type = 1
    inputs = 'pin2'
    pattern = '0 0;
               0 0'
  []

  [amg2]
    type = AssemblyMeshGenerator
    assembly_type = 2
    inputs = 'pin3 pin1 pin2'
    pattern = '0 1;
               1 2'
  []

  [cmg]
    type = CoreMeshGenerator
    inputs = 'amg2 amg1 empty'
    dummy_assembly_name = empty
    pattern = '1 0;
               0 1'
    extrude = true
  []
  [test_rgmb]
    type = TestReactorGeometryMeshBuilderMeshGenerator
    input = cmg
  []
  [transform]
    type = TransformGenerator
    input = test_rgmb
    transform = scale
    vector_value = '1 1 1'
  []
  data_driven_generator = test_rgmb
[]

[Problem]
  solve = false
[]

[Executioner]
  type = Steady
[]

[Reporters/metadata]
  type = MeshMetaDataReporter
[]

[Outputs]
  [out]
    type = JSON
    execute_on = FINAL
    execute_system_information_on = none
  []
[]
