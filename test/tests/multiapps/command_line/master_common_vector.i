[Mesh]
  type = GeneratedMesh
  dim = 1
[]

[Problem]
  type = FEProblem
  solve = false
[]

[Executioner]
  type = Steady
[]

[MultiApps]
  [sub]
    type = FullSolveMultiApp
    positions = '0 0 0
                 1 1 1'
    input_files = 'sub.i'
    cli_args = 'Mesh/mesh/type=CartesianMeshGenerator;Mesh/mesh/dim=1;Mesh/mesh/dx="1 2 1";Mesh/mesh/ix="4 3 1"'
  []
[]
