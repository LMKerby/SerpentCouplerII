[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 2
  ny = 2
  nz = 2
  xmin = -5
  xmax = 5
  ymin = -5
  ymax = 5
  zmin = -5
  zmax = 5
  uniform_refine = 3
  elem_type = PRISM6
[]

[Variables]
  [./u]
    initial_condition = 600
  [../]
[]

[Kernels]
  [./diff]
    type = Diffusion
    variable = u
  [../]
  [./elemsrc]
    type = ElementHeatSource
    to_moose_object = heatin
    variable = u
  [../]
[]

[BCs]
  [./top]
    type = DirichletBC
    variable = u
    boundary = top
    value = 600
  [../]
  [./bot]
    type = DirichletBC
    variable = u
    boundary = bottom
    value = 600
  [../]
  [./left]
    type = DirichletBC
    variable = u
    boundary = left
    value = 600
  [../]
  [./right]
    type = DirichletBC
    variable = u
    boundary = right
    value = 600
  [../]
  [./front]
    type = DirichletBC
    variable = u
    boundary = front
    value = 600
  [../]
  [./back]
    type = DirichletBC
    variable = u
    boundary = back
    value = 600
  [../]
[]

[UserObjects]
  [./elemtrans]
    type = ElementTransfer
    variable = u
    block = 0
  [../]
  [./runserpent]
    type = RunSerpent
    execute_on = timestep_end
    transfer_user_object = elemtrans
  [../]
  [./heatin]
    type = HeatToMoose
    execute_on = timestep_begin
  [../]
[]

[Executioner]
  # Preconditioned JFNK (default)
  type = Transient
  num_steps = 5
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

